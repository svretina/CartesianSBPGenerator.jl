using CartesianSBPGenerator
using CartesianSBPGenerator: Experimental
using LinearAlgebra
using SparseArrays

const DESIGNS = ((2, 5, 17), (4, 9, 25), (6, 15, 37))

function seed(p, depth, N)
    system = Experimental.assemble_stencil_experiment(p; N, period=2,
                 closure_depth=depth, stencil_radius=p)
    result = Experimental.search_stencil(system)
    result.status === :success || error("could not reproduce order $p seed")
    return system, result.operator
end

function operator_metrics(op, depth)
    h = (last(op.x) - first(op.x)) / (length(op.x) - 1)
    D = op.D
    rows = vcat(1:depth, (size(D,1)-depth+1):size(D,1))
    return (; depth, width=maximum(maximum(findall(!iszero,D[i,:])) -
                   minimum(findall(!iszero,D[i,:])) + 1 for i in rows),
            weight_extrema=Float64.(extrema(diag(op.H) / h)),
            norm_condition=cond(Float64.(op.H)),
            maximum_coefficient=Float64(h * maximum(abs, D)),
            nonzeros_per_row=extrema([count(!iszero,D[i,:]) for i in axes(D,1)]))
end

function dispersion(data)
    W = Diagonal(Float64.(collect(data.weights)))
    sqrtW = sqrt(W)
    invsqrtW = inv(sqrtW)
    thetas = range(1e-4, π/2; length=2049)
    rows = NamedTuple[]
    physical = 0
    previous = nothing
    maximum_skew_residual = 0.0
    zero_crossings = zeros(Int, 2)
    for (i, theta) in enumerate(thetas)
        M = Experimental.bloch_symbol(data, theta)
        Mp = Experimental.bloch_symbol(data, theta; derivative=true)
        maximum_skew_residual = max(maximum_skew_residual, norm(W * M + M' * W, Inf))
        eig = eigen(Hermitian(-im * sqrtW * M * invsqrtW))
        smooth = sqrtW * ones(2)
        energy = sum(abs2, smooth)
        fractions = [abs2(dot(eig.vectors[:, j], smooth)) / energy for j in 1:2]
        i == 1 && (physical = argmax(fractions))
        indices = (physical, 3 - physical)
        frequencies = Tuple(eig.values[j] for j in indices)
        group = Tuple(real(dot(eig.vectors[:, j], (-im * sqrtW * Mp * invsqrtW) *
                               eig.vectors[:, j])) for j in indices)
        if !isnothing(previous)
            for j in 1:2
                previous[j] * frequencies[j] < 0 && (zero_crossings[j] += 1)
            end
        end
        previous = frequencies
        push!(rows, (; theta, frequencies, phase=frequencies ./ theta, group,
                     folded_smooth_energy=fractions[3 - physical]))
    end
    conventional(s, theta) =
        sum(Float64(c) * sin(j * theta) for (c,j) in zip(coefficients(s),s.offsets))
    stencils = map(order -> centered_stencil(derivative_order=1, accuracy_order=order),
                   (data.order, data.order + 2))
    i = argmin(abs.([r.theta for r in rows] .- 0.1))
    low = rows[i]
    centered_radius(s) = maximum(abs(conventional(s, theta)) for theta in range(0,π;length=4097))
    return (; rows, physical, maximum_skew_residual, zero_crossings,
            zero_at_origin=data.parity_compatible,
            radius=maximum(maximum(abs, r.frequencies) for r in rows), low,
            edge=last(rows), max_folded_smooth_energy=maximum(r.folded_smooth_energy for r in rows),
            centered_low=map(s -> conventional(s, low.theta), stencils),
            centered_radius=map(centered_radius, stencils))
end

function truncation(op, depth)
    return setprecision(192) do
        D, x = sparse(BigFloat.(op.D)), BigFloat.(op.x)
        interior = (depth + 1):(length(x) - depth)
        edge = 1:min(2, length(interior))
        transition = unique(vcat(interior[edge], reverse(interior)[edge]))
        groups = (even=filter(isodd, interior), odd=filter(iseven, interior),
                  left=1:depth, right=(length(x)-depth+1):length(x), transition=transition)
        e = exp.(x)
        s, exact_s = sin.(2 .* BigFloat(π) .* x), 2 .* BigFloat(π) .* cos.(2 .* BigFloat(π) .* x)
        return map(residual -> map(indices -> Float64(maximum(abs, residual[indices])), groups),
                   (D * e - e, D * s - exact_s))
    end
end

rk4_amplification(z) = 1 + z * (1 + z * (1/2 + z * (1/6 + z/24)))

function rk4_limit(L)
    vals = eigvals(Matrix(L))
    rho = maximum(abs, vals)
    maximum(real, vals) <= 1e-10 * rho || error("SAT spectrum has positive real part")
    stable(dt) = maximum(x -> abs(rk4_amplification(dt*x)), vals) <= 1 + 1e-12
    lo, hi = 0.0, 4 / rho
    for _ in 1:48
        mid = (lo + hi) / 2
        if stable(mid)
            lo = mid
        else
            hi = mid
        end
    end
    return lo
end

function sat_matrices(op)
    Q = weighted_derivative_matrix(op)
    Q + Q' == boundary_matrix(op) || error("SAT benchmark requires an exactly certified SBP operator")
    x, H, D = Float64.(op.x), Float64.(op.H), Float64.(op.D)
    n = length(x)
    e0, e1 = zeros(n), zeros(n)
    e0[1], e1[end] = 1, 1
    left, right = H \ e0, H \ e1
    plus = -D - left * e0'
    minus = D - right * e1'
    wave = [plus zeros(n,n); right * e1' minus]
    wave_norm = [H zeros(n,n); zeros(n,n) H]
    for (mass, generator) in ((H, plus), (wave_norm,wave))
        residual = mass * generator
        allowance = 256eps(Float64) * max(1, norm(residual, Inf))
        maximum(eigvals(Symmetric(residual + residual'))) <= allowance ||
            error("SAT generator fails the semidiscrete energy inequality")
    end
    max_step = min(rk4_limit(plus), rk4_limit(wave))
    return (; x, H, plus, wave, left, max_step, wave_norm)
end

# For constant generators the powered RK4 update avoids millions of
# tiny explicit steps in the stiffest boundary closures.
function time_errors(generator, initial, tf, step)
    steps = ceil(Int, tf/step)
    dt = tf / steps
    Z = dt * generator
    base = I + Z + Z^2/2 + Z^3/6 + Z^4/24
    Zhalf = Z/2
    half = I + Zhalf + Zhalf^2/2 + Zhalf^3/6 + Zhalf^4/24
    u = exp(tf * generator) * initial
    half_exp = exp((tf/2) * generator)
    semigroup_error = maximum(abs, half_exp * (half_exp * initial) - u)
    coarse = base^steps * initial
    fine = half^(2steps) * initial
    return u, maximum(abs, coarse-u), maximum(abs, fine-u), semigroup_error
end

function pde_benchmark(sat; step)
    (; x, H, plus, wave, left, max_step, wave_norm) = sat
    n, h = length(x), x[2]-x[1]
    step < max_step || error("time step exceeds the spectral RK4 bound")
    # The extra oscillator appends the exact, time-dependent inflow.
    incoming = [plus -left zeros(n); zeros(2,n) [0.0 2π; -2π 0.0]]
    initial = vcat(sinpi.(2 .* x),0,1)
    pulse(x) = abs((x-0.3)/0.18) < 1 ? exp(1 - 1/(1-((x-0.3)/0.18)^2)) : 0.0
    wave_initial = vcat(pulse.(x),pulse.(2 .- x))
    ta, tw = 8.0, 1.0
    au, atime, afine, asemigroup = time_errors(incoming, initial, ta, step)
    v, wtime, wfine, wsemigroup = time_errors(wave, wave_initial, tw, step)
    u = au[1:n]
    ua = sinpi.(2 .* (x .- ta))
    vw = vcat(pulse.(x .- tw), pulse.(2 .- x .- tw))
    mode = cis.(2π .* (x .- ta))
    gridscale = abs(dot((-1.0) .^ (0:n-1), H * (u - ua)))
    function energy_path(G, initial, H, tf, nphysical)
        i = 1:nphysical
        e0 = dot(initial[i], H * initial[i])
        return [let u = (exp(t * G) * initial)[i]; dot(u,H*u)/e0; end
                for t in range(0,tf;length=9)]
    end
    return (; max_cfl=max_step/h, step, advection_error=maximum(abs,u-ua), atime, afine, asemigroup,
             aen=energy_path(incoming,initial,H,ta,n),
             advection_phase=angle(dot(mode,H*u)/dot(mode,H*ua)), gridscale,
             wave_error=maximum(abs,v-vw), wtime, wfine, wsemigroup,
             wen=energy_path(wave,wave_initial,wave_norm,tw,2n),
             reflection_error=maximum(abs,view(v-vw,(n+1):2n)))
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
for (p, depth, N) in DESIGNS
    system, first_op = seed(p, depth, N)
    data = Experimental.alternating_characteristics(system, first_op)
    println((; task=:exact, p, data...))
    report = dispersion(data)
    if "--dispersion-table" in ARGS
        foreach(r -> println((; task=:dispersion_row, p, r...)), report.rows)
    end
    println((; task=:dispersion, p, Base.structdiff(report, (rows=report.rows,))...))
    println((; task=:metrics, p, family=:alternating, operator_metrics(first_op, depth)...))
    comparisons = map(((:diagonal, 2p), (:restricted_full, p+2))) do (norm, interior)
        family = design_sbp(interior,p; norm)
        op = instantiate(family; N=max(N, 2family.closure_width+1))
        println((; task=:metrics, p, family=norm, operator_metrics(op, family.closure_depth)...))
        family
    end
    for k in 2:6
        n = div(k * (N-1),2) + 1
        op = n == N ? first_op : Experimental.refine_stencil(system,first_op;N=n)
        println((; task=:truncation, p,n, errors=truncation(op,depth),
                 op.diagnostics.nullity, extrema=extrema(op.diagnostics.row_orders)))
    end
    for delta in (2, 4, 6, 8, 10, 12)
        op = Experimental.refine_stencil(system,first_op;N=N+delta,domain=(-1,1))
        op.diagnostics.nullity == 1 && all(==(p),op.diagnostics.row_orders) &&
            op.D == -reverse(op.D; dims=(1,2)) || error("parity grid validation failed")
    end
    println((; task=:phase_aligned_grid_scan, p, tested=collect(N:2:N+12)))
    if "--pde" in ARGS
        minimum_N = max(N, maximum(2f.closure_width+1 for f in comparisons))
        first_N = isodd(minimum_N) ? minimum_N : minimum_N + 1
        for n in (first_N, 2first_N - 1, 3first_N - 2)
            op = n == N ? first_op : Experimental.refine_stencil(system,first_op;N=n)
            candidates = (op, (instantiate(f; N=n) for f in comparisons)...)
            sats = map(sat_matrices, candidates)
            # Common step on each grid, checked against all SAT spectra.
            dt = min(0.01/(n-1), minimum(sat.max_step for sat in sats)/8)
            for (name, sat) in zip((:alternating,:diagonal,:restricted_full), sats)
                println((; task=:pde, p,n,family=name, pde_benchmark(sat; step=dt)...))
            end
        end
    end
end
end
