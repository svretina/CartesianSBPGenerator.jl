function x = lsqminnorm(A, b)
    % Octave compatibility for the two-argument MATLAB least-norm solve.
    x = pinv(full(A)) * b;
end
