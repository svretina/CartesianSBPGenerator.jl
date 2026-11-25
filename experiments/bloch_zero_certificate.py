"""Read the six coefficient table rows from stdin and count Bloch zeros exactly."""
import re
import sys

import sympy as sp

z = sp.symbols("z")
rows = {}
for line in sys.stdin:
    match = re.match(r"\| (\d+) \| ([AB]) \| `[^`]+` \| `([^`]+)` \|", line)
    if match:
        order, phase, taps = match.groups()
        rows[int(order), phase] = [sp.Rational(tap) for tap in taps.split(", ")]

for order in (2, 4, 6):
    def part(phase, parity):
        return sum(tap * sp.chebyshevu(j - 1, z)
                   for j, tap in enumerate(rows[order, phase], 1) if j % 2 == parity)

    polynomial = sp.Poly(part("A", 0) * part("B", 0) -
                         part("A", 1) * part("B", 1), z)
    roots = polynomial.count_roots(0, 1)
    assert roots == 0 and polynomial.eval(0) < 0 and polynomial.eval(1) < 0
    print(f"order={order}, degree={polynomial.degree()}, roots_on_[0,1]={roots}")
