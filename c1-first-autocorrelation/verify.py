#!/usr/bin/env python3
"""Exact verifier for an upper bound on the first autocorrelation constant C1 (Tao's optimization constants, 1a).

f = sum_i a_i 1[-1/4 + i h, -1/4 + (i+1) h), h = 1/(2n), integer a_i >= 0.  f*f is piecewise linear with corners
at the cell boundaries, so  max_t (f*f)(t) / (int f)^2 = 2n * max_k c_k / S^2,  c_k = sum_{i+j=k} a_i a_j,  S = sum a_i.
All c_k come from one big-integer square: with X = sum_i a_i B^i, the base-B digit k of X^2 is c_k exactly when B
exceeds every coefficient (no carries).  A sample of c_k is also summed term by term.  Standard library only.

Usage: python3 verify_c1.py heights.txt [--beat NUM/DEN]     (heights: whitespace-separated nonnegative integers)
"""
import sys, random
from fractions import Fraction

def main():
    a = [int(t) for t in open(sys.argv[1]).read().split() if not t.startswith('#')]
    assert all(x >= 0 for x in a) and any(a)
    n = len(a); S = sum(a)
    bits = 2 * max(a).bit_length() + n.bit_length() + 2
    width = (bits + 3) // 4
    X = int("".join(format(x, "0%dx" % width) for x in reversed(a)), 16)
    h = format(X * X, "x")
    h = "0" * (2 * n * width - len(h)) + h
    c = [int(h[len(h) - (k + 1) * width:len(h) - k * width], 16) for k in range(2 * n - 1)]
    rng = random.Random(1)
    for k in [0, n - 1, 2 * n - 2] + rng.sample(range(2 * n - 1), 20):
        lo, hi = max(0, k - n + 1), min(k, n - 1)
        assert c[k] == sum(a[i] * a[k - i] for i in range(lo, hi + 1)), f"carry check failed at k={k}"
    m = max(c); km = c.index(m)
    v = Fraction(2 * n * m, S * S)
    print(f"n = {n} pieces, maximum of f*f at corner k = {km}")
    print(f"C1 <= {v.numerator}/{v.denominator}")
    print(f"   = {-(-v.numerator * 10**18 // v.denominator) / 10**18:.18f} (rounded up)")
    if len(sys.argv) > 3 and sys.argv[2] == "--beat":
        p, q = sys.argv[3].split("/"); prev = Fraction(int(p), int(q))
        assert v < prev; print(f"below {float(prev):.15f} by {float(prev - v):.6e}")

if __name__ == "__main__":
    main()
