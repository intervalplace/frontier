#!/usr/bin/env python3
"""Exact check of a C1d upper-bound certificate (Python standard library only).

Input: n signed integer heights v_0..v_{n-1} of a step function on n equal cells of [-1/4, 1/4]
(a text file, whitespace separated). f*f is linear between the knots, where it is proportional to
c = v*v, so  max|f*f| / (int f)^2 = 2n max_m |c_m| / S^2  with S = sum(v).

Every c_m is computed exactly: c = |v|*|v| - 4 (P*Q), where v = P - Q splits v into its positive and
negative parts, and each product is one big-integer multiplication (Kronecker substitution). The whole
sequence is then confirmed by evaluating both sides of (sum v_i x^i)^2 = sum c_m x^m at random points
modulo the prime 2^61 - 1, and the deciding knot is summed again directly.

usage: python3 verify.py heights.txt [--beat NUM/DEN] [--claim DECIMAL]
"""
import sys, random
from fractions import Fraction
from decimal import Decimal, getcontext

def kron(a, b, nbytes):
    A = int.from_bytes(b"".join(x.to_bytes(nbytes, "little") for x in a), "little")
    B = A if b is a else int.from_bytes(b"".join(x.to_bytes(nbytes, "little") for x in b), "little")
    L = len(a) + len(b) - 1
    C = (A * A if b is a else A * B).to_bytes(L * nbytes, "little")
    return [int.from_bytes(C[k * nbytes:(k + 1) * nbytes], "little") for k in range(L)]

def main():
    args = sys.argv[1:]; fn = args[0]
    v = [int(t) for t in open(fn).read().split()]
    n = len(v); S = sum(v); mx = max(abs(x) for x in v)
    assert n > 0 and S > 0
    nbytes = ((n * mx * mx).bit_length() + 1 + 7) // 8          # every coefficient fits in a slot
    absv = [abs(x) for x in v]; P = [x if x > 0 else 0 for x in v]; Q = [-x if x < 0 else 0 for x in v]
    AA = kron(absv, absv, nbytes); PQ = kron(P, Q, nbytes)
    c = [AA[m] - 4 * PQ[m] for m in range(2 * n - 1)]
    # confirm the whole sequence modulo a prime at random points
    p = (1 << 61) - 1; rng = random.Random(2026)
    for _ in range(4):
        x = rng.randrange(2, p - 1)
        lhs = 0
        for t in reversed(v): lhs = (lhs * x + t) % p
        rhs = 0
        for t in reversed(c): rhs = (rhs * x + t) % p
        assert lhs * lhs % p == rhs, "polynomial identity failed"
    M = max(abs(t) for t in c); arg = [m for m in range(2 * n - 1) if abs(c[m]) == M]
    m0 = arg[0]
    direct = sum(v[i] * v[m0 - i] for i in range(max(0, m0 - n + 1), min(m0, n - 1) + 1))
    assert abs(direct) == M, "direct re-summation at the deciding knot disagrees"
    val = Fraction(2 * n * M, S * S)
    getcontext().prec = 40
    print(f"n = {n}, max |v_i| = {mx}, S = {S}")
    print(f"max_m |c_m| = {M} at knot(s) m = {arg[:10]}{' ...' if len(arg) > 10 else ''} ({len(arg)} knot(s))")
    print(f"exact bound = {val.numerator}/{val.denominator}")
    print(f"            = {Decimal(val.numerator) / Decimal(val.denominator)}")
    ok = True
    if "--claim" in args:
        cl = Fraction(Decimal(args[args.index("--claim") + 1]))
        r = val <= cl; ok &= r; print(f"claim C_1d <= {args[args.index('--claim') + 1]}: {'implied' if r else 'NOT implied'}")
    for k, a in enumerate(args):
        if a == "--beat":
            rec = Fraction(args[k + 1]); r = val < rec; ok &= r
            print(f"beats {args[k + 1]} strictly: {r} (margin {float(rec - val):.6e})")
    print("RESULT:", "PASS" if ok else "FAIL")

main()
