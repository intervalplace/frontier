#!/usr/bin/env python3
"""Exact verifier for an upper bound on the signed third autocorrelation constant C_1d (Tao's optimization constants, 1d).

f has n equal cells on [-1/4, 1/4] with signed integer values v_i. f*f is piecewise linear with knots at multiples of
h = 1/(2n), so  max_t |f*f(t)| / (int f)^2 = 2n * max_m |c_m| / (sum v)^2,  c = v*v the full discrete autoconvolution.
Standard library only, all integer arithmetic. Writing v = p - q with p, q >= 0, every c_m = (p*p + q*q - 2 p*q)_m comes
from three big-integer products (no-carry base-B digits); 30 values of m are also summed term by term as a check.

Usage: python3 verify.py certificate.json [--beat NUM/DEN]
"""
import json, sys, random
from fractions import Fraction

def conv_digits(a, b, n, width):
    X = int("".join(format(x, "0%dx" % width) for x in reversed(a)), 16)
    Y = int("".join(format(x, "0%dx" % width) for x in reversed(b)), 16)
    h = format(X * Y, "x").rjust((2 * n) * width, "0")
    return [int(h[len(h) - (t + 1) * width: len(h) - t * width] or "0", 16) for t in range(2 * n - 1)]

def main():
    c = json.load(open(sys.argv[1])); n, v = c["n"], c["v_int"]
    assert isinstance(n, int) and len(v) == n and all(isinstance(x, int) for x in v)
    S = sum(v); assert S != 0
    p = [x if x > 0 else 0 for x in v]; q = [-x if x < 0 else 0 for x in v]
    width = (2 * max(abs(x) for x in v).bit_length() + n.bit_length() + 2 + 3) // 4
    pp, qq, pq = conv_digits(p, p, n, width), conv_digits(q, q, n, width), conv_digits(p, q, n, width)
    cv = [pp[t] + qq[t] - 2 * pq[t] for t in range(2 * n - 1)]
    rng = random.Random(1)
    for m in [0, n - 1, 2 * n - 2] + rng.sample(range(2 * n - 1), 27):
        direct = sum(v[i] * v[m - i] for i in range(max(0, m - n + 1), min(m, n - 1) + 1))
        assert cv[m] == direct, f"carry check failed at m={m}"
    best, arg = max((abs(x), m) for m, x in enumerate(cv))
    val = Fraction(2 * n * best, S * S)
    print(f"n = {n}, maximum |c_m| at m = {arg}")
    print(f"C_1d <= {val.numerator}/{val.denominator}")
    d = str(val.numerator * 10**20 // val.denominator).rjust(21, "0")
    print(f"     = {d[0]}.{d[1:]} (exact digits, truncated)")
    if "score_num" in c:
        assert val == Fraction(c["score_num"], c["score_den"]), "stated score does not match"
        print("matches the stated score_num/score_den")
    if len(sys.argv) > 3 and sys.argv[2] == "--beat":
        a, b = sys.argv[3].split("/"); prev = Fraction(int(a), int(b))
        assert val < prev; print(f"beats {prev} by {float(prev - val):.6e}")

if __name__ == "__main__":
    main()
