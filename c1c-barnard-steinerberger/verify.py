#!/usr/bin/env python3
"""Exact verifier for a lower bound on the Barnard-Steinerberger constant C_1c (Tao's optimization constants, 1c).

For f = sum_i w_i 1[i/K, (i+1)/K) with integer w_i >= 0 the autocorrelation is linear in t between multiples of
1/K, so  min_{0<=t<=1} int f(x) f(x+t) dx / (int f)^2 = K * min_{0<=k<=K} A_k / S^2,  A_k = sum_i w_i w_{i+k},
S = sum_i w_i  (Kravitz, arXiv:2004.06611, Thm 3.1).  Standard library only, all integer arithmetic.

All A_k come from one big-integer product: with X = sum_i w_i B^i and Y = sum_i w_{n-1-i} B^i, the base-B digit
n-1+k of X*Y is A_k, exactly, when B exceeds every possible coefficient (no carries).  As a cross-check, a sample
of lags is also summed term by term.

Usage: python3 verify.py certificate.json [--beat NUM/DEN]
"""
import json, sys, random
from fractions import Fraction

def score(K, w):
    n = len(w); S = sum(w)
    bits = (max(w).bit_length() * 2 + n.bit_length() + 2)
    # pack: build hex strings for speed
    width = (bits + 3) // 4
    X = int("".join(format(x, "0%dx" % width) for x in reversed(w)), 16)
    Y = int("".join(format(x, "0%dx" % width) for x in w), 16)
    P = X * Y
    h = format(P, "x")
    h = "0" * ((2 * n) * width - len(h)) + h if len(h) < 2 * n * width else h
    def digit(j):  # base-B digit j of P, B = 16**width
        end = len(h) - j * width
        return int(h[max(0, end - width):end] or "0", 16)
    A = [digit(n - 1 + k) for k in range(n)]
    return K, n, S, A

def main():
    c = json.load(open(sys.argv[1])); K = c["K"]; w = c["weights_int"]
    assert isinstance(K, int) and K >= 1
    assert all(isinstance(x, int) and x >= 0 for x in w) and any(w)
    K, n, S, A = score(K, w)
    rng = random.Random(1)
    for k in [0, min(K, n - 1)] + rng.sample(range(min(K, n - 1) + 1), 25):
        assert A[k] == sum(w[i] * w[i + k] for i in range(n - k)), f"carry check failed at k={k}"
    vals = [(A[k] if k < n else 0, k) for k in range(K + 1)]
    m, kmin = min(vals)
    v = Fraction(K * m, S * S)
    print(f"K = {K}, {n} weights, minimum at k = {kmin}")
    print(f"C_1c >= {v.numerator}/{v.denominator}")
    print(f"     = {v.numerator * 10**18 // v.denominator / 10**18:.18f} (truncated)")
    if "score_num" in c:
        assert v == Fraction(c["score_num"], c["score_den"]), "stated score does not match"
        print("matches the stated score_num/score_den")
    if len(sys.argv) > 3 and sys.argv[2] == "--beat":
        a, b = sys.argv[3].split("/"); prev = Fraction(int(a), int(b))
        assert v > prev; print(f"beats {prev} by {float(v - prev):.6e}")

if __name__ == "__main__":
    main()
