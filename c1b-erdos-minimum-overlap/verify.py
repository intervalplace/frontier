#!/usr/bin/env python3
"""Exact verifier for an upper bound on the Erdos minimum overlap constant C_1b (Tao's optimization constants, 1b).

h is a step function on [0,2] with n equal steps of width w = 2/n and heights h_i = k_i/Q, 0 <= k_i <= Q, mass
w * sum h_i = 1 exactly (sum k_i = n*Q/2). With 1 - h extended by zero outside [0,2], the overlap
O(x) = int h(t) (1 - h(t+x)) dt is linear in x between multiples of w, so
    C_1b <= max_x O(x) = w * max_m sum_i h_i (1 - h_{i+m}) = 2 * max_m (Q * S_m - X_m) / (n * Q^2),
where S_m is the sum of k_i over the overlapping range and X_m = sum_i k_i k_{i+m}, -(n-1) <= m <= n-1.
Standard library only, all integer arithmetic. Every X_m comes from one big-integer product (base-B digits of
X*Y with X = sum k_i B^i, Y = sum k_{n-1-i} B^i, B larger than any X_m so no digit carries); a sample of shifts is
also summed term by term as a cross-check.

Usage: python3 verify.py certificate.json [--beat NUM/DEN]
"""
import json, sys, random
from fractions import Fraction

def main():
    c = json.load(open(sys.argv[1])); n, Q, k = c["n"], c["Q"], c["heights_int"]
    assert isinstance(n, int) and len(k) == n and n % 2 == 0
    assert all(isinstance(x, int) and 0 <= x <= Q for x in k), "heights must lie in [0, 1]"
    assert 2 * sum(k) == n * Q, "mass must be exactly 1"
    width = (2 * Q.bit_length() + n.bit_length() + 2 + 3) // 4          # hex digits per base-B digit
    X = int("".join(format(x, "0%dx" % width) for x in reversed(k)), 16)
    Y = int("".join(format(x, "0%dx" % width) for x in k), 16)
    h = format(X * Y, "x").rjust(2 * n * width, "0")
    def digit(j):
        end = len(h) - j * width
        return int(h[max(0, end - width):end] or "0", 16)
    pre = [0]
    for x in k: pre.append(pre[-1] + x)
    def overlap(m, cross):
        lo, hi = max(0, -m), min(n - 1, n - 1 - m)
        return Q * (pre[hi + 1] - pre[lo]) - cross
    vals = []
    for m in range(-(n - 1), n):
        vals.append((overlap(m, digit((n - 1) - m)), m))
    rng = random.Random(1)
    for m in [-(n - 1), 0, n - 1] + rng.sample(range(-(n - 1), n), 25):
        direct = sum(k[i] * k[i + m] for i in range(max(0, -m), min(n - 1, n - 1 - m) + 1))
        assert digit((n - 1) - m) == direct, f"carry check failed at m={m}"
    best, marg = max(vals)
    v = Fraction(2 * best, n * Q * Q)
    print(f"n = {n}, Q = {Q}, maximum at shift m = {marg}")
    print(f"C_1b <= {v.numerator}/{v.denominator}")
    d = str(v.numerator * 10**20 // v.denominator).rjust(21, "0")
    print(f"     = {d[0]}.{d[1:]} (exact digits, truncated)")
    if "score_num" in c:
        assert v == Fraction(c["score_num"], c["score_den"]), "stated score does not match"
        print("matches the stated score_num/score_den")
    if len(sys.argv) > 3 and sys.argv[2] == "--beat":
        a, b = sys.argv[3].split("/"); prev = Fraction(int(a), int(b))
        assert v < prev; print(f"beats {prev} by {float(prev - v):.6e}")

if __name__ == "__main__":
    main()
