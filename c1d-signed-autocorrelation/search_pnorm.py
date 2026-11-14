#!/usr/bin/env python3
# Search code for the C1d certificate in this folder (needs numpy and scipy).
# p-norm continuation: minimise log ||v*v||_p - 2 log(sum v) with L-BFGS, p doubled every stage
# (`iters` iterations a stage below p = 10^6, twice that above), keeping the vector with the lowest
# true ratio 2n max|v*v| / (sum v)^2 seen at the end of any stage.  Writes <tag>.log and <tag>_best.npy.
# usage: python3 search_pnorm.py start.npy tag p0 pmax iters [split]   (split=2 repeats every cell: the same function)
import sys, time, json, os
import numpy as np, scipy.fft as sf, scipy.optimize as so
start, tag, p0, pmax, iters = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
split = int(sys.argv[6]) if len(sys.argv) > 6 else 1
v = np.load(start).astype(np.float64)
if split > 1: v = np.repeat(v, split)
n = len(v); N = sf.next_fast_len(2 * n - 1, real=True)   # N > 2n-2: no wrap-around in the correlation
def ratio(v):
    V = sf.rfft(v, N); c = sf.irfft(V * V, N)[:2 * n - 1]
    return 2 * n * np.abs(c).max() / v.sum() ** 2
def fg(v, p):
    V = sf.rfft(v, N); c = sf.irfft(V * V, N)[:2 * n - 1]
    S = v.sum(); a = np.abs(c); M = a.max()
    if not (S > 0 and np.isfinite(S) and M > 0): return 1e6, -np.ones(n)
    with np.errstate(divide='ignore'):
        lr = np.log(a / M)
    t1 = np.exp((p - 1) * lr); T = (t1 * (a / M)).sum()
    J = np.log(M) + np.log(T) / p - 2 * np.log(S)
    w = np.sign(c) * t1 / (M * T)
    g = 2 * sf.irfft(sf.rfft(w, N) * np.conj(V), N)[:n] - 2.0 / S
    return J, g
v *= n / v.sum()
best = ratio(v); log = open(f"{tag}.log", "a")
print(time.strftime("%H:%M:%S"), f"n={n} start ratio {best:.13f} p0={p0:g}", file=log, flush=True)
np.save(f"{tag}_best.npy", v)
p = p0; t0 = time.time()
while p <= pmax:
    r = so.minimize(fg, v, args=(p,), jac=True, method="L-BFGS-B",
                    options=dict(maxiter=iters if p < 1e6 else 2 * iters, maxcor=10, ftol=0, gtol=0, maxls=40))
    v = r.x * (n / r.x.sum()); R = ratio(v)
    flag = ""
    if R < best: best = R; np.save(f"{tag}_best.npy", v); flag = " *"
    print(time.strftime("%H:%M:%S"), f"p={p:.4g} it={r.nit} {r.message[:12]} ratio {R:.13f} best {best:.13f} ({time.time()-t0:.0f}s){flag}", file=log, flush=True)
    p *= 2
print(time.strftime("%H:%M:%S"), "done best", f"{best:.13f}", file=log, flush=True)
