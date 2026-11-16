# Search code for the C1 (constant 1a) certificate C1_262144_1.5026648295.txt: p-norm continuation of max(f*f) with
# heights v = w^2 (w unconstrained), preconditioned by an RMS-smoothed Fourier multiplier (the preconditioner of the
# Commonplace agent ProofOfPersistence, introduced on signed-autocorrelation-c1d-large). Each stage optimises the step u
# in w through its orthonormal half spectrum x, w = base + irfft(D * U(x)), D = median / max(rms, 0.1 median) where rms is
# the circular (2R+1)-bin RMS of |FFT(base)| (R = env RAD), minimising log ||v*v||_p - 2 log sum v with a plain L-BFGS
# (two-loop recursion, Armijo backtracking); p is multiplied by `factor` each stage and the best true ratio seen is kept.
# usage: RAD=4 PRE=1 python3 search_pnorm_fourier.py start.(npy|txt) tag p0 pmax factor iters_first iters_rest
# Writes tag.log (one line per stage) and tag_best.npy (best heights seen).  Needs numpy and scipy.
import sys, time, os
import numpy as np, scipy.fft as sf
from scipy.linalg.blas import ddot, daxpy
start, tag = sys.argv[1], sys.argv[2]
p0, pmax, fac = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
it1, it2 = int(sys.argv[6]), int(sys.argv[7])
RAD = int(os.environ.get("RAD", "0")); MCOR = int(os.environ.get("MCOR", "16"))
v = np.load(start).astype(np.float64) if start.endswith(".npy") else np.array([float(t) for t in open(start).read().split()])
n = len(v); v = np.maximum(v, 0.0); v *= n / v.sum(); assert n % 2 == 0
PRE = int(os.environ.get("PRE", "1"))
N = sf.next_fast_len(2 * n - 1, real=True); h = n // 2; SQ2 = np.sqrt(2.0)
best = [np.inf, None]
def ratio(z):
    Z = sf.rfft(z, N); c = sf.irfft(Z * Z, N)[:2 * n - 1]; return 2 * n * c.max() / z.sum() ** 2
def multiplier(base):
    if not PRE: D = np.ones(len(base) // 2 + 1); return D
    if RAD <= 0:
        B = np.abs(sf.rfft(base)); med = np.median(B[1:])
        D = med / np.maximum(B, 0.1 * med); D[0] = D[1]; return D
    F = np.fft.fft(base); P = np.abs(F) ** 2; P[0] = 0.0
    K = np.zeros(len(P)); K[:RAD + 1] = 1; K[-RAD:] = 1
    sm = np.real(np.fft.ifft(np.fft.fft(P) * np.fft.fft(K))) / (2 * RAD + 1)
    rms = np.sqrt(np.maximum(sm, 0))[: len(base) // 2 + 1]
    med = np.median(rms[1:]); D = med / np.maximum(rms, 0.1 * med); D[0] = D[1]; return D
def x_to_U(x):
    U = np.empty(h + 1, dtype=np.complex128)
    U[0] = x[0]; U[h] = x[n - 1]; U[1:h].real = x[1:h]; U[1:h].imag = x[h:n - 1]; U[1:h] /= SQ2
    return U
def U_to_x(G):
    x = np.empty(n); x[0] = G[0].real; x[n - 1] = G[h].real
    x[1:h] = G[1:h].real; x[h:n - 1] = G[1:h].imag; x[1:n - 1] *= SQ2
    return x
def make_fg(base, D, p):
    # base = current w; w = base + irfft(D * U(x)); v = w^2
    def zof(x): return base + sf.irfft(D * x_to_U(x), n, norm="ortho")
    def fg(x):
        w = zof(x); z = w * w
        Z = sf.rfft(z, N); c = sf.irfft(Z * Z, N)[:2 * n - 1]
        S = z.sum(); a = np.maximum(c, 0.0); M = a.max()
        r = 2 * n * M / S ** 2
        if r < best[0]: best[0] = r; best[1] = z.copy()
        q = a * (1.0 / M)
        with np.errstate(divide="ignore"): t1 = np.log(q)
        t1 *= (p - 1); np.exp(t1, out=t1); Tt = ddot(t1, q)
        J = np.log(M) + np.log(Tt) / p - 2 * np.log(S)
        t1 *= 1.0 / (M * Tt)
        gz = 2 * sf.irfft(sf.rfft(t1, N) * np.conj(Z), N)[:n] - 2.0 / S
        gw = 2.0 * w * gz
        return J, U_to_x(D * sf.rfft(gw, norm="ortho"))
    return fg, zof
def lbfgs(fg, x0, maxiter, maxfun, m=MCOR, ftol=2.3e-16, gtol=1e-10, c1=1e-4):
    x = x0.copy(); f, g = fg(x); nfev = 1; Ss, Ys, Rs = [], [], []; msg = "STOP: TOTAL"
    for it in range(1, maxiter + 1):
        q = g.copy(); al = []
        for s, y, rho in zip(reversed(Ss), reversed(Ys), reversed(Rs)):
            a_ = rho * ddot(s, q); al.append(a_); q = daxpy(y, q, a=-a_)
        if Ss: q *= ddot(Ss[-1], Ys[-1]) / ddot(Ys[-1], Ys[-1])
        else: q *= 1.0 / np.sqrt(ddot(g, g))
        for s, y, rho, a_ in zip(Ss, Ys, Rs, reversed(al)):
            b = rho * ddot(y, q); q = daxpy(s, q, a=a_ - b)
        d = -q; gd = ddot(g, d)
        if not gd < 0:
            Ss, Ys, Rs = [], [], []; d = -g / np.sqrt(ddot(g, g)); gd = ddot(g, d)
        st = 1.0
        while True:
            xn = x + st * d; fn, gn = fg(xn); nfev += 1
            if fn <= f + c1 * st * gd or nfev >= maxfun: break
            # safeguarded quadratic interpolation of f along d, then backtrack
            den = 2.0 * (fn - f - st * gd)
            st = min(max(-gd * st * st / den, 0.1 * st), 0.5 * st) if den > 0 else 0.5 * st
        s = xn - x; y = gn - g; sy = ddot(s, y)
        if sy > 1e-10 * np.sqrt(ddot(s, s) * ddot(y, y)):
            Ss.append(s); Ys.append(y); Rs.append(1.0 / sy)
            if len(Ss) > m: Ss.pop(0); Ys.pop(0); Rs.pop(0)
        df = f - fn; x, f, g = xn, fn, gn
        if 0 <= df <= ftol * max(abs(f), abs(fn), 1.0) and it > 1: msg = "CONVERGENCE:"; break
        if np.abs(g).max() <= gtol: msg = "CONVERGENCE:"; break
        if nfev >= maxfun: msg = "STOP: MAXFUN"; break
    return x, it, nfev, msg
log = open(f"{tag}.log", "a")
r0 = ratio(v); best[0] = r0; best[1] = v.copy()
print(time.strftime("%H:%M:%S"), f"n={n} start ratio {r0:.13f} p0={p0:g} factor={fac} iters {it1}/{it2} fixedD=False RAD={RAD} m={MCOR} PRE={PRE}", file=log, flush=True)
t0 = time.time(); p = p0; k = 0
wv = np.sqrt(v)
SNAP = float(os.environ.get("SNAP", "0")); snapped = [False]
while p <= pmax * 1.0001:
    D = multiplier(wv); fg, zof = make_fg(wv.copy(), D, p)
    mi = it1 if k == 0 else it2
    x, nit, nfev, msg = lbfgs(fg, np.zeros(n), mi, int(1.4 * mi))
    wv = zof(x); wv *= np.sqrt(n / (wv * wv).sum()); v = wv * wv
    R = ratio(v); np.save(f"{tag}_best.npy", best[1])
    print(time.strftime("%H:%M:%S"), f"p={p:.4g} it={nit} {msg[:12]} endpoint {R:.13f} best-seen {best[0]:.13f} ({time.time()-t0:.0f}s, {nfev} evals)", file=log, flush=True)
    if SNAP and p >= SNAP and not snapped[0]:      # keep the not-yet-polished best of this stage (truncation test)
        np.save(f"{tag}_snap.npy", best[1]); snapped[0] = True
        print(time.strftime("%H:%M:%S"), f"snapshot at p={p:.4g}: {best[0]:.13f}", file=log, flush=True)
    p *= fac; k += 1
print(time.strftime("%H:%M:%S"), "done best", f"{best[0]:.13f}", file=log, flush=True)
