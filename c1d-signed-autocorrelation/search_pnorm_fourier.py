# Fourier-preconditioned p-norm continuation for C1d (after the Commonplace agent ProofOfPersistence's
# "Fourier-scaled p-norm continuation", notebook pages 0xdf9c10b9 / 0x58a7bb1f on signed-autocorrelation-c1d-large).
# At each stage: base = current point (mean 1), z = base + T(u) with T a fixed real symmetric n-periodic Fourier multiplier
# D_k = median_{k>0}|B_k| / max(|B_k|, 0.1 median), D_0 = 0 (B = FFT of base), so the mean stays 1;
# minimise log||z*z||_p - 2 log(sum z) over u with L-BFGS-B (maxcor 16), gradient T(grad_z).
# The best true ratio 2n max|z*z| / (sum z)^2 over ALL evaluated points is kept.
# usage: python3 pnF.py start.(npy|txt) tag p0 pmax factor iters_first iters_rest [fixedD]
import sys, time
import numpy as np, scipy.fft as sf, scipy.optimize as so
start, tag = sys.argv[1], sys.argv[2]
p0, pmax, fac = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
it1, it2 = int(sys.argv[6]), int(sys.argv[7]); fixedD = len(sys.argv) > 8 and sys.argv[8] == "fixed"
v = np.load(start).astype(np.float64) if start.endswith(".npy") else np.array([float(t) for t in open(start).read().split()])
n = len(v); v *= n / v.sum()
N = sf.next_fast_len(2 * n - 1, real=True)
best = [np.inf, None]
def conv(z):
    Z = sf.rfft(z, N); return sf.irfft(Z * Z, N)[:2 * n - 1], Z
def ratio(z):
    c, _ = conv(z); return 2 * n * np.abs(c).max() / z.sum() ** 2
def multiplier(base):
    B = np.abs(sf.rfft(base)); med = np.median(B[1:])
    D = med / np.maximum(B, 0.1 * med); D[0] = 0.0; return D
def make_fg(base, D, p):
    def T(x): return sf.irfft(D * sf.rfft(x), n)
    def fg(u):
        z = base + T(u)
        c, Z = conv(z); S = z.sum(); a = np.abs(c); M = a.max()
        r = 2 * n * M / S ** 2
        if r < best[0]: best[0] = r; best[1] = z.copy()
        with np.errstate(divide="ignore"): lr = np.log(a / M)
        t1 = np.exp((p - 1) * lr); Tt = (t1 * (a / M)).sum()
        J = np.log(M) + np.log(Tt) / p - 2 * np.log(S)
        w = np.sign(c) * t1 / (M * Tt)
        gz = 2 * sf.irfft(sf.rfft(w, N) * np.conj(Z), N)[:n] - 2.0 / S
        return J, T(gz)
    return fg, T
log = open(f"{tag}.log", "a")
r0 = ratio(v); best[0] = r0; best[1] = v.copy()
print(time.strftime("%H:%M:%S"), f"n={n} start ratio {r0:.13f} p0={p0:g} factor={fac} iters {it1}/{it2} fixedD={fixedD}", file=log, flush=True)
t0 = time.time(); p = p0; k = 0; D = multiplier(v)
while p <= pmax * 1.0001:
    if not fixedD: D = multiplier(v)
    fg, T = make_fg(v.copy(), D, p)
    r = so.minimize(fg, np.zeros(n), jac=True, method="L-BFGS-B",
                    options=dict(maxiter=it1 if k == 0 else it2, maxfun=int(1.4 * (it1 if k == 0 else it2)), maxcor=16, ftol=2.3e-16, gtol=1e-10, maxls=30))
    v = v + T(r.x); v *= n / v.sum()
    R = ratio(v); np.save(f"{tag}_best.npy", best[1])
    print(time.strftime("%H:%M:%S"), f"p={p:.4g} it={r.nit} {r.message[:12]} endpoint {R:.13f} best-seen {best[0]:.13f} ({time.time()-t0:.0f}s)", file=log, flush=True)
    p *= fac; k += 1
print(time.strftime("%H:%M:%S"), "done best", f"{best[0]:.13f}", file=log, flush=True)
