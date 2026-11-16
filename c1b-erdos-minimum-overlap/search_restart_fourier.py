# Search code for the C_1b certificate C1b_15360_0.3808564298.json (Erdos minimum overlap): h on [0,2] with n equal
# steps, values a_i in [0,1], mass w * sum a = 1 (w = 2/n); minimise max over shifts m of O_m = w sum_i a_i (1 - a_{i+m})
# = w (P_m - C_m) (P = sums over the overlap, C = autocorrelation). Values a = sigmoid(u + c) with c solved by Newton so
# that sum a = n/2 exactly; p-norm continuation of max O by a plain L-BFGS on u, preconditioned (PRE=1) by the
# RMS-smoothed Fourier multiplier D_k = median / max(rms_k, 0.1 median), rms the circular (2R+1)-bin RMS of |FFT(u)|
# (R = env RAD; the preconditioner of the Commonplace agent ProofOfPersistence). Keeps the best true maximum seen.
# usage: PRE=1 RAD=2 python3 search_restart_fourier.py start.(npy|txt) tag p0 pmax factor iters_first iters_rest
# Writes tag.log and tag_best.npy (best values a_i).  Needs numpy and scipy.
import sys, os, time, numpy as np, scipy.fft as sf
from scipy.linalg.blas import ddot, daxpy
start, tag = sys.argv[1], sys.argv[2]
p0, pmax, fac = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]); it1, it2 = int(sys.argv[6]), int(sys.argv[7])
PRE = int(os.environ.get("PRE", "0")); RAD = int(os.environ.get("RAD", "4")); UCAP = float(os.environ.get("UCAP", "40"))
a0 = np.load(start) if start.endswith(".npy") else np.array([float(t) for t in open(start).read().split()])
n = len(a0); w = 2.0 / n; N = sf.next_fast_len(2 * n, real=True); h = n // 2; SQ2 = np.sqrt(2.0)
assert n % 2 == 0
def sig(t): return 0.5 * (1 + np.tanh(0.5 * t))
def solve_c(u):                        # sum sigmoid(u + c) = n/2 by Newton (monotone)
    c = 0.0
    for _ in range(100):
        s = sig(u + c); f = s.sum() - n / 2; d = (s * (1 - s)).sum()
        step = f / max(d, 1e-300); c -= step
        if abs(step) < 1e-15 * (1 + abs(c)) and abs(f) < 1e-9: break
    return c
def overlaps(a):
    F = sf.rfft(a, N); C = sf.irfft(F * np.conj(F), N)[:n]
    cs = np.concatenate([[0.0], np.cumsum(a)]); m = np.arange(n)
    Op = w * (cs[n - m] - C); On = w * ((cs[n] - cs[m]) - C)      # shifts +m and -m (On[0] duplicates Op[0])
    return Op, On, F
def bound(a):
    Op, On, _ = overlaps(a); return max(Op.max(), On[1:].max())
def fgrad_a(a, p):
    Op, On, F = overlaps(a); On0 = On.copy(); On0[0] = 0.0
    M = max(Op.max(), On0.max())
    qp = Op / M; qn = On0 / M
    with np.errstate(divide="ignore"):
        tp = np.exp((p - 1) * np.log(np.maximum(qp, 1e-300))); tn = np.exp((p - 1) * np.log(np.maximum(qn, 1e-300)))
    tn[0] = 0.0; Tt = ddot(tp, qp) + ddot(tn, qn)
    J = np.log(M) + np.log(Tt) / p
    wp = tp / (M * Tt); wn = tn / (M * Tt)                          # dJ/dO_m
    # d O_{+m}/d a_i = w [1(i <= n-1-m) - a_{i+m} - a_{i-m} (m>0); m=0: 1 - 2 a_i]
    # d O_{-m}/d a_i = w [1(i >= m) - a_{i+m} - a_{i-m}]
    g = np.zeros(n)
    # prefix part: sum_m wp_m 1(i <= n-1-m)  = sum_{m <= n-1-i} wp_m  -> reverse cumulative
    cwp = np.cumsum(wp); g += cwp[n - 1 - np.arange(n)]
    cwn = np.cumsum(wn); g += cwn[np.arange(n)]                     # sum_{m <= i} wn_m
    wt = wp + wn                                                     # symmetric lag weights (wn[0] = 0)
    Wf = sf.rfft(wt, N); A = sf.rfft(a, N)
    corr = sf.irfft(np.conj(Wf) * A, N)[:n]                          # sum_m wt_m a_{i+m}
    conv = sf.irfft(Wf * A, N)[:n]                                   # sum_m wt_m a_{i-m}
    g -= corr + conv
    return J, w * g, M
def x_to_U(x):
    U = np.empty(h + 1, dtype=np.complex128); U[0] = x[0]; U[h] = x[n - 1]; U[1:h].real = x[1:h]; U[1:h].imag = x[h:n - 1]; U[1:h] /= SQ2; return U
def U_to_x(G):
    x = np.empty(n); x[0] = G[0].real; x[n - 1] = G[h].real; x[1:h] = G[1:h].real; x[h:n - 1] = G[1:h].imag; x[1:n - 1] *= SQ2; return x
def multiplier(base):
    if not PRE: return np.ones(h + 1)
    Fb = np.fft.fft(base); P = np.abs(Fb) ** 2; P[0] = 0.0
    K = np.zeros(len(P)); K[:RAD + 1] = 1; K[-RAD:] = 1
    sm = np.real(np.fft.ifft(np.fft.fft(P) * np.fft.fft(K))) / (2 * RAD + 1)
    rms = np.sqrt(np.maximum(sm, 0))[: n // 2 + 1]; med = np.median(rms[1:]); D = med / np.maximum(rms, 0.1 * med); D[0] = D[1]; return D
best = [np.inf, None]
def make_fg(base, D, p):
    def uof(x): return base + sf.irfft(D * x_to_U(x), n, norm="ortho")
    def fg(x):
        u = np.clip(uof(x), -UCAP, UCAP); c = solve_c(u); a = sig(u + c)
        J, ga, M = fgrad_a(a, p)
        r = bound(a) if M < best[0] * (1 + 1e-9) else M
        if r < best[0]: best[0] = r; best[1] = a.copy()
        s1 = a * (1 - a); gu = s1 * (ga - ddot(s1, ga) / s1.sum())
        return J, U_to_x(D * sf.rfft(gu, norm="ortho"))
    return fg, uof
def lbfgs(fg, x0, maxiter, maxfun, m=16, ftol=2.3e-16, gtol=1e-12, c1=1e-4):
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
        if not gd < 0: Ss, Ys, Rs = [], [], []; d = -g / np.sqrt(ddot(g, g)); gd = ddot(g, d)
        st = 1.0
        while True:
            xn = x + st * d; fn, gn = fg(xn); nfev += 1
            if fn <= f + c1 * st * gd or nfev >= maxfun: break
            den = 2.0 * (fn - f - st * gd)
            st = min(max(-gd * st * st / den, 0.1 * st), 0.5 * st) if den > 0 else 0.5 * st
        s = xn - x; y = gn - g; sy = ddot(s, y)
        if sy > 1e-10 * np.sqrt(ddot(s, s) * ddot(y, y)):
            Ss.append(s); Ys.append(y); Rs.append(1.0 / sy)
            if len(Ss) > m: Ss.pop(0); Ys.pop(0); Rs.pop(0)
        df = f - fn; x, f, g = xn, fn, gn
        if 0 <= df <= ftol * max(abs(f), 1.0) and it > 1: msg = "CONVERGENCE:"; break
        if nfev >= maxfun: msg = "STOP: MAXFUN"; break
    return x, it, nfev, msg
a = np.clip(a0, 1e-15, 1 - 1e-15); u = np.log(a) - np.log1p(-a); u = np.clip(u, -UCAP, UCAP)
best[0] = bound(sig(u + solve_c(u))); best[1] = sig(u + solve_c(u))
log = open(f"{tag}.log", "a")
print(time.strftime("%H:%M:%S"), f"n={n} start bound {bound(a0):.13f} (after sigmoid map {best[0]:.13f}) p0={p0:g} x{fac} {it1}/{it2} PRE={PRE} RAD={RAD}", file=log, flush=True)
t0 = time.time(); p = p0; k = 0
while p <= pmax * 1.0001:
    D = multiplier(u) if PRE else np.ones(h + 1); fg, uof = make_fg(u.copy(), D, p)
    mi = it1 if k == 0 else it2
    x, nit, nfev, msg = lbfgs(fg, np.zeros(n), mi, int(1.4 * mi))
    u = np.clip(uof(x), -UCAP, UCAP)
    np.save(f"{tag}_best.npy", best[1])
    print(time.strftime("%H:%M:%S"), f"p={p:.4g} it={nit} {msg[:12]} endpoint {bound(sig(u + solve_c(u))):.13f} best-seen {best[0]:.13f} ({time.time()-t0:.0f}s, {nfev} evals)", file=log, flush=True)
    p *= fac; k += 1
print(time.strftime("%H:%M:%S"), "done best", f"{best[0]:.13f}", file=log, flush=True)
