# Restart chain used for C1b_15360_0.3808564298.json: each round runs 3 configurations of search_restart_fourier.py in
# parallel from the best so far
# (optionally best + m (best - previous best)), keeps the best result. Logs to chain.log.  usage: python3 chain_restart.py start.npy
import os, sys, re, time, subprocess, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); os.chdir(D)
LOG = open("chain.log", "a")
def log(*a): print(time.strftime("%H:%M:%S"), *a, file=LOG, flush=True)
src = sys.argv[1]; best = np.load(src); prev = None
def bnd(path):
    t = open(path).read(); m = re.search(r"done best ([0-9.]+)", t); return float(m.group(1)) if m else None
import importlib.util
spec_src = open("search_restart_fourier.py").read(); pre = spec_src[:spec_src.index("a = np.clip(a0, 1e-15")]
g = {}; sys_argv = sys.argv; sys.argv = ["x", src, "x", "1", "1", "2", "1", "1"]; exec(pre, g); sys.argv = sys_argv
bval = g["bound"](best); log("chain start", f"{bval:.13f}")
cfgs = [("2", "3e6"), ("2", "1e6"), ("3", "3e6"), ("2", "1e7"), ("1", "3e6"), ("3", "1e6")]
k = 0
while True:
    procs = []
    for j in range(3):
        rad, p0 = cfgs[(3 * k + j) % len(cfgs)]
        st = best if prev is None or j == 0 else np.clip(best + 0.5 * (best - prev), 0, 1)
        if prev is not None and j > 0:
            st = st * (len(st) / 2) / st.sum()            # keep w * sum = 1 approximately (exact sum is restored inside)
            st = np.clip(st, 0, 1)
        np.save(f"cs_{j}.npy", st); tag = f"C{k}_{j}"
        env = dict(os.environ, PRE="1", RAD=rad, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
        procs.append((tag, rad, p0, subprocess.Popen(["python3", "search_restart_fourier.py", f"cs_{j}.npy", tag, p0, "1e11", "1.5", "2000", "1000"],
                                                       env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)))
    for tag, rad, p0, p in procs: p.wait()
    res = []
    for tag, rad, p0, p in procs:
        v = bnd(f"{tag}.log")
        if v is not None: res.append((v, tag, rad, p0))
    res.sort()
    if res and res[0][0] < bval:
        v, tag, rad, p0 = res[0]; prev, best = best, np.load(f"{tag}_best.npy");
        log("round", k, "best", tag, "RAD", rad, "p0", p0, f"{v:.13f}", "improves", f"{bval - v:.3e}", "| others", " ".join(f"{r[0]:.13f}" for r in res[1:]))
        bval = v; np.save("chain_best.npy", best)
    else:
        log("round", k, "no gain:", " ".join(f"{r[0]:.13f}" for r in res), "| best", f"{bval:.13f}"); prev = None
    k += 1
