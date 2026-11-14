# polish signed autocorrelation C1d. usage: acd.py in out budget_s beta0 [noise seed]
import numpy as np, sys, time
from scipy.optimize import minimize
from scipy.signal import fftconvolve
src,out,budget,beta0=sys.argv[1],sys.argv[2],float(sys.argv[3]),float(sys.argv[4])
noise=float(sys.argv[5]) if len(sys.argv)>5 else 0; seed=int(sys.argv[6]) if len(sys.argv)>6 else 0
v=np.array(open(src).read().replace(',',' ').split(),dtype=float); n=len(v)
def score(v): c=fftconvolve(v,v); return 2*n*np.abs(c).max()/v.sum()**2
print('seed n',n,'score',score(v),flush=True)
v=v/np.abs(v).mean(); best=v.copy(); bs=score(v)
def f(v,beta):
    s=v.sum(); K=2*n/s**2; c=fftconvolve(v,v)*K; cc=np.concatenate([c,-c]); m=cc.max()
    w=np.exp(beta*(cc-m)); Z=w.sum(); val=m+np.log(Z)/beta; w/=Z
    wd=w[:len(c)]-w[len(c):]            # d val / d c_k
    g=2*K*fftconvolve(wd,v[::-1])[n-1:2*n-1] - 2*(wd*c).sum()/s
    return val,g
rng=np.random.default_rng(seed); t0=time.time(); rnd=0
while time.time()-t0<budget:
    x=best*(1+noise*rng.normal(size=n)) if rnd>0 and noise>0 else best.copy()
    beta=beta0 if rnd>0 or noise==0 else beta0
    while beta<2e10 and time.time()-t0<budget:
        r=minimize(f,x,args=(beta,),jac=True,method='L-BFGS-B',options={'maxiter':1500,'maxcor':30,'ftol':1e-16,'gtol':1e-14}); x=r.x; beta*=1.6
    s=score(x); rnd+=1
    print(f'round {rnd} score {s:.12f} best {bs:.12f} {"BEST" if s<bs else ""} t {time.time()-t0:.0f}',flush=True)
    if s<bs: bs=s; best=x.copy(); np.savetxt(out,best/np.abs(best).max()*1e4)
    if noise==0: break
print('final',bs,flush=True)
