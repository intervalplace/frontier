# polish Erdos minimum overlap. usage: eo.py in out rep budget beta0
import numpy as np, sys, time
from scipy.optimize import minimize
from scipy.signal import fftconvolve
src,out,rep,budget,beta0=sys.argv[1],sys.argv[2],int(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5])
a=np.repeat(np.array(open(src).read().replace(',',' ').split(),dtype=float),rep); n=len(a); w=2/n
def overl(a):  # C_m = sum_i a_i (1-a_{i+m}), all m in -(n-1)..n-1 : corr(a, 1) - corr(a,a)
    one=np.ones(n); return fftconvolve(a,one[::-1]) - fftconvolve(a,a[::-1])
def score(a):
    return w*overl(a).max()
def proj(a):  # fix mass to exactly n/2 while staying in [0,1]
    for _ in range(60):
        d=(n/2-a.sum()); free=(a>0)&(a<1) if d!=0 else None
        if abs(d)<1e-13*n: break
        a=np.clip(a+d/max(free.sum(),1)*free,0,1)
    return a
print('seed n',n,'score',score(a),'mass',w*a.sum(),flush=True)
LAM=1e4; MU=[0.0]
def f(a,beta):
    C=w*overl(a); m=C.max(); e=np.exp(beta*(C-m)); Z=e.sum(); val=m+np.log(Z)/beta; p=e/Z
    # dC_m/da_j = w*[ (1-a_{j+m}) - a_{j-m} ]  ; sum_m p_m (.)
    one=np.ones(n)
    g1=fftconvolve(p,one)[n-1:2*n-1]            # sum_m p_m * 1[j+m in range]
    g2=fftconvolve(p,a)[n-1:2*n-1]              # sum_m p_m a_{j+m}
    g3=fftconvolve(a,p[::-1])[n-1:2*n-1]        # sum_m p_m a_{j-m}
    g=w*(g1-g2-g3)
    r_=w*a.sum()-1; pen=MU[0]*r_+LAM*r_**2; g=g+(MU[0]+2*LAM*r_)*w
    return val+pen, g
best=proj(a.copy()); bs=score(best); t0=time.time(); x=best.copy(); beta=beta0
while beta<3e10 and time.time()-t0<budget:
    for _al in range(8):
        r=minimize(f,x,args=(beta,),jac=True,method='L-BFGS-B',bounds=[(0,1)]*n,options={'maxiter':3000,'maxcor':30,'ftol':1e-16,'gtol':1e-14})
        x=r.x; MU[0]+=2*LAM*(w*x.sum()-1)
        if abs(w*x.sum()-1)<1e-13: break
    y=proj(x.copy()); s=score(y)
    if s<bs: bs=s; best=y.copy(); np.savetxt(out,best,fmt='%.17g')
    print(f'beta {beta:.1e} score {s:.12f} best {bs:.12f} t {time.time()-t0:.0f}',flush=True); beta*=1.6
print('final',bs,'mass',w*best.sum(),flush=True)
