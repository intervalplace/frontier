# Sharp LSE ramp on a large C1 step function (a = x^2). Args: src r b0 b1 fac its out [sig seed]
import os; os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np, sys, time
from scipy.signal import fftconvolve
from scipy.optimize import minimize
src,r,b0,b1,fac,its,out=sys.argv[1],int(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]),int(sys.argv[6]),sys.argv[7]
a=np.repeat(np.load(src) if src.endswith(".npy") else np.loadtxt(src),r); n=len(a)
if len(sys.argv)>8: a=a*np.exp(float(sys.argv[8])*np.random.default_rng(int(sys.argv[9])).normal(size=n))
a*=n/a.sum()
sc=lambda a: 2*len(a)*fftconvolve(a,a).max()/a.sum()**2
def obj(x,beta):
    a=x*x; S=a.sum(); c=fftconvolve(a,a)*(2*n/S**2); m=c.max(); w=np.exp(beta*(c-m)); Z=w.sum()
    p=w/Z; g=(2*n/S**2)*2*fftconvolve(p,a[::-1])[n-1:2*n-1]-2*np.dot(p,c)/S
    return m+np.log(Z)/beta, g*2*x
x=np.sqrt(a); best=sc(a); t0=time.time(); beta=b0
print('start',n,repr(best),flush=True)
while beta<=b1:
    x=minimize(obj,x,args=(beta,),jac=True,method='L-BFGS-B',options={'maxiter':its,'gtol':1e-18,'ftol':1e-20,'maxcor':30}).x
    a=x*x; v=sc(a); tag=''
    if v<best: best=v; np.save(out,a*n/a.sum()); tag='BEST'
    print(f'{time.time()-t0:.0f}s beta {beta:.3g} score {v:.13f} best {best:.13f} {tag}',flush=True)
    beta*=fac
print('final',repr(best))
