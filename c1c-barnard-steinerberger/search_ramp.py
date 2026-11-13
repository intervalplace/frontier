# C1c: maximise K min_{0<=k<=K} A_k / S^2, A_k = sum_i w_i w_{i+k}, via a slow soft-min ramp.
# Parametrise w = x^2; objective: -softmin_beta(K A_k / S^2) over k = 0..K (k < n).
import numpy as np, json, sys, time
from scipy.signal import fftconvolve
from scipy.optimize import minimize
src=sys.argv[1]; Kt=int(sys.argv[2]); b0=float(sys.argv[3]); b1=float(sys.argv[4]); fac=float(sys.argv[5]); its=int(sys.argv[6]); out=sys.argv[7]
d=json.load(open(src)) if src.endswith('.json') else None
if d: K0=d['K']; w=np.array(d['weights_int'],float)
else:
    t=open(src).read().split(); K0=int(t[0]); w=np.array(t[1:],float)
r=K0//Kt; assert K0==r*Kt, 'grid must nest'
n=len(w)//r*r if len(w)%r==0 else (len(w)+r-1)//r*r
w=np.concatenate([w,np.zeros(n-len(w))]).reshape(-1,r).sum(1); K=Kt; n=len(w)
def corr(w): return fftconvolve(w,w[::-1])[n-1:n-1+K+1] if K+1<=n else np.concatenate([fftconvolve(w,w[::-1])[n-1:],np.zeros(K+1-n)])
def score(w): return K*corr(w).min()/w.sum()**2
def obj(x,beta):
    w=x*x; S=w.sum(); A=corr(w); v=K*A/S**2; m=v.min(); e=np.exp(-beta*(v-m)); Z=e.sum(); f=m-np.log(Z)/beta
    p=e/Z
    P=np.zeros(n); m_=min(K+1,n); P[:m_]=p[:m_]
    t1=fftconvolve(w[::-1],P)[:n][::-1]; t2=fftconvolve(w,P)[:n]
    g=K/S**2*(t1+t2) - 2*np.dot(p,v)/S
    return -f, -g*2*x
x=np.sqrt(w); best=score(w); print('start n',n,'K',K,f'{best:.13f}',flush=True); t0=time.time(); beta=b0
np.savetxt(out,w,fmt='%.17g')
while beta<=b1:
    x=minimize(obj,x,args=(beta,),jac=True,method='L-BFGS-B',options={'maxiter':its,'gtol':1e-16,'ftol':1e-18}).x
    v=score(x*x)
    if v>best: best=v; np.savetxt(out,x*x,fmt='%.17g')
    print(f'{time.time()-t0:.0f}s beta {beta:.3g} score {v:.13f} best {best:.13f}',flush=True); beta*=fac
print('final',f'{best:.13f}')
