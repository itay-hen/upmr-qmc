# Exact-diagonalization references for the UPMR benchmarks (QMCofEverything_v17).
import numpy as np, sys
from scipy.linalg import eigh
def sb_chain(Ns,eps,om,g,J,nmax,h=0.,kap=0.,h3=False):
    d=nmax+1; a=np.diag(np.sqrt(np.arange(1,d)),1); n=a.T@a; x=a+a.T
    sx=np.array([[0.,1],[1,0]]); sz=np.diag([1.,-1]); I2=np.eye(2); Ib=np.eye(d)
    loc=[I2,Ib]*Ns
    def op(ops):  # ops: dict comp->matrix
        M=np.array([[1.]])
        for c in range(2*Ns): M=np.kron(M,ops.get(c,loc[c]))
        return M
    H=0; Hoff=0; SZ=0; N=0; dSZ=0
    for i in range(Ns):
        H=H+eps/2*op({2*i:sz})+om*op({2*i+1:n}); SZ=SZ+op({2*i:sz})/Ns; N=N+op({2*i+1:n})/Ns
        c=(-g if h3 else g)*op({2*i:sx,2*i+1:x})-h*op({2*i:sx})-kap*op({2*i+1:x})
        H=H+c; Hoff=Hoff+c
    for i in range(Ns-1): H=H+J*op({2*i:sz,2*i+2:sz})
    return H,Hoff,SZ,N
def thermo(H,Hoff,SZ,N,beta):
    w,v=eigh(H); p=np.exp(-beta*(w-w.min())); p/=p.sum()
    ex=lambda O: float(np.sum(p*(v*(O@v)).sum(0)))
    return dict(E=float(np.sum(p*w)),sz=ex(SZ),n=ex(N),q=-beta*ex(Hoff))
def chi(Ns,eps,om,g,J,nmax,beta,d=1e-4):
    a=thermo(*sb_chain(Ns,eps+d,om,g,J,nmax),beta)['sz']; b=thermo(*sb_chain(Ns,eps-d,om,g,J,nmax),beta)['sz']
    return (a-b)/(2*d)
if __name__=="__main__":
    
    for b in [0.1,0.5,1,2,5,10]:
        r=thermo(*sb_chain(1,1,1,0.6,0,8),b); print(("A",b,r,chi(1,1,1,0.6,0,8,b)), flush=True)
    for g in [0,0.5,0.8,1,2]:
        r=thermo(*sb_chain(1,1,1,g,0,10),5); print(("B",g,r,chi(1,1,1,g,0,10,5)), flush=True)
    r=thermo(*sb_chain(1,1,1,0.6,0,40,h=0.5,kap=0.3,h3=True),2); print(("C",2,r,None), flush=True)
    for Ns in [1,2,3]:
        r=thermo(*sb_chain(Ns,1,1,0.6,0.5,10),5); print(("D",Ns,r,None), flush=True)
    pass
