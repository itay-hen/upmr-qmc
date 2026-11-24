# Lanczos thermal reference for the spin-boson chain (Table D, QMCofEverything_v17).
# Uses the conserved local parities P_i = sigma^z_i (-1)^{n_i} (16 sectors for Ns=4), takes the
# lowest k eigenpairs per sector with ARPACK (eigsh), and forms Boltzmann averages. The truncation
# is checked by doubling k. Usage: python ed_lanczos.py Ns k
import numpy as np, scipy.sparse as sp, sys, itertools
from scipy.sparse.linalg import eigsh
def build(Ns,eps=1.,om=1.,g=0.6,J=0.5,nmax=10):
    d=nmax+1; a=sp.diags(np.sqrt(np.arange(1,d)),1); n=sp.diags(np.arange(d)*1.); x=(a+a.T).tocsr()
    sx=sp.csr_matrix([[0.,1],[1,0]]); sz=sp.diags([1.,-1]); loc=[sp.identity(2),sp.identity(d)]*Ns
    def op(o):
        M=sp.identity(1,format='csr')
        for c in range(2*Ns): M=sp.kron(M,o.get(c,loc[c]),format='csr')
        return M
    H=0; Hoff=0; SZ=0; N=0
    for i in range(Ns):
        H=H+eps/2*op({2*i:sz})+om*op({2*i+1:n}); SZ=SZ+op({2*i:sz})/Ns; N=N+op({2*i+1:n})/Ns
        c=g*op({2*i:sx,2*i+1:x}); H=H+c; Hoff=Hoff+c
    for i in range(Ns-1): H=H+J*op({2*i:sz,2*i+2:sz})
    # local parity labels of each basis state (component order s0,n0,s1,n1,...)
    idx=np.indices([2,d]*Ns).reshape(2*Ns,-1)
    par=np.zeros(idx.shape[1],dtype=int)
    for i in range(Ns): par|=((idx[2*i]+idx[2*i+1])%2)<<i
    return H.tocsr(),Hoff.tocsr(),SZ.diagonal(),N.diagonal(),par
def spectrum(Ns,k):
    H,Hoff,sz,nn,par=build(Ns); out=[]
    for s in range(2**Ns):
        I=np.where(par==s)[0]; Hs=H[I][:,I]; Os=Hoff[I][:,I]
        kk=min(k,len(I)-2)
        w,v=eigsh(Hs,k=kk,which='SA',tol=1e-12)
        out.append((w,(v*(sz[I][:,None]*v)).sum(0),(v*(nn[I][:,None]*v)).sum(0),(v*(Os@v)).sum(0)))
    W=np.concatenate([o[0] for o in out]); cols=[np.concatenate([o[j] for o in out]) for j in (1,2,3)]
    return W,cols,[o[0].max() for o in out]
def thermo(W,cols,beta):
    p=np.exp(-beta*(W-W.min())); p/=p.sum()
    return dict(E=float(p@W),sz=float(p@cols[0]),n=float(p@cols[1]),q=float(-beta*(p@cols[2])))
if __name__=="__main__":
    Ns=int(sys.argv[1]); k=int(sys.argv[2]); beta=5.
    W,cols,tops=spectrum(Ns,k)
    # weight of the highest retained level in each sector bounds the truncation
    tail=max(np.exp(-beta*(t-W.min())) for t in tops)/np.exp(-beta*(W-W.min())).sum()
    print(("D",Ns,thermo(W,cols,beta),None),flush=True)
    print(f"# Ns={Ns} k={k} E0={W.min():.10f} max_tail_weight={tail:.2e}",flush=True)
