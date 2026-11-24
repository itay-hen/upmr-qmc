# Accuracy test of log_edd (via ./upmr_qmc --test) against 2000-digit divided differences.
# Validated range: q <= 160, beta <= 10 (max abs error in ln|f| ~1e-13). Beyond q ~ 180 at
# beta >= 5 the double-precision Taylor scheme runs out of range and the samplers abort.
import mpmath as mp, random, subprocess
mp.mp.dps=2000
def dd(xs,beta):
    # divided differences of exp(-beta x) by recursion in high precision
    f=[mp.e**(-beta*x) for x in xs]; n=len(xs)
    for k in range(1,n):
        for i in range(n-1,k-1,-1):
            den=xs[i]-xs[i-k]
            f[i]=(f[i]-f[i-1])/den if den!=0 else None
    return f[-1]
random.seed(1); worst=0
for q in [10,40,80,120,160]:
    for beta in [1,5,10]:
        for rep in range(3):
            xs=[mp.mpf(round(random.uniform(-3,3),6))+mp.mpf(i)*mp.mpf('1e-7') for i in range(q+1)]  # distinct
            ref=mp.log(abs(dd(xs,beta)))
            out=subprocess.run(['./upmr_qmc','--test',' '.join([str(beta)]+[mp.nstr(x,20) for x in xs])],capture_output=True,text=True).stdout
            err=abs(float(out)-float(ref)); worst=max(worst,err)
    print(q,'max abs err in ln|f| so far',worst,flush=True)
