from pathlib import Path
import re,json,numpy as np,sys
from scipy.linalg import eigh
R=Path(__file__).parent;orig=R.parent/'runs'
checks={}
def read(f):
 s=f.read_text();a=np.array([[float(x) for x in l.split()[2:]] for l in s.splitlines() if l.startswith('RAW_BIN')]);c=a[:,0];raw=a[:,1:];tot=raw.sum(0);jk=(tot-raw)/(tot[0]-raw[:,0])[:,None]
 est=tot/tot[0];err=np.sqrt((len(a)-1)/len(a)*((jk-jk.mean(0))**2).sum(0));return s,c,raw,est,err
for fn,of in [('I_new_seed41.txt','W_sb.out'),('I_h3_worm_corrected.txt','C_worm.out'),('I_h3_all_corrected.txt','C_all.out')]:
 s,c,r,m,e=read(R/fn);old=(orig/of).read_text();accmatch=re.findall(r'  acc .*',s)==re.findall(r'  acc .*',old);assert accmatch
 fields=re.search('RESULT.*',s)[0].split()[1:];res={fields[i]:(float(fields[i+1]),float(fields[i+2])) for i in range(0,len(fields),3)}
 for name,j in [('E',1),('szroot',2),('szw',3),('n',4),('q',5)]:
  assert abs(res[name][0]-m[j])<6e-6 and abs(res[name][1]-e[j])<6e-6
 checks[fn]={'matches_original_acceptance_counts':accmatch,'raw_bins':len(c),'closed_total':int(sum(c)),'observables':{n:{'mean':float(m[j]),'jackknife_SEM':float(e[j])} for n,j in [('E',1),('szroot',2),('szw',3),('n',4),('q',5)]}}
# independent beta=1 dense 18-state reference for a new run with fully recorded settings.
a=np.diag(np.sqrt(np.arange(1,9)),1);z=np.diag([1.,-1]);x=np.array([[0.,1],[1,0]]);off=.6*np.kron(x,a+a.T);H=np.kron(z/2,np.eye(9))+np.kron(np.eye(2),np.diag(np.arange(9)))+off
w,v=eigh(H);p=np.exp(-(w-w.min()));p/=p.sum();rho=(v*p)@v.T
ref={'E':float(w@p),'szw':float(np.trace(np.kron(z,np.eye(9))@rho)),'n':float(np.trace(np.kron(np.eye(2),np.diag(np.arange(9)))@rho)),'q':float(-np.trace(off@rho))}
s,c,r,m,e=read(R/'I_beta1_worm_corrected.txt');checks['new_beta1_seed45']={}
for n,j in [('E',1),('szw',3),('n',4),('q',5)]:checks['new_beta1_seed45'][n]={'mean':float(m[j]),'SEM':float(e[j]),'ED':ref[n],'z':float((m[j]-ref[n])/e[j])}
# generic complex-phase driver: phase means and errors must use all closed visits.
s,c,r,_,_=read(R/'I_generic_rawsum_check.txt');tot=r.sum(0);B=len(c)
vals=[]
for j in (0,1):
 mean=tot[j]/sum(c);jk=(tot[j]-r[:,j])/(sum(c)-c);err=np.sqrt((B-1)/B*sum((jk-jk.mean())**2));vals.extend([mean,err])
reported=[float(x) for x in re.search(r'PHASE re (\S+) (\S+) im (\S+) (\S+)',s).groups()]
assert np.max(abs(np.array(vals)-reported))<6e-7
checks['generic_phase_rawsum_check']={'computed':vals,'reported':reported,'PASS':True}
(R/'I_new_code_checks.json').write_text(json.dumps(checks,indent=2))
# correct every table entry affected by the bin-weight bug, not only the long worm paragraph.
pth=R/'QMC of everything/QMCofEverything_v18.tex';t=pth.read_text()
t=t.replace('local + worm & $-0.35620(78)$ & $0.8958(43)$ & $-1.1927(16)$ & $3.8167(79)$', 'local + worm & $-0.35630(78)$ & $0.8948(43)$ & $-1.1927(15)$ & $3.8150(79)$')
t=t.replace('all & $-0.35709(77)$ & $0.8860(35)$ & $-1.1952(19)$ & $3.8039(80)$','all & $-0.35717(77)$ & $0.8853(34)$ & $-1.1952(19)$ & $3.8027(79)$')
t=t.replace('for different update sets ($10^7$ sweeps each).','for different update sets ($10^7$ sweeps each). Worm-containing rows are recomputed from raw numerator/denominator sums in replays with identical acceptance counts.')
t=t.replace('Per-run uncertainties use the raw-sum jackknife,','For this controlled comparison, per-run uncertainties use 100 equal-sweep replay blocks and the raw-sum jackknife,')
zmax=max(abs(u['z']) for u in checks['new_beta1_seed45'].values())
t=t.replace(r'At $\beta=1$ the worm run agrees with ED within $1.2\sigma$.',f'A new $\\beta=1$ worm run (seed 45, $10^7$ sweeps, $10^5$ thermalization sweeps, 100 raw-sum bins) gives $\\langle H\\rangle=0.04696(250)$ and $\\langle\\sigma^z\\rangle=-0.37570(70)$; all four tested observables agree with independent ED within ${zmax:.2f}$ standard errors.')
pth.write_text(t);print(json.dumps({k:v for k,v in checks.items() if k!='generic_phase_rawsum_check'},indent=2))
