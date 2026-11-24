# Builds the LaTeX benchmark tables for QMCofEverything_v19 from runs/*.out and ed_ref.out.
# The worm-containing rows of the three-label table (C_worm, C_all) are taken from the raw-sum
# replays in revision_checks/ (same seeds and parameters, identical acceptance counts): the
# original runs/ outputs averaged bin means with equal weights, which misweights bins with
# different numbers of closed visits. Runs without the worm have no open sector, so every
# measurement sweep ends closed and the two estimators coincide.
import re, ast, math
OVERRIDE={'C_worm':'revision_checks/I_h3_worm_corrected.txt','C_all':'revision_checks/I_h3_all_corrected.txt'}
def res(f):
    t=open(OVERRIDE.get(f,'runs/'+f+'.out')).read(); m=re.search(r"RESULT (.*)",t).group(1).split()
    d={m[i]:(float(m[i+1]),float(m[i+2])) for i in range(0,len(m),3)}
    s=re.search(r"steps_per_sweep (\d+) closed_frac ([\d.]+)\s+ms_per_sweep ([\d.]+)\s+tau_int_E ([\d.]+)",t)
    d['meta']=dict(steps=int(s.group(1)),closed=float(s.group(2)),ms=float(s.group(3)),tau=float(s.group(4)))
    rb=[[float(x) for x in l.split()[3:]] for l in t.splitlines() if l.startswith('RAW_BIN')]
    if rb:
        # full-precision ratio estimates and delete-one-bin jackknife from the raw sums
        # (columns: sum w, E, szroot, szw, n, q, ...), avoiding double rounding of RESULT
        B=len(rb); tot=[sum(r[i] for r in rb) for i in range(len(rb[0]))]
        for j,k in enumerate(['E','szroot','szw','n','q'],start=1):
            jk=[(tot[j]-r[j])/(tot[0]-r[0]) for r in rb]; mj=sum(jk)/B
            d[k]=(tot[j]/tot[0], math.sqrt((B-1)/B*sum((x-mj)**2 for x in jk)))
    return d
def cpu(n):
    # single-core CPU time per sweep from the short timing reruns (timing/), same parameters
    return float(re.search(r"cpu_ms_per_sweep ([\d.]+)", open('timing/'+n+'.out').read()).group(1))
ED={}
for line in open('ed_ref.out'):
    k,x,r,c=ast.literal_eval(line); ED[(k,x)]=dict(r,chi=c)
def pe(v,e):
    # value(error) with 2 significant digits of error
    if e==0: return f"{v:.5f}"
    nd=max(0,1-int(math.floor(math.log10(e))))
    ee=round(e*10**nd)
    return f"{v:.{nd}f}({ee})"
def f5(x): return f"{x:.5f}" if abs(x)<10 else f"{x:.3f}"
L=[]
# Table A
L.append(r"""\begin{table*}[t]
\caption{UPMR-QMC versus ED for a single spin--boson pair, Eq.~\eqref{eq:H_SB} with $\epsilon=\omega_0=1$, $g=0.6$, $n_{\max}=8$. QMC: $2\times10^7$ sweeps with local and nonlocal moves; $\langle\sigma^z\rangle$ from the weighted estimator~\eqref{eq:diag_estimator_avg}; errors (one standard deviation, jackknife over 100 bins) in the last digits.}
\label{tab:sb_ed}
\small
\begin{tabular}{r|ll|ll|ll|ll}
\hline\hline
$\beta$ & $\langle\sigma^z\rangle_{\rm QMC}$ & $\langle\sigma^z\rangle_{\rm ED}$ & $\langle n\rangle_{\rm QMC}$ & $\langle n\rangle_{\rm ED}$ & $\langle H\rangle_{\rm QMC}$ & $\langle H\rangle_{\rm ED}$ & $\langle q\rangle_{\rm QMC}$ & $-\beta\langle H_{\rm off}\rangle_{\rm ED}$\\
\hline""")
for b,f in [(0.1,'A_b0.1'),(0.5,'A_b0.5'),(1,'A_b1'),(2,'A_b2'),(5,'A_b5'),(10,'A_b10')]:
    r=res(f); e=ED[('A',b)]
    L.append(f"{b} & ${pe(*r['szw'])}$ & ${f5(e['sz'])}$ & ${pe(*r['n'])}$ & ${f5(e['n'])}$ & ${pe(*r['E'])}$ & ${f5(e['E'])}$ & ${pe(*r['q'])}$ & ${f5(e['q'])}$ \\\\")
L.append(r"""\hline\hline
\end{tabular}
\end{table*}
""")
# Table B
L.append(r"""\begin{table*}[t]
\caption{UPMR-QMC versus ED for a single spin--boson pair at $\beta=5$, $\epsilon=\omega_0=1$, $n_{\max}=10$, as a function of $g$. The susceptibility $\chi_z=\partial\langle\sigma^z\rangle/\partial\epsilon$ is measured with the covariance estimator described in the text; ED values are from a finite difference. Values for $n_{\max}=40$ agree with the ED column except at $g=2$ (see text). QMC: $2\times10^7$ sweeps ($5\times10^6$ at $g=2$, $2\times10^6$ at $g=0$).}
\label{tab:sb_g}
\small
\begin{tabular}{r|ll|ll|ll|ll}
\hline\hline
$g$ & $\langle\sigma^z\rangle_{\rm QMC}$ & $\langle\sigma^z\rangle_{\rm ED}$ & $\langle n\rangle_{\rm QMC}$ & $\langle n\rangle_{\rm ED}$ & $\langle q\rangle_{\rm QMC}$ & $-\beta\langle H_{\rm off}\rangle_{\rm ED}$ & $\chi_{z,\rm QMC}$ & $\chi_{z,\rm ED}$\\
\hline""")
for g,f in [(0,'B_g0'),(0.5,'B_g0.5'),(0.8,'B_g0.8'),(1,'B_g1'),(2,'B_g2')]:
    r=res(f); e=ED[('B',g)]
    qq = "$0$" if g==0 else f"${pe(*r['q'])}$"
    L.append(f"{g} & ${pe(*r['szw'])}$ & ${f5(e['sz'])}$ & ${pe(*r['n'])}$ & ${f5(e['n'])}$ & {qq} & ${f5(abs(e['q']))}$ & ${pe(*r['chi'])}$ & ${f5(e['chi'])}$ \\\\")
L.append(r"""\hline\hline
\end{tabular}
\end{table*}
""")
# Table C
e=ED[('C',2)]
L.append(r"""\begin{table*}[t]
\caption{Three-label benchmark, Eq.~\eqref{eq:H3} with $\epsilon=\omega_0=1$, $h=0.5$, $\kappa=0.3$, $g=0.6$, $\beta=2$, $n_{\max}=40$, for different update sets ($10^7$ sweeps each). Worm-containing rows are recomputed from raw numerator/denominator sums in replays with identical acceptance counts. ``Local'': basis-state, swap, and pair moves. ``All'': local, subcycle, re-rooting, reversal, segment exchange, and worm.}
\label{tab:h3}
\begin{tabular}{lllll}
\hline\hline
updates & $\langle\sigma^z\rangle$ & $\langle n\rangle$ & $\langle H\rangle$ & $\langle q\rangle$\\
\hline""")
for lab,f in [("local","C_base"),("local + subcycle","C_sub"),("local + worm","C_worm"),("all","C_all")]:
    r=res(f)
    L.append(f"{lab} & ${pe(*r['szw'])}$ & ${pe(*r['n'])}$ & ${pe(*r['E'])}$ & ${pe(*r['q'])}$ \\\\")
L.append(f"ED & ${f5(e['sz'])}$ & ${f5(e['n'])}$ & ${f5(e['E'])}$ & ${f5(e['q'])}$ \\\\")
L.append(r"""\hline\hline
\end{tabular}
\end{table*}
""")
# Table D
L.append(r"""\begin{table*}[t]
\caption{Chains of $N_s$ spin--boson pairs, Eq.~\eqref{eq:H_SB} with $\epsilon=\omega_0=1$, $g=0.6$, $J=0.5$ (open chain), $n_{\max}=10$, $\beta=5$; $4\times10^6$ sweeps per run. ``loc'': basis-state, swap, and pair moves; ``nl'': in addition re-rooting, inner-loop reversal, and segment exchange. ``steps'' is the number of elementary steps per sweep (for the runs with nonlocal moves); $\sigma(\cdot)$ denotes the statistical error of the root estimator of $\langle\sigma^z\rangle$ and of $\langle n\rangle$ (per site). Timings are single-core CPU (process) times per sweep, measured in separate short runs with the same parameters. ED values for $N_s=4$ ($\dim\mathcal H=234{,}256$) are from Lanczos diagonalization in the 16 sectors of the conserved local parities $\sigma^z_i(-1)^{n_i}$, keeping the lowest 160 levels per sector.}
\label{tab:scaling}
\footnotesize
\begin{tabular}{r|lll|ll|r|c|c|c|c}
\hline\hline
$N_s$ & $\langle H\rangle_{\rm loc}$ & $\langle H\rangle_{\rm nl}$ & $\langle H\rangle_{\rm ED}$ & $\langle q\rangle_{\rm nl}$ & $-\beta\langle H_{\rm off}\rangle_{\rm ED}$ & steps & ms/sweep & $\tau_{\rm int}(H)$ & $\sigma(\sigma^z_{\rm root})$ & $\sigma(n_{\rm root})$\\
 & & & & & & & loc / nl & loc / nl & loc / nl & loc / nl\\
\hline""")
for N in [1,2,3,4]:
    a=res(f'D_N{N}_loc'); b=res(f'D_N{N}_nl'); e=ED.get(('D',N))
    edE = f5(e['E']) if e else "{-}"; edq = f5(e['q']) if e else "{-}"
    L.append(f"{N} & ${pe(*a['E'])}$ & ${pe(*b['E'])}$ & ${edE}$ & ${pe(*b['q'])}$ & ${edq}$ & {b['meta']['steps']} & {cpu(f'D_N{N}_loc'):.3f} / {cpu(f'D_N{N}_nl'):.3f} & {a['meta']['tau']:.1f} / {b['meta']['tau']:.1f} & {a['szroot'][1]:.4f} / {b['szroot'][1]:.4f} & {a['n'][1]:.4f} / {b['n'][1]:.4f}\\\\")
L.append(r"""\hline\hline
\end{tabular}
\end{table*}
""")
open('tables.tex','w').write("\n".join(L)); print("\n".join(L))
for N in [2,3]:
    r=res(f'D_N{N}_nl'); print(N, 'sz', r['szw'], ED[('D',N)]['sz'], 'n', r['n'], ED[('D',N)]['n'])
