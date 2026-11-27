# Update-set comparison from 8 independent runs per set (ab/runs/QS_*, LS_*).
# sigma_jk: mean per-run jackknife error; sigma_seed: scatter between runs / sqrt(1) (per-run).
# R = sigma_seed / sigma_jk (about 1 when the jackknife is reliable); efficiency 1/(sigma^2 t_cpu)
# relative to re-rooting ("rr"), using sigma = max(sigma_seed, sigma_jk) per run.
import re, glob, numpy as np, json, sys
SYS = {'QS': ('chiral', {'E': -12.8202, 'chiral': -0.1999}), 'LS': ('Scdw', {'E': -8.694, 'Scdw': 0.5984})}
VAR = [('base', 'local + subcycle'), ('rr', '+ re-rooting'), ('rev', '+ re-rooting + reversal'), ('exch', '+ re-rooting + exchange'), ('rrsep', '+ re-rooting + separated pairs'), ('all', '+ all of the above')]
out = {}
for P, (ob, ref) in SYS.items():
    for v, _ in VAR:
        fs = sorted(glob.glob(f'ab/runs/{P}_{v}_*.out')); val = {k: [] for k in ref}; jk = {k: [] for k in ref}; cpu = []
        for f in fs:
            t = open(f).read(); m = re.search(r"RESULT (.*)", t).group(1).split(); d = {m[i]: (float(m[i+1]), float(m[i+2])) for i in range(0, len(m), 3)}
            cpu.append(float(re.search(r"cpu_ms_per_sweep ([\d.]+)", t).group(1)))
            for k in ref: val[k].append(d[k][0]); jk[k].append(d[k][1])
        r = {'n': len(fs), 'cpu': np.mean(cpu)}
        for k in ref:
            x = np.array(val[k]); sj = np.sqrt(np.mean(np.square(jk[k]))); ss = x.std(ddof=1)
            r[k] = dict(mean=x.mean(), se=ss / np.sqrt(len(x)), sj=sj, ss=ss, R=ss / sj, z=(x.mean() - ref[k]) / (ss / np.sqrt(len(x))))
        out[(P, v)] = r
for P, (ob, ref) in SYS.items():
    rr = out[(P, 'rr')]
    print(P)
    for v, lab in VAR:
        r = out[(P, v)]; s = f"  {v:6s} cpu {r['cpu']:.3f}"
        for k in ref:
            e = lambda q: 1 / (max(q[k]['ss'], q[k]['sj']) ** 2 * q['cpu'])
            s += f" | {k}: eff {e(r)/e(rr):6.2f} R {r[k]['R']:.2f} mean {r[k]['mean']:.5f}({r[k]['se']:.5f}) z {r[k]['z']:+.1f}"
        print(s)
# ---- LaTeX table
L = [r"""\begin{table}[t]
\caption{Update sets compared on the quantum link model ($L=16$, $\beta=4$, Table~\ref{tab:qlm}) and the long-range chain ($L=8$, Table~\ref{tab:lr}), from 8 independent runs of $5\times10^5$ sweeps per set. ``Eff.'' is the efficiency $1/(\sigma^2t_{\rm CPU})$ relative to the set with re-rooting, with $\sigma$ the larger of the mean jackknife error and the scatter between runs; the ratios are uncertain by about $\pm30\%$. $R$ is the scatter between runs divided by the jackknife error of the $\omega_k$-weighted estimator of $O$, the chiral condensate (QLM) and $S_{\rm CDW}$ (LR); $R\gg1$ signals that the jackknife underestimates the error. Every set includes basis shifts, swaps, pair moves, and subcycle replacement; reversal uses the targeted proposal.}
\label{tab:updates}
\begin{tabular}{l|ccc|ccc}
\hline\hline
 & \multicolumn{3}{c|}{QLM} & \multicolumn{3}{c}{LR}\\
update set & $H$ & $O$ & $R$ & $H$ & $O$ & $R$\\
\hline"""]
for v, lab in VAR:
    row = lab.replace('local + subcycle', 'without re-rooting')
    for P, (ob, ref) in SYS.items():
        r = out[(P, v)]; rr = out[(P, 'rr')]; e = lambda q, k: 1 / (max(q[k]['ss'], q[k]['sj']) ** 2 * q['cpu'])
        row += f" & {e(r,'E')/e(rr,'E'):.1f} & {e(r,ob)/e(rr,ob):.{2 if e(r,ob)/e(rr,ob) < 0.1 else 1}f} & {r[ob]['R']:.1f}"
    L.append(row + r" \\")
L.append(r"""\hline\hline
\end{tabular}
\end{table}""")
open('ab/tab_updates.tex', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
