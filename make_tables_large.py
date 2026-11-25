# LaTeX tables for the large-system section of QMCofEverything_v17:
#   tab:lr  -- fermion-spin-boson chain with 1/r interactions (lr/), references in lr/ref.out, lr/ref256.out
#   tab:qlm -- U(1) quantum link model (qlm/), references in qlm/ref.out, qlm/ref256.out
#   tab:nl  -- efficiency of nonlocal moves per CPU second (nl/)
import re, json, math, os
def res(path):
    if not os.path.exists(path): return None
    t = open(path).read(); m = re.search(r"RESULT (.*)", t)
    if not m: return None
    m = m.group(1).split(); d = {m[i]: (float(m[i+1]), float(m[i+2])) for i in range(0, len(m), 3)}
    for k, pat in [('cpu', r"cpu_ms_per_sweep ([\d.]+)"), ('qmax', r"q_max (\d+)"), ('tau', r"tau_int_E ([\d.]+)")]:
        d[k] = float(re.search(pat, t).group(1))
    return d
def pe(v, e):
    if e <= 0: return f"{v:.4f}"
    nd = max(0, 1 - int(math.floor(math.log10(e)))); return f"{v:.{nd}f}({round(e*10**nd)})"
def refs(files, key):
    out = {}
    for fn in files:
        if not os.path.exists(fn): continue
        for line in open(fn):
            if not line.startswith('{'): continue
            r = json.loads(line); d = r.get(key) or r.get('ed') or r.get('typ') or r.get('ref')
            out[(r['kw']['L'], r['beta'])] = (d, 'ed' in r or all((not isinstance(v, list)) or v[1] == 0 for v in d.values()))
    return out
def fmt_ref(d, k):
    v = d[k]
    if isinstance(v, list): return f"${v[0]:.4f}$" if v[1] == 0 else f"${pe(v[0], v[1])}$"
    return f"${v:.4f}$"
L_ = []; dev = []
# ---------------- long-range chain
R = refs(['lr/ref.out', 'lr/ref256.out'], 'typ')   # later files override earlier ones
L_.append(r"""\begin{table*}[t]
\caption{Fermion--spin--boson chain with long-range interactions, Eq.~\eqref{eq:H_mix} plus Eq.~\eqref{eq:H_V} with $V=1$, open boundaries, $g=0.6/\sqrt L$, $\beta=2$, other parameters as in Table~\ref{tab:mix_open}. Reference values: full ED for $L=4$, quantum typicality with 256 random vectors for $L=6,8$ (errors in parentheses); for $L=12,16$ no reference is available. $\dim$ is the dimension of the $N=L$ sector. QMC: all moves, $2\times10^6$ sweeps; $q_{\max}$ is the largest expansion order reached, and the last column the single-core CPU time per sweep.}
\label{tab:lr}
\small
\begin{tabular}{r r|ll|ll|ll|ll|r|r}
\hline\hline
$L$ & $\dim$ & $\langle H\rangle_{\rm QMC}$ & ref. & $\langle n_d\rangle_{\rm QMC}$ & ref. & $\langle a^\dagger a\rangle_{\rm QMC}$ & ref. & $S_{\rm CDW}^{\rm QMC}$ & ref. & $q_{\max}$ & ms\\
\hline""")
dims = {4: 1890, 6: 24948, 8: 347490, 12: 2704156*27, 16: 601080390*27}
for L in [4, 6, 8, 12, 16]:
    r = res(f'lr/runs/LR_L{L}.out')
    if r is None: L_.append(f"{L} & & \\multicolumn{{10}}{{c}}{{(running)}} \\\\"); continue
    ref = R.get((L, 2.0)); D = dims[L]
    row = f"{L} & " + (f"{D:,}".replace(',', '{,}') if D < 10**7 else f"${D/10**int(math.log10(D)):.1f}\\times10^{{{int(math.log10(D))}}}$")
    for k in ['E', 'docc', 'nph', 'Scdw']:
        row += f" & ${pe(*r[k])}$ & " + (fmt_ref(ref[0], k) if ref else "--")
        if ref:
            v = ref[0][k]; rv, re_ = (v if isinstance(v, list) else (v, 0.0))
            dev.append((abs(r[k][0] - rv) / math.hypot(r[k][1], re_), f'LR{L}', k))
    row += f" & {int(r['qmax'])} & {r['cpu']:.2f} \\\\"
    L_.append(row)
L_.append(r"""\hline\hline
\end{tabular}
\end{table*}
""")
# ---------------- quantum link model
Rq = refs(['qlm/ref.out', 'qlm/ref256.out'], 'ref')
def dimq(L, S=1):
    from functools import lru_cache
    @lru_cache(None)
    def f(x, E):
        if x == L: return 1 if E == 0 else 0
        tot = 0
        for n in (0, 1):
            Ex = E + n - (1 - (-1) ** x) // 2
            if x < L - 1 and abs(Ex) > S: continue
            if x == L - 1 and Ex != 0: continue
            tot += f(x + 1, Ex) if x < L - 1 else 1
        return tot
    return f(0, 0)
L_.append(r"""\begin{table*}[t]
\caption{U(1) quantum link model with staggered fermions and spin-$1$ links, Eq.~\eqref{eq:H_qlm}, $w=1$, $m=0.3$, $J=0.5$, $\theta=0$, open chain with zero boundary fields. $\dim$ is the number of gauge-invariant states. Reference values: ED in the gauge-invariant sector for $L\le16$, quantum typicality for $L=20,24$; none for $L\ge32$. QMC: all moves, $2\times10^6$ sweeps ($10^6$ for $L\ge32$ and at $\beta=4$; for $L=16$, $\beta=4$ two independent runs of $10^6$ and $2\times10^6$ sweeps are combined). $\langle\bar\psi\psi\rangle=L^{-1}\sum_x(-1)^x\langle n_x\rangle$ is the chiral condensate and $\langle E^2\rangle$ the mean squared electric field per link.}
\label{tab:qlm}
\small
\begin{tabular}{r l r|ll|ll|ll|r|r}
\hline\hline
$L$ & $\beta$ & $\dim$ & $\langle H\rangle_{\rm QMC}$ & ref. & $\langle\bar\psi\psi\rangle_{\rm QMC}$ & ref. & $\langle E^2\rangle_{\rm QMC}$ & ref. & $q_{\max}$ & ms\\
\hline""")
for L, b in [(8, 2), (12, 2), (16, 2), (20, 2), (24, 2), (32, 2), (48, 2), (64, 2), (16, 4), (24, 4)]:
    r = res(f'qlm/runs/Q_L{L}_b{b}.out'); D = dimq(L)
    r2 = res(f'qlm/runs/Q_L{L}_b{b}_s2.out')          # independent rerun (other seed): combine by inverse variance
    if r and r2:
        for k in r:
            if isinstance(r[k], tuple) and r[k][1] > 0 and r2[k][1] > 0:
                w1, w2 = 1 / r[k][1] ** 2, 1 / r2[k][1] ** 2
                r[k] = ((w1 * r[k][0] + w2 * r2[k][0]) / (w1 + w2), (w1 + w2) ** -0.5)
        r['qmax'] = max(r['qmax'], r2['qmax'])
    Ds = f"{D:,}".replace(',', '{,}') if D < 10**7 else f"${D/10**int(math.log10(D)):.1f}\\times10^{{{int(math.log10(D))}}}$"
    if r is None: L_.append(f"{L} & {b} & {Ds} & \\multicolumn{{8}}{{c}}{{(running)}} \\\\"); continue
    ref = Rq.get((L, float(b))); row = f"{L} & {b} & {Ds}"
    for k in ['E', 'chiral', 'E2avg']:
        row += f" & ${pe(*r[k])}$ & " + (fmt_ref(ref[0], k) if ref else "--")
        if ref:
            v = ref[0][k]; rv, re_ = v
            dev.append((abs(r[k][0] - rv) / math.hypot(r[k][1], re_), f'Q{L}b{b}', k))
    row += f" & {int(r['qmax'])} & {r['cpu']:.2f} \\\\"
    L_.append(row)
L_.append(r"""\hline\hline
\end{tabular}
\end{table*}
""")
# ---------------- nonlocal moves
def nlres(n): return res(f'nl/runs/{n}.out')
L_.append(r"""\begin{table}[t]
\caption{Statistical efficiency $1/(\sigma^2 t_{\rm CPU})$ of update sets relative to the baseline (basis shifts, swaps, pair moves, and subcycle replacement), from runs of $2\times10^5$ sweeps. QLM: Table~\ref{tab:qlm} model at $L=32$, $\beta=4$ ($q_{\max}\approx150$); LR: Table~\ref{tab:lr} model at $L=12$. ``root'' and ``$\omega$'' denote the root estimator $A(z_0)$ and the $\omega_k$-weighted estimator, Eq.~\eqref{eq:diag_estimator_avg}, of a diagonal observable: the chiral condensate (QLM) and the CDW structure factor (LR). Ratios carry an uncertainty of about $\pm20\%$.}
\label{tab:nl}
\begin{tabular}{l|ccc|ccc}
\hline\hline
 & \multicolumn{3}{c|}{QLM} & \multicolumn{3}{c}{LR}\\
added moves & $H$ & root & $\omega$ & $H$ & root & $\omega$\\
\hline""")
for v, lab in [('reroot', 're-rooting'), ('rev', 'reversal (random pairs)'), ('revT', 'reversal (targeted)'), ('exch', 'segment exchange'), ('allT', 'all three')]:
    row = lab
    for sysn, ob in [('Q32', 'chiral'), ('LR12', 'Scdw')]:
        b0 = nlres(f'{sysn}_base'); r = nlres(f'{sysn}_{v}')
        if not (b0 and r): row += " & -- & -- & --"; continue
        cpu0 = b0['cpu']; cpu = r['cpu']
        eff = lambda d, k, c: 1 / (d[k][1] ** 2 * c)
        row += ''.join(f" & {eff(r, k, cpu) / eff(b0, k, cpu0):.1f}" for k in ['E', ob + '_root', ob])
    L_.append(row + r" \\")
L_.append(r"""\hline\hline
\end{tabular}
\end{table}
""")
open('tables_large.tex', 'w').write('\n'.join(L_)); print('\n'.join(L_))
dev.sort(reverse=True); print('largest deviations (sigma):', [(round(a, 2), f, k) for a, f, k in dev[:8]], 'n>2:', sum(a > 2 for a, _, _ in dev), 'of', len(dev))
