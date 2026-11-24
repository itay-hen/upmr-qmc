# LaTeX tables for the fermion-spin-boson benchmark (QMCofEverything_v17) from gen/runs and gen/ed_mixed.out
import re, json, math, sys
def res(f):
    t = open('gen/runs/' + f + '.out').read()
    m = re.search(r"RESULT (.*)", t).group(1).split()
    d = {m[i]: (float(m[i+1]), float(m[i+2])) for i in range(0, len(m), 3)}
    p = re.search(r"PHASE re (\S+) (\S+) im (\S+) (\S+)", t); d['phase'] = (float(p.group(1)), float(p.group(2)))
    s = re.search(r"steps_per_sweep (\d+) closed_frac ([\d.]+) ms_per_sweep ([\d.]+) tau_int_E ([\d.]+) weight_check_fail (\d+)", t)
    d['meta'] = dict(steps=int(s.group(1)), ms=float(s.group(3)), tau=float(s.group(4)), wfail=int(s.group(5)))
    h = re.search(r"labels=(\d+) diag_strings=(\d+) templates=(\d+)", t); d['hdr'] = tuple(map(int, h.groups()))
    return d
ED = {}
for line in open('gen/ed_mixed.out'):
    if not line.startswith('{'): continue
    r = json.loads(line); kw = r['kw']
    key = ('open' if kw.get('pbc') is False else 'ring%.4f' % kw.get('phi', 0), r['beta'])
    ED[key] = r['ed']
def pe(v, e):
    if e == 0: return f"{v:.4f}"
    nd = max(0, 1 - int(math.floor(math.log10(e)))); return f"{v:.{nd}f}({round(e*10**nd)})"
def f4(x): return f"{x:.4f}"
OBS = [('E', 'E', r'$\langle H\rangle$'), ('q', 'q', r'$\langle q\rangle$'), ('H_kin', 'H_kin', r'$\langle H_{\rm hop}\rangle$'),
       ('H_ssh', 'H_ssh', r'$\langle H_g\rangle$'), ('H_Kflip', 'H_Kflip', r'$\langle H_{K,\perp}\rangle$'),
       ('H_Sx', 'H_Sx', r'$\langle H_\gamma\rangle$'), ('docc', 'docc', r'$\langle n_d\rangle$'),
       ('nph', 'nph', r'$\langle a^\dagger a\rangle$'), ('Sz2', 'Sz2', r'$\langle (S^z)^2\rangle$'), ('SzSz0', 'SzSz0', r'$\langle S^zs^z_0\rangle$')]
def zs(r, e, k):
    v, s = r[k]; return abs(v - e) / s if s > 0 else 0
out = []; maxdev = []
# Table 1: open chain, beta = 1, 2, 4
cols = [('O_b1', 1.0), ('O_b2', 2.0), ('O_b4', 4.0)]
R = {f: res(f) for f, _ in cols}
out.append(r"""\begin{table*}[t]
\caption{Fermion--spin--boson benchmark, Eq.~\eqref{eq:H_mix}, open chain: QMC versus ED (dimension of the $N=4$ sector: 1890). $L=4$, $t=1$, $U=4$, $\omega_0=1$, $g=0.3$, $J_K=1$, $D=0.5$, $\gamma=0.3$, $S=1$, $n_{\max}=8$. All moves (basis shifts, swap, pair, subcycle replacement with automatically generated templates, re-rooting, reversal, segment exchange), $2\times10^6$ sweeps per $\beta$. $H_{\rm hop}$, $H_g$, $H_{K,\perp}=\tfrac{J_K}{2}(S^+s_0^-+S^-s_0^+)$, and $H_\gamma$ are the off-diagonal parts of $H$ (hopping, boson-assisted hopping, transverse Kondo, and impurity--boson terms), measured from vertex counts; $n_d=\sum_i n_{i\uparrow}n_{i\downarrow}$ is the total double occupancy; the remaining diagonal observables use the $\omega_k$-weighted estimator. Every sampled weight was nonnegative.}
\label{tab:mix_open}
\begin{tabular}{l|ll|ll|ll}
\hline\hline
 & \multicolumn{2}{c|}{$\beta=1$} & \multicolumn{2}{c|}{$\beta=2$} & \multicolumn{2}{c}{$\beta=4$}\\
 & QMC & ED & QMC & ED & QMC & ED\\
\hline""")
for k, ek, lab in OBS:
    row = lab
    for f, b in cols:
        e = ED[('open', b)][ek]; row += f" & ${pe(*R[f][k])}$ & ${f4(e)}$"; maxdev.append((zs(R[f], e, k), f, k))
    out.append(row + r" \\")
out.append(r"""\hline\hline
\end{tabular}
\end{table*}
""")
# Table 2: ring with flux
rcols = [('R0_b1', 0.0, 1.0, r'$0$', 1), ('Rh_b1', math.pi/2, 1.0, r'$\pi/2$', 1), ('Rh_b2', math.pi/2, 2.0, r'$\pi/2$', 2),
         ('Rp_b1', math.pi, 1.0, r'$\pi$', 1), ('Rp_b2', math.pi, 2.0, r'$\pi$', 2), ('Rp_b4', math.pi, 4.0, r'$\pi$', 4)]
RR = {c[0]: res(c[0]) for c in rcols}
out.append(r"""\begin{table*}[t]
\caption{Fermion--spin--boson benchmark on the ring (periodic boundary, flux $\phi$), same parameters as Table~\ref{tab:mix_open}. Weights now carry fermionic exchange signs ($\phi=0,\pi$) and complex phases ($\phi=\pi/2$); $\langle e^{i\theta}\rangle$ is the average phase of $w$ under $|w|$ sampling (its imaginary part vanishes within errors). $4\times10^6$ sweeps ($2\times10^6$ at $\beta=4$).}
\label{tab:mix_ring}
\small
\begin{tabular}{ll|l|ll|ll|ll|ll}
\hline\hline
$\phi$ & $\beta$ & $\langle e^{i\theta}\rangle$ & $\langle H\rangle_{\rm QMC}$ & $\langle H\rangle_{\rm ED}$ & $\langle n_d\rangle_{\rm QMC}$ & ED & $\langle H_{K,\perp}\rangle_{\rm QMC}$ & ED & $\langle q\rangle_{\rm QMC}$ & ED\\
\hline""")
for f, phi, b, plab, bl in rcols:
    e = ED[('ring%.4f' % phi, b)]; r = RR[f]
    row = f"{plab} & {bl} & ${pe(*r['phase'])}$"
    for k in ['E', 'docc', 'H_Kflip', 'q']:
        row += f" & ${pe(*r[k])}$ & ${f4(e[k])}$"; maxdev.append((zs(r, e[k], k), f, k))
    out.append(row + r" \\")
out.append(r"""\hline\hline
\end{tabular}
\end{table*}
""")
open('tables_mixed.tex', 'w').write('\n'.join(out)); print('\n'.join(out))
maxdev.sort(reverse=True); print('largest deviations (sigma):', [(round(a, 2), f, k) for a, f, k in maxdev[:6]])
import collections
print('n>2sigma:', sum(a > 2 for a, _, _ in maxdev), 'of', len(maxdev))
for f in list(R) + list(RR) + ['O_b2_loc', 'O_b2_sub']:
    try:
        r = res(f); print(f, r['hdr'], r['meta'], 'phase', r['phase'])
    except Exception as ex: print(f, 'missing', ex)
for f in ['O_b2_loc', 'O_b2_sub', 'O_b2']:
    try:
        r = res(f); e = ED[('open', 2.0)]
        print(f, ' '.join(f"{k}={pe(*r[k])}[{e[k]:.4f}]" for k, _, _ in OBS))
    except Exception as ex: print(f, 'missing')
