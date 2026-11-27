# Q16 beta=4 multi-seed comparison: errors from seed-to-seed scatter, bias vs ED, efficiency 1/(var*cpu)
import re, glob, numpy as np
R = {'E': -12.8202, 'chiral': -0.1999, 'E2avg': 0.1655}
base_eff = {}
for v in ['base', 'sep', 'rr', 'rrsep']:
    V = {k: [] for k in R}; S = {k: [] for k in R}; cpu = []
    fs = sorted(glob.glob(f'ab/runs/QS_{v}_*.out'))
    for f in fs:
        t = open(f).read(); m = re.search(r"RESULT (.*)", t).group(1).split(); d = {m[i]: (float(m[i+1]), float(m[i+2])) for i in range(0, len(m), 3)}
        cpu.append(float(re.search(r"cpu_ms_per_sweep ([\d.]+)", t).group(1)))
        for k in R: V[k].append(d[k][0]); S[k].append(d[k][1])
    n = len(fs); print(f"{v:6s} n={n} cpu={np.mean(cpu):.3f} ms/sweep")
    for k in R:
        x = np.array(V[k]); sd = x.std(ddof=1); se = sd / np.sqrt(n); eff = 1 / (sd ** 2 * np.mean(cpu))
        if v == 'base': base_eff[k] = eff
        print(f"   {k:7s} mean {x.mean():.5f}±{se:.5f} (z vs ED {(x.mean()-R[k])/se:+.2f})  seed-sd/jackknife {sd/np.mean(S[k]):.1f}  eff/base {eff/base_eff[k]:.2f}")
