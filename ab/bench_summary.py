# efficiency 1/(sigma^2 t_cpu) of update sets relative to the baseline (local + subcycle), ab/runs/
import re, glob, os
def res(f):
    if not os.path.exists(f): return None
    t = open(f).read(); m = re.search(r"RESULT (.*)", t)
    if not m: return None
    m = m.group(1).split(); d = {m[i]: (float(m[i+1]), float(m[i+2])) for i in range(0, len(m), 3)}
    for k, p in [('cpu', r"cpu_ms_per_sweep ([\d.]+)"), ('tau', r"tau_int_E ([\d.]+)"), ('cf', r"closed_frac ([\d.]+)")]: d[k] = float(re.search(p, t).group(1))
    d['acc'] = {a: round(int(b) / max(1, int(c)), 3) for a, b, c in re.findall(r"acc (\w+)\s+(\d+)/(\d+)", t) if a.startswith(('sep', 'rot', 'worm'))}
    return d
for sysn, ob in [('Q32', 'chiral'), ('LR12', 'Scdw')]:
    b = res(f'ab/runs/{sysn}_base.out'); eff = lambda d, k: 1 / (d[k][1] ** 2 * d['cpu'])
    print(f"{sysn} base: tau_E {b['tau']} cpu {b['cpu']} ms  E {b['E']}")
    names = [os.path.basename(f)[len(sysn)+1:-4] for f in sorted(glob.glob(f'ab/runs/{sysn}_*.out'))]
    for v in [n for n in names if n != 'base']:
        r = res(f'ab/runs/{sysn}_{v}.out')
        if r is None: print('  ', v, '(no result)'); continue
        print(f"  {v:18s} H {eff(r,'E')/eff(b,'E'):.2f}  root {eff(r,ob+'_root')/eff(b,ob+'_root'):.2f}  w {eff(r,ob)/eff(b,ob):.2f}  q {eff(r,'q')/eff(b,'q'):.2f} | tau_E {r['tau']:6.1f} cpu {r['cpu']:.2f} closed {r['cf']:.2f} E {r['E'][0]:.3f}({r['E'][1]:.3f}) {r['acc']}")
