# Efficiency of update sets: eff = 1/(sigma^2 * CPU seconds) for several observables.
import re, sys, glob
def parse(f):
    t = open(f).read()
    m = re.search(r"RESULT (.*)", t).group(1).split(); d = {m[i]: (float(m[i+1]), float(m[i+2])) for i in range(0, len(m), 3)}
    cpu = float(re.search(r"cpu_ms_per_sweep ([\d.]+)", t).group(1)); sw = 2e5
    tau = float(re.search(r"tau_int_E ([\d.]+)", t).group(1)); qmax = int(re.search(r"q_max (\d+)", t).group(1))
    return d, cpu * 1e-3 * sw, tau, qmax
for sysn in sys.argv[1:]:
    print(f"== {sysn}")
    base = None
    for v in ['base', 'reroot', 'rev', 'exch', 'all', 'revT', 'allT']:
        f = f'nl/runs/{sysn}_{v}.out'
        try: d, cpu, tau, qmax = parse(f)
        except Exception: print(v, 'missing'); continue
        obs = [k for k in d if k in ('E', 'q') or k.endswith('_root') or k in ('docc', 'nph', 'Scdw', 'chiral', 'E2avg')]
        obs = [k for k in obs if d[k][1] > 0]; eff = {k: 1.0 / (d[k][1] ** 2 * cpu) for k in obs}
        if base is None: base = eff
        print(f"{v:7s} cpu={cpu:7.0f}s tauE={tau:6.1f} qmax={qmax:4d} " + ' '.join(f"{k}:{eff[k]/base[k]:.2f}" for k in obs))
