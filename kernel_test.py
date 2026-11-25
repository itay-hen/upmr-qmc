# Accuracy test of log_edd (edd.h, via ./upmr_qmc --test) against 2000-digit divided differences
# of exp(-beta x). Energies: random in [-W, W], made distinct by tiny offsets (the reference uses
# the plain recursion). Prints the largest absolute error in ln|f| for each q.
import mpmath as mp, random, subprocess
mp.mp.dps = 2000
def dd(xs, beta):
    f = [mp.e**(-beta * x) for x in xs]; n = len(xs)
    for k in range(1, n):
        for i in range(n - 1, k - 1, -1):
            f[i] = (f[i] - f[i - 1]) / (xs[i] - xs[i - k])
    return f[-1]
random.seed(1)
for q in [10, 40, 80, 160, 240, 320, 400]:
    worst = 0
    for W in [3, 10]:
        for beta in [1, 5, 10]:
            for rep in range(2):
                xs = [mp.mpf(round(random.uniform(-W, W), 6)) + mp.mpf(i) * mp.mpf('1e-7') for i in range(q + 1)]
                ref = mp.log(abs(dd(xs, beta)))
                out = subprocess.run(['./upmr_qmc', '--test', ' '.join([str(beta)] + [mp.nstr(x, 20) for x in xs])],
                                     capture_output=True, text=True)
                val = float(out.stdout) if out.returncode == 0 else float('nan')
                err = abs(val - float(ref)); worst = max(worst, err) if err == err else float('nan')
    print(q, 'max abs err in ln|f|:', worst, flush=True)
