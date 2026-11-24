"""Seed-41 deviation and first-half vs second-half drift for the four beta=5 worm runs.
Usage: python3 revision_checks/I_worm_drift.py  (reads the --rawbins 1 outputs in revision_checks/)"""
import sys, json, numpy as np
from pathlib import Path
from scipy.linalg import eigh
D = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent
FN = {41: "I_new_seed41.txt", 42: "I_replay_W2_sb_42.txt", 43: "I_replay_W2_sb_43.txt", 44: "I_replay_W2_sb_44.txt"}
b = np.diag(np.sqrt(np.arange(1, 9)), 1); z = np.diag([1., -1.]); x = np.array([[0., 1], [1, 0]])
off = .6*np.kron(x, b+b.T); H = np.kron(z/2, np.eye(9)) + np.kron(np.eye(2), np.diag(np.arange(9.))) + off
w, V = eigh(H); p = np.exp(-5*(w-w.min())); p /= p.sum(); rho = (V*p) @ V.T
ED = {'E': w@p, 'szw': np.trace(np.kron(z, np.eye(9))@rho), 'n': np.trace(np.kron(np.eye(2), np.diag(np.arange(9.)))@rho), 'q': -5*np.trace(off@rho)}
COL = {'E': 1, 'szw': 3, 'n': 4, 'q': 5}
def load(f):
    a = np.array([[float(t) for t in l.split()[2:]] for l in open(f) if l.startswith('RAW_BIN')])
    return a[:, 1:]          # col 0 = sum of weights (= closed count), then E, szroot, szw, n, q, ...
def est(raw):
    tot = raw.sum(0); B = len(raw); jk = (tot-raw)/(tot[0]-raw[:, 0])[:, None]
    return tot/tot[0], np.sqrt((B-1)/B*((jk-jk.mean(0))**2).sum(0))
out = {'ED': {k: float(v) for k, v in ED.items()}, 'runs': {}}
diffs = {k: [] for k in COL}
for s in (41, 42, 43, 44):
    raw = load(D/f'{FN[s]}'); m, e = est(raw); h = len(raw)//2
    m1, e1 = est(raw[:h]); m2, e2 = est(raw[h:])
    r = {'bins': len(raw)}
    for k, j in COL.items():
        dz = (m2[j]-m1[j])/np.hypot(e1[j], e2[j]); diffs[k].append((m2[j]-m1[j], np.hypot(e1[j], e2[j])))
        r[k] = dict(mean=float(m[j]), sem=float(e[j]), z_vs_ED=float((m[j]-ED[k])/e[j]), half1=float(m1[j]), half2=float(m2[j]), z_second_minus_first=float(dz))
    out['runs'][s] = r
pooled = {}
for k, L in diffs.items():
    d = np.array([a for a, _ in L]); s = np.array([c for _, c in L])
    wts = 1/s**2
    pooled[k] = dict(stouffer=float((d/s).sum()/2), invvar=float((wts@d/wts.sum())/np.sqrt(1/wts.sum())))
out['pooled_second_minus_first'] = pooled
print(json.dumps(out, indent=1))
