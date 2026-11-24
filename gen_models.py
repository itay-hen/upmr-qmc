# Model front end for the generic UPMR sampler (upmr_gen.cpp).
# A Hamiltonian is written as a sum of products of elementary operators
#   fermions:  cd c n      (Jordan-Wigner order = order of the F components)
#   spin S:    Sp Sm Sz    (local value v = 0..2S, S^z = v - S)
#   bosons:    bd b nb     (local value n = 0..nmax)
# The script writes the operator strings to a spec file read by upmr_gen and computes an
# independent ED reference from explicit sparse matrices (JW strings built as Kronecker
# products of Z's), so the two codes share no amplitude logic.
import numpy as np, scipy.sparse as sp, sys, json, argparse, cmath
from collections import defaultdict
from scipy.linalg import eigh

class O:
    def __init__(s, d=None): s.d = defaultdict(complex); [s.d.__setitem__(k, v) for k, v in (d or {}).items()]
    def __add__(s, o):
        r = O(s.d)
        for k, v in o.d.items(): r.d[k] += v
        return r
    def __sub__(s, o): return s + (-1) * o
    def __mul__(s, o):
        if not isinstance(o, O): return O({k: v * o for k, v in s.d.items()})
        r = O()
        for k1, v1 in s.d.items():
            for k2, v2 in o.d.items(): r.d[k1 + k2] += v1 * v2
        return r
    __rmul__ = lambda s, c: s * c
    def dag(s):
        sw = {'cd': 'c', 'c': 'cd', 'Sp': 'Sm', 'Sm': 'Sp', 'bd': 'b', 'b': 'bd'}
        return O({tuple((sw.get(o, o), c) for o, c in reversed(k)): np.conj(v) for k, v in s.d.items()})
def el(op, c): return O({((op, c),): 1.0})
ZERO = O()
def osum(xs):
    r = O()
    for x in xs: r = r + x
    return r

# ------------------------------------------------------------------ ED from matrices
def local_mats(comp):
    k = comp[0]
    if k == 'F':
        c = sp.csr_matrix(np.array([[0., 1.], [0., 0.]]))
        return dict(c=c, cd=c.T.tocsr(), n=sp.diags([0., 1.]), Z=sp.diags([1., -1.]), I=sp.identity(2))
    if k == 'S':
        d = comp[1]; S = (d - 1) / 2; m = np.arange(d) - S
        Sp = sp.diags(np.sqrt(S * (S + 1) - m[:-1] * (m[:-1] + 1)), -1)   # |v> -> |v+1>
        return dict(Sp=Sp.tocsr(), Sm=Sp.T.tocsr(), Sz=sp.diags(m), I=sp.identity(d))
    nmax = comp[1]; b = sp.diags(np.sqrt(np.arange(1, nmax + 1.)), 1)
    return dict(b=b.tocsr(), bd=b.T.tocsr(), nb=sp.diags(np.arange(nmax + 1.)), I=sp.identity(nmax + 1))

def build_matrices(comps, ops_needed):
    L = [local_mats(c) for c in comps]
    cache = {}
    def elmat(op, c):
        if (op, c) in cache: return cache[(op, c)]
        M = sp.identity(1, format='csr')
        for j, comp in enumerate(comps):
            if j == c: f = L[j][op]
            elif comp[0] == 'F' and j < c and op in ('c', 'cd'): f = L[j]['Z']
            else: f = L[j]['I']
            M = sp.kron(M, f, format='csr')
        cache[(op, c)] = M; return M
    def mat(o):
        dim = int(np.prod([2 if c[0] == 'F' else (c[1] if c[0] == 'S' else c[1] + 1) for c in comps]))
        H = sp.csr_matrix((dim, dim), dtype=complex)
        for k, v in o.d.items():
            if abs(v) < 1e-15: continue
            M = sp.identity(dim, format='csr', dtype=complex)
            for op, c in k: M = M @ elmat(op, c)
            H = H + v * M
        return H
    return mat

def basis_values(comps):
    dims = [2 if c[0] == 'F' else (c[1] if c[0] == 'S' else c[1] + 1) for c in comps]
    return np.indices(dims).reshape(len(comps), -1)

def ed(model, betas):
    comps, terms, obs = model['comps'], model['terms'], model['obs']
    mat = build_matrices(comps, None)
    H = osum(o for _, o in terms); Hm = mat(H)
    assert abs(Hm - Hm.getH()).max() < 1e-12, "H not Hermitian"
    # anticommutation check of the JW fermions used here
    fi = [i for i, c in enumerate(comps) if c[0] == 'F']
    for i in fi[:3]:
        for j in fi[:3]:
            A = mat(el('c', i)); B = mat(el('cd', j))
            ac = (A @ B + B @ A).toarray(); assert np.allclose(ac, np.eye(ac.shape[0]) * (i == j))
    vals = basis_values(comps)
    sel = np.ones(vals.shape[1], bool)
    if model.get('N') is not None: sel = vals[fi].sum(0) == model['N']
    I = np.where(sel)[0]
    def sub(M): return M[I][:, I].toarray()
    w, v = eigh(sub(Hm))
    tagmats = defaultdict(lambda: O())
    for tag, o in terms: tagmats[tag] = tagmats[tag] + o
    exs = {('H_' + t): sub(mat(o)) for t, o in tagmats.items()}
    for name, o in obs: exs[name] = sub(mat(o))
    exv = {k: np.real((v.conj() * (M @ v)).sum(0)) for k, M in exs.items()}
    out = {}
    for beta in betas:
        p = np.exp(-beta * (w - w.min())); p /= p.sum()
        r = dict(E=float(p @ w), dim=len(I))
        for k, x in exv.items(): r[k] = float(p @ x)
        Hoff = sum(r['H_' + t] for t in tagmats if t in model['offdiag_tags'])
        r['q'] = -beta * Hoff
        out[beta] = r
    return out

# ------------------------------------------------------------------ spec writer
def write_spec(model, fn):
    with open(fn, 'w') as f:
        for c in model['comps']:
            f.write('comp ' + ' '.join(map(str, c)) + '\n')
        f.write(f"canonical {1 if model.get('N') is not None else 0}\n")
        f.write('root ' + ' '.join(map(str, model['root'])) + '\n')
        for tag, o in model['terms']:
            for k, v in o.d.items():
                if abs(v) < 1e-15: continue
                f.write(f"term {tag} {v.real:.17g} {v.imag:.17g} {len(k)} " + ' '.join(f'{op} {c}' for op, c in k) + '\n')
        for name, o in model['obs']:
            for k, v in o.d.items():
                if abs(v) < 1e-15: continue
                assert abs(v.imag) < 1e-14
                f.write(f"obs {name} {v.real:.17g} {len(k)} " + ' '.join(f'{op} {c}' for op, c in k) + '\n')

# ------------------------------------------------------------------ models
def model_sb(Ns=1, eps=1., om=1., g=0.6, J=0., nmax=8):
    comps = []; S = []; B = []
    for i in range(Ns): S.append(len(comps)); comps.append(('S', 2)); B.append(len(comps)); comps.append(('B', nmax))
    sx = lambda i: el('Sp', S[i]) + el('Sm', S[i]); sz = lambda i: 2 * el('Sz', S[i])
    x = lambda i: el('b', B[i]) + el('bd', B[i])
    terms = [('spin', osum(eps / 2 * sz(i) for i in range(Ns))), ('bos', osum(om * el('nb', B[i]) for i in range(Ns))),
             ('coup', osum(g * sx(i) * x(i) for i in range(Ns)))]
    if Ns > 1: terms.append(('ising', osum(J * sz(i) * sz(i + 1) for i in range(Ns - 1))))
    obs = [('sz', osum(sz(i) * (1 / Ns) for i in range(Ns))), ('n', osum(el('nb', B[i]) * (1 / Ns) for i in range(Ns)))]
    root = [1, 0] * Ns
    return dict(comps=comps, terms=terms, obs=obs, root=root, N=None, offdiag_tags=['coup'])

def model_h3(eps=1., om=1., g=0.6, h=0.5, kap=0.3, nmax=40):
    comps = [('S', 2), ('B', nmax)]
    sx = el('Sp', 0) + el('Sm', 0); x = el('b', 1) + el('bd', 1)
    terms = [('diag', eps * el('Sz', 0) + om * el('nb', 1)), ('h', -h * sx), ('kap', -kap * x), ('coup', -g * sx * x)]
    obs = [('sz', 2 * el('Sz', 0)), ('n', el('nb', 1))]
    return dict(comps=comps, terms=terms, obs=obs, root=[1, 0], N=None, offdiag_tags=['h', 'kap', 'coup'])

def model_mixed(L=4, t=1., U=4., om=1., g=0.3, phi=0., JK=1., Dan=0.5, gam=0.3, nmax=8, N=None, Simp=1, pbc=True, coup='ssh'):
    """Hubbard ring (L sites, flux phi) whose paramagnetic current couples to a cavity mode,
    with a spin-Simp Kondo impurity on site 0 that is also coupled to the cavity.
      H = -t sum_{i,s} (e^{i phi/L} c+_{i+1,s} c_{i,s} + h.c.) + U sum_i n_iu n_id + om a+a
          + g (a + a+) sum_{i,s} (i e^{i phi/L} c+_{i+1,s} c_{i,s} + h.c.)
          + JK S.s_0 + D (S^z)^2 + gam S^x (a + a+)                                       """
    comps = []; f = {}
    for i in range(L):
        for s in (0, 1): f[(i, s)] = len(comps); comps.append(('F',))
    si = len(comps); comps.append(('S', int(round(2 * Simp + 1)))); bi = len(comps); comps.append(('B', nmax))
    cd = lambda i, s: el('cd', f[(i % L, s)]); c = lambda i, s: el('c', f[(i % L, s)]); n = lambda i, s: el('n', f[(i % L, s)])
    ph = cmath.exp(1j * phi / L)
    hopf = osum(ph * cd(i + 1, s) * c(i, s) for i in range(L if pbc else L - 1) for s in (0, 1))       # forward hops
    Hkin = -t * (hopf + hopf.dag())
    if coup == 'cur': cur = 1j * hopf; cur = cur + cur.dag()     # paramagnetic current: i e^{i phi/L} c+c + h.c.
    else: cur = hopf + hopf.dag()                                  # bond modulation (SSH-type): e^{i phi/L} c+c + h.c.
    x = el('b', bi) + el('bd', bi)
    Sz, Sp, Sm = el('Sz', si), el('Sp', si), el('Sm', si)
    sz0 = 0.5 * (n(0, 0) - n(0, 1)); sp0 = cd(0, 0) * c(0, 1); sm0 = cd(0, 1) * c(0, 0)
    terms = [('kin', Hkin), ('U', osum(U * n(i, 0) * n(i, 1) for i in range(L))), ('cav', om * el('nb', bi)),
             ('ssh' if coup == 'ssh' else 'cur', g * x * cur), ('Kzz', JK * Sz * sz0), ('Kflip', 0.5 * JK * (Sp * sm0 + Sm * sp0)),
             ('Dan', Dan * Sz * Sz), ('Sx', 0.5 * gam * (Sp + Sm) * x)]
    obs = [('docc', osum(n(i, 0) * n(i, 1) for i in range(L))), ('nph', el('nb', bi)),
           ('Sz2', Sz * Sz), ('SzSz0', Sz * sz0), ('Nf', osum(n(i, s) for i in range(L) for s in (0, 1)))]
    if N is None: N = L
    root = []
    for i in range(L):                                    # Neel-like start, N fermions
        for s in (0, 1): root.append(0)
    k = 0
    for i in range(L):
        if k < N: root[f[(i, i % 2)]] = 1; k += 1
    for i in range(L):
        if k < N: root[f[(i, 1 - i % 2)]] = 1; k += 1
    root += [int(round(Simp)), 0]
    return dict(comps=comps, terms=terms, obs=obs, root=root, N=N, offdiag_tags=['kin', 'cur', 'ssh', 'Kflip', 'Sx'])

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('model'); ap.add_argument('--spec', default=None); ap.add_argument('--betas', default='1')
    ap.add_argument('--kw', default='{}', help='json dict of model parameters')
    a = ap.parse_args()
    kw = json.loads(a.kw)
    m = dict(sb=model_sb, h3=model_h3, mixed=model_mixed)[a.model](**kw)
    if a.spec: write_spec(m, a.spec)
    betas = [float(b) for b in a.betas.split(',')]
    if betas[0] > 0:
        r = ed(m, betas)
        for b in betas: print(json.dumps(dict(model=a.model, kw=kw, beta=b, ed=r[b])), flush=True)
