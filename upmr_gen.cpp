// upmr_gen.cpp -- model-agnostic UPMR closed-walk QMC (QMCofEverything_v17).
//
// Input: a spec file listing the components (fermion modes F, spin-S qudits "S d", bosons
// "B nmax") and the Hamiltonian as a sum of operator strings (cd c n | Sp Sm Sz | bd b nb),
// see gen_models.py. Everything else is derived automatically:
//   * each string's displacement vector (mod 2 for F, mod d for S, integer for B);
//   * D_0 = all zero-displacement strings; labels t = classes of strings with equal
//     displacement, d_t(z) = sum of their amplitudes on the departure state z (complex,
//     with Jordan--Wigner signs evaluated on the fly);
//   * the inverse label of each t (Hermitian pairing, checked numerically);
//   * subcycle templates = all minimal relations of length 3..tmax among the labels.
// Moves (as in upmr_qmc.cpp): basis (single component, grand-canonical only), basis shift
// by a label displacement, swap, pair, subcycle replacement, re-rooting, inner-loop
// reversal, segment exchange, worm. Sampling |w|, measuring with the complex phase of w.
//
// Build: clang++ -O3 -std=c++17 -o upmr_gen upmr_gen.cpp

#include <algorithm>
#include <chrono>
#include <ctime>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using std::vector; using std::string;
using cplx = std::complex<double>;

// ---------------------------------------------------------------- kernel (as upmr_qmc.cpp)
#include "edd.h"
#include <unordered_map>
static const int MAXC = 512;   // max number of components

// ---------------------------------------------------------------- generic model
enum { CD, C, NF, SP, SM, SZ, BD, B, NB };
struct Elem { int op, c; };
struct Str { cplx coef; vector<Elem> ops; int tag; };
struct Comp { char kind; int d; int nmax; double S; };
struct Label { vector<std::pair<int,int>> disp; int inv = -1; vector<int> strs; int tag = -1; };

struct Model {
    vector<Comp> comps; int nc = 0; bool canonical = false; vector<int> root0;
    vector<string> tags; vector<Str> diag; vector<Label> labels;
    vector<string> obs_names; vector<vector<Str>> obs;
    vector<string> cons_names; vector<vector<Str>> cons;   // conserved diagonal quantities (root moves keep them)
    mutable std::unordered_map<string, double> Ecache;     // diagonal energy cache
    vector<int> offtag;   // per tag: 1 if it owns labels

    int tag_id(const string& s) {
        for (size_t i = 0; i < tags.size(); i++) if (tags[i] == s) return i;
        tags.push_back(s); return tags.size() - 1;
    }
    static int opcode(const string& s) {
        static const char* nm[] = {"cd","c","n","Sp","Sm","Sz","bd","b","nb"};
        for (int i = 0; i < 9; i++) if (s == nm[i]) return i;
        fprintf(stderr, "unknown op %s\n", s.c_str()); exit(1);
    }
    static int shift(int op) { return (op == CD || op == SP || op == BD) ? 1 : (op == C || op == SM || op == B) ? -1 : 0; }
    int norm(int c, int x) const {
        const Comp& K = comps[c];
        if (K.kind == 'F') return ((x % 2) + 2) % 2;
        if (K.kind == 'S') return ((x % K.d) + K.d) % K.d;
        return x;
    }
    vector<std::pair<int,int>> disp_of(const Str& s) const {
        std::map<int,int> m; for (auto& e : s.ops) m[e.c] += shift(e.op);
        vector<std::pair<int,int>> r;
        for (auto& kv : m) { int x = norm(kv.first, kv.second); if (x) r.push_back({kv.first, x}); }
        return r;
    }
    vector<std::pair<int,int>> neg(const vector<std::pair<int,int>>& d) const {
        vector<std::pair<int,int>> r; for (auto& p : d) r.push_back({p.first, norm(p.first, -p.second)}); return r;
    }
    void load(const string& fn) {
        std::ifstream in(fn); if (!in) { fprintf(stderr, "cannot open %s\n", fn.c_str()); exit(1); }
        string line; std::map<vector<std::pair<int,int>>, int> idx;
        vector<Str> off;
        while (std::getline(in, line)) {
            std::istringstream is(line); string kw; is >> kw;
            if (kw == "comp") {
                string k; is >> k; Comp K{k[0], 2, 0, 0.5};
                if (k == "S") { is >> K.d; K.S = (K.d - 1) / 2.0; }
                if (k == "B") { is >> K.nmax; K.d = 0; }
                comps.push_back(K);
            } else if (kw == "canonical") { int x; is >> x; canonical = x; }
            else if (kw == "root") { int x; while (is >> x) root0.push_back(x); }
            else if (kw == "term" || kw == "obs" || kw == "conserve") {
                string tg; double re, im = 0; int k; is >> tg >> re; if (kw == "term") is >> im; is >> k;
                Str s; s.coef = cplx(re, im);
                for (int i = 0; i < k; i++) { string o; int c; is >> o >> c; s.ops.push_back({opcode(o), c}); }
                if (kw == "conserve") {
                    size_t j = 0; for (; j < cons_names.size(); j++) if (cons_names[j] == tg) break;
                    if (j == cons_names.size()) { cons_names.push_back(tg); cons.emplace_back(); }
                    cons[j].push_back(s);
                } else if (kw == "obs") {
                    size_t j = 0; for (; j < obs_names.size(); j++) if (obs_names[j] == tg) break;
                    if (j == obs_names.size()) { obs_names.push_back(tg); obs.emplace_back(); }
                    obs[j].push_back(s);
                } else { s.tag = tag_id(tg); off.push_back(s); }
            }
        }
        nc = comps.size();
        if (nc > MAXC) { fprintf(stderr, "too many components (%d > %d)\n", nc, MAXC); exit(1); }
        for (auto& s : off) {
            auto d = disp_of(s);
            if (d.empty()) { diag.push_back(s); continue; }
            auto it = idx.find(d); int t;
            if (it == idx.end()) { t = labels.size(); idx[d] = t; labels.push_back(Label{d}); labels[t].tag = s.tag; }
            else t = it->second;
            labels[t].strs.push_back(labels.size() ? 0 : 0);
            labels[t].strs.back() = (int)strs.size(); strs.push_back(s);
            if (labels[t].tag != s.tag) labels[t].tag = tag_id("mixed");
        }
        for (auto& L : labels) {
            auto it = idx.find(neg(L.disp));
            if (it == idx.end()) { fprintf(stderr, "label without inverse: H not Hermitian?\n"); exit(1); }
            L.inv = it->second;
        }
        offtag.assign(tags.size(), 0); for (auto& L : labels) offtag[L.tag] = 1;
    }
    vector<Str> strs;
    bool fpar(const int* z, int c) const { int p = 0; for (int j = 0; j < c; j++) if (comps[j].kind == 'F') p ^= z[j]; return p; }
    // amplitude <z'|s|z> of one string; tmp receives z' (operators act right to left)
    cplx act(const Str& s, const int* z, int* tmp) const {
        std::copy(z, z + nc, tmp); cplx a = s.coef;
        for (int k = (int)s.ops.size() - 1; k >= 0; k--) {
            const Elem& e = s.ops[k]; int& v = tmp[e.c]; const Comp& K = comps[e.c];
            switch (e.op) {
                case CD: if (v == 1) return 0; if (fpar(tmp, e.c)) a = -a; v = 1; break;
                case C:  if (v == 0) return 0; if (fpar(tmp, e.c)) a = -a; v = 0; break;
                case NF: if (v == 0) return 0; break;
                case SP: { if (v == K.d - 1) return 0; double m = v - K.S; a *= std::sqrt(K.S*(K.S+1) - m*(m+1)); v++; break; }
                case SM: { if (v == 0) return 0; double m = v - K.S; a *= std::sqrt(K.S*(K.S+1) - m*(m-1)); v--; break; }
                case SZ: a *= (v - K.S); break;
                case BD: if (v >= K.nmax) return 0; a *= std::sqrt(v + 1.0); v++; break;
                case B:  if (v <= 0) return 0; a *= std::sqrt((double)v); v--; break;
                case NB: a *= (double)v; break;
            }
            if (a == 0.0) return 0;
        }
        return a;
    }
    double E(const int* z) const {
        string key(nc, '\0'); for (int c = 0; c < nc; c++) key[c] = (char)z[c];
        auto it = Ecache.find(key); if (it != Ecache.end()) return it->second;
        int tmp[MAXC]; double e = 0;
        for (auto& s : diag) e += act(s, z, tmp).real();
        if (Ecache.size() > 500000) Ecache.clear();   // bounded memory (~60 MB)
        Ecache.emplace(std::move(key), e);
        return e;
    }
    double consval(int o, const int* z) const {
        int tmp[MAXC]; double e = 0;
        for (auto& s : cons[o]) e += act(s, z, tmp).real();
        return e;
    }
    // root moves must keep every conserved quantity (and, if canonical, the fermion number)
    bool same_sector(const int* z) const {
        if (canonical && nferm(z) != nferm(&root0[0])) return false;
        for (size_t o = 0; o < cons.size(); o++) if (std::fabs(consval(o, z) - consval(o, &root0[0])) > 1e-9) return false;
        return true;
    }
    double obsval(int o, const int* z) const {
        int tmp[MAXC]; double e = 0;
        for (auto& s : obs[o]) e += act(s, z, tmp).real();
        return e;
    }
    cplx amp(int t, const int* z) const {
        int tmp[MAXC]; cplx a = 0;
        for (int i : labels[t].strs) a += act(strs[i], z, tmp);
        return a;
    }
    void apply(int t, int* z) const { for (auto& d : labels[t].disp) z[d.first] = norm(d.first, z[d.first] + d.second); }
    int nferm(const int* z) const { int n = 0; for (int c = 0; c < nc; c++) if (comps[c].kind == 'F') n += z[c]; return n; }
    bool in_domain(const int* z) const {
        for (int c = 0; c < nc; c++) if (comps[c].kind == 'B' && (z[c] < 0 || z[c] > comps[c].nmax)) return false;
        return true;
    }
    // numerical check of the Hermitian pairing d_tbar(pi_t z) = conj d_t(z) on random states
    long selftest(std::mt19937_64& rng, int ntrial) const {
        long bad = 0; vector<int> z(nc), z2(nc);
        for (int n = 0; n < ntrial; n++) {
            for (int c = 0; c < nc; c++) {
                int r = comps[c].kind == 'F' ? 2 : comps[c].kind == 'S' ? comps[c].d : comps[c].nmax + 1;
                z[c] = std::uniform_int_distribution<int>(0, r - 1)(rng);
            }
            for (size_t t = 0; t < labels.size(); t++) {
                cplx a = amp(t, &z[0]); z2 = z; apply(t, &z2[0]);
                if (!in_domain(&z2[0])) { if (std::abs(a) > 1e-14) bad++; continue; }
                cplx b = amp(labels[t].inv, &z2[0]);
                if (std::abs(b - std::conj(a)) > 1e-12) bad++;
            }
        }
        return bad;
    }
};

// ---------------------------------------------------------------- sampler
struct Eval { bool ok; double logw; cplx ph; vector<double> En; vector<int> states; };

struct Sampler {
    Model& m; double beta; std::mt19937_64 rng; int nc;
    vector<int> root, word; double logw = 0; cplx ph = 1; vector<double> En; vector<int> st;
    bool use_sub = false, use_nonlocal = false, use_worm = false; double mu = 1.0;
    bool use_reroot = false, use_rev = false, use_exch = false;   // individual nonlocal moves
    vector<vector<int>> templ; vector<int> templ_inv;
    std::map<string, std::pair<long,long>> stats; long weight_check_fail = 0;

    Sampler(Model& mm, double b, unsigned long seed) : m(mm), beta(b), rng(seed) {
        nc = m.nc; root = m.root0;
        auto e = eval(root, word); logw = e.logw; ph = e.ph; En = e.En; st = e.states;
    }
    double U() { return std::uniform_real_distribution<double>(0, 1)(rng); }
    int Ri(int n) { return std::uniform_int_distribution<int>(0, n - 1)(rng); }
    bool closed_of(const vector<int>& s, int q) const { for (int c = 0; c < nc; c++) if (s[(size_t)q*nc+c] != s[c]) return false; return true; }
    Eval eval(const vector<int>& r, const vector<int>& w) const {
        Eval e; e.ok = true; e.ph = 1; double la = 0; int q = w.size();
        e.states.resize((size_t)(q + 1) * nc); std::copy(r.begin(), r.end(), e.states.begin());
        e.En.resize(q + 1); e.En[0] = m.E(&e.states[0]);
        for (int k = 0; k < q; k++) {
            int* cur = &e.states[(size_t)k*nc]; int* nxt = &e.states[(size_t)(k+1)*nc];
            cplx d = m.amp(w[k], cur); double ad = std::abs(d);
            if (ad < 1e-14) { e.ok = false; e.logw = -INFINITY; return e; }
            e.ph *= d / ad; la += std::log(ad);
            std::copy(cur, cur + nc, nxt); m.apply(w[k], nxt);
            e.En[k+1] = m.E(nxt);
        }
        if (q % 2) e.ph = -e.ph;
        e.logw = la + log_edd(e.En, beta);
        if (!closed_of(e.states, q)) e.logw += std::log(mu);
        return e;
    }
    bool is_closed() const { return closed_of(st, (int)word.size()); }
    bool metro(const vector<int>& r2, const vector<int>& w2, double logR, const char* nm, bool expect_equal = false) {
        auto& s = stats[nm]; s.first++;
        Eval e = eval(r2, w2); if (!e.ok) return false;
        if (expect_equal && std::fabs(e.logw - logw) > 1e-8) weight_check_fail++;
        double a = e.logw - logw + logR;
        if (a >= 0 || U() < std::exp(a)) { root = r2; word = w2; logw = e.logw; ph = e.ph; En = e.En; st = e.states; s.second++; return true; }
        return false;
    }
    // ---- moves
    void basis_comp() {
        vector<int> r2 = root; int c = Ri(nc); int d = U() < 0.5 ? 1 : -1;
        r2[c] = m.norm(c, r2[c] + d);
        if (!m.in_domain(&r2[0]) || !m.same_sector(&r2[0])) { stats["basis"].first++; return; }
        metro(r2, word, 0, "basis");
    }
    void basis_label() {                 // shift every state of the walk by a label displacement
        vector<int> r2 = root; m.apply(Ri(m.labels.size()), &r2[0]);
        // a Z_2 displacement does not fix the direction of a hop: stay in the fermion-number sector
        if (!m.in_domain(&r2[0]) || !m.same_sector(&r2[0])) { stats["basis_lab"].first++; return; }
        metro(r2, word, 0, "basis_lab");
    }
    void swap_move() {
        int q = word.size(); if (q < 2) return;
        int j = Ri(q - 1); if (word[j] == word[j+1]) return;
        vector<int> w2 = word; std::swap(w2[j], w2[j+1]); metro(root, w2, 0, "swap");
    }
    void pair_move() {
        const double pdel = 0.5; int M = m.labels.size(); int q = word.size();
        if (U() < 1 - pdel) {
            int t = Ri(M); int g = Ri(q + 1); vector<int> w2 = word; w2.insert(w2.begin() + g, {t, m.labels[t].inv});
            metro(root, w2, std::log(pdel / (1 - pdel) * M), "pair_ins");
        } else {
            if (q < 2) { stats["pair_del"].first++; return; }
            int j = Ri(q - 1); if (word[j+1] != m.labels[word[j]].inv) { stats["pair_del"].first++; return; }
            vector<int> w2 = word; w2.erase(w2.begin() + j, w2.begin() + j + 2);
            metro(root, w2, std::log((1 - pdel) / pdel / M), "pair_del");
        }
    }
    static double log_norderings(const vector<int>& cnt) {
        int s = 0; double l = 0; for (int c : cnt) { s += c; l -= std::lgamma(c + 1.0); } return l + std::lgamma(s + 1.0);
    }
    int count_windows(const vector<int>& w, int r, const vector<int>& tc, vector<int>* pos) const {
        int q = w.size(), N = 0;
        if (r == 0) { if (pos) for (int g = 0; g <= q; g++) pos->push_back(g); return q + 1; }
        vector<int> c(tc.size());
        for (int j = 0; j + r <= q; j++) {
            std::fill(c.begin(), c.end(), 0); bool ok = true;
            for (int k = 0; k < r && ok; k++) if (++c[w[j+k]] > tc[w[j+k]]) ok = false;
            if (ok) { N++; if (pos) pos->push_back(j); }
        }
        return N;
    }
    void subcycle_move() {
        if (templ.empty()) return;
        auto& s = stats["subcycle"];
        int a = Ri(templ.size()); const vector<int>& tc = templ[a]; const vector<int>& tci = templ[templ_inv[a]];
        int ell = 0; for (int c : tc) ell += c;
        int r = Ri(ell + 1);
        vector<int> pos; int N = count_windows(word, r, tc, &pos);
        if (N == 0) { s.first++; return; }
        int j = pos[Ri(pos.size())];
        vector<int> wc(tc.size(), 0); for (int k = 0; k < r; k++) wc[word[j+k]]++;
        vector<int> ins, Tbar(tc.size(), 0);
        for (size_t t = 0; t < tc.size(); t++) for (int k = 0; k < tc[t] - wc[t]; k++) { ins.push_back(m.labels[t].inv); Tbar[m.labels[t].inv]++; }
        std::shuffle(ins.begin(), ins.end(), rng);
        vector<int> w2(word.begin(), word.begin() + j);
        w2.insert(w2.end(), ins.begin(), ins.end()); w2.insert(w2.end(), word.begin() + j + r, word.end());
        int N2 = count_windows(w2, ell - r, tci, nullptr);
        double logR = std::log((double)N) + log_norderings(Tbar) - std::log((double)N2) - log_norderings(wc);
        metro(root, w2, logR, "subcycle");
    }
    void reroot_move() {
        int q = word.size(); if (q < 2) return;
        int k = 1 + Ri(q - 1);
        vector<int> r2(st.begin() + (size_t)k*nc, st.begin() + (size_t)(k+1)*nc);
        vector<int> w2(word.begin() + k, word.end()); w2.insert(w2.end(), word.begin(), word.begin() + k);
        metro(r2, w2, 0, "reroot");
    }
    bool same_state(int a, int b) const { for (int c = 0; c < nc; c++) if (st[(size_t)a*nc+c] != st[(size_t)b*nc+c]) return false; return true; }
    // Inner-loop reversal with a targeted proposal: choose a position a uniformly in 0..q, then b
    // uniformly among the other positions whose state equals z_a. The pair {a,b} is proposed with
    // probability 2/((q+1)(n-1)), n = multiplicity of z_a; reversal preserves the multiset of visited
    // states and the states at a and b, so the reverse proposal has the same probability.
    bool targeted_reversal = true;
    void reversal_move() {
        int q = word.size(); if (q < 2) return;
        int a = Ri(q + 1), b;
        if (targeted_reversal) {
            vector<int> P;
            for (int k = 0; k <= q; k++) if (k != a && same_state(a, k)) P.push_back(k);
            if (P.empty()) { stats["reversal"].first++; return; }
            b = P[Ri(P.size())];
        } else { b = Ri(q + 1); if (a == b) return; }
        if (a > b) std::swap(a, b);
        if (!same_state(a, b)) { stats["reversal"].first++; return; }
        if (a == 0 && b == q && q >= 2) { /* whole word: allowed */ }
        vector<int> w2 = word; for (int k = 0; k < b - a; k++) w2[a + k] = m.labels[word[b - 1 - k]].inv;
        metro(root, w2, 0, "reversal", true);
    }
    void exchange_move() {
        int q = word.size(); if (q < 2) return;
        vector<int> slots(q + 2); for (int i = 0; i < q + 2; i++) slots[i] = i;
        for (int i = 0; i < 4; i++) std::swap(slots[i], slots[i + Ri(q + 2 - i)]);
        vector<int> bars(slots.begin(), slots.begin() + 4); std::sort(bars.begin(), bars.end());
        int pA = bars[0], pS1 = bars[1]-bars[0], pB = bars[2]-bars[1]-1, pS2 = bars[3]-bars[2];
        int s1a = pA, s1b = pA + pS1, s2a = s1b + pB, s2b = s2a + pS2;
        bool opposite = U() < 0.5; auto& s = stats[opposite ? "exch_opp" : "exch_eq"];
        vector<int> A(word.begin(), word.begin()+s1a), S1(word.begin()+s1a, word.begin()+s1b),
                    Bm(word.begin()+s1b, word.begin()+s2a), S2(word.begin()+s2a, word.begin()+s2b), Cc(word.begin()+s2b, word.end());
        vector<int> w2;
        if (!opposite) {
            if (!(same_state(s1a, s2a) && same_state(s1b, s2b)) || S1 == S2) { s.first++; return; }
            w2 = A; w2.insert(w2.end(), S2.begin(), S2.end()); w2.insert(w2.end(), Bm.begin(), Bm.end());
            w2.insert(w2.end(), S1.begin(), S1.end()); w2.insert(w2.end(), Cc.begin(), Cc.end());
        } else {
            if (!(same_state(s1a, s2b) && same_state(s1b, s2a))) { s.first++; return; }
            auto inv = [&](const vector<int>& S) { vector<int> R; for (int k = S.size() - 1; k >= 0; k--) R.push_back(m.labels[S[k]].inv); return R; };
            vector<int> S2i = inv(S2), S1i = inv(S1);
            w2 = A; w2.insert(w2.end(), S2i.begin(), S2i.end()); w2.insert(w2.end(), Bm.begin(), Bm.end());
            w2.insert(w2.end(), S1i.begin(), S1i.end()); w2.insert(w2.end(), Cc.begin(), Cc.end());
        }
        metro(root, w2, 0, opposite ? "exch_opp" : "exch_eq", true);
    }
    void worm_move() {
        int M = m.labels.size();
        if (U() < 0.5) { vector<int> w2 = word; w2.push_back(Ri(M)); metro(root, w2, std::log((double)M), "worm_app"); }
        else {
            if (word.empty()) { stats["worm_ret"].first++; return; }
            vector<int> w2 = word; w2.pop_back(); metro(root, w2, -std::log((double)M), "worm_ret");
        }
    }
    void step() {
        vector<int> menu = {8, 1, 2};
        if (!m.canonical && m.cons.empty()) menu.push_back(0);   // single-component moves break most constraints
        if (use_sub) menu.push_back(3);
        if (use_reroot) menu.push_back(4);
        if (use_rev) menu.push_back(5);
        if (use_exch) menu.push_back(6);
        if (use_worm) { menu.push_back(7); menu.push_back(7); }
        int mv = menu[Ri(menu.size())];
        if (!is_closed() && mv != 7) return;
        switch (mv) {
            case 0: basis_comp(); break; case 8: basis_label(); break; case 1: swap_move(); break; case 2: pair_move(); break;
            case 3: subcycle_move(); break; case 4: reroot_move(); break; case 5: reversal_move(); break;
            case 6: exchange_move(); break; case 7: worm_move(); break;
        }
    }
    // all minimal relations of length 3..tmax (multisets with zero total displacement and no
    // pair {t, tbar}); for tmax <= 4 this is exactly minimality.
    void make_templates(int tmax) {
        int M = m.labels.size(); vector<int> cur;
        std::function<void(int, vector<int>&)> rec = [&](int start, vector<int>& z) {
            int len = cur.size();
            if (len >= 3) {
                bool zero = true; for (int c = 0; c < nc; c++) if (z[c]) { zero = false; break; }
                if (zero) { vector<int> cnt(M, 0); for (int t : cur) cnt[t]++; templ.push_back(cnt); }
            }
            if (len == tmax) return;
            for (int t = start; t < M; t++) {
                bool pair = false; for (int u : cur) if (u == m.labels[t].inv) { pair = true; break; }
                if (pair) continue;
                vector<int> z2 = z; for (auto& d : m.labels[t].disp) z2[d.first] = m.comps[d.first].kind == 'B' ? z2[d.first] + d.second : m.norm(d.first, z2[d.first] + d.second);
                cur.push_back(t); rec(t, z2); cur.pop_back();
            }
        };
        vector<int> z(nc, 0); rec(0, z);
        for (auto& T : templ) {
            vector<int> Ti(M, 0); for (int t = 0; t < M; t++) Ti[m.labels[t].inv] += T[t];
            int j = std::find(templ.begin(), templ.end(), Ti) - templ.begin();
            if (j == (int)templ.size()) { fprintf(stderr, "template set not closed under inversion\n"); exit(1); }
            templ_inv.push_back(j);
        }
    }
};

// ---------------------------------------------------------------- driver
static double argd(int argc, char** argv, const char* k, double def) { for (int i = 1; i + 1 < argc; i++) if (string(argv[i]) == k) return atof(argv[i+1]); return def; }
static string args(int argc, char** argv, const char* k, const char* def) { for (int i = 1; i + 1 < argc; i++) if (string(argv[i]) == k) return argv[i+1]; return def; }

int main(int argc, char** argv) {
    Model m; m.load(args(argc, argv, "--spec", "model.spec"));
    double beta = argd(argc, argv, "--beta", 1);
    long therm = (long)argd(argc, argv, "--therm", 1e4), prod = (long)argd(argc, argv, "--sweeps", 1e5);
    int nbins = (int)argd(argc, argv, "--bins", 100);
    unsigned long seed = (unsigned long)argd(argc, argv, "--seed", 1);
    Sampler S(m, beta, seed);
    S.use_sub = argd(argc, argv, "--sub", 0) > 0; S.use_nonlocal = argd(argc, argv, "--nonlocal", 0) > 0;
    // --nonlocal 1 enables all three; --reroot/--reversal/--exchange override individually
    S.use_reroot = argd(argc, argv, "--reroot", S.use_nonlocal) > 0;
    S.use_rev = argd(argc, argv, "--reversal", S.use_nonlocal) > 0;
    S.use_exch = argd(argc, argv, "--exchange", S.use_nonlocal) > 0;
    S.targeted_reversal = argd(argc, argv, "--targeted", 1) > 0;
    S.use_worm = argd(argc, argv, "--worm", 0) > 0; S.mu = argd(argc, argv, "--mu", 1.0);
    if (S.use_sub) S.make_templates((int)argd(argc, argv, "--tmax", 4));
    long bad = m.selftest(S.rng, 2000);
    printf("spec=%s comps=%d labels=%zu diag_strings=%zu templates=%zu hermitian_pairing_failures=%ld\n",
           args(argc, argv, "--spec", "").c_str(), m.nc, m.labels.size(), m.diag.size(), S.templ.size(), bad);
    if (bad) return 1;
    // fixed number of steps per sweep, set during thermalization
    double qsum = 0; long qn = 0;
    for (long s = 0; s < therm; s++) {
        for (long k = (long)S.word.size() + S.nc + 1; k > 0; k--) S.step();
        if (s >= therm / 2) { qsum += S.word.size(); qn++; }
    }
    long nsteps = (long)std::lround(qsum / std::max(1L, qn)) + S.nc + 1;
    if (argd(argc, argv, "--steps", 0) > 0) nsteps = (long)argd(argc, argv, "--steps", 0);
    // observables: 0 Re ph, 1 Im ph, 2 E, 3 q, then per tag n_tag, then obs root, obs weighted
    int T = m.tags.size(), NO = m.obs.size(), NV = 4 + T + 2 * NO;
    long qmax = 0; auto t0 = std::chrono::steady_clock::now(); std::clock_t c0 = std::clock();
    long per = prod / nbins; vector<vector<double>> bins; vector<long> closed_counts; long nmeas = 0, nsweep = 0;
    vector<double> Eser;
    for (int b = 0; b < nbins; b++) {
        vector<double> sum(NV, 0); long cnt = 0;
        for (long s = 0; s < per; s++) {
            for (long k = nsteps; k > 0; k--) S.step();
            nsweep++; qmax = std::max(qmax, (long)S.word.size());
            if (!S.is_closed()) continue;
            int q = S.word.size(); double p = S.ph.real();
            double lf0 = log_edd(S.En, beta), Eest = S.En[0];
            if (q > 0) { vector<double> E1(S.En.begin() + 1, S.En.end()); Eest = S.En[0] - std::exp(log_edd(E1, beta) - lf0); }
            vector<double> v(NV, 0); v[0] = p; v[1] = S.ph.imag(); v[2] = p * Eest; v[3] = p * q;
            for (int t : S.word) v[4 + m.labels[t].tag] += p;
            vector<double> om(q + 1);
            for (int k = 0; k <= q; k++) { vector<double> E2 = S.En; E2.push_back(S.En[k]); om[k] = std::exp(log_edd(E2, beta) - lf0) / beta; }
            for (int o = 0; o < NO; o++) {
                v[4 + T + o] = p * m.obsval(o, &S.st[0]);
                double a = 0; for (int k = 0; k <= q; k++) a += om[k] * m.obsval(o, &S.st[(size_t)k*S.nc]);
                v[4 + T + NO + o] = p * a;
            }
            for (int i = 0; i < NV; i++) sum[i] += v[i];
            cnt++; Eser.push_back(Eest);
        }
        // Preserve raw sums; all closed visits have equal weight.
        closed_counts.push_back(cnt);
        if (argd(argc, argv, "--rawbins", 0) > 0) {
            printf("RAW_BIN %d %ld", b, cnt); for (double v : sum) printf(" %.17g", v); printf("\n");
        }
        bins.push_back(sum); nmeas += cnt;
    }
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); double cpusecs = double(std::clock() - c0) / CLOCKS_PER_SEC;
    int Bn = bins.size();
    auto stat = [&](auto f) {
        if (Bn < 2 || nmeas == 0) throw std::runtime_error("Need at least two bins and closed measurements");
        for (long c : closed_counts) if (nmeas == c) throw std::runtime_error("Jackknife needs closed visits in more than one bin");
        vector<double> tot(NV, 0); for (auto& b : bins) for (int i = 0; i < NV; i++) tot[i] += b[i];
        vector<double> full(NV); for (int i = 0; i < NV; i++) full[i] = tot[i] / nmeas;
        double est = f(full); vector<double> jk(Bn); double mj = 0, m2 = 0;
        for (int b = 0; b < Bn; b++) { vector<double> x(NV); for (int i = 0; i < NV; i++) x[i] = (tot[i] - bins[b][i]) / (nmeas - closed_counts[b]); jk[b] = f(x); mj += jk[b]; }
        mj /= Bn; for (double x : jk) m2 += (x - mj) * (x - mj);
        return std::make_pair(est, std::sqrt(m2 * (Bn - 1) / Bn));
    };
    auto sre = stat([](const vector<double>& x) { return x[0]; });
    auto sim = stat([](const vector<double>& x) { return x[1]; });
    printf("beta=%g sub=%d reroot=%d reversal=%d exchange=%d worm=%d mu=%g seed=%lu\n", beta, (int)S.use_sub, (int)S.use_reroot, (int)S.use_rev, (int)S.use_exch, (int)S.use_worm, S.mu, seed);
    printf("PHASE re %.6f %.6f im %.6f %.6f\n", sre.first, sre.second, sim.first, sim.second);
    printf("RESULT");
    auto pr = [&](const string& nm, int i, double fac) {
        auto r = stat([&](const vector<double>& x) { return fac * x[i] / x[0]; });
        printf(" %s %.6f %.6f", nm.c_str(), r.first, r.second);
    };
    pr("E", 2, 1); pr("q", 3, 1);
    for (int t = 0; t < T; t++) if (m.offtag[t]) pr("H_" + m.tags[t], 4 + t, -1.0 / beta);
    for (int o = 0; o < NO; o++) { pr(m.obs_names[o] + "_root", 4 + T + o, 1); pr(m.obs_names[o], 4 + T + NO + o, 1); }
    printf("\n");
    double tau = 0;
    { int Nn = Eser.size(); double mu_ = 0; for (double x : Eser) mu_ += x; mu_ /= Nn;
      double var = 0; for (double x : Eser) var += (x - mu_) * (x - mu_); var /= Nn;
      int bs = std::max(1, Nn / 50), nb = Nn / bs; double vb = 0;
      for (int b = 0; b < nb; b++) { double s = 0; for (int i = 0; i < bs; i++) s += Eser[b*bs+i]; s /= bs; vb += (s - mu_) * (s - mu_); }
      vb /= nb; tau = 0.5 * vb * bs / var; }
    printf("steps_per_sweep %ld closed_frac %.4f ms_per_sweep %.4f tau_int_E %.2f weight_check_fail %ld cpu_ms_per_sweep %.4f q_max %ld\n",
           nsteps, (double)nmeas / nsweep, 1e3 * secs / nsweep, tau, S.weight_check_fail, 1e3 * cpusecs / nsweep, qmax);
    for (auto& kv : S.stats) printf("  acc %-9s %ld/%ld = %.3f\n", kv.first.c_str(), kv.second.second, kv.second.first,
                                    kv.second.first ? (double)kv.second.second / kv.second.first : 0.0);
    return 0;
}
