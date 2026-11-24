// upmr_qmc.cpp -- compact UPMR closed-walk QMC sampler used for the benchmarks
// in QMCofEverything_v17 (Sec. "Benchmark application").
//
// Implements the algorithm as specified in v17:
//   H = D_0 + sum_t P_t D_t, amplitudes d_t(z) evaluated on the departure state,
//   weights w = prod_k d_{i_k}(z_{k-1}) * e^{-beta[E_{z_0},...,E_{z_q}]},
//   sampling |w|; moves: basis-state, swap, pair insertion/deletion,
//   subcycle replacement (inverse-complement templates, App. B),
//   re-rooting, inner-loop reversal, segment exchange, worm (open-sector
//   multiplier mu only on open words).
//
// Models (single spin = Z_2 component, single boson = 0..nmax component):
//   sb : H = sum_i eps/2 sz_i + g sx_i (a_i + a_i^dag) + om n_i + J sum_<i,i+1> sz_i sz_{i+1}
//        (open chain of Ns pairs)
//   h3 : H = eps/2 sz + om n - h sx - kap (a + a^dag) - g sx (a + a^dag)
//
// Divided differences of exp(-beta x) are computed from the Opitz matrix
// with all diagonal entries shifted to be nonnegative, so that the Taylor
// series has only nonnegative terms (no cancellation; high relative accuracy).
//
// Build: clang++ -O3 -std=c++17 -o upmr_qmc upmr_qmc.cpp

#include <algorithm>
#include <chrono>
#include <ctime>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <random>
#include <string>
#include <vector>

using std::vector;

// ---------------------------------------------------------------- kernel
// returns log|f[x_0..x_L]| for f(x)=exp(-beta x); sign is (-1)^L.
static double log_edd(const vector<double>& E, double beta) {
    const int L = (int)E.size() - 1;
    double ymax = -1e300;
    for (double e : E) ymax = std::max(ymax, beta * e);
    if (L == 0) return -beta * E[0];
    vector<double> a(L + 1), v(L + 1, 0.0), t(L + 1);
    double amax = 0;
    for (int i = 0; i <= L; i++) { a[i] = ymax - beta * E[i]; amax = std::max(amax, a[i]); }
    const double N = amax + 1.0;
    v[L] = 1.0;                       // term_0 = e_L
    double logscale = 0.0, S = 0.0;   // true term = v * exp(logscale)
    for (int k = 1;; k++) {
        // t = Z' v / k  (upper bidiagonal, diag a, superdiag 1)
        double mx = 0;
        for (int i = 0; i < L; i++) { t[i] = (a[i] * v[i] + v[i + 1]) / k; mx = std::max(mx, t[i]); }
        t[L] = a[L] * v[L] / k; mx = std::max(mx, t[L]);
        v.swap(t);
        if (mx > 1e150 || (k < L && mx > 0 && mx < 1e-150)) { for (double& x : v) x /= mx; logscale += std::log(mx); S /= mx; }
        if (k >= L) {
            S += v[0];
            if (k > 2 * N && mx < 1e-18 * S) break;
            if (k > 100000) { fprintf(stderr, "log_edd: no convergence (q=%d); aborting\n", L); exit(2); }
        }
    }
    return L * std::log(beta) - ymax + std::log(S) + logscale;
}

// ---------------------------------------------------------------- model
struct Label { vector<std::pair<int,int>> disp; int inv; int kind; int site; int dir; };
// kind: 0 = sx (a+a^dag) branch [sb/h3 coupling], 1 = plain sx (h3), 2 = plain boson (h3)

struct Model {
    std::string name;
    int Ns = 1, nmax = 8;
    double eps = 1, om = 1, g = 0.6, J = 0, h = 0, kap = 0;
    int ncomp() const { return 2 * Ns; }      // comps 2i: spin i, 2i+1: boson i
    vector<Label> labels;
    void build() {
        labels.clear();
        for (int i = 0; i < Ns; i++) {
            int base = (int)labels.size();
            labels.push_back({{{2*i,1},{2*i+1,+1}}, base+1, 0, i, +1});   // sx a^dag
            labels.push_back({{{2*i,1},{2*i+1,-1}}, base,   0, i, -1});   // sx a
            if (name == "h3") {
                int b2 = (int)labels.size();
                labels.push_back({{{2*i,1}}, b2, 1, i, 0});               // sx
                labels.push_back({{{2*i+1,+1}}, b2+2, 2, i, +1});         // a^dag
                labels.push_back({{{2*i+1,-1}}, b2+1, 2, i, -1});         // a
            }
        }
    }
    static int sz(int s) { return s == 0 ? 1 : -1; }
    double E(const int* z) const {
        double e = 0;
        for (int i = 0; i < Ns; i++) e += 0.5 * eps * sz(z[2*i]) + om * z[2*i+1];
        for (int i = 0; i + 1 < Ns; i++) e += J * sz(z[2*i]) * sz(z[2*i+2]);
        return e;
    }
    // amplitude on departure state z; 0 if the hop leaves the physical domain
    double amp(int t, const int* z) const {
        const Label& L = labels[t];
        int n = z[2*L.site+1];
        double bamp = 1.0;
        if (L.dir == +1) { if (n >= nmax) return 0; bamp = std::sqrt(n + 1.0); }
        if (L.dir == -1) { if (n <= 0) return 0; bamp = std::sqrt((double)n); }
        if (L.kind == 0) return (name == "h3" ? -g : g) * bamp;
        if (L.kind == 1) return -h;
        return -kap * bamp;
    }
    void apply(int t, int* z) const {
        for (auto& d : labels[t].disp) {
            if (d.first % 2 == 0) z[d.first] ^= 1; else z[d.first] += d.second;
        }
    }
    bool in_domain(const int* z) const {
        for (int i = 0; i < Ns; i++) if (z[2*i+1] < 0 || z[2*i+1] > nmax) return false;
        return true;
    }
};

// ---------------------------------------------------------------- sampler
struct Eval { bool ok; double logw; int sign; vector<double> En; vector<int> states; };

struct Sampler {
    Model& m; double beta; std::mt19937_64 rng;
    int nc;
    vector<int> root; vector<int> word;   // configuration (word may be open in worm mode)
    double logw = 0; int sgn = 1; vector<double> En; vector<int> st;
    bool use_sub = false, use_nonlocal = false, use_worm = false;
    double mu = 1.0;
    vector<vector<int>> templ; vector<int> templ_inv;   // template count vectors
    std::map<std::string, std::pair<long,long>> stats;  // attempts, accepts
    long weight_check_fail = 0;

    Sampler(Model& mm, double b, unsigned long seed) : m(mm), beta(b), rng(seed) {
        nc = m.ncomp(); root.assign(nc, 0);
        auto e = eval(root, word); logw = e.logw; sgn = e.sign; En = e.En; st = e.states;
    }
    double U() { return std::uniform_real_distribution<double>(0, 1)(rng); }
    int Ri(int n) { return std::uniform_int_distribution<int>(0, n - 1)(rng); }

    bool closed_of(const vector<int>& s, int q) const {
        for (int c = 0; c < nc; c++) if (s[q*nc+c] != s[c]) return false;
        return true;
    }
    Eval eval(const vector<int>& r, const vector<int>& w) const {
        Eval e; e.ok = true; e.sign = 1; double la = 0;
        int q = (int)w.size();
        e.states.resize((size_t)(q + 1) * nc);
        std::copy(r.begin(), r.end(), e.states.begin());
        e.En.resize(q + 1); e.En[0] = m.E(&e.states[0]);
        for (int k = 0; k < q; k++) {
            int* cur = &e.states[(size_t)k*nc]; int* nxt = &e.states[(size_t)(k+1)*nc];
            double d = m.amp(w[k], cur);
            if (d == 0) { e.ok = false; e.logw = -INFINITY; return e; }
            if (d < 0) e.sign = -e.sign;
            la += std::log(std::fabs(d));
            std::copy(cur, cur + nc, nxt); m.apply(w[k], nxt);
            e.En[k+1] = m.E(nxt);
        }
        double lf = log_edd(e.En, beta);
        if (q % 2) e.sign = -e.sign;
        e.logw = la + lf;
        bool closed = closed_of(e.states, q);
        if (!closed) e.logw += std::log(mu);   // open-sector multiplier
        return e;
    }
    bool is_closed() const { return closed_of(st, (int)word.size()); }
    bool metro(const vector<int>& r2, const vector<int>& w2, double logR, const char* nm, bool expect_equal = false) {
        auto& s = stats[nm]; s.first++;
        Eval e = eval(r2, w2);
        if (!e.ok) return false;
        if (expect_equal && std::fabs(e.logw - logw) > 1e-8) weight_check_fail++;
        double a = e.logw - logw + logR;
        if (a >= 0 || U() < std::exp(a)) {
            root = r2; word = w2; logw = e.logw; sgn = e.sign; En = e.En; st = e.states; s.second++;
            return true;
        }
        return false;
    }
    // ---- moves
    void basis_move() {
        vector<int> r2 = root; int c = Ri(nc); int d = (U() < 0.5) ? 1 : -1;
        if (c % 2 == 0) r2[c] ^= 1; else r2[c] += d;
        if (!m.in_domain(&r2[0])) { stats["basis"].first++; return; }
        metro(r2, word, 0, "basis");
    }
    void swap_move() {
        int q = word.size(); if (q < 2) return;
        int j = Ri(q - 1); if (word[j] == word[j+1]) return;
        vector<int> w2 = word; std::swap(w2[j], w2[j+1]);
        metro(root, w2, 0, "swap");
    }
    void pair_move() {
        const double pdel = 0.5; int M = m.labels.size(); int q = word.size();
        if (U() < 1 - pdel) {
            int t = Ri(M); int g = Ri(q + 1);
            vector<int> w2 = word; w2.insert(w2.begin() + g, {t, m.labels[t].inv});
            metro(root, w2, std::log(pdel / (1 - pdel) * M), "pair_ins");
        } else {
            if (q < 2) { stats["pair_del"].first++; return; }
            int j = Ri(q - 1);
            if (word[j+1] != m.labels[word[j]].inv) { stats["pair_del"].first++; return; }
            vector<int> w2 = word; w2.erase(w2.begin() + j, w2.begin() + j + 2);
            metro(root, w2, std::log((1 - pdel) / pdel / M), "pair_del");
        }
    }
    static double log_norderings(const vector<int>& cnt) {
        int s = 0; double l = 0; for (int c : cnt) { s += c; l -= std::lgamma(c + 1.0); }
        return l + std::lgamma(s + 1.0);
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
        vector<int> ins;                                   // inverse complement
        vector<int> Tbar(tc.size(), 0);
        for (size_t t = 0; t < tc.size(); t++) for (int k = 0; k < tc[t] - wc[t]; k++) { ins.push_back(m.labels[t].inv); Tbar[m.labels[t].inv]++; }
        std::shuffle(ins.begin(), ins.end(), rng);
        vector<int> w2(word.begin(), word.begin() + j);
        w2.insert(w2.end(), ins.begin(), ins.end());
        w2.insert(w2.end(), word.begin() + j + r, word.end());
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
    bool same_state(int a, int b) const {
        for (int c = 0; c < nc; c++) if (st[(size_t)a*nc+c] != st[(size_t)b*nc+c]) return false; return true;
    }
    void reversal_move() {
        int q = word.size(); if (q < 2) return;
        int a = Ri(q + 1), b = Ri(q + 1); if (a == b) return; if (a > b) std::swap(a, b);
        if (!same_state(a, b)) { stats["reversal"].first++; return; }
        vector<int> w2 = word;
        for (int k = 0; k < b - a; k++) w2[a + k] = m.labels[word[b - 1 - k]].inv;
        metro(root, w2, 0, "reversal", true);
    }
    void exchange_move() {
        int q = word.size(); if (q < 2) return;
        // uniform composition (A,S1,B,S2,C) of q, |S1|,|S2|>=1: 4 bars among q+2 slots
        vector<int> slots(q + 2); for (int i = 0; i < q + 2; i++) slots[i] = i;
        for (int i = 0; i < 4; i++) std::swap(slots[i], slots[i + Ri(q + 2 - i)]);
        vector<int> bars(slots.begin(), slots.begin() + 4); std::sort(bars.begin(), bars.end());
        int pA = bars[0], pS1 = bars[1]-bars[0]-1+1, pB = bars[2]-bars[1]-1, pS2 = bars[3]-bars[2]-1+1, pC = q+1-bars[3];
        // positions (state indices)
        int s1a = pA, s1b = pA + pS1, s2a = s1b + pB, s2b = s2a + pS2;
        if (s2b + pC != q) { fprintf(stderr, "composition bug\n"); exit(1); }
        bool opposite = U() < 0.5;
        auto& s = stats[opposite ? "exch_opp" : "exch_eq"];
        vector<int> w2;
        vector<int> A(word.begin(), word.begin()+s1a), S1(word.begin()+s1a, word.begin()+s1b),
                    B(word.begin()+s1b, word.begin()+s2a), S2(word.begin()+s2a, word.begin()+s2b), C(word.begin()+s2b, word.end());
        if (!opposite) {
            if (!(same_state(s1a, s2a) && same_state(s1b, s2b))) { s.first++; return; }
            if (S1 == S2) { s.first++; return; }
            w2 = A; w2.insert(w2.end(), S2.begin(), S2.end()); w2.insert(w2.end(), B.begin(), B.end());
            w2.insert(w2.end(), S1.begin(), S1.end()); w2.insert(w2.end(), C.begin(), C.end());
        } else {
            if (!(same_state(s1a, s2b) && same_state(s1b, s2a))) { s.first++; return; }
            auto inv = [&](const vector<int>& S) { vector<int> R; for (int k = S.size() - 1; k >= 0; k--) R.push_back(m.labels[S[k]].inv); return R; };
            vector<int> S2i = inv(S2), S1i = inv(S1);
            w2 = A; w2.insert(w2.end(), S2i.begin(), S2i.end()); w2.insert(w2.end(), B.begin(), B.end());
            w2.insert(w2.end(), S1i.begin(), S1i.end()); w2.insert(w2.end(), C.begin(), C.end());
        }
        metro(root, w2, 0, opposite ? "exch_opp" : "exch_eq", true);
    }
    void worm_move() {
        int M = m.labels.size();
        if (U() < 0.5) {
            int t = Ri(M); vector<int> w2 = word; w2.push_back(t);
            metro(root, w2, std::log((double)M), "worm_app");
        } else {
            if (word.empty()) { stats["worm_ret"].first++; return; }
            vector<int> w2 = word; w2.pop_back();
            metro(root, w2, -std::log((double)M), "worm_ret");
        }
    }
    // one elementary step; move-type probabilities do not depend on the state
    void step() {
        // move menu
        vector<int> menu = {0, 1, 2};               // basis, swap, pair
        if (use_sub) menu.push_back(3);
        if (use_nonlocal) { menu.push_back(4); menu.push_back(5); menu.push_back(6); }
        if (use_worm) { menu.push_back(7); menu.push_back(7); }
        int mv = menu[Ri(menu.size())];
        bool closed = is_closed();
        if (!closed && mv != 7) return;             // self-loop in open sector
        switch (mv) {
            case 0: basis_move(); break; case 1: swap_move(); break; case 2: pair_move(); break;
            case 3: subcycle_move(); break; case 4: reroot_move(); break; case 5: reversal_move(); break;
            case 6: exchange_move(); break; case 7: worm_move(); break;
        }
    }
    // ---- estimators (closed configurations only)
    struct Meas { double E, szroot, szw, n, q, chiA, sgn; };
    Meas measure() const {
        Meas r; int q = word.size(); r.sgn = sgn; r.q = q;
        double sz0 = 0, n0 = 0;
        for (int i = 0; i < m.Ns; i++) { sz0 += Model::sz(st[2*i]); n0 += st[2*i+1]; }
        r.szroot = sz0 / m.Ns; r.n = n0 / m.Ns;
        double lf0 = log_edd(En, beta);
        if (q == 0) r.E = En[0];
        else {
            vector<double> E1(En.begin() + 1, En.end());
            r.E = En[0] - std::exp(log_edd(E1, beta) - lf0);
        }
        // omega_k-weighted sz:  omega_k = -f[E,E_k]/(beta f[E]) (positive)
        double szw = 0;
        for (int k = 0; k <= q; k++) {
            vector<double> E2 = En; E2.push_back(En[k]);
            double om = std::exp(log_edd(E2, beta) - lf0) / beta;
            double s = 0; for (int i = 0; i < m.Ns; i++) s += Model::sz(st[(size_t)k*m.ncomp()+2*i]);
            szw += om * s / m.Ns;
        }
        r.szw = szw;
        r.chiA = szw;   // d ln w / d eps = -beta * Ns/2 * szw  (per-site normalisation handled in post)
        return r;
    }
};

// ---------------------------------------------------------------- driver
static double argd(int argc, char** argv, const char* k, double def) {
    for (int i = 1; i + 1 < argc; i++) if (std::string(argv[i]) == k) return atof(argv[i+1]); return def; }
static std::string args(int argc, char** argv, const char* k, const char* def) {
    for (int i = 1; i + 1 < argc; i++) if (std::string(argv[i]) == k) return argv[i+1]; return def; }

int main(int argc, char** argv) {
    if (args(argc, argv, "--test", "")[0]) {   // kernel test: --test "b x0 x1 ..."
        std::string t = args(argc, argv, "--test", ""); vector<double> v; char* p = &t[0]; char* e;
        double b = strtod(p, &e); p = e; while (true) { double x = strtod(p, &e); if (e == p) break; v.push_back(x); p = e; }
        printf("%.15e\n", log_edd(v, b)); return 0;
    }
    Model m; m.name = args(argc, argv, "--model", "sb");
    m.Ns = (int)argd(argc, argv, "--Ns", 1); m.nmax = (int)argd(argc, argv, "--nmax", 8);
    m.eps = argd(argc, argv, "--eps", 1); m.om = argd(argc, argv, "--om", 1); m.g = argd(argc, argv, "--g", 0.6);
    m.J = argd(argc, argv, "--J", 0); m.h = argd(argc, argv, "--h", 0); m.kap = argd(argc, argv, "--kap", 0);
    m.build();
    double beta = argd(argc, argv, "--beta", 1);
    long therm = (long)argd(argc, argv, "--therm", 1e4), prod = (long)argd(argc, argv, "--sweeps", 1e5);
    int nbins = (int)argd(argc, argv, "--bins", 100);
    Sampler S(m, beta, (unsigned long)argd(argc, argv, "--seed", 1));
    S.use_sub = argd(argc, argv, "--sub", 0) > 0; S.use_nonlocal = argd(argc, argv, "--nonlocal", 0) > 0;
    S.use_worm = argd(argc, argv, "--worm", 0) > 0; S.mu = argd(argc, argv, "--mu", 1.0);
    if (m.name == "h3" && S.use_sub) {
        // labels per site: 0 sx a^dag, 1 sx a, 2 sx, 3 a^dag, 4 a.  template {sx a^dag, a, sx} and its inverse
        vector<int> t1(5, 0), t2(5, 0);
        t1[0] = 1; t1[4] = 1; t1[2] = 1;  t2[1] = 1; t2[3] = 1; t2[2] = 1;
        S.templ = {t1, t2}; S.templ_inv = {1, 0};
    }
    // a sweep has a FIXED number of steps (state-independent measurement times):
    // set from the mean word length over the second half of thermalization.
    long nsteps = 0; double qsum = 0; long qn = 0;
    for (long s = 0; s < therm; s++) {
        for (long k = (long)S.word.size() + S.nc + 1; k > 0; k--) S.step();
        if (s >= therm / 2) { qsum += S.word.size(); qn++; }
    }
    nsteps = (long)std::lround(qsum / std::max(1L, qn)) + S.nc + 1;
    if (argd(argc, argv, "--steps", 0) > 0) nsteps = (long)argd(argc, argv, "--steps", 0);
    auto steps_in_sweep = [&]() { return nsteps; };
    long qmax = 0; auto t0 = std::chrono::steady_clock::now(); std::clock_t c0 = std::clock();
    long per = prod / nbins; vector<vector<double>> bins;   // per bin: sums
    vector<double> acc(9, 0); long nmeas = 0, nsweep = 0;
    // per-sweep series of E for autocorrelation (closed measurements only)
    vector<double> Eseries;
    for (int b = 0; b < nbins; b++) {
        vector<double> sum(9, 0); long cnt = 0;
        for (long s = 0; s < per; s++) {
            for (long k = steps_in_sweep(); k > 0; k--) S.step();
            nsweep++; qmax = std::max(qmax, (long)S.word.size());
            if (!S.is_closed()) continue;
            auto r = S.measure(); cnt++;
            double sg = r.sgn;
            // sums: sgn, sgn*E, sgn*szroot, sgn*szw, sgn*n, sgn*q, sgn*szroot*szw, sgn*E^2-ish unused
            double v[9] = {sg, sg*r.E, sg*r.szroot, sg*r.szw, sg*r.n, sg*r.q, sg*r.szroot*r.szw, 0, 0};
            for (int i = 0; i < 9; i++) sum[i] += v[i];
            Eseries.push_back(r.E);
        }
        for (int i = 0; i < 9; i++) sum[i] /= std::max(1L, cnt);
        bins.push_back(sum); nmeas += cnt;
    }
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); double cpusecs = double(std::clock() - c0) / CLOCKS_PER_SEC;
    // jackknife over bins
    int B = bins.size();
    auto stat = [&](auto f) {
        vector<double> tot(9, 0); for (auto& b : bins) for (int i = 0; i < 9; i++) tot[i] += b[i];
        vector<double> full(9); for (int i = 0; i < 9; i++) full[i] = tot[i] / B;
        double est = f(full), m2 = 0; vector<double> jk(B);
        for (int b = 0; b < B; b++) { vector<double> x(9); for (int i = 0; i < 9; i++) x[i] = (tot[i] - bins[b][i]) / (B - 1); jk[b] = f(x); }
        double mj = 0; for (double x : jk) mj += x; mj /= B;
        for (double x : jk) m2 += (x - mj) * (x - mj);
        return std::make_pair(est, std::sqrt(m2 * (B - 1) / B));
    };
    auto E = stat([](const vector<double>& x) { return x[1] / x[0]; });
    auto szr = stat([](const vector<double>& x) { return x[2] / x[0]; });
    auto szw = stat([](const vector<double>& x) { return x[3] / x[0]; });
    auto n = stat([](const vector<double>& x) { return x[4] / x[0]; });
    auto q = stat([](const vector<double>& x) { return x[5] / x[0]; });
    // chi_z = d<sz>/d eps = -(beta Ns/2) [ <sz_root * szw> - <sz_root><szw> ]   (sz per site)
    double Nsd = m.Ns;
    auto chi = stat([&](const vector<double>& x) { return -(beta * Nsd / 2) * (x[6]/x[0] - (x[2]/x[0]) * (x[3]/x[0])); });
    // integrated autocorrelation time of E (in sweeps, between closed measurements) via binning
    double tau = 0;
    {
        int Nn = Eseries.size(); double mu_ = 0; for (double x : Eseries) mu_ += x; mu_ /= Nn;
        double var = 0; for (double x : Eseries) var += (x - mu_) * (x - mu_); var /= Nn;
        int bs = std::max(1, Nn / 50); int nb = Nn / bs; double vb = 0;
        for (int b = 0; b < nb; b++) { double s = 0; for (int i = 0; i < bs; i++) s += Eseries[b*bs+i]; s /= bs; vb += (s - mu_) * (s - mu_); }
        vb /= nb; tau = 0.5 * vb * bs / var;
    }
    printf("model=%s Ns=%d nmax=%d eps=%g om=%g g=%g J=%g h=%g kap=%g beta=%g sub=%d nonlocal=%d worm=%d mu=%g\n",
           m.name.c_str(), m.Ns, m.nmax, m.eps, m.om, m.g, m.J, m.h, m.kap, beta, (int)S.use_sub, (int)S.use_nonlocal, (int)S.use_worm, S.mu);
    auto sgnav = stat([](const vector<double>& x) { return x[0]; });
    printf("SIGN %.6f %.6f\n", sgnav.first, sgnav.second);
    printf("RESULT E %.6f %.6f szroot %.6f %.6f szw %.6f %.6f n %.6f %.6f q %.5f %.5f chi %.5f %.5f\n",
           E.first, E.second, szr.first, szr.second, szw.first, szw.second, n.first, n.second, q.first, q.second, chi.first, chi.second);
    printf("steps_per_sweep %ld closed_frac %.4f  ms_per_sweep %.4f  tau_int_E %.2f  weight_check_fail %ld cpu_ms_per_sweep %.4f q_max %ld\n",
           nsteps, (double)nmeas / nsweep, 1e3 * secs / nsweep, tau, S.weight_check_fail, 1e3 * cpusecs / nsweep, qmax);
    for (auto& kv : S.stats) printf("  acc %-9s %ld/%ld = %.3f\n", kv.first.c_str(), kv.second.second, kv.second.first,
                                    kv.second.first ? (double)kv.second.second / kv.second.first : 0.0);
    return 0;
}
