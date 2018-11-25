// Various markov-chain computations in O(N^3)

template<typename T>
vector<T> matmul(vector<vector<T> > const&a, vector<T> const&b){
    const int n = a.size(), m = b.size();
    assert(a.empty() || (int)a[0].size() == m);
    vector<T> ret(n, T{0});
    for(int i=0;i<n;++i){
        for(int j=0;j<m;++j){
            ret[i]+= a[i][j] * b[j];
        }
    }
    return ret;
}

template<typename T = long double>
struct LUP{
private:
    const int n;
    const T eps;
    vector<vector<T> > data;
    vector<int> p;

    void decompose(){
        for(int i=0;i<n;++i){
            if(1){
                int k = i;
                for(int j=i+1;j<n;++j){
                    if(abs(data[j][i]) > abs(data[k][i])){
                        k=j;
                    }
                }
                data[k].swap(data[i]);
                swap(p[i], p[k]);
            }
            if(abs(data[i][i]) < eps){
                cerr << "Warning: Matrix singular\n";
            } else {
                {
                    T factor = T{1} / data[i][i];
                    for(int j=i+1;j<n;++j){
                        data[j][i]*=factor;
                    }
                }
                for(int j=i+1;j<n;++j){
                    T factor = data[j][i];
                    for(int k=i+1;k<n;++k){
                        data[j][k]-= factor*data[i][k];
                    }
                }
            }
        }
    }

    void solve_L(vector<T> &a) const {
        for(int i=0;i<n;++i){
            for(int j=0;j<i;++j){
                a[i]-= a[j] * data[i][j];
            }
        }
    }
    void solve_R(vector<T> &a) const {
        for(int i = n-1;i>=0;--i){
            for(int j=i+1;j<n;++j){
                a[i] -= a[j] * data[i][j];
            }
            a[i] /= data[i][i];
        }
    }
    void solve_P(vector<T> &a) const {
        static vector<T> tmp;
        tmp = a;
        for(int i=0;i<n;++i){
            a[i] = tmp[p[i]];
        }
    }

    T scal(vector<T> const&u, vector<T> const&v) const {
        assert(u.size() == v.size());
        T ret{0};
        for(size_t i=0;i<u.size();++i){
            ret+=u[i]*v[i];
        }
        return ret;
    }


public:
    LUP(vector<vector<T> > data_, const T eps_ = 1e-9):n(data_.size()), eps(eps_), data(move(data_)), p(n){
        iota(p.begin(), p.end(), 0);
        decompose();
    }
    /// A^{-1} b
    vector<T> solve(vector<T> b) const {
        solve_P(b);
        solve_L(b);
        solve_R(b);
        return b;
    }

    /**
     *  computes: (A + u v^T)^{-1} b
     *
     *  Sherman–Morrison formula
     *  (A + u v^T)^{-1} = A^{-1} - (A^{-1} u v^T A^{-1}) / (1 + v^T A^{-1} u)
     */
    vector<T> solve_rank_one_update(vector<T> b, vector<T> const&u, vector<T> const&v) const {
        assert((int)b.size() == n && b.size() == u.size() && b.size() == v.size());
        // A^{-1} u
        vector<T> Ai_u = solve(u);
        // 1 + v^T A^{-1} u
        T denom = T{1} + scal(v, Ai_u);
        if(abs(denom) < eps) cerr << "Warning: Matrix singular\n";
        // A^{-1} b
        b = solve(move(b));
        // v^T A^{-1}
        T numer = scal(v, b);
        // (v^T A^{-1}) / (1 + v^T A^{-1} u)
        T factor = numer / denom;
        for(size_t i=0;i<b.size();++i){
            b[i] -= Ai_u[i] * factor;
        }
        return b;
    }
    /**
     *  computes: (A + U V^T)^{-1} b
     *
     *  Woodbury matrix identity
     *  (A + U V^T)^{-1} = A^{-1} - A^{-1} U (I + V^T A^{-1} U)^{-1} V^T A^{-1}
     */
    template<size_t m>
    vector<T> solve_rank_update(vector<T> b, array<vector<T>, m> const&U, array<vector<T>, m> const&V) const {
        if(n == 0) return solve(b);
        assert((int)b.size() == n && b.size() == U[0].size() && b.size() == V[0].size());
        // A^{-1} U
        array<vector<T>, m> Ai_U;
        for(size_t i=0;i<m;++i){
            Ai_U[i] = solve(U[i]);
        }
        // 1 + V^T A^{-1} U
        vector<vector<T> > denom(m, vector<T>(m, 0));
        for(size_t i=0;i<m;++i){
            for(size_t j=0;j<m;++j){
                denom[i][j] = scal(V[i], Ai_U[j]) + (i==j);
            }
        }
        LUP denom_lup(denom, eps);

        // A^{-1} b
        b = solve(move(b));
        // V^T A^{-1} b
        vector<T> numer(n);
        for(size_t i=0;i<m;++i){
            numer[i] = scal(V[i], b);
        }
        // (V^T A^{-1}) / (1 + V^T A^{-1} U)
        numer = denom_lup.solve(move(numer));
        // compute answer
        for(size_t i=0;i<m;++i){
            for(int j=0;j<n;++j){
                b[j] -= Ai_U[i][j] * numer[i];
            }
        }
        return b;
    }
};


template<typename T = long double>
class Markov_Chain{
private:
    const int n;
    const T eps;
    vector<vector<pair<int, T> > > g;
    /// solve over-determined system A*x = b
    static vector<T> linsolve(vector<vector<T>> A, vector<T> b, const T eps){
        assert(A.size() == b.size());
        if(A.empty()) return vector<T>();
        const int n = A.size(), m = A[0].size();
        assert(n >= m);
        for(int i=0;i<m;++i){
            for(int j=i;j<n;++j){
                if(abs(A[j][i]) > abs(A[i][i])){
                    A[i].swap(A[j]);
                    swap(b[i], b[j]);
                }
            }
            if(abs(A[i][i]) > eps)
            for(int j=0;j<n;++j) if(j!=i){
                T factor = A[j][i]/A[i][i];
                for(int k=i;k<m;++k){
                    A[j][k]-=factor*A[i][k];
                }
                b[j]-=factor*b[i];
            }
        }
        for(int i=0;i<m;++i){
            if(abs(b[i]) + abs(A[i][i]) > eps)
                b[i]/=A[i][i];
        }
        b.resize(m);
        return b;
    }

    vector<vector<T> > build_matrix() const {
        vector<vector<T> > A(n, vector<T>(n, 0));
        for(int i=0;i<n;++i){
            for(auto const&e:g[i]){
                A[e.first][i] = e.second;
            }
        }
        return A;
    }
    vector<vector<T> > build_matrix_transpose() const {
        vector<vector<T> > A(n, vector<T>(n, 0));
        for(int i=0;i<n;++i){
            for(auto const&e:g[i]){
                A[i][e.first] = e.second;
            }
        }
        return A;
    }

public:
    Markov_Chain(int n_, T  eps_ = 1e-9):n(n_), eps(eps_), g(n){}
    void add_edge(int from, int to, T prob){
        g[from].emplace_back(to, prob);
    }
    void assert_probabilities() const {
        for(int i=0;i<n;++i){
            T sum = 0;
            for(auto const&e:g[i]){
                sum+=e.second;
            }
            assert(abs(T{1}-sum) < eps);
        }
    }
    vector<T> reach_probability(int t) const {
        vector<vector<T> > A(n, typename decltype(A)::value_type(n, T{0}));
        vector<T> b(n, 0);
        for(int i=0;i<n;++i){
            for(auto const&e:g[i]){
                A[i][e.first]-= e.second;
            }
            swap(b[i], A[i][t]);
            b[i]*= -1;
            A[i][i]+= 1;
        }
        auto ans = linsolve(move(A), move(b), eps);
        return ans;
    }
    /// expected time to reach t
    vector<T> reach_time(int t) const {
        vector<vector<T> > A(n, typename decltype(A)::value_type(n, T{0}));
        vector<T> b(n, 1);
        for(int i=0;i<n;++i){
            for(auto const&e:g[i]){
                A[i][e.first]-= e.second;
            }
            A[i][i]+= 1;
            A[i][t] = 0;
        }
        auto ans = linsolve(move(A), move(b), eps);
        return ans;
    }

    vector<T> stationary_distribution() const {
        vector<vector<T> > A = build_matrix();
        for(int i=0;i<n;++i) --A[i][i];
        A.emplace_back(n, 1);
        vector<T> b(n+1, 0);
        b.back() = 1;
        auto ans = linsolve(move(A), move(b), eps);
        return ans;
    }
    /** 
	 *  ans[i][j] = time for j -> i
	 *  chain has to be irreducible
	 */
    vector<vector<T> > reach_time_all_pairs() const {
        vector<T> s = stationary_distribution();
        vector<vector<T> > ret(n, vector<T>(n, 0));
        vector<vector<T> > A = build_matrix_transpose();
        for(int i=0;i<n;++i){
            A[i][i] -= 1;
        }
        array<vector<T>, 2> U, V;
        for(auto &e:U) e.resize(n, 0);
        for(auto &e:V) e.resize(n, 0);
        V[0][0] = 1;
        for(int i=0;i<n;++i){
            U[0][i] = A[i][0];
            A[i][0] = 0;
        }
        U[0][0]+= 1;
        A[0][0] = -1;
        vector<T> b(n, -1);
        b[0] = 0;
        LUP<T> lup(A, eps);
        ret[0] = lup.solve(b);
        auto res = matmul(A, ret[0]);
        for(int i=1;i<n;++i){
            b[i-1] = -1;
            b[i] = 0;
            V[1][i-1] = 0;
            V[1][i] = -1;
            for(int j=0;j<n;++j){
                U[1][j] = A[j][i];
            }
            U[1][i]-= 1;

            ret[i] = lup.solve_rank_update(b, U, V);
        }
        for(int i=0;i<n;++i){
            ret[i][i] = T{1} / s[i];
        }
        return ret;
    }
    int N() const {
        return n;
    }
};
