#ifndef TEST_HPP
#define TEST_HPP

#include <iostream>
#include <fstream>
#include <functional>
#include "utility.hpp"
#include "bignum.hpp"
#include "lp_seidel.hpp"
#include "lp_clarkson.hpp"
#include "lp_seidel_barrierless.hpp"
#include "lp_clarkson_barrierless.hpp"

using namespace std;

template<typename Big_Int>
void test_factorial(unsigned int n, bool debug = false){
    cout << "Computing " << n << "!\n";
    Big_Int val(1);
    for(unsigned int i=1;i<=n;++i){
        val*=i;
        if(debug)
            cerr << i << " : " << val << "\n";
    }
    cerr << "Printing\n";
    cout << val << endl;
}

template<typename Big_Int>
void test_addsub(int n){
    Big_Int a(5432123456789ll);
    Big_Int b(9876543210123ll);
    Big_Int c;
    for(int it=0;it<4*n;++it){
        cout << a << " + " << b << " = " << a+b << "\n";
        cout << b << " + " << a << " = " << b+a << "\n";
        cout << a << " - " << b << " = " << a-b << "\n";
        cout << b << " - " << a << " = " << b-a << "\n";
        c = a; c+=b;
        cout << a << " + " << b << " = " << c << "\n";
        c = b; c+=a;
        cout << b << " + " << a << " = " << c << "\n";
        c = a; c-=b;
        cout << a << " - " << b << " = " << c << "\n";
        c = b; c-=a;
        cout << b << " - " << a << " = " << c << "\n";
        a=-a;
        if(it%4 == 3){
            a*=it;
        } else {
            swap(a, b);
        }
    }
}

template<typename Big_Int>
void test_addmul(int n, int k){
    cout << "Addmul test " << n << " " << k << "\n";;
    Rng::set_seed(635241);
    Big_Int ret(0);
    for(int i=0;i<n;++i){
        Big_Int cur(1);
        for(int j=0;j<k;++j){
            Big_Int mu(Rng::uniform(1ll, (long long)1e15));
            cur*=mu;
        }
        //cerr << cur << "\n";
        ret+=cur;
        //ret = ret + cur;
    }
    cout << ret << "\n";
}

void test_timer(){
    function<int64_t(int64_t)> fib = [&](int64_t x){return x<2?x:fib(x-1)+fib(x-2);};
    Timer::execute_timed<int64_t>(fib, "Timer test", 35);
    Timer::execute_timed<int64_t>(fib, "Timer test", 40);
}


template<typename Big_Int, typename Solver>
void test_lp(string const&path){
    cerr << "LP Test: " << path << "\n";
    ifstream in;
    in.open(path.c_str());
    int variables, constraints;
    in >> variables >> constraints;
    vector<vector<Big_Int> > A(constraints, vector<Big_Int>(variables));
    vector<Big_Int> b(constraints), c(variables);
    // read objective function
    for(auto &e:c) in >> e;
    // read constraints
    for(int i=0;i<constraints;++i){
        for(auto &e:A[i]) in >> e;
        in >> b[i];
    }
    function<vector<Big_Int>()> solve = [&](){
        Solver solver;
        return solver.solve(A, b, c);
    };
    vector<Big_Int> solution = Timer::execute_timed<vector<Big_Int>>(solve, "solve LP");
    cerr << "Solution:\n";
    for(auto &e:solution){
        cerr << e << " ";
    }
    cerr << "\n";
}

template<typename Big_Int, typename Solver>
void test_enclosing_annulus(int n){
    cerr << "Enclosing annulus test (" << n << ")\n";
    cerr << "Int: " << typeid(Big_Int).name()<< ", Solver: " << typeid(Solver).name() << "\n";
    // fixed seed for consistency
    Rng::set_seed(1234);
    vector<pair<int, int> > pts(n);
    for(auto &e:pts){
        e.first = Rng::uniform(-10000, 10000);
        e.second = Rng::uniform(-10000, 10000);
    }
    //Rng::set_seed(24680);
    Rng::timebased_seed();
    vector<vector<Big_Int> > A;
    vector<Big_Int> b, c;
    c = {Big_Int(0), Big_Int(0), -1, 1};
    // 2x a + 2y b + (-a^2-b^2+r^2) <= x^2 y^2
    // -2x a + -2y b - (-a^2-b^2+R^2) <= -x^2 -y^2
    for(auto const&e:pts){
        A.push_back({Big_Int(-2*e.first), Big_Int(-2*e.second), -1, 0});
        b.emplace_back(-e.first*e.first-e.second*e.second);
        A.push_back({Big_Int(2*e.first), Big_Int(2*e.second), 0, 1});
        b.emplace_back(e.first*e.first+e.second*e.second);
    }
    function<vector<Big_Int>()> solve = [&](){
        Solver solver;
        return solver.solve(A, b, c);
    };
    vector<Big_Int> solution = Timer::execute_timed<vector<Big_Int>>(solve, "solve enclosing circle LP");
    cerr << "Solution:\n";
    for(auto &e:solution){
        cerr << e << " ";
    }
    cerr << "\n";
    if(solution.empty()) cerr << "infeasible\n";
    if(!solution.empty()){
		// BUG? should it be [0] and [1] instead??
        Big_Int r_sq = solution[2]*solution[4] - solution[1]*solution[1] + solution[2]*solution[2];
        Big_Int R_sq = solution[3]*solution[4] - solution[1]*solution[1] + solution[2]*solution[2];
        Big_Int area_numerator = solution[3] - solution[2];
        //cerr << "area^2: " << area_numerator << "/" << solution[4] << "\n";
        //cerr << "r^2: " << r_sq << "/" << solution[4]*solution[4] << "\n";
        //cerr << "R^2: " << R_sq << "/" << solution[4]*solution[4] << "\n";
        cerr << "area^2: " << (long double)area_numerator / (long double)solution[4] << "\n";
        cerr << "r^2: " << (long double)r_sq / (long double)(solution[4]*solution[4]) << "\n";
        cerr << "R^2: " << (long double)R_sq / (long double)(solution[4]*solution[4]) << "\n";
    }
    cerr << "\n";
}
template<typename Big_Int, typename Solver, typename Ret_Int>
void test_enclosing_annulus_2(int n){
    cerr << "Enclosing annulus test (" << n << ")\n";
    cerr << "Int: " << typeid(Big_Int).name()<< ", Solver: " << typeid(Solver).name() << "\n";
    // fixed seed for consistency
    Rng::set_seed(8125345);
    vector<pair<int, int> > pts(n);
    for(auto &e:pts){
        e.first = Rng::uniform(-1000, 1000);
        e.second = Rng::uniform(-1000, 1000);
    }
    Rng::set_seed(24680);
    //Rng::timebased_seed();
    vector<vector<Big_Int> > A;
    vector<Big_Int> b, c;
    c = {Big_Int(0), Big_Int(0), -1, 1};
    // 2x a + 2y b + (-a^2-b^2+r^2) <= x^2 y^2
    // -2x a + -2y b - (-a^2-b^2+R^2) <= -x^2 -y^2
    for(auto const&e:pts){
        A.push_back({Big_Int(-2*e.first), Big_Int(-2*e.second), -1, 0});
        b.emplace_back(-e.first*e.first-e.second*e.second);
        A.push_back({Big_Int(2*e.first), Big_Int(2*e.second), 0, 1});
        b.emplace_back(e.first*e.first+e.second*e.second);
    }
    function<vector<Ret_Int>()> solve = [&](){
        Solver solver;
        return solver.solve(A, b, c);
    };
    vector<Ret_Int> solution = Timer::execute_timed<vector<Ret_Int>>(solve, "solve enclosing circle LP");
    cerr << "Solution:\n";
    for(auto &e:solution){
        cerr << e << " ";
    }
    cerr << "\n";
    if(solution.empty()) cerr << "infeasible\n";
    if(!solution.empty()){
		// BUG? should it be [0] and [1] instead??
        long double r_sq = (long double)solution[2]*(long double)solution[4] - (long double)solution[1]*(long double)solution[1] + (long double)solution[2]*(long double)solution[2];
        long double R_sq = (long double)solution[3]*(long double)solution[4] - (long double)solution[1]*(long double)solution[1] + (long double)solution[2]*(long double)solution[2];
        long double area_numerator = solution[3] - solution[2];
        //cerr << "area^2: " << area_numerator << "/" << solution[4] << "\n";
        //cerr << "r^2: " << r_sq << "/" << solution[4]*solution[4] << "\n";
        //cerr << "R^2: " << R_sq << "/" << solution[4]*solution[4] << "\n";
        cerr << "area^2: " << area_numerator / (long double)solution[4] << "\n";
        cerr << "r^2: " << r_sq / ((long double)solution[4]*(long double)solution[4]) << "\n";
        cerr << "R^2: " << R_sq / ((long double)solution[4]*(long double)solution[4]) << "\n";
    }
    cerr << "\n";
}


void test_codeforces_549E(){
    // Sasha circle
    auto codeforces_549E = [&](){
        // 10 might be enough too, didn't check
        using DOUBLE = Bigint_Fixedsize_Fast<12>;
        int N, M;
        cin >> N >> M;
        vector<vector<DOUBLE> > A;
        vector<DOUBLE> b, c(4);
        A.reserve(N+M);
        vector<pair<int, int> > p2(N), p(M);
        for(int i=0;i<N;++i){
            cin >> p2[i].first >> p2[i].second;
        }
        for(int i=0;i<M;++i){
            cin >> p[i].first >> p[i].second;
        }
        for(int it=0;it<2;++it){
            N=p.size();
            M=p2.size();
            A.clear();
            b.clear();
            for(int i=0;i<N;++i){
                int x, y, z;
                tie(x, y) = p[i];
                z=x*x+y*y;
                //A.push_back({x, y, 1, 0, z});
                A.push_back({x, y, 1, 1});
                b.push_back(z);
            }
            for(int i=0;i<M;++i){
                int x, y, z;
                tie(x, y) = p2[i];
                z=x*x+y*y;
                //A.push_back({-x, -y, 0, 1, -z});
                A.push_back({-x, -y, -1, 1});
                b.push_back(-z);
            }
            //c[2] = c[3] = 1;
            c[3] = 1;
            //Lp_Seidel<DOUBLE> solver;
            Lp_Clarkson_Barrierless<DOUBLE> solver;
            //Lp_Clarkson<DOUBLE, true> solver;
            vector<Lp_Clarkson_Barrierless<DOUBLE>::Solution_Int> x = solver.solve(A, b, c);
            //cerr << "len:\n";
            //for(size_t i=0;i<x.size();++i) cerr << x[i].data.size() << "\n";
            //cerr << "\n";
            if(!x.empty()){
                Lp_Clarkson_Barrierless<DOUBLE>::Solution_Int res(0);
                for(size_t i=0;i<c.size();++i) res+=x[i]*c[i];
                Lp_Clarkson_Barrierless<DOUBLE>::Solution_Int dist = res;
                if(dist>Lp_Clarkson_Barrierless<DOUBLE>::Solution_Int(0)){
                    cout << "YES\n";
                    return;
                }
            } else {
                //cerr << "nope.\n";
            }
            p.swap(p2);
        }
        cout << "NO\n";
    };
    freopen("codeforces_549E.in", "r", stdin);
    int TTT; cin >> TTT;
    while(TTT--){
        codeforces_549E();
    }
}


signed test_codejam_2017_4_D(){
    using DOUBLE = Bigint_Fixedsize_Fast<25>;

    freopen("codejam_d.in", "r", stdin);
    freopen("codejam_d.out", "w", stdout);
    cin.tie(0); ios_base::sync_with_stdio(false);
    int T; cin >> T;
    for(int cas=1;cas<=T;++cas){
        cout << "Case #" << cas << ": ";
        //cerr << "Case #" << cas << ": ";
        int n; cin >> n;
        vector<vector<DOUBLE> > A;
        vector<DOUBLE> b, c(4);
        A.reserve(n);
        lp_debug("n: " << n);
        for(int i=0;i<n;++i){
            int x, y, z;
            cin >> x >> y >> z;
            A.push_back({x, y, z, 1});
            b.push_back(0);
        }
        c[3] = 1;
        Lp_Clarkson_Barrierless<DOUBLE> solver;

		vector<Lp_Clarkson_Barrierless<DOUBLE>::Solution_Int> x = solver.solve(A, b, c);
		bool did = false;
		if(!x.empty()){
            if(x[3] > Lp_Clarkson_Barrierless<DOUBLE>::Solution_Int(0)){
                did = true;
            }
		}
        if(did) cout << "NO\n";
        else cout << "YES\n";
        //cerr << (did?"NO\n":"YES\n");
    }
    return 0;
}
void run_tests(){
    //test_factorial<Bigint_Fixedsize<100>>(200);
    //test_factorial<Bigint_Fixedsize<500>>(1000);
    //test_factorial<Bigint_Fixedsize_Fast<500>>(1000);

    //test_addsub<Bigint_Fixedsize_Fast<50>>(1);
    //test_addmul<Bigint_Fixedsize_Fast<50>>(2, 2);

    //test_timer();

    //test_lp<Bigint_Fixedsize<3>, Lp_Seidel<Bigint_Fixedsize<3> >>("small.lp");
    //test_lp<Bigint_Fixedsize<3>, Lp_Clarkson<Bigint_Fixedsize<3> >>("small.lp");
    //test_lp<Bigint_Fixedsize_Fast<5>, Lp_Seidel<Bigint_Fixedsize_Fast<5> >>("small.lp");


    constexpr int annulus_n = 100000;
    for(int it=0;it<1;++it){
        //test_enclosing_annulus<Bigint_Fixedsize_Fast<15>, Lp_Seidel<Bigint_Fixedsize_Fast<15>>>(annulus_n);
        //test_enclosing_annulus<Bigint_Fixedsize_Fast<15>, Lp_Clarkson<Bigint_Fixedsize_Fast<15>>>(annulus_n);
        //test_enclosing_annulus<Bigint_Fixedsize_Fast<15>, Lp_Clarkson<Bigint_Fixedsize_Fast<15>, true>>(annulus_n);
        test_enclosing_annulus_2<Bigint_Fixedsize_Fast<15>, Lp_Seidel<Bigint_Fixedsize_Fast<15>>, Bigint_Fixedsize_Fast<15>>(annulus_n);
        test_enclosing_annulus_2<Bigint_Fixedsize_Fast<15>, Lp_Clarkson_Barrierless<Bigint_Fixedsize_Fast<15>, true>, Barrier_Int<Bigint_Fixedsize_Fast<15>>>(annulus_n);
        //test_enclosing_annulus<Bigint_Fixedsize_Fast<25>, Lp_Clarkson<Bigint_Fixedsize_Fast<25>, true>>(annulus_n);
        //test_enclosing_annulus<Bigint_Fixedsize_Fast<13>, Lp_Seidel<Bigint_Fixedsize_Fast<13>>>(annulus_n);
        cerr << "\n\n";
    }

    //test_enclosing_annulus<Bigint_Fixedsize<10>, Lp_Clarkson<Bigint_Fixedsize<10>>>(450);
}


#endif // TEST_HPP
