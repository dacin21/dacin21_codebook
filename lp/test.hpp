#ifndef TEST_HPP
#define TEST_HPP

#include <iostream>
#include <fstream>
#include <functional>
#include "utility.hpp"
#include "bignum.hpp"
#include "lp_seidel.hpp"
#include "lp_clarkson.hpp"
#include <boost/multiprecision/cpp_int.hpp>

template<typename Big_Int>
void test_factorial(unsigned int n, bool debug = false){
    std::cout << "Computing " << n << "!\n";
    Big_Int val(1);
    for(unsigned int i=1;i<=n;++i){
        val*=i;
        if(debug)
            std::cerr << i << " : " << val << "\n";
    }
    std::cerr << "Printing\n";
    std::cout << val << endl;
}

template<typename Big_Int>
void test_addsub(int n){
    Big_Int a(5432123456789ll);
    Big_Int b(9876543210123ll);
    Big_Int c;
    for(int it=0;it<4*n;++it){
        std::cout << a << " + " << b << " = " << a+b << "\n";
        std::cout << b << " + " << a << " = " << b+a << "\n";
        std::cout << a << " - " << b << " = " << a-b << "\n";
        std::cout << b << " - " << a << " = " << b-a << "\n";
        c = a; c+=b;
        std::cout << a << " + " << b << " = " << c << "\n";
        c = b; c+=a;
        std::cout << b << " + " << a << " = " << c << "\n";
        c = a; c-=b;
        std::cout << a << " - " << b << " = " << c << "\n";
        c = b; c-=a;
        std::cout << b << " - " << a << " = " << c << "\n";
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
    std::cout << "Addmul test " << n << " " << k << "\n";;
    Rng::set_seed(635241);
    Big_Int ret(0);
    for(int i=0;i<n;++i){
        Big_Int cur(1);
        for(int j=0;j<k;++j){
            Big_Int mu(Rng::uniform(1ll, (long long)1e15));
            cur*=mu;
        }
        //std::cerr << cur << "\n";
        ret+=cur;
        //ret = ret + cur;
    }
    std::cout << ret << "\n";
}

void test_timer(){
    std::function<int64_t(int64_t)> fib = [&](int64_t x){return x<2?x:fib(x-1)+fib(x-2);};
    Timer::execute_timed<int64_t>(fib, "Timer test", 35);
    Timer::execute_timed<int64_t>(fib, "Timer test", 40);
}


template<typename Big_Int, typename Solver>
void test_lp(std::string const&path){
    std::cerr << "LP Test: " << path << "\n";
    std::ifstream in;
    in.open(path.c_str());
    int variables, constraints;
    in >> variables >> constraints;
    std::vector<std::vector<Big_Int> > A(constraints, std::vector<Big_Int>(variables));
    std::vector<Big_Int> b(constraints), c(variables);
    // read objective function
    for(auto &e:c) in >> e;
    // read constraints
    for(int i=0;i<constraints;++i){
        for(auto &e:A[i]) in >> e;
        in >> b[i];
    }
    std::function<std::vector<Big_Int>()> solve = [&](){
        Solver solver;
        return solver.solve(A, b, c);
    };
    std::vector<Big_Int> solution = Timer::execute_timed<vector<Big_Int>>(solve, "solve LP");
    std::cerr << "Solution:\n";
    for(auto &e:solution){
        std::cerr << e << " ";
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
    std::function<std::vector<Big_Int>()> solve = [&](){
        Solver solver;
        return solver.solve(A, b, c);
    };
    std::vector<Big_Int> solution = Timer::execute_timed<vector<Big_Int>>(solve, "solve enclosing circle LP");
    std::cerr << "Solution:\n";
    for(auto &e:solution){
        std::cerr << e << " ";
    }
    cerr << "\n";
    if(!solution.empty()){
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
using boost_int = boost::multiprecision::cpp_int;
using boost_int1024 = boost::multiprecision::int1024_t;
using boost_int512 = boost::multiprecision::int512_t;

void run_tests(){
    //test_factorial<Bigint_Fixedsize<100>>(200);
    //test_factorial<Bigint_Fixedsize<500>>(1000);
    //test_factorial<Bigint_Fixedsize_Fast<500>>(1000);

    //test_addsub<Bigint_Fixedsize_Fast<50>>(1);
    //test_addmul<boost_int>(2, 2);
    //test_addmul<Bigint_Fixedsize_Fast<50>>(2, 2);

    //test_timer();

    //test_lp<Bigint_Fixedsize<3>, Lp_Seidel<Bigint_Fixedsize<3> >>("small.lp");
    //test_lp<Bigint_Fixedsize<3>, Lp_Clarkson<Bigint_Fixedsize<3> >>("small.lp");
    //test_lp<boost::multiprecision::cpp_int, Lp_Clarkson<boost::multiprecision::cpp_int>>("small.lp");
    //test_lp<Bigint_Fixedsize_Fast<5>, Lp_Seidel<Bigint_Fixedsize_Fast<5> >>("small.lp");


    constexpr int annulus_n = 50000;
    for(int it=0;it<5;++it){
        //test_enclosing_annulus<boost_int, Lp_Seidel<boost_int>>(annulus_n);
        //test_enclosing_annulus<boost_int, Lp_Clarkson<boost_int>>(annulus_n);
        test_enclosing_annulus<Bigint_Fixedsize_Fast<13>, Lp_Clarkson<Bigint_Fixedsize_Fast<13>>>(annulus_n);
        //test_enclosing_annulus<Bigint_Fixedsize_Fast<13>, Lp_Seidel<Bigint_Fixedsize_Fast<13>>>(annulus_n);
        //test_enclosing_annulus<boost_int1024, Lp_Clarkson<boost_int1024>>(annulus_n);
        test_enclosing_annulus<boost_int512, Lp_Clarkson<boost_int512>>(annulus_n);
        cerr << "\n\n";
    }
    //test_enclosing_annulus<boost_int, Lp_Seidel<boost_int>>(4500);

    //test_enclosing_annulus<Bigint_Fixedsize<10>, Lp_Clarkson<Bigint_Fixedsize<10>>>(450);
}


#endif // TEST_HPP
