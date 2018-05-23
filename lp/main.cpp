
#include "test.hpp"

}

template<size_t len>
using Big_Int = lowdim_lp::Bigint_Fixedsize_Fast<len>;
template<size_t len>
using Lp_Solver = lowdim_lp::Lp_Clarkson_Barrierless<Big_Int<len>>;

#include <bits/stdc++.h>
using namespace std;
using lowdim_lp::Rng;
using lowdim_lp::Timer;


int main()
{
    //test_codeforces_549E();
    //Timer::execute_timed<signed>(test_codejam_2017_4_D, "Codejam test");

    lowdim_lp::run_tests();
    return 0;
}
