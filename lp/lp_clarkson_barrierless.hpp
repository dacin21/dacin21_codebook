#ifndef LP_CLARKSON_BARRIERLESS_HPP
#define LP_CLARKSON_BARRIERLESS_HPP

#include <vector>
#include <algorithm>
#include "lp_seidel_barrierless.hpp"

/**
 *  Randomized LP, runtime is expected
 *  O(d^2 n + d^4 sqrt{n} log n + d^{O(d)} log(n)) | two phase
 *  O(d n log n + d^{O(d)} log n)                  | one phase
 *  Does exact calculations.
 */
template<typename Big_Int, bool use_two_phase = true>
class Lp_Clarkson_Barrierless{
public:
    using Solution_Int = Barrier_Int<Big_Int>;
private:
    /**
     *  Returns a sub-multiset of size siz uniformly at random
     *  out of the set where i is present weight[i] times.
     *
     *  Runs in O(|weight| + siz^2) expected time.
     *  Could be optimized
     */
    vector<int> sample_subset(vector<int64_t> const&weight, unsigned int siz){
        int64_t total_weight = accumulate(weight.begin(), weight.end(), 0ll);
        vector<int64_t> samples;
        while(samples.size() < siz){
            int64_t new_sample = Rng::uniform<int64_t>(0, total_weight-1);
            if(find(samples.begin(), samples.end(), new_sample) == samples.end()){
                samples.push_back(new_sample);
            }
        }
        sort(samples.begin(), samples.end());
        vector<int> ret;
        int64_t left_weight = 0;
        for(unsigned int i=0, j=0;i<weight.size() && j<samples.size();){
            if(samples[j] < left_weight + weight[i]){
                ret.push_back(i);
                ++j;
            } else {
                left_weight+=weight[i];
                ++i;
            }
        }
        return ret;
    }
    /// violation check
    bool is_satisfied(vector<Solution_Int> const&x, vector<Big_Int> const&a){
        assert(x.size() == a.size());
        Solution_Int ret(0);
        for(size_t i=0;i<x.size();++i){
            ret+=x[i]*a[i];
        }
        return ret <= Solution_Int(0);
    }
    vector<Solution_Int> solve_two(vector<vector<Big_Int> > const&A, vector<Big_Int> const&c){
        const unsigned int sample_size = c.size()*c.size()*4;
        Lp_Seidel_Barierless<Big_Int> sub_lp;
        // to few constrains -> use other solver
        if(A.size() < sample_size){
            return sub_lp.solve(A, c);
        } else {
            int constraints = A.size();
            int variables = c.size();
            vector<int64_t> weight(constraints, 1);
            vector<Solution_Int> x;
            vector<vector<Big_Int> > subproblem_A;
            vector<char> is_violated(constraints, 0);
            for(unsigned int iteration=1;;++iteration){
                subproblem_A.clear();
                vector<int> subspace = sample_subset(weight, sample_size);
                for(int const&e:subspace){
                    subproblem_A.push_back(A[e]);
                }

                x = sub_lp.solve(subproblem_A, c);
                // infeasible case
                if(x.empty()){
                    return x;
                }

                int64_t total_violated = 0;
                for(int i=0;i<constraints;++i){
                    is_violated[i] = !is_satisfied(x, A[i]);
                    if(is_violated[i]){
                        total_violated+=weight[i];
                    }
                }
                if(total_violated == 0){
                    cerr << "Iterations: " << iteration;
                    cerr << ", max weight: " << *max_element(weight.begin(), weight.end()) << "\n";
                    break;
                }
                if(total_violated*3*variables <= accumulate(weight.begin(), weight.end(), 0ll)){
                    for(int i=0;i<constraints;++i){
                        if(is_violated[i]){
                            weight[i]*=2;
                        }
                    }
                    assert_msg(accumulate(weight.begin(), weight.end(), 0ll) < (1ll<<62), "Weight overflow");
                }
            }
            return x;
        }
    }
    vector<Solution_Int> solve_one(vector<vector<Big_Int> > const&A, vector<Big_Int> const&c){
        const unsigned int constraints = A.size(), variables = c.size();
        if(constraints <= variables*variables*6){
            return solve_two(A, c);
        } else {
            const unsigned int sqrt_constraints = sqrt(constraints);
            const unsigned int sample_size = variables * sqrt(constraints);
            vector<Solution_Int> x;
            vector<vector<Big_Int> > subproblem_A;
            vector<int> violations;
            for(unsigned int iteration=1;;++iteration){
                //function<vector<int>()> samp = [&, this](){return sample_subset(vector<int64_t>(constraints, 1), sample_size);};
                //vector<int> subspace = Timer::execute_timed<vector<int>>(samp, "Sampling");
                vector<int> subspace = sample_subset(vector<int64_t>(constraints, 1), sample_size);
                for(int const&e:subspace){
                    subproblem_A.push_back(A[e]);
                }

                x = solve_two(subproblem_A, c);
                // infeasible case
                if(x.empty()){
                    return x;
                }
                violations.clear();
                for(unsigned int i=0;i<constraints;++i){
                    if(!is_satisfied(x, A[i])){
                        violations.push_back(i);
                    }
                }
                cerr << "Violations: " << violations.size() << " / " << 2*sqrt_constraints << "\n";
                if(violations.empty()){
                    cerr << "Iterations: " << iteration;
                    cerr << ", used constraints:" << subproblem_A.size() << "\n";
                    break;
                }
                subproblem_A.erase(subproblem_A.end()-sample_size, subproblem_A.end());
                if(violations.size() <= 2*sqrt_constraints){
                    for(int const&e : violations){
                        subproblem_A.push_back(A[e]);
                    }
                }
            }
            return x;
        }
    }



public:
    vector<Solution_Int> solve(vector<vector<Big_Int> > const&A, vector<Big_Int> const&c){
        if(use_two_phase){
            return solve_one(A, c);
        } else {
            return solve_two(A, c);
        }
    }

    /**
     *  Maximize c^T x
     *  Subject to Ax <= b
     *
     *  Returns empty vector if infeasible
     */
    vector<Solution_Int> solve(vector<vector<Big_Int> > A, vector<Big_Int> const&b, vector<Big_Int> const&c){
        assert(A.size() == b.size());
        for(unsigned int i=0;i<A.size();++i){
            A[i].push_back(-b[i]);
        }
        return solve(A, c);
    }
};

#endif // LP_CLARKSON_HPP
