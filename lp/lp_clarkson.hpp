#ifndef LP_CLARKSON_HPP
#define LP_CLARKSON_HPP

/**
 *  randomized LP in expected O(d n log n + ~ d^4 2^d d! log n)
 */
template<typename Big_Int>
class Lp_Clarkson{
private:
    /**
     *  Returns a sub-multiset of size siz uniformly at random
     *  out of the set where i is present weight[i] times.
     *
     *  Runs in O(|weight| + siz^2) expected time.
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

    bool is_satisfied(vector<Big_Int> const&x, vector<Big_Int> const&a){
        assert(x.size() == a.size());
        Big_Int ret=0;
        for(size_t i=0;i<x.size();++i){
            ret+=x[i]*a[i];
        }
        return ret <= Big_Int(0);
    }
    vector<Big_Int> solve(vector<vector<Big_Int> > const&A, vector<Big_Int> const&c, const unsigned int sample_size){
        Lp_Seidel<Big_Int> sub_lp;
        // to few constrains -> use other solver
        if(A.size() < sample_size){
            return sub_lp.solve(A, c);
        } else {
            int equations = A.size();
            int variables = c.size();
            vector<int64_t> weight(equations, 1);
            vector<Big_Int> x;
            vector<vector<Big_Int> > subproblem_A;
            vector<char> is_violated(equations, 0);
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
                for(int i=0;i<equations;++i){
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
                    for(int i=0;i<equations;++i){
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



public:
    vector<Big_Int> solve(vector<vector<Big_Int> > const&A, vector<Big_Int> const&c){
        return solve(A, c, c.size()*c.size()*4);
    }

    /**
     *  Maximize c^T x
     *  Subject to Ax <= b
     */
    vector<Big_Int> solve(vector<vector<Big_Int> > A, vector<Big_Int> const&b, vector<Big_Int> const&c){
        assert(A.size() == b.size());
        for(unsigned int i=0;i<A.size();++i){
            A[i].push_back(-b[i]);
        }
        return solve(A, c);
    }


};

#endif // LP_CLARKSON_HPP
