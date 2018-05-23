#ifndef LP_SEIDEL_HPP
#define LP_SEIDEL_HPP

#include <vector>
#include <iostream>
#include "utility.hpp"
using namespace std;

//#define lp_debug(x) do{cerr << x; }while(0)
#define lp_debug(x) do{}while(0)

/**
 *  Randomized LP in expected
 *  O(d! 4^d n)
 *  Does exact calculations.
 */
template<typename FLOAT>
class Lp_Seidel{
private:


    // orthogonal projection of 'vec' into 'plane'
    vector<FLOAT> proj_down(vector<FLOAT> const&vec, vector<FLOAT> const&plane, size_t j){
        assert(vec.size() <= plane.size() && plane.size()<=vec.size()+1);
        assert(j+1 < plane.size());
        assert(!!plane[j]);
        vector<FLOAT> ret (vec.size()-1);
        //FLOAT tmp;
        if(plane[j] < FLOAT(0)){
            for(size_t i=0;i+1<vec.size();++i){
                ret[i] = vec[j] * plane[i+(i>=j)] - vec[i+(i>=j)]* plane[j];
            }
        } else {
            for(size_t i=0;i+1<vec.size();++i){
                ret[i] = vec[i+(i>=j)]*plane[j] - vec[j]*plane[i+(i>=j)];
            }
        }
        return ret;
    }

    // orthogonal projection of 'vec' out of 'plane'
    vector<FLOAT> proj_up(vector<FLOAT> const&vec, vector<FLOAT> const&plane, size_t j){
        assert(vec.size()+1 == plane.size());
        assert(j+1 < plane.size());
        assert(!!plane[j]);
        vector<FLOAT> ret(vec.size()+1);
        copy(vec.begin(), vec.begin()+j, ret.begin());
        copy(vec.begin()+j, vec.end(), ret.begin()+j+1);
        for(size_t i=0;i<vec.size();++i){
            ret[j]+=vec[i]*plane[i+(i>=j)];
        }
        FLOAT denom = plane[j];
        if(denom < FLOAT(0)){
            denom = -denom;
        }
        for(size_t i=0;i<vec.size();++i){
            ret[i+(i>=j)]*=denom;
        }
        if(plane[j] >= FLOAT(0)){
            ret[j] = -ret[j];
        }
        return ret;
    }
    FLOAT planescal(vector<FLOAT> const&x, vector<FLOAT> const&a){
        assert(x.size() == a.size());
        FLOAT ret=0;
        for(size_t i=0;i<x.size();++i){
            ret+=x[i]*a[i];
        }
        return ret;
    }

    // solve lp recursively
    vector<FLOAT> solve(vector<vector<FLOAT> > const &A, vector<FLOAT> const&c, int d, FLOAT const& barier_loc){
        int n=A.size();
        if(d==1){ // base case: single dimension
            vector<FLOAT> ret(2);
            ret[0] = (c[0]<FLOAT(0) ? -barier_loc : barier_loc);
            ret[1] = 1ull;
            for(int i=0;i<n;++i){
                if(ret[0]*A[i][0]+ret[1]*A[i].back()>FLOAT(0)){
                    if(!A[i][0]){
                        lp_debug("infeasible single\n");
                        return vector<FLOAT>();
                    }
                    ret[0] = -A[i].back();
                    ret[1] = A[i][0];
                    /*if(ret[1].is_negative()){
                        ret[1].negate();
                        ret[0].negate();
                    }*/
                    if(ret[1] < FLOAT(0)){
                        ret[1] = -ret[1];
                        ret[0] = -ret[0];
                    }
                    lp_debug(" -> " << i << "\n");
                }
            }
            for(int i=0;i<n;++i){
                if(ret[0]*A[i][0]+ret[1]*A[i].back()>FLOAT(0)){
                    lp_debug("infeasible\n");
                    return vector<FLOAT>();
                }
            }
            return ret;
        }
        FLOAT barier_next = barier_loc * barier_loc;
        // initial solution
        vector<FLOAT> x(d+1);
        for(int i=0;i<d;++i){
            x[i] = (c[i]<FLOAT(0)?-barier_loc:barier_loc);
        }
        x.back() = FLOAT(1);
        for(size_t i=0;i<A.size();++i){
            if(planescal(x, A[i])>FLOAT(0)){
                int k = 0;
                while(k<d && !A[i][k]) ++k;
                // recurse
                if(k==d) {lp_debug("what?\n"); return vector<FLOAT>();} // degenerate failing plane??????
                vector<vector<FLOAT> > A2(i);
                for(size_t j=0;j<A2.size();++j){
                    A2[j] = proj_down(A[j], A[i], k);
                }
                shuffle(A2.begin(), A2.end(), Rng::get_engine());
                lp_debug(string(2*d, ' ') << i << "\n");
                vector<FLOAT> c2 = proj_down(c, A[i], k);
                vector<FLOAT> x2 = solve(A2, c2, d-1, barier_next);
                if(x2.empty()) return x2; // infeasible
                x = proj_up(x2, A[i], k);
                lp_debug(string(2*d, ' ') << ":");
                lp_debug("";for(auto const&e:x) lp_debug(" " << e));
                lp_debug("\n");
            }
        }
        return x;
    }

public:
    vector<FLOAT> solve(vector<vector<FLOAT> > const &A, vector<FLOAT> const&c, const FLOAT& barier = FLOAT((long long)1e9)){
        assert(A.empty() || A[0].size() == c.size()+1);
        return solve(A, c, c.size(), barier);
    }
    /**
     *  Maximize c^T x
     *  subject to Ax <= b
     *
     *  Returns empty vector if infeasible
     */
    vector<FLOAT> solve(vector<vector<FLOAT> > A, vector<FLOAT> const&b, vector<FLOAT> const&c){
        assert(A.size() == b.size());
        for(unsigned int i=0;i<A.size();++i){
            A[i].push_back(-b[i]);
        }
        return solve(A, c);
    }
};

#endif // LP_SEIDEL_HPP
