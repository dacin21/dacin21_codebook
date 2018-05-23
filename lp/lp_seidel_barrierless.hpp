// TODO: finish implementation
#ifndef LP_SEIDEL_BARIERLESS_HPP
#define LP_SEIDEL_BARIERLESS_HPP

#include <vector>
#include <iostream>
#include "utility.hpp"
using namespace std;

//#define lp_debug(x) do{cerr << x; }while(0)
#define lp_debug(x) do{}while(0)

/**
 *  Represents an integer of the form
 *  <val> + <inf_part> * inf
 *  where inf is a symbol bigger than any integer
 */
template<typename INT>
class Barrier_Int{
private:
    Barrier_Int(INT const&_val, INT const&_inf_part):val(_val), inf_part(_inf_part){}
public:
    Barrier_Int():val(0), inf_part(0){}
    explicit Barrier_Int(INT const&_val):val(_val), inf_part(0){}
    static Barrier_Int infinity(){
        return Barrier_Int(0, 1);
    }
    static Barrier_Int negative_infinity(){
        return Barrier_Int(0, -1);
    }

    Barrier_Int operator-()const{
        return Barrier_Int(-val, -inf_part);
    }
    Barrier_Int& operator+=(Barrier_Int const&o){
        val+=o.val;
        inf_part+=o.inf_part;
        return *this;
    }
    Barrier_Int operator+(Barrier_Int const&o)const{
        return Barrier_Int(val+o.val, inf_part+o.inf_part);
    }
    Barrier_Int& operator-=(Barrier_Int const&o){
        val-=o.val;
        inf_part-=o.inf_part;
        return *this;
    }
    Barrier_Int operator-(Barrier_Int const&o)const{
        return Barrier_Int(val-o.val, inf_part-o.inf_part);
    }
    Barrier_Int& operator*=(INT const&o){
        val*=o;
        inf_part*=o;
        return *this;
    }
    Barrier_Int operator*(INT const&o)const{
        return Barrier_Int(val*o, inf_part*o);
    }
    bool operator<(Barrier_Int const&o)const{
        if(inf_part != o.inf_part) return inf_part < o.inf_part;
        return val < o.val;
    }
    bool operator>(Barrier_Int const&o)const{
        return o < *this;
    }
    bool operator>=(Barrier_Int const&o)const{
        return !operator<(o);
    }
    bool operator<=(Barrier_Int const&o)const{
        return !operator>(o);
    }
    bool operator==(Barrier_Int const&o)const{
        return val == o.val && inf_part == o.inf_part;
    }
    bool operator!=(Barrier_Int const&o)const{
        return val!=o.val || inf_part != o.inf_part;
    }
    friend ostream& operator<<(ostream&o, Barrier_Int const&b){
        if(!b.inf_part){
            return o << b.val;
        }
        o << b.inf_part << numeric_limits<double>::infinity() << "+" << b.val;
        return o;
    }
    explicit operator long double()const{
        if(inf_part != INT(0)){
            return inf_part<INT(0) ? -numeric_limits<long double>::infinity() : numeric_limits<long double>::infinity();
        }
        return (long double)val;
    }
public:
    INT val;
    INT inf_part;
};

/**
 *  Randomized LP in expected
 *  O(d! 4^d n)
 *  Does exact calculations.
 *  Does not need a barrier bound for dealing with unbounded parts
 */
template<typename FLOAT>
class Lp_Seidel_Barierless{
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
    vector<Barrier_Int<FLOAT> > proj_up(vector<Barrier_Int<FLOAT> > const&vec, vector<FLOAT> const&plane, size_t j){
        assert(vec.size()+1 == plane.size());
        assert(j+1 < plane.size());
        assert(!!plane[j]);
        vector<Barrier_Int<FLOAT> > ret(vec.size()+1);
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
    Barrier_Int<FLOAT> planescal(vector<Barrier_Int<FLOAT> > const&x, vector<FLOAT> const&a){
        assert(x.size() == a.size());
        Barrier_Int<FLOAT> ret(0);
        for(size_t i=0;i<x.size();++i){
            ret+=x[i]*a[i];
        }
        return ret;
    }

    // solve lp recursively
    vector<Barrier_Int<FLOAT>> solve(vector<vector<FLOAT> > const &A, vector<FLOAT> const&c, int d){
        int n=A.size();
        if(d==1){ // base case: single dimension
            vector<Barrier_Int<FLOAT> > ret(2);
            if(c[0] != FLOAT(0)){
                ret[0] = (c[0]<FLOAT(0) ? Barrier_Int<FLOAT>::negative_infinity() : Barrier_Int<FLOAT>::infinity());
            }
            ret[1].val = FLOAT(1ull);
            for(int i=0;i<n;++i){
                if(ret[0]*A[i][0]+ret[1]*A[i].back()>Barrier_Int<FLOAT>(0)){
                    if(!A[i][0]){
                        lp_debug("infeasible single\n");
                        return vector<Barrier_Int<FLOAT>>();
                    }
                    ret[0] = Barrier_Int<FLOAT>(-A[i].back());
                    ret[1] = Barrier_Int<FLOAT>(A[i][0]);
                    if(ret[1] < Barrier_Int<FLOAT>(0)){
                        ret[1] = -ret[1];
                        ret[0] = -ret[0];
                    }
                    lp_debug(" -> " << i << " " << ret[0] << " " << ret[1] << "\n");
                }
            }
            for(int i=0;i<n;++i){
                if(ret[0]*A[i][0]+ret[1]*A[i].back()>Barrier_Int<FLOAT>(0)){
                    lp_debug("infeasible\n");
                    return vector<Barrier_Int<FLOAT>>();
                }
            }
            return ret;
        }
        // initial solution
        vector<Barrier_Int<FLOAT>> x(d+1);
        for(int i=0;i<d;++i){
            if(c[i] != FLOAT(0)){
                x[i] = (c[i]<FLOAT(0) ? Barrier_Int<FLOAT>::negative_infinity() : Barrier_Int<FLOAT>::infinity());
            }
        }
        x.back() = Barrier_Int<FLOAT>(1);
        for(size_t i=0;i<A.size();++i){
            if(planescal(x, A[i])>Barrier_Int<FLOAT>(0)){
                int k = 0;
                while(k<d && !A[i][k]) ++k;
                // recurse
                if(k==d) {lp_debug("what?\n"); return vector<Barrier_Int<FLOAT>>();} // degenerate failing plane??????
                vector<vector<FLOAT> > A2(i);
                for(size_t j=0;j<A2.size();++j){
                    A2[j] = proj_down(A[j], A[i], k);
                }
                shuffle(A2.begin(), A2.end(), Rng::get_engine());
                lp_debug(string(2*d, ' ') << i << "\n");
                vector<FLOAT> c2 = proj_down(c, A[i], k);
                vector<Barrier_Int<FLOAT>> x2 = solve(A2, c2, d-1);
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
    vector<Barrier_Int<FLOAT>> solve(vector<vector<FLOAT> > const &A, vector<FLOAT> const&c){
        assert(A.empty() || A[0].size() == c.size()+1);
        return solve(A, c, c.size());
    }
    /**
     *  Maximize c^T x
     *  subject to Ax <= b
     *
     *  Returns empty vector if infeasible
     */
    vector<Barrier_Int<FLOAT>> solve(vector<vector<FLOAT> > A, vector<FLOAT> const&b, vector<FLOAT> const&c){
        assert(A.size() == b.size());
        for(unsigned int i=0;i<A.size();++i){
            A[i].push_back(-b[i]);
        }
        return solve(A, c);
    }
};

#endif // LP_SEIDEL_BARIERLESS_HPP
