
// Linear programming in low dimension
// Uses Bigintegers for exact calculations

#include <random>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <bitset>

namespace lowdim_lp{
#define assert_msg(x,y) assert((x && y))
//#define lp_debug(x) do{cerr << x; }while(0)
#define lp_debug(x) do{}while(0)
//#define lp_debug_2(x) do{cerr << x; }while(0)
#define lp_debug_2(x) do{}while(0)


class Rng{
private:
    static std::mt19937 engine;

public:
    static std::mt19937& get_engine(){
        return engine;
    }
    template<typename T>
    static void set_seed(T const& seed){
        engine = std::mt19937(seed);
    }
    static void timebased_seed(){
        engine = std::mt19937(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count());
    }
    template<typename T>
    static typename std::enable_if<std::is_integral<T>::value, T>::type uniform(T l, T r){
        return std::uniform_int_distribution<T>(l, r)(engine);
    }
    template<typename T>
    static typename std::enable_if<std::is_floating_point<T>::value, T>::type uniform(T l, T r){
        return std::uniform_real_distribution<T>(l, r)(engine);
    }

};
std::mt19937 Rng::engine(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count());


class Timer{
public:
    template<typename S, typename T, typename... U>
    static S execute_timed(T func, std::string const& name, U&&... u){
        auto time_begin = std::chrono::high_resolution_clock::now();
        S ret = func(std::forward<U>(u)...);
        auto time_end = std::chrono::high_resolution_clock::now();
        auto timespan = std::chrono::duration_cast<std::chrono::nanoseconds>(time_end - time_begin);
        std::cerr << "Execution of '" << name << "' took " << std::fixed << std::setprecision(6) << timespan.count()*1e-9 << "\n";
        return ret;
    }
};


/**
 *  signed Bigint in base 10^18, used for Input / Output
 *  don't use for computations, doesn't support most operations
 *
 */
class Bigint_base10{
public:
    static constexpr int64_t BASE = 1e18;
    static constexpr int DIGITS = 18;
private:
    bool is_neg;
    std::vector<int64_t> data;

public:
    Bigint_base10():is_neg(false), data(1, 0){}
    Bigint_base10(int64_t const&val):is_neg(val<0){
        int64_t abs_val = abs(val);
        if(abs_val < BASE){
            data = {abs_val};
        } else {
            data = {abs_val%BASE, abs_val/BASE};
        }
    }

    Bigint_base10 operator+(Bigint_base10 const&o)const{
        assert_msg(is_neg == o.is_neg, "Addition operands need to have equal sign");
        Bigint_base10 ret;
        ret.is_neg = is_neg;
        ret.data.assign(1+std::max(data.size(), o.data.size()), 0);
        std::copy(data.begin(), data.end(), ret.data.begin());
        int64_t carry = 0;
        for(unsigned int i=0;i<o.data.size();++i){
            ret.data[i]+=o.data[i] + carry;
            carry = 0;
            if(ret.data[i] >= BASE){
                carry = 1;
                ret.data[i]-=BASE;
            }
        }
        for(unsigned int i=o.data.size();carry;++i){
            ret.data[i]+=carry;
            carry = 0;
            if(ret.data[i] >= BASE){
                carry = 1;
                ret.data[i]-=BASE;
            }
        }
        return ret.trim();
    }

    Bigint_base10 operator*(int64_t const&o)const{
        if(o == 0){
            return Bigint_base10(0);
        }
        if(o<0){
            return operator*(-o).negate();
        }
        if(o&1){
            return operator+(operator*(o-1));
        }
        return operator+(*this)*(o/2);
    }
    Bigint_base10& operator+=(Bigint_base10 const&o){
        *this = operator+(o);
        return *this;
    }
    Bigint_base10& operator*=(int64_t const&o){
        *this = operator*(o);
        return *this;
    }

    Bigint_base10& trim(){
        while(data.size()>1 && data.back() == 0){
            data.pop_back();
        }
        return *this;
    }

    bool is_zero()const{
        for(auto const&e:data) if(e) return false;
        return true;
    }

    Bigint_base10& negate(){
        is_neg = !is_neg;
        if(is_zero()) is_neg = false;
        return *this;
    }

    friend std::ostream& operator<<(std::ostream&o, Bigint_base10 const&b){
        if(b.is_neg) o << '-';
        o << b.data.back();
        o << std::setfill('0');
        for(auto it = next(b.data.rbegin());it != b.data.rend();++it){
            o << std::setw(9) << *it;
        }
        o << std::setw(0) << std::setfill(' ');
        return o;
    }
    friend std::istream& operator>>(std::istream&in, Bigint_base10 &b){
        static std::string tmp;
        in >> tmp;
        assert_msg(in, "input should be readable as a string");
        if(tmp[0] == '-'){
            b.is_neg = true;
            tmp = tmp.substr(1, -1);
        } else {
            b.is_neg = false;
        }
        assert_msg(std::all_of(tmp.begin(), tmp.end(), [](char const&c){return '0'<=c && c<='9';}), "Input should consist of digits and possibly a '-'");
        assert_msg(!tmp.empty(), "Input should contain at least one digit");

        b.data.resize((tmp.size()+DIGITS-1)/DIGITS);
        unsigned int i, j;
        for(i=tmp.size()-DIGITS, j=0;i>0;i-=DIGITS, ++j){
            b.data[j] = stoll(tmp.substr(i, DIGITS));
        }
        b.data[j] = stoll(tmp.substr(0, i+DIGITS));
        return in;
    }
};

/**
 *  Biginteger with fixed precision
 *  Has 32*len bits, is signed
 *
 *  Speed is ok.
 */
template<size_t len>
struct Bigint_Fixedsize_Fast{
    unsigned int data[len];
    uint16_t siz;
    bool sign;
    static constexpr unsigned int bits = 32;
    Bigint_Fixedsize_Fast(){
        data[0] = 0;
        siz = 1;
        sign = false;
    }
    Bigint_Fixedsize_Fast(long long a){
        sign = false;
        if(a<0){
            sign = true;
            a=-a;
        }
        siz = 0;
        do{
            long long b = a>>bits;
            data[siz] = a - (b<<bits);
            a = b;
            ++siz;
        } while(a);
    }
    void trim(){
        while(siz>1 && !data[siz-1]) --siz;
    }
    int comp_unsigned(Bigint_Fixedsize_Fast const&o)const{
        uint16_t lim = std::min(siz, o.siz);
        for(unsigned int i=lim;i<siz;++i){
            if(data[i]){
                return 1;
            }
        }
        for(unsigned int i=lim;i<o.siz;++i){
            if(o.data[i]){
                return -1;
            }
        }
        for(unsigned int i=lim-1;i+1;--i){
            if(data[i]!=o.data[i]){
                return data[i] < o.data[i]?-1:1;
            }
        }
        return 0;
    }
    int comp(Bigint_Fixedsize_Fast const&o)const{
        int sign_mul = 1-2*sign;
        if(sign != o.sign){
            return sign_mul;
        }
        return sign_mul * comp_unsigned(o);
    }
    bool operator<(Bigint_Fixedsize_Fast const&o)const{
        return comp(o)<0;
    }
    bool operator>(Bigint_Fixedsize_Fast const&o)const{
        return comp(o)>0;
    }
    bool operator<=(Bigint_Fixedsize_Fast const&o)const{
        return comp(o)<=0;
    }
    bool operator>=(Bigint_Fixedsize_Fast const&o)const{
        return comp(o)>=0;
    }
    bool operator==(Bigint_Fixedsize_Fast const&o)const{
        return comp(o)==0;
    }
    bool operator!=(Bigint_Fixedsize_Fast const&o)const{
        return comp(o)!=0;
    }
    bool operator!()const{
        return operator==(ZERO);
    }
    Bigint_Fixedsize_Fast operator-()const{
        Bigint_Fixedsize_Fast ret(*this);
        if(!!ret){
            ret.sign ^=1;
        }
        return ret;
    }
    Bigint_Fixedsize_Fast operator*(Bigint_Fixedsize_Fast const&o)const{
        Bigint_Fixedsize_Fast ret;
        ret.siz = std::min(siz+o.siz, (int)len);
        ret.sign = (sign!=o.sign);
        std::fill(ret.data, ret.data+ret.siz, 0);
        for(unsigned int i=0;i<siz;++i){
            unsigned long long carry = 0, carry_2;
            for(unsigned int j=0;j<o.siz;++j){
                carry+= data[i]*(unsigned long long)o.data[j] + ret.data[i+j];
                carry_2 = carry >> bits;
                ret.data[i+j] = carry - (carry_2<<bits);
                carry = carry_2;
            }
            for(unsigned int j=i+o.siz;carry;++j){
                carry+= ret.data[j];
                carry_2 = carry >> bits;
                ret.data[j] = carry - (carry_2<<bits);
                carry = carry_2;
            }
        }
        ret.trim();
        return ret;
    }
    Bigint_Fixedsize_Fast& operator*=(Bigint_Fixedsize_Fast const&o){
        *this = operator*(o);
        return *this;
    }
    static void unsigned_add(Bigint_Fixedsize_Fast &ret, Bigint_Fixedsize_Fast const&A, Bigint_Fixedsize_Fast const&B){
        const Bigint_Fixedsize_Fast *a = &A, *b = &B;
        if(a->siz < b->siz) std::swap(a, b);
        ret.sign = A.sign;
        unsigned long long carry = 0, carry_2;
        unsigned int j;
        for(j=0;j<b->siz;++j){
            carry+=(unsigned long long)a->data[j] + (unsigned long long)b->data[j];
            carry_2 = carry>>bits;
            ret.data[j] = carry - (carry_2<<bits);
            carry = carry_2;
        }
        for(;j<a->siz;++j){
            carry+=a->data[j];
            carry_2 = carry>>bits;
            ret.data[j] = carry - (carry_2<<bits);
            carry = carry_2;
        }
        if(carry){
            ret.data[j++] = carry;
        }
        ret.siz = j;
        ret.trim();
    }
    static void unsigned_subtract(Bigint_Fixedsize_Fast &ret, Bigint_Fixedsize_Fast const&A, Bigint_Fixedsize_Fast const&B){
        int com = A.comp_unsigned(B);
        if(com == 0){
            ret.sign = false;
            ret.siz = 1;
            ret.data[0] = 0;
            return;
        }
        ret.sign = A.sign;
        const Bigint_Fixedsize_Fast *a = &A, *b = &B;
        if(com < 0){
            ret.sign ^= true;
            std::swap(a, b);
        }
        // deal with case then o is not trimed.
        unsigned int min_siz = std::min(A.siz, B.siz);
        unsigned long long carry = 0, carry_2;
        unsigned int j;
        for(j=0;j<min_siz;++j){
            carry+=(unsigned long long)a->data[j] - (unsigned long long)b->data[j];
            carry_2 = carry>>bits;
            ret.data[j] = carry - (carry_2<<bits);
            carry = -(carry_2 & 1u);
        }
        for(;carry;++j){
            assert(j < a->siz);
            carry+=a->data[j];
            carry_2 = carry>>bits;
            ret.data[j] = carry - (carry_2<<bits);
            carry = -(carry_2 & 1u);
        }
        std::copy(a->data+j, a->data+a->siz, ret.data+j);
        ret.siz = a->siz;
        ret.trim();
    }
    static void add(Bigint_Fixedsize_Fast &ret, Bigint_Fixedsize_Fast const&A, Bigint_Fixedsize_Fast const&B){
        if(A.sign == B.sign){
            unsigned_add(ret, A, B);
        } else {
            unsigned_subtract(ret, A, B);
        }
    }
    static void sub(Bigint_Fixedsize_Fast &ret, Bigint_Fixedsize_Fast const&A, Bigint_Fixedsize_Fast const&B){
        if(A.sign != B.sign){
            unsigned_add(ret, A, B);
        } else {
            unsigned_subtract(ret, A, B);
        }
    }
    Bigint_Fixedsize_Fast operator+(Bigint_Fixedsize_Fast const&o)const{
        Bigint_Fixedsize_Fast ret;
        add(ret, *this, o);
        return ret;
    }
    Bigint_Fixedsize_Fast& operator+=(Bigint_Fixedsize_Fast const&o){
        add(*this, *this, o);
        return *this;
    }
    Bigint_Fixedsize_Fast operator-(Bigint_Fixedsize_Fast const&o)const{
        Bigint_Fixedsize_Fast ret;
        sub(ret, *this, o);
        return ret;
    }
    Bigint_Fixedsize_Fast operator-=(Bigint_Fixedsize_Fast const&o){
        sub(*this, *this, o);
        return *this;
    }
    /// TODO: more operators
    void print_binary(std::ostream&o, Bigint_Fixedsize_Fast const&b){
        o << "[";
        o << sign << "/" << len << "/";
        for(size_t i=siz;i>0;--i){
            o << std::bitset<bits>(b.data[i-1]);
        }
        o << "]";
    }
    friend std::ostream& operator<<(std::ostream&o, Bigint_Fixedsize_Fast const&b){
        if(b.sign){
            return o << '-' << -b;
        }
        Bigint_base10 ret(0);
        int64_t base = 1ll<<bits;
        for(int i = b.siz-1;i>=0;--i){
            ret*=base;
            ret+=Bigint_base10(b.data[i]);
        }
        o << ret;

        return o;
    }
    explicit operator long double()const{
        if(sign){
            return -(long double)operator-();
        }
        long double ret = 0.0;
        long double base = 1ll<<bits;
        for(int i = siz-1;i>=0;--i){
            ret = ret * base + data[i];
        }
        return ret;
    }
    /// TODO: implement for larger inputs
    friend std::istream& operator>>(std::istream&i, Bigint_Fixedsize_Fast &b){
        int64_t tmp;
        i >> tmp;
        b = Bigint_Fixedsize_Fast(tmp);
        return i;
    }
    static const Bigint_Fixedsize_Fast ZERO;
};
template<size_t len>
const Bigint_Fixedsize_Fast<len> Bigint_Fixedsize_Fast<len>::ZERO(0);

/**
 *  Randomized LP in expected
 *  O(d! 4^d n)
 *  10x slower that 2-phase Clarkson
 *  Used a backend for small LPs in Clarkson
 *  Does exact calculations.
 */
template<typename FLOAT>
class Lp_Seidel{
private:


    // orthogonal projection of 'vec' into 'plane'
    std::vector<FLOAT> proj_down(std::vector<FLOAT> const&vec, std::vector<FLOAT> const&plane, size_t j){
        assert(vec.size() <= plane.size() && plane.size()<=vec.size()+1);
        assert(j+1 < plane.size());
        assert(!!plane[j]);
        std::vector<FLOAT> ret (vec.size()-1);
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
    std::vector<FLOAT> proj_up(std::vector<FLOAT> const&vec, std::vector<FLOAT> const&plane, size_t j){
        assert(vec.size()+1 == plane.size());
        assert(j+1 < plane.size());
        assert(!!plane[j]);
        std::vector<FLOAT> ret(vec.size()+1);
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
    FLOAT planescal(std::vector<FLOAT> const&x, std::vector<FLOAT> const&a){
        assert(x.size() == a.size());
        FLOAT ret=0;
        for(size_t i=0;i<x.size();++i){
            ret+=x[i]*a[i];
        }
        return ret;
    }

    // solve lp recursively
    std::vector<FLOAT> solve(std::vector<std::vector<FLOAT> > const &A, std::vector<FLOAT> const&c, int d, FLOAT const& barier_loc){
        int n=A.size();
        if(d==1){ // base case: single dimension
            std::vector<FLOAT> ret(2);
            ret[0] = (c[0]<FLOAT(0) ? -barier_loc : barier_loc);
            ret[1] = 1ull;
            for(int i=0;i<n;++i){
                if(ret[0]*A[i][0]+ret[1]*A[i].back()>FLOAT(0)){
                    if(!A[i][0]){
                        lp_debug("infeasible single\n");
                        return std::vector<FLOAT>();
                    }
                    ret[0] = -A[i].back();
                    ret[1] = A[i][0];
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
                    return std::vector<FLOAT>();
                }
            }
            return ret;
        }
        FLOAT barier_next = barier_loc * barier_loc;
        // initial solution
        std::vector<FLOAT> x(d+1);
        for(int i=0;i<d;++i){
            x[i] = (c[i]<FLOAT(0)?-barier_loc:barier_loc);
        }
        x.back() = FLOAT(1);
        for(size_t i=0;i<A.size();++i){
            if(planescal(x, A[i])>FLOAT(0)){
                int k = 0;
                while(k<d && !A[i][k]) ++k;
                // recurse
                if(k==d) {lp_debug("what?\n"); return std::vector<FLOAT>();} // degenerate failing plane??????
                std::vector<std::vector<FLOAT> > A2(i);
                for(size_t j=0;j<A2.size();++j){
                    A2[j] = proj_down(A[j], A[i], k);
                }
                shuffle(A2.begin(), A2.end(), Rng::get_engine());
                lp_debug(string(2*d, ' ') << i << "\n");
                std::vector<FLOAT> c2 = proj_down(c, A[i], k);
                std::vector<FLOAT> x2 = solve(A2, c2, d-1, barier_next);
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
    std::vector<FLOAT> solve(std::vector<std::vector<FLOAT> > const &A, std::vector<FLOAT> const&c, const FLOAT& barier = FLOAT((long long)1e9)){
        assert(A.empty() || A[0].size() == c.size()+1);
        return solve(A, c, c.size(), barier);
    }
    /**
     *  Maximize c^T x
     *  subject to Ax <= b
     *
     *  Returns empty vector if infeasible
     */
    std::vector<FLOAT> solve(std::vector<std::vector<FLOAT> > A, std::vector<FLOAT> const&b, std::vector<FLOAT> const&c){
        assert(A.size() == b.size());
        for(unsigned int i=0;i<A.size();++i){
            A[i].push_back(-b[i]);
        }
        return solve(A, c);
    }
};

/**
 *  Randomized LP, runtime is expected
 *  O(d^2 n + d^4 sqrt{n} log n + d^{O(d)} log(n)) | two phase
 *  O(d n log n + d^{O(d)} log n)                  | one phase
 *  Quite fast for low d.
 *  Does exact calculations.
 */
template<typename Big_Int, bool use_two_phase = true>
class Lp_Clarkson{
private:
    /**
     *  Returns a sub-multiset of size siz uniformly at random
     *  out of the set where i is present weight[i] times.
     *
     *  Runs in O(|weight| + siz^2) expected time.
     *  Could be optimized
     */
    std::vector<int> sample_subset(std::vector<int64_t> const&weight, unsigned int siz){
        int64_t total_weight = accumulate(weight.begin(), weight.end(), 0ll);
        std::vector<int64_t> samples;
        while(samples.size() < siz){
            int64_t new_sample = Rng::uniform<int64_t>(0, total_weight-1);
            if(find(samples.begin(), samples.end(), new_sample) == samples.end()){
                samples.push_back(new_sample);
            }
        }
        sort(samples.begin(), samples.end());
        std::vector<int> ret;
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
    bool is_satisfied(std::vector<Big_Int> const&x, std::vector<Big_Int> const&a){
        assert(x.size() == a.size());
        Big_Int ret=0;
        for(size_t i=0;i<x.size();++i){
            ret+=x[i]*a[i];
        }
        return ret <= Big_Int(0);
    }
    std::vector<Big_Int> solve_two(std::vector<std::vector<Big_Int> > const&A, std::vector<Big_Int> const&c){
        const unsigned int sample_size = c.size()*c.size()*4;
        Lp_Seidel<Big_Int> sub_lp;
        // to few constrains -> use other solver
        if(A.size() < sample_size){
            return sub_lp.solve(A, c);
        } else {
            int constraints = A.size();
            int variables = c.size();
            std::vector<int64_t> weight(constraints, 1);
            std::vector<Big_Int> x;
            std::vector<std::vector<Big_Int> > subproblem_A;
            std::vector<char> is_violated(constraints, 0);
            for(unsigned int iteration=1;;++iteration){
                subproblem_A.clear();
                std::vector<int> subspace = sample_subset(weight, sample_size);
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
                    lp_debug_2(cerr << "Iterations: " << iteration);
                    lp_debug_2(cerr << ", max weight: " << *max_element(weight.begin(), weight.end()) << "\n");
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
    std::vector<Big_Int> solve_one(std::vector<std::vector<Big_Int> > const&A, std::vector<Big_Int> const&c){
        const unsigned int constraints = A.size(), variables = c.size();
        if(constraints <= variables*variables*6){
            return solve_two(A, c);
        } else {
            const unsigned int sqrt_constraints = sqrt(constraints);
            const unsigned int sample_size = variables * sqrt(constraints);
            std::vector<Big_Int> x;
            std::vector<std::vector<Big_Int> > subproblem_A;
            std::vector<int> violations;
            for(unsigned int iteration=1;;++iteration){
                std::vector<int> subspace = sample_subset(std::vector<int64_t>(constraints, 1), sample_size);
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
                lp_debug_2("Violations: " << violations.size() << " / " << 2*sqrt_constraints << "\n");
                if(violations.empty()){
                    lp_debug_2(cerr << "Iterations: " << iteration);
                    lp_debug_2(cerr << ", used constraints:" << subproblem_A.size() << "\n");
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
    std::vector<Big_Int> solve(std::vector<std::vector<Big_Int> > const&A, std::vector<Big_Int> const&c){
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
    std::vector<Big_Int> solve(std::vector<std::vector<Big_Int> > A, std::vector<Big_Int> const&b, std::vector<Big_Int> const&c){
        assert(A.size() == b.size());
        for(unsigned int i=0;i<A.size();++i){
            A[i].push_back(-b[i]);
        }
        return solve(A, c);
    }
};

}
template<size_t len>
using Big_Int = lowdim_lp::Bigint_Fixedsize_Fast<len>;
template<size_t len>
using Lp_Solver = lowdim_lp::Lp_Clarkson<Big_Int<len>>;

#include <bits/stdc++.h>
using namespace std;
using lowdim_lp::Rng;
using lowdim_lp::Timer;
