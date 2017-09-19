// linear programming in low dimension
// runs in O(d! * n)
// optimized for 32 bit architectures (i.e. codeforces)
// takes ~2 seconds for d=4, n=20'000

// more speed:
//#undef assert
//#define assert(x) do{}while(0)

struct Bignum{
    static const unsigned int base;
    static const unsigned long long half_base;
    vector<unsigned int> data;
    bool sign;
    Bignum():data(1, 0), sign(false){}
    Bignum(long long const&a):data(1+(abs(a)>=base)), sign(a<0ll){
        unsigned long long b = abs(a);
        if(b>=base){
            data[1] = b/base;
            data[0] = b%base;
        } else {
            data[0] = b;
        }
    }
    void normalize(){
        while(data.size()>1 && data.back()==0) data.pop_back();
        if(data.size()==1 && data.back()==0) sign = false;
    }
    void abs_add(Bignum const&o){
        if(o.data.size()>=data.size()) data.resize(o.data.size()+1, 0ull);
        else data.push_back(0);
        unsigned int carry=0;
        for(unsigned int i=0;i<o.data.size();++i){
            data[i]+=o.data[i]+carry;
            carry = 0;
            if(data[i]>=base){
                data[i]-=base;
                carry = 1;
            }
        }
        for(unsigned int i=o.data.size();carry;++i){
            data[i]+=carry;
            carry=0;
            if(data[i]>=base){
                data[i]-=base;
                ++carry;
            }
        }
        normalize();
    }
    void abs_sub(Bignum const&o){
        unsigned int carry=0;
        for(unsigned int i=0;i<o.data.size();++i){
            data[i]-=o.data[i]+carry;
            carry = 0;
            if(data[i]&(1ull<<31)){
                data[i]+=base;
                carry = 1;
            }
        }
        for(unsigned int i=o.data.size();carry;++i){
            data[i]-=carry;
            carry=0;
            if(data[i]&(1ull<<31)){
                data[i]+=base;
                carry = 1;
            }
        }
        normalize();
    }
    void abs_isub(Bignum const&o){
        data.resize(o.data.size(), 0ull);
        unsigned int carry=0;
        for(unsigned int i=0;i<o.data.size();++i){
            data[i] = o.data[i]-data[i]-carry;
            carry = 0;
            if(data[i]&(1ull<<31)){
                data[i]+=base;
                carry = 1;
            }
        }
        assert(!carry);
        normalize();
    }
    bool abs_lt(Bignum const&o)const{
        if(o.data.size()!=data.size()) return data.size()<o.data.size();
        for(int i = data.size()-1;i>=0;--i){
            if(data[i]!=o.data[i]) return data[i]<o.data[i];
        }
        return false;
    }
    Bignum& operator+=(Bignum const&o){
        if(sign == o.sign) abs_add(o);
        else if(abs_lt(o)){
            sign = !sign;
            abs_isub(o);
        } else {
            abs_sub(o);
        }
        return *this;
    }
    Bignum& operator-=(Bignum const&o){
        if(sign == !o.sign) abs_add(o);
        else if(abs_lt(o)){
            sign = !sign;
            abs_isub(o);
        } else {
            abs_sub(o);
        }
        return *this;
    }
    Bignum operator+(Bignum const&o)const{
        if(data.size()<o.data.size()) return o+*this;
        Bignum ret(*this);
        ret+=o;
        return ret;
    }
    Bignum operator-(Bignum const&o)const{
        if(data.size()<o.data.size()){
            Bignum ret = o-*this;
            ret.sign = !ret.sign;
            return ret;
        }
        Bignum ret(*this);
        ret-=o;
        return ret;
    }
    void multo(Bignum const&o, Bignum&ret)const{
        ret.sign = (sign!=o.sign);
        static vector<unsigned long long> tmp;
        tmp.resize(data.size()+o.data.size()+1);
        fill(tmp.begin(), tmp.end(), 0ull);
        for(size_t i=0;i<data.size();++i){
            for(size_t j=0;j<o.data.size();++j){
                tmp[i+j]+=data[i]*(unsigned long long)o.data[j];
                if(tmp[i+j]>half_base){
                    ++tmp[i+j+1];
                    tmp[i+j]-=half_base;
                }
            }
        }
        ret.data.resize(data.size()+o.data.size());
        for(size_t i=0;i<ret.data.size();++i){
            unsigned int low = tmp[i]%base;
            unsigned long long high = tmp[i]/base;
            ret.data[i]=low;
            tmp[i+1]+=high;
        }
        assert(ret.data.back()<base);
        ret.normalize();
    }
    Bignum operator*(Bignum const&o)const{
        Bignum ret;
        multo(o, ret);
        return ret;
    }
    Bignum operator*=(Bignum const&o){
        multo(o, *this);
        return *this;
    }
    Bignum operator-()const{
        Bignum ret = *this;
        ret.sign = !ret.sign;
        return ret;
    }
    Bignum& negate(){
        if(!!*this) sign = !sign;
        return *this;
    }
    bool operator<(Bignum const&o)const{
        if(sign!=o.sign) return sign;
        return sign!=abs_lt(o);
    }
    bool operator>(Bignum const&o)const{
        return o<*this;
    }
    bool operator==(Bignum const&o)const{
        if(sign!=o.sign) return false;
        if(data.size()!=o.data.size()) return false;
        for(size_t i=0;i<data.size();++i){
            if(data[i]!=o.data[i]) return false;
        }
        return true;
    }
    bool operator!=(Bignum const&o)const{
        return !operator==(o);
    }
    friend ostream& operator<<(ostream&o, Bignum const&b){
        if(b.sign) o << "-";
        o << b.data.back() << setfill('0') ;
        for(int i=(int)b.data.size()-2;i>=0;--i){
            o << setw(18) << b.data[i];
        }
        return o<<setw(0);
    }
    bool operator!()const{
        return data.size()==1 && data[0]==0;
    }
    bool is_negative()const{
        return sign;
    }
};
unsigned int const Bignum::base = 1e9;
unsigned long long const Bignum::half_base = 1e18;

struct Fsbignum{
    static constexpr size_t len = 8;
    unsigned int data[len];
    static constexpr unsigned int bits = 31;
    Fsbignum(){memset(data, 0, sizeof(data));}
    Fsbignum(long long const&_a){
        memset(data, 0, sizeof(data));
        unsigned long long a = _a;
        data[0] = a&((1u<<bits)-1);
        data[1] = a>>bits;
        data[1]&=~(1u<<bits);
        if(a>~a){ // negative number, use complement
            for(size_t i=2;i<len;++i){
                data[i] = (1u<<bits)-1;
            }
        }
    }
    __attribute__((optimize("unroll-loops")))
    Fsbignum& operator+=(Fsbignum const&o){
        unsigned int carry = 0;
        for(size_t i=0;i<len;++i){
            data[i]+=o.data[i] + carry;
            carry = data[i]>>bits;
            data[i]&=~(1<<bits);
        }
        return *this;
    }
    __attribute__((optimize("unroll-loops")))
    Fsbignum& operator-=(Fsbignum const&o){
        unsigned int carry = 0;
        for(size_t i=0;i<len;++i){
            data[i]-=o.data[i] + carry;
            carry = data[i]>>bits;
            data[i]&=~(1<<bits);
        }
        return *this;
    }
    Fsbignum operator+(Fsbignum const&o)const{
        Fsbignum ret(*this);
        ret+=o;
        return ret;
    }
    Fsbignum operator-(Fsbignum const&o)const{
        Fsbignum ret(*this);
        ret-=o;
        return ret;
    }
    __attribute__((optimize("unroll-loops")))
    void multo(Fsbignum const&o, Fsbignum &ret)const{
        static unsigned int tmp[len+1];
        memset(tmp, 0, sizeof(tmp));
        for(size_t i=0;i<len;++i){
            unsigned long long val = 0;
            for(size_t j=0;j<len-i;++j){
                val+= data[i]*(unsigned long long)o.data[j]+tmp[i+j];
                tmp[i+j] = val&((1u<<bits)-1);
                val>>=bits;
            }
        }
        memcpy(ret.data, tmp, sizeof(ret.data));
    }
    Fsbignum& operator*=(Fsbignum const&o){
        multo(o, *this);
        return *this;
    }
    Fsbignum operator*(Fsbignum const&o)const{
        Fsbignum ret;
        multo(o, ret);
        return ret;
    }
    Fsbignum& negate(){
        unsigned int carry = 0;
        for(size_t i=0;i<len;++i){
            data[i] = ~data[i] + !carry;
            carry = (data[i]>>bits);
            data[i]&=~(1<<bits);
        }
        return *this;
    }

    Fsbignum operator-()const{
        Fsbignum ret(*this);
        ret.negate();
        return ret;
    }
    bool operator<(Fsbignum const&o)const{
        if(data[len-1]>>(bits-1)!=o.data[len-1]>>(bits-1)){
            return data[len-1]>>(bits-1);
        }
        for(size_t i=len-1;~i;--i){
            if(data[i]!=o.data[i]) return data[i]<o.data[i];
        }
        return false;
    }
    bool operator>(Fsbignum const&o)const{
        return o<*this;
    }
    bool operator==(Fsbignum const&o)const{
        for(size_t i=0;i<len;++i){
            if(data[i] !=o.data[i]) return false;
        }
        return true;
    }
    bool operator!=(Fsbignum const&o)const{
        return !operator==(o);
    }
    bool operator!()const{
        for(size_t i=0;i<len;++i){
            if(data[i]) return false;
        }
        return true;
    }
    bool is_negative()const{
        return data[len-1]>>(bits-1);
    }
    friend ostream& operator<<(ostream&o, Fsbignum const&b){
        o << "[";
        for(size_t i=len;i>0;--i){
            o << bitset<bits>(b.data[i-1]);
        }
        o << "]";
        return o;
    }
};

// adjust barier and Fsbignum::len for bigger inputs
// 1e12, len=8 worked fine for enclosing circle with coordinates <= 1e4
using FLOAT = Bignum; // adaptive length
//using FLOAT = Fsbignum; // fixed size, way faster
const FLOAT barier = FLOAT((long long)1e9);
const FLOAT barier2 = barier*barier;
const FLOAT eps = FLOAT(0);
const bool square_barier = false; // squares the barrier value at each recursion step

//#define lp_debug(x) do{cerr << x; }while(0)
#define lp_debug(x) do{}while(0)

// orthogonal projection of 'vec' into 'plane'
vector<FLOAT> proj_down(vector<FLOAT> const&vec, vector<FLOAT> const&plane, size_t j){
    assert(vec.size() <= plane.size() && plane.size()<=vec.size()+1);
    assert(j+1 < plane.size());
    assert(!!plane[j]);
    vector<FLOAT> ret (vec.size()-1);
    FLOAT tmp;
    if(plane[j].is_negative()){
        for(size_t i=0;i+1<vec.size();++i){
            vec[j].multo(plane[i+(i>=j)], ret[i]);
            vec[i+(i>=j)].multo(plane[j], tmp);
            ret[i] -=  tmp;
        }
    } else {
        for(size_t i=0;i+1<vec.size();++i){
            vec[j].multo(plane[i+(i>=j)], tmp);
            vec[i+(i>=j)].multo(plane[j], ret[i]);
            ret[i] -=  tmp;
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
    if(denom.is_negative()){
        denom.negate();
    }
    for(size_t i=0;i<vec.size();++i){
        ret[i+(i>=j)]*=denom;
    }
    if(!plane[j].is_negative()){
        ret[j].negate();
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


//mt19937 rng(4679784);
mt19937 rng(std::chrono::duration_cast<std::chrono::nanoseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count());
vector<FLOAT> solve(vector<vector<FLOAT> > const &A, vector<FLOAT> const&c, int d, FLOAT const& barier_loc){
    int n=A.size();
    if(d==1){ // base case: single dimension
        vector<FLOAT> ret(2);
        if(square_barier){
            ret[0] = (c[0]<FLOAT(0) ? -barier_loc : barier_loc);
        } else {
            ret[0] = (c[0]<FLOAT(0) ? -barier2 : barier2);
        }
        ret[1] = 1ull;
        for(int i=0;i<n;++i){
            if(ret[0]*A[i][0]+ret[1]*A[i].back()>eps){
                if(!A[i][0]) return vector<FLOAT>();
                ret[0] = -A[i].back();
                ret[1] = A[i][0];
                if(ret[1].is_negative()){
                    ret[1].negate();
                    ret[0].negate();
                }
                lp_debug(" -> " << i << "\n");
            }
        }
        for(int i=0;i<n;++i){
            if(ret[0]*A[i][0]+ret[1]*A[i].back()>eps){
                lp_debug("unfeasible\n");
                return vector<FLOAT>();
            }
        }
        return ret;
    }
    FLOAT barier_next = (square_barier ? barier_loc * barier_loc: barier_loc);
    // initial solution
    vector<FLOAT> x(d+1);
    for(int i=0;i<d;++i){
        x[i] = (c[i]<FLOAT(0)?-barier:barier);
    }
    x.back() = FLOAT(1);
    for(size_t i=0;i<A.size();++i){
        if(planescal(x, A[i])>eps){
            int k = 0;
            while(k<d && !A[i][k]) ++k;
            // recurse
            if(k==d) {lp_debug("what?\n"); return vector<FLOAT>();} // degenerate failing plane??????
            vector<vector<FLOAT> > A2(i);
            for(size_t j=0;j<A2.size();++j){
                A2[j] = proj_down(A[j], A[i], k);
            }
            shuffle(A2.begin(), A2.end(), rng);
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

vector<FLOAT> solve(vector<vector<FLOAT> > const &A, vector<FLOAT> const&c){
	return solve(A, c, c.size(), barier);
}
