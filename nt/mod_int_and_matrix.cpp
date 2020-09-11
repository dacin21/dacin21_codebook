template<typename traits>
class Mod_Int{
public:
    using int_t = typename traits::int_t;
    using long_t = typename traits::long_t;
    static constexpr int_t mod(){ return traits::get_mod(); };

    struct Summer{
    public:
        static constexpr long_t modmod(){ return traits::get_mod()*(long_t)traits::get_mod(); };
        static long_t modmod_step(long_t const& val){
            return val >= modmod() ? val-modmod() : val;
        }

        Summer() : val{} {}
        explicit Summer(Mod_Int const&o) : val(o.get_value()){}
        operator Mod_Int() const {
            return Mod_Int(mod_full(val));
        }

        Summer operator-(Summer const&o){
            return Summer(modmod() - modmod_step(o.val));
        }
        Summer& operator+=(Summer const&o){
            val = modmod_step(val+o.val);
            return *this;
        }
        Summer& operator-=(Summer const&o){
            return operator-=(-o);
        }
        Summer& addmul(Mod_Int const&a, Mod_Int const&b){
            val = modmod_step(val + a.get_value()*(long_t)b.get_value());
            return *this;
        }

    private:
        long_t val;
    };


    Mod_Int() : value(0) {}
    Mod_Int(int_t value_) : value(value_) {}

    friend ostream& operator<<(ostream&o, Mod_Int const&val){
        return o << val.value;
    }
    friend istream& operator>>(istream&i, Mod_Int &val){
         i >> val.value;
         val.value = mod_full(val.value);
         return i;
    }
    int_t const& get_value() const {
        return value;
    }

    Mod_Int operator+(Mod_Int const&o) const {
        Mod_Int ret(mod_step(value + o.value));
        return ret;
    }
    Mod_Int& operator+=(Mod_Int const&o){
        return *this = *this + o;
    }
    Mod_Int& operator++(){
        return operator+=(int_t{1});
    }
    Mod_Int operator-() const{
        Mod_Int ret(mod_step(mod() - value));
        return ret;
    }
    Mod_Int operator-(Mod_Int const&o) const {
        return operator+(-o);
    }
    Mod_Int& operator-=(Mod_Int const&o){
        return operator+=(-o);
    }
    Mod_Int& operator--(){
        return operator-=(int_t{1});
    }
    Mod_Int operator*(Mod_Int const&o) const {
        Mod_Int ret(mod_full(value * static_cast<long_t>(o.value)));
        return ret;
    }
    Mod_Int& operator*=(Mod_Int const&o){
        return *this = *this * o;
    }
    Mod_Int inv() const {
        return inv_impl(value);
    }
    Mod_Int operator/(Mod_Int const&o) const {
        return operator*(o.inv());
    }
    Mod_Int& operator/=(Mod_Int const&o) {
        return *this = *this / o;
    }

    bool operator==(Mod_Int const&o) const {
        return value == o.value;
    }
    bool operator!=(Mod_Int const&o) const {
        return !(*this == o);
    }
    bool operator!() const {
        return !value;
    }

private:
    static int_t mod_step(int_t const& val){
        return val >= mod() ? val-mod() : val;
    }
    static int_t mod_full(long_t const&val){
        return mod() ? val%mod() : val;
    }
    static Mod_Int inv_impl(Mod_Int const& val){
        if(mod() == 0){
            assert(val*val == 1);
            return val;
        }
        int_t value = val.value;
        constexpr size_t cache_size = traits::inv_cache_size()+1;
        static_assert(cache_size > 1);
        static array<Mod_Int, cache_size> cache = [](){
            array<Mod_Int, cache_size> ret;
            ret[1] = 1;
            for(int_t i=2;i<cache_size && i < mod();++i){
                ret[i] = -ret[mod()%i] * (mod()/i);
            }
            return ret;
        } ();
        assert(value != 0);
        Mod_Int factor = 1;
        while(value >= cache_size){
            factor *= - Mod_Int(mod() / value);
            value = mod() % value;
        }
        assert(value != 0);
        assert(factor != 0  && value != 0);
        return factor * cache[value];
    }

    int_t value;
};

template<uint32_t mod>
struct fixed_mod{
    using int_t = uint32_t;
    using long_t = uint64_t;
    static_assert(mod != 0, "Negative numbers won't work.");
    static_assert(numeric_limits<int_t>::max()/2 >= mod, "Addition might overflow.");
    static constexpr size_t inv_cache_size(){ return 30000; }
    static constexpr int_t get_mod(){ return mod; }
};
#ifdef __SIZEOF_INT128__
template<uint64_t mod>
struct fixed_mod_long{
    using int_t = uint64_t;
    using long_t = __int128;
    static_assert(mod != 0, "Negative numbers won't work.");
    static_assert(numeric_limits<int_t>::max()/2 >= mod, "Addition might overflow.");
    static constexpr size_t inv_cache_size(){ return 30000; }
    static constexpr int_t get_mod(){ return mod; }
};
#endif // __SIZEOF_INT128__
struct no_mod{
    using int_t  = int64_t;
    using long_t = int_t;
    static constexpr size_t inv_cache_size(){ return 1; }
    static constexpr int_t get_mod(){ return 0; }
};
struct mutable_mod{
    using int_t = uint32_t;
    using long_t = uint64_t;
    static int_t mod;
    // can't use cache if mod is changing
    static constexpr size_t inv_cache_size(){ return 1; }
    static int_t get_mod(){ return mod; }
};
mutable_mod::int_t mutable_mod::mod = 1000000007;

template<typename T>
struct Matrix{
    Matrix() {}
    Matrix(int R_, int C_) : R(R_), C(C_), data(R*C) {}
    // column vector
    Matrix(std::vector<T> data_) : R(data_.size()), C(1), data(move(data_)) {}
    // matrix from vector of rows
    Matrix(std::vector<std::vector<T> > const&data_) : R(data_.size()), C(data_.empty() ? 0 : data_[0].size()), data(R*C) {
        for(int i=0;i<R;++i){
            copy(data_[i].begin(), data_[i].end(), data.begin()+i*C);
        }
    }

    static Matrix eye(int n, const T diag_val = 1){
        Matrix ret(n, n);
        for(int i=0;i<n;++i){
            ret.at(i, i) = diag_val;
        }
        return ret;
    }

    Matrix operator*(Matrix const&other) const {
        assert(C == other.R);
        Matrix ret(R, other.C);
        const Matrix o = other.transposed();
        for(int i=0;i<R;++i){
            for(int j=0;j<other.C;++j){
                typename T::Summer su{};
                for(int k=0;k<C;++k){
                    //ret.at(i, j) += at(i, k) * o.at(j, k);
                    su.addmul(data[i*C+k], o.data[j*o.C+k]);
                }
                ret.data[i*ret.C+j] = T(su);
            }
        }
        return ret;
    }

    Matrix transposed() const {
        Matrix ret(C, R);
        for(int i=0;i<C;++i){
            for(int j=0;j<R;++j){
                ret.at(j, i) = at(i, j);
            }
        }
        return ret;
    }

    T& at(int i, int j){
        assert(0 <= i && i < R);
        assert(0 <= j && j < C);
        return data[i*C+j];
    }
    T const& at(int i, int j) const {
        assert(0 <= i && i < R);
        assert(0 <= j && j < C);
        return data[i*C+j];
    }

    int R, C;
    vector<T> data;
};

// computes base^exp * start
template<typename T>
T powmul(T base, T start, uint64_t exp){
    for(;exp;exp>>=1){
        if(exp&1){
            start = base * start;
        }
        base = base * base;
    }
    return start;
}


template<typename T>
struct Slow_Poly : vector<T>{
    using Summer = typename T::Summer;

    Slow_Poly(){}
    Slow_Poly(vector<T> const&a) : vector<T>(a) {}

    Slow_Poly& operator*=(Slow_Poly const&o){
        static Slow_Poly tmp;
        tmp.assign(this->size()+o.size()-1, T{});
        for(int i=0; i<(int)tmp.size(); ++i){
            Summer su;
            for(int j = max(0, i-(int)o.size()+1), j_end = min((int)this->size(), i+1); j<j_end; ++j){
                su.addmul((*this)[j], o[i-j]);
            }
            tmp[i] = su;
        }
        this->resize(tmp.size());
        //debug << "mul " << (vector<T>)*this << " " << (vector<T>)o << " : " << (vector<T>)(tmp) << "\n";
        copy(tmp.begin(), tmp.end(), this->begin());
        return *this;
    }
    Slow_Poly operator%=(Slow_Poly const&o){
        auto inv = -T{1}/o.back();
        static vector<Summer> tmp;
        tmp.resize(this->size());
        for(int i=0; i<(int)this->size(); ++i){
            tmp[i] = Summer((*this)[i]);
        }
        for(int i=tmp.size()-1; i>=(int)o.size()-1; --i){
            T fa = T(tmp[i])*inv;
            for(int j=0;j<(int)o.size(); ++j){
                tmp[i-j].addmul(fa, o.rbegin()[j]);
            }
        }
        const int k = min(tmp.size(), o.size()-1);
        this->resize(k, T{});
        for(int i=0; i<k; ++i){
            (*this)[i] = tmp[i];
        }
        return *this;
    }
};

template<typename T>
vector<T> berlekamp_massey(vector<T> const&vals, bool make_monic=true){
    const int n = vals.size();
    int e=1, p=1;
    T delta_old{1};
    vector<T> rho_old{1}, rho{1}, rho_new;
    for(int i=0; i<n; ++i){
        typename T::Summer x;
        for(int j=0; j<(int)rho.size(); ++j){
            x.addmul(rho[j], vals[i-j]);
        }
        T delta(x);
        if(!!delta){
            rho_new.assign(max(rho.size(), p+rho_old.size()), T{});
            for(int j=0; j<(int)rho.size(); ++j){
                rho_new[j] = delta_old*rho[j];
            }
            for(int j=0; j<(int)rho_old.size(); ++j){
                rho_new[p+j] -= delta*rho_old[j];
            }
            rho.swap(rho_new);
            if(e > 0){
                rho_old.swap(rho_new);
                delta_old = delta;
                p = 0;
                e = -e;
            }
        }
        ++p;
        ++e;
    }
    if(make_monic){
        T fac = T{1}/rho[0];
        for(auto &f:rho) f*=fac;
    }
    reverse(rho.begin(), rho.end());
    return rho;

}

template<typename T>
T extend_linear(vector<T> vals, ll pos){
    Slow_Poly<T> coeff(berlekamp_massey<T>(vals));
    Slow_Poly<T> X(vector<T>{0, 1}), cur(vector<T>{1});
    for(int j = __lg(pos);j>=0;--j){
        cur *= cur;
        cur %= coeff;
        if((pos>>j)&1){
            cur *= X;
            cur %= coeff;
        }
    }
    typename T::Summer su;
    for(int i=0; i<(int)cur.size(); ++i){
        su.addmul(cur[i], vals[i]);
    }
    return T(su);
}

using num = Mod_Int<fixed_mod<1000000007> >;
using mut_num = Mod_Int<mutable_mod>;
using ntt_num = Mod_Int<fixed_mod<998244353> >;
