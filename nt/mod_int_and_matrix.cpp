template<typename traits>
class Mod_Int{
public:
    using int_t = typename traits::int_t;
    using long_t = typename traits::long_t;
    static constexpr int_t mod(){ return traits::get_mod(); };


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
    int_t const& get_value (){
        return value;
    }

    Mod_Int operator+(Mod_Int const&o) const {
        Mod_Int ret(mod_step(value + o.value));
        return ret;
    }
    Mod_Int& operator+=(Mod_Int const&o){
        return *this = *this + o;
    }
    Mod_Int operator-() const{
        Mod_Int ret(mod_step(mod() - value));
        return ret;
    }
    Mod_Int operator-(Mod_Int const&o) const {
        return operator+(-o);
    }
    Mod_Int operator-=(Mod_Int const&o){
        return operator+=(-o);
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

using num = Mod_Int<fixed_mod_long<1000000007>>;
using mut_num = Mod_Int<mutable_mod>;


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

    Matrix operator*(Matrix const&o) const {
        assert(C == o.R);
        Matrix ret(R, o.C);
        for(int i=0;i<R;++i){
            for(int k=0;k<C;++k){
                for(int j=0;j<o.C;++j){
                    //ret.at(i, j) += at(i, k) * o.at(k, j);
                    ret.data[i*ret.C+j] += data[i*C+k] * o.data[k*o.C+j];
                }
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
