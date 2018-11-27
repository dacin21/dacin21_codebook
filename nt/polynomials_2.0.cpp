/*
 *  fft with doubles
 *  fft with doubles
 *  operations in prime fields
 *  operations on polynomials and power series
 */


#define USE_PRAGMA 0

#if USE_PRAGMA
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")
#pragma GCC target("fma")
#include <smmintrin.h>
#include <immintrin.h>
// vectorized complex multiplication
// deals with stack alignment problems by forced inlining :-)
__attribute__((always_inline))
__m256d complex_mul(__m256d const&a, __m256d const&b){
    __m256d c = _mm256_movedup_pd(a);
    __m256d d = _mm256_shuffle_pd(a, a, 15);
    __m256d cb = _mm256_mul_pd(c, b);
    __m256d db = _mm256_mul_pd(d, b);
    __m256d e = _mm256_shuffle_pd(db, db, 5);
    __m256d r = _mm256_addsub_pd(cb, e);
    return r;
}
#endif // USE_PRAGMA



#include <bits/stdc++.h>

template<typename complex_t>
class Fft{
    static constexpr double PI(){return 3.1415926535897932384626;};
    static constexpr int SPLIT_BITS = 15;
    using cvec = std::vector<complex_t>;
    using ivec = std::vector<int64_t>;

    static int lg(int x){
        return __builtin_clz(1) - __builtin_clz(x);
    }
    static int next_2pow(int x){
        return x>1 ? 1<<lg(2*x-1): 1;
    }

    static std::vector<complex_t> roots;
    static void gen_roots(int const&depth){
        if(roots.size() < (1u<<depth)){
            if(roots.empty()){
                roots.emplace_back(std::polar(0.0, 0.0));
                roots.emplace_back(std::polar(1.0, 0.0));
            }
            int cur_d = lg(roots.size())-1;
            roots.resize(1<<(depth+1));
            for(int i = cur_d;i<depth;++i){
                int n = 1 << i;
                complex_t ome = std::polar(1.0, PI()/2.0/n);
                for(int j=n;j<2*n;++j){
                    roots[2*j] = roots[j];
                    roots[2*j+1] = roots[j]*ome;
                }
            }
        }
    }
    static void bitflip(cvec &v){
        int n = v.size();
        for(int i=1, j=0;i<n;++i){
            int m = n >> 1;
            for(;j>=m;m >>= 1)
                j -= m;
            j += m;
            if(i < j) std::swap(v[i], v[j]);
        }
    }
    template<bool is_rev>
    static void fft_regular(cvec& vec){
        bitflip(vec);
        int n = vec.size();
        gen_roots(lg(n));
        for(int iter=1, sh=0;iter<n;iter*=2, ++sh){
            for(int x=0;x<n;x+=2*iter){
                for(int y=0;y<iter;++y){
                    complex_t ome = roots[(1 << sh) + y];
                    if (is_rev) ome = conj(ome);
                    complex_t v = vec[x+y], w=vec[x+y+iter];
                    vec[x+y] = v+ome*w;
                    vec[x+y+iter] = v-ome*w;
                }
            }
        }
    }
    template<bool is_rev>
    static void fft_inplace(cvec&v){
        fft_regular<is_rev>(v);
    }

    static void complexify(cvec& vout, ivec const&vin, size_t const&min_size){
        vout.resize(next_2pow(min_size));
        for(size_t i=0;i<vin.size();++i){
            vout[i] = complex_t(static_cast<double>(vin[i]), 0.0);
        }
        std::fill(vout.begin()+vin.size(), vout.end(), complex_t(0.0, 0.0));
    }
    static void complexify_pair(cvec& vout, ivec const& vin1, ivec const& vin2, size_t const&min_size){
        vout.resize(next_2pow(min_size));
        size_t min_siz = std::min(vin1.size(), vin2.size());
        for(size_t i=0;i<min_siz;++i){
            vout[i] = complex_t(static_cast<double>(vin1[i]), static_cast<double>(vin2[i]));
        }
        for(size_t i=min_siz;i<vin1.size();++i){
            vout[i] = complex_t(static_cast<double>(vin1[i]), 0.0);
        }
        for(size_t i=min_siz;i<vin2.size();++i){
            vout[i] = complex_t(0.0, static_cast<double>(vin2[i]));
        }
        std::fill(vout.begin()+min_siz, vout.end(), complex_t(0.0, 0.0));
    }
    static void uncomplexify(ivec &vout, cvec const& vin, size_t const&result_size){
        vout.resize(result_size);
        for(size_t i=0;i<result_size;++i){
            vout[i] = llround(vin[i].real());
        }
    }
    static void complexify_split(cvec& vout, ivec const& vin, size_t const&min_size){
        vout.resize(next_2pow(min_size));
        constexpr int HALF_BITS = 15;
        for(size_t i=0;i<vin.size();++i){
            int64_t high = vin[i] >> HALF_BITS, low = vin[i]&((1<<HALF_BITS)-1);
            vout[i] = complex_t(static_cast<double>(low), static_cast<double>(high));
        }
        std::fill(vout.begin()+vin.size(), vout.end(), complex_t(0.0, 0.0));
    }

public:
    static void poly_mul(ivec& vout, ivec const& v1, ivec const& v2, const size_t mod = 0){
        static cvec c1, c2;
        const size_t out_size = v1.size() + v2.size() - 1;
        complexify(c1, v1, out_size);
        complexify(c2, v2, out_size);
        fft_inplace<false>(c1);
        fft_inplace<false>(c2);
        const double factor = 1 / static_cast<double>(c1.size());
        for(size_t i=0;i<c1.size();++i){
            c1[i]*=c2[i]*factor;
        }
        fft_inplace<true>(c1);
        uncomplexify(vout, c1, out_size);
        if(mod){
            for(auto &e:vout) e%=mod;
        }
    }
    static void poly_mul_faster(ivec& vout, ivec const& v1, ivec const& v2, const size_t mod = 0){
        static cvec c1;
        const size_t out_size = v1.size() + v2.size() - 1;
        complexify_pair(c1, v1, v2, out_size);
        fft_inplace<false>(c1);
        const double factor = 0.25 / static_cast<double>(c1.size());
        const int n = c1.size();
        for(int i=0;i<=n/2;++i){
            int j = (n-i)&(n-1); //reverse index
            complex_t rx = (c1[i] + conj(c1[j]));
            complex_t ix = (c1[i] - conj(c1[j]))*complex_t(0, -1);
            c1[i] = (rx*ix) * factor;
            c1[j] = conj(rx*ix) * factor;
        }
        fft_inplace<true>(c1);
        uncomplexify(vout, c1, out_size);
        if(mod){
            for(auto &e:vout) e%=mod;
        }
    }
    static void poly_mul_split(ivec& vout, ivec const&v1, ivec const& v2, const size_t mod = 0){
        static cvec c1, c2;
        const size_t out_size = v1.size() + v2.size() - 1;
        complexify_split(c1, v1, out_size);
        complexify_split(c2, v2, out_size);
        fft_inplace<false>(c1);
        fft_inplace<false>(c2);

        const double factor = 0.25 / static_cast<double>(c1.size());
        const complex_t IMAG (0.0, 1.0);
        const complex_t NIMAG = conj(IMAG);
        // use that fft(conj(x)) = reverse(conj(fft(x)))
        // to recover fft(real(x)) and fft(imag(x))
        const int n = c1.size();
        for(int i=0;i<=n/2;++i){
            int j = (n-i)&(n-1); //reverse index
            complex_t rx = (c1[i] + conj(c1[j]));
            complex_t ix = (c1[i] - conj(c1[j])); // * -i
            complex_t ry = (c2[i] + conj(c2[j]));
            complex_t iy = (c2[i] - conj(c2[j])); // * -i
            c1[i] = (rx*ry - ix*iy * IMAG) * factor;
            c2[i] = (rx*iy + ix*ry) * factor;
            c1[j] = conj(rx*ry - ix*iy * NIMAG) * factor;
            c2[j] = -conj(rx*iy + ix*ry) * factor;
        }
        fft_inplace<true>(c1);
        fft_inplace<true>(c2);
        vout.resize(out_size);
        for(size_t i=0;i<out_size;++i){
            int64_t l = llround(c1[i].real()), m = llround(c2[i].imag()), r=llround(c1[i].imag());
            if(mod) vout[i] = (l + (m%mod<<SPLIT_BITS) + (r%mod<<(2*SPLIT_BITS)))%mod;
            else vout[i] = l + (m<<SPLIT_BITS) + (r<<(2*SPLIT_BITS));
        }
    }
        static void poly_mul_third(ivec& vout, ivec const&v1, ivec const&v2, const size_t mod = 0){
        if(v1.size() > v2.size()) return poly_mul_third(vout, v2, v1, mod);
        constexpr int BITS = 17, MASK = (1<<BITS) - 1;
        const size_t out_size = v1.size() + v2.size() - 1;
        const size_t n = next_2pow(out_size);

        static cvec c1, c2, c3;
        c1.resize(n); c2.resize(n), c3.resize(n);
        for(size_t i=0;i<v1.size();++i){
            int64_t a = v1[i] >> (2*BITS), b = (v1[i] >> BITS) & MASK, c = v1[i] & MASK;
            int64_t d = v2[i] >> (2*BITS), e = (v2[i] >> BITS) & MASK, f = v2[i] & MASK;
            if(MASK-c < c){ c-= 1<<BITS; ++b; }
            if(MASK-f < f){ f-= 1<<BITS; ++e; }
            if(MASK-b < b){ b-= 1<<BITS; ++a; }
            if(MASK-e < e){ e-= 1<<BITS; ++d; }
            c1[i] = complex_t(a, b);
            c2[i] = complex_t(c, d);
            c3[i] = complex_t(e, f);
        }
        for(size_t i=v1.size();i<v2.size();++i){
            int64_t d = v2[i] >> (2*BITS), e = (v2[i] >> BITS) & MASK, f = v2[i] & MASK;
            if(MASK-f < f){ f-= 1<<BITS; ++e; }
            if(MASK-e < e){ e-= 1<<BITS; ++d; }
            c2[i] = complex_t(0, d);
            c3[i] = complex_t(e, f);
        }
        fill(c1.begin()+v1.size(), c1.end(), complex_t(0, 0));
        fill(c2.begin()+v2.size(), c2.end(), complex_t(0, 0));
        fill(c3.begin()+v2.size(), c3.end(), complex_t(0, 0));
        fft_inplace<false>(c1);
        fft_inplace<false>(c2);
        fft_inplace<false>(c3);

        const double factor = 0.25 / static_cast<double>(c1.size());
        const complex_t IMAG (0.0, 1.0);
        const complex_t NIMAG = conj(IMAG);
        const complex_t NIMAG_factor = NIMAG * factor;

        /*for(int i=0;i<=n/2;++i){
            int j = (n-i)&(n-1); //reverse index
            complex_t rx = (c1[i] + conj(c1[j]));
            complex_t ix = (c1[i] - conj(c1[j])) * NIMAG;
            complex_t ry = (c2[i] + conj(c2[j]));
            complex_t iy = (c2[i] - conj(c2[j])) * NIMAG_factor;
            complex_t rz = (c3[i] + conj(c3[j])) * factor;
            complex_t iz = (c3[i] - conj(c3[j])) * NIMAG_factor;

            c1[i] = (ry*iz + (ix*iz + ry*rz) * IMAG);
            c2[i] = (rx*iz + ix*rz + ry*iy + (ix*iy + rx*rz) * IMAG);
            c3[i] = (rx*iy);
            c1[j] = conj(ry*iz + (ix*iz + ry*rz) * NIMAG);
            c2[j] = conj(rx*iz + ix*rz + ry*iy + (ix*iy + rx*rz) * NIMAG);
            c3[j] = conj(rx*iy);
        }*/
        for(int i=0;i<=n/2;++i){
            int j = (n-i)&(n-1); //reverse index
            complex_t rx = (c1[i] + conj(c1[j]));
            complex_t ix = (c1[i] - conj(c1[j]));
            complex_t ry = (c2[i] + conj(c2[j]));
            complex_t iy = (c2[i] - conj(c2[j]));
            complex_t rz = (c3[i] + conj(c3[j]));
            complex_t iz = (c3[i] - conj(c3[j]));

            c1[i] = (ry*iz + (ix*iz - ry*rz)) * NIMAG_factor;
            c2[i] = (rx*iz + ix*rz + ry*iy + (ix*iy - rx*rz)) * NIMAG_factor;
            c3[i] = (rx*iy) * NIMAG_factor;
            c1[j] = conj((ry*iz + (ry*rz - ix*iz)) * NIMAG_factor);
            c2[j] = conj((rx*iz + ix*rz + ry*iy + (rx*rz - ix*iy)) * NIMAG_factor);
            c3[j] = conj(rx*iy * NIMAG_factor);
        }
        fft_inplace<true>(c1);
        fft_inplace<true>(c2);
        fft_inplace<true>(c3);

        vout.resize(out_size);
        //double err = 0.0;
        //double maxval = 0.0;
        for(int i=0;i<out_size;++i){
            int64_t a = llround(c1[i].real()), b = llround(c1[i].imag()), c = llround(c2[i].real());
            int64_t d = llround(c2[i].imag()), e = llround(c3[i].real()), f = llround(c3[i].imag());
            //double real_error = max(abs(a-c1[i].real()), max(abs(b-c1[i].imag()), max(abs(c-c2[i].real()), max(abs(d-c2[i].imag()), max(abs(e-c3[i].real()), abs(f-c3[i].imag())))))); err = max(err, real_error);
            //maxval = max(maxval, max(abs(c1[i].real()), max(abs(c1[i].imag()), max(abs(c2[i].real()), max(abs(c2[i].imag()), max(abs(c3[i].real()), abs(c3[i].imag())))))));
            vout[i] = (((((((((((f << BITS) % mod + e) << BITS) % mod + d) << BITS) % mod + c) << BITS) % mod) + b) << BITS) % mod + a)% mod;
            if(vout[i] < 0) vout[i] += mod;
        }
        //cerr << "=> " << err << "  " << maxval << "\n";
    }
    static void poly_square_split(ivec& vout, ivec const&v1, const size_t mod = 0){
        static cvec c1, c2;
        const size_t out_size = 2*v1.size() - 1;
        complexify_split(c1, v1, out_size);
        fft_inplace<false>(c1);

        const double factor = 0.25 / static_cast<double>(c1.size());
        const complex_t IMAG (0.0, 1.0);
        const complex_t NIMAG = conj(IMAG);
        // use that fft(conj(x)) = reverse(conj(fft(x)))
        // to recover fft(real(x)) and fft(imag(x))
        const int n = c1.size();
        c2.resize(n);
        for(int i=0;i<=n/2;++i){
            int j = (n-i)&(n-1); //reverse index
            complex_t rx = (c1[i] + conj(c1[j]));
            complex_t ix = (c1[i] - conj(c1[j])); // * -i
            c1[i] = (rx*rx - ix*ix * IMAG) * factor;
            c2[i] = (rx*ix) * 2.0 * factor;
            c1[j] = conj(rx*rx - ix*ix * NIMAG) * factor;
            c2[j] = -2.0 * conj(rx*ix) * factor;
        }
        fft_inplace<true>(c1);
        fft_inplace<true>(c2);
        vout.resize(out_size);
        for(size_t i=0;i<out_size;++i){
            int64_t l = llround(c1[i].real()), m = llround(c2[i].imag()), r=llround(c1[i].imag());
            if(mod) vout[i] = (l + (m%mod<<SPLIT_BITS) + (r%mod<<(2*SPLIT_BITS)))%mod;
            else vout[i] = l + (m<<SPLIT_BITS) + (r<<(2*SPLIT_BITS));
        }
    }
};

template <typename floating_t>
struct Complex
{
    floating_t re, im;
    Complex() : re(0), im(0) {}
    Complex (const floating_t& a, const floating_t& b) : re(a), im(b) {}
    Complex(const std::complex<floating_t>& c) : re(c.real()), im(c.imag()) {}
    Complex& operator = (const std::complex<floating_t>& c){
        this->re = c.real();
        this->im = c.imag();
        return *this;
    }
    floating_t real() const{return this->re;}
    floating_t imag() const{return this->im;}
    void real(floating_t const& r){this->re = r;}
    void imag(floating_t const& i){this->im = i;}
    void polar(const floating_t& rho, const floating_t& theta = 0){
        this->re = rho * cos(theta);
        this->im = rho * sin(theta);
    }
    Complex conj() const{return Complex(this->re, -this->im);}
    floating_t norm() const{return sqrt(this->re * this->re + this->im * this->im);}
    floating_t norm_squared() const{return this->re * this->re + this->im * this->im;}
    Complex inverse() const{return this->conj() / this->norm_squared();}

    Complex operator+ (const Complex<floating_t>& c) const{return Complex(this->re + c.re, this->im + c.im);}
    Complex& operator+= (const Complex<floating_t>& c){return *this = *this + c;}
    Complex operator - (const Complex<floating_t>& c) const{return Complex(this->re - c.re, this->im - c.im);}
    Complex operator- () const{return Complex(-this->re , -this->im);}
    Complex& operator-= (const Complex<floating_t>& c){return *this = *this + c;}

    Complex operator* (const Complex<floating_t>& c) const{return Complex(this->re * c.real() - this->im * c.imag(), this->re * c.imag() + this->im * c.real());}
    Complex& operator*= (const Complex<floating_t>& c){return *this = *this * c;}
    Complex operator* (const floating_t& c) const{return Complex(this->re * c, this->im * c);}
    Complex& operator*= (const floating_t& c){return *this = *this * c;}
    Complex operator/ (const Complex<floating_t>& c) const{return *this * c.inverse();}
    Complex& operator/= (const Complex<floating_t>& c){return *this = *this / c;}
    Complex operator/ (const floating_t& c) const{ return Complex(this->re / c, this->im / c);}
    Complex& operator/= (const floating_t& c){ return *this = *this / c;}

    friend std::istream& operator >> (std::istream& stream, Complex& C){
        return stream >> C.re >> C.im;
    }
    friend std::ostream& operator << (std::ostream& stream, const Complex& C){
        return stream << "(" << C.re << "," << C.im << ")";
    }
};
template<typename T>
Complex<T> conj(Complex<T> const&x){return x.conj();}
template<typename complex_t>
std::vector<complex_t> Fft<complex_t>::roots;



template<int64_t mod>
struct NT64 {
    static_assert(mod >= 0, "mod should be non-negative.");
    static int64_t add(int64_t const&a, int64_t const&b){
        int64_t ret = a+b;
        if(ret>=mod) ret-=mod;
        return ret;
    }
    static int64_t& xadd(int64_t& a, int64_t const&b){
        a+=b;
        if(a>=mod) a-=mod;
        return a;
    }
    static int64_t sub(int64_t const&a, int64_t const&b){
        return add(a, mod-b);
    }
    static int64_t& xsub(int64_t& a, int64_t const&b){
        return xadd(a, mod-b);
    }
    static int64_t mul(int64_t const&a, int64_t const&b){
        return a*static_cast<int64_t>(b)%mod;
    }
    static int64_t& xmul(int64_t &a, int64_t const&b){
        return a=mul(a, b);
    }
    static int64_t inv_rec(int64_t const&a, int64_t const&m){
        assert(a!=0);
        if(a==1) return 1;
        int64_t ret = m+(1-inv_rec(m%a, a)*static_cast<int64_t>(m))/a;
        return ret;
    }
    // this is soooo great, can even be used for a sieve
    static int64_t inv_rec_2(int64_t const&a, int64_t const&m){
        assert(a!=0);
        if(a==1) return 1;
        int64_t ret = m-mul((m/a), inv_rec_2(m%a, m));
        return ret;
    }
    static int64_t inv(int64_t const&a){
        return inv_rec_2(a, mod);
    }
    static int64_t div(int64_t const&a, int64_t const&b){
        return mul(a, inv(b));
    }
};

template<>
struct NT64 <int64_t{0}>{
    static int64_t add(int64_t const&a, int64_t const&b){
        return a+b;
    }
    static int64_t& xadd(int64_t& a, int64_t const&b){
        return a+=b;
    }
    static int64_t sub(int64_t const&a, int64_t const&b){
        return a-b;;
    }
    static int64_t& xsub(int64_t& a, int64_t const&b){
        return a-=b;
    }
    static int64_t mul(int64_t const&a, int64_t const&b){
        return a*b;
    }
    static int64_t& xmul(int64_t &a, int64_t const&b){
        return a=mul(a, b);
    }
    static int64_t inv_rec(int64_t const&a, int64_t const&m) = delete;
    static int64_t inv_rec_2(int64_t const&a, int64_t const&m) = delete;
    static int64_t inv(int64_t const&a){
        assert(a==1 || a==-1); // There are only 2 units in Z
        return a;
    }
    static int64_t div(int64_t const&a, int64_t const&b){
        return mul(a, inv(b));
    }
};

// Not sure if the move semantics for Polynomials work all the time.
template<int64_t mod>
class Polynomial : public std::vector<int64_t> {
public:
    using fft = Fft<std::complex<double>>;

    using std::vector<int64_t>::vector;
    Polynomial() = default;
    ~Polynomial() = default;
    Polynomial(Polynomial const& o) = default;
    Polynomial(Polynomial && o) = default;
    Polynomial& operator=(Polynomial const& o) = default;
    Polynomial& operator=(Polynomial && o) = default;
    explicit Polynomial(std::vector<int64_t> const& v):vector<int64_t>(v){};
    explicit Polynomial(std::vector<int64_t> && v):vector<int64_t>(std::move(v)){};

private:
    static void poly_mul(std::vector<int64_t>& vout, std::vector<int64_t> const& v1, std::vector<int64_t> const& v2){
        if(mod!=0 && mod <= 5000){
            return fft::poly_mul_faster(vout, v1, v2, mod);
        } else {
            return fft::poly_mul_split(vout, v1, v2, mod);
        }
    }
    static void poly_square(std::vector<int64_t>& vout, std::vector<int64_t> const& v1){
        if(mod!=0 && mod <= 5000){
            return fft::poly_mul_faster(vout, v1, v1, mod);
        } else {
            return fft::poly_square_split(vout, v1, mod);
        }
    }
public:
    Polynomial& normalize() & {
        while(size()>1 && back()==0) pop_back();
        return *this;
    }
    Polynomial normalize() && {
        while(size()>1 && back()==0) pop_back();
        return std::move(*this);
    }
    Polynomial reversed() const {
        return Polynomial (rbegin(), rend());
    }
    Polynomial operator+(Polynomial const&o) const{
        Polynomial ret(std::max(size(), o.size()), 0);
        std::copy(begin(), end(), ret.begin());
        for(size_t i=0;i<o.size();++i){
            NT64<mod>::xadd(ret[i], o[i]);
        }
        ret.normalize();
        return ret;
    }
    Polynomial operator-(Polynomial const&o) const{
        Polynomial ret(std::max(size(), o.size()), 0);
        std::copy(begin(), end(), ret.begin());
        for(size_t i=0;i<o.size();++i){
            NT64<mod>::xsub(ret[i], o[i]);
        }
        ret.normalize();
        return ret;
    }

private:
    void add_inplace_impl(Polynomial const&o){
        if(o.size() > size()) resize(o.size(), 0);
        for(size_t i=0;i<o.size();++i){
            NT64<mod>::xadd((*this)[i], o[i]);
        }
        normalize();
    }
    void sub_inplace_impl(Polynomial const&o){
        if(o.size() > size()) resize(o.size(), 0);
        for(size_t i=0;i<o.size();++i){
            NT64<mod>::xsub((*this)[i], o[i]);
        }
        normalize();
    }
public:
    Polynomial& operator+=(Polynomial const&o) & {
        add_inplace_impl(o);
        return *this;
    }
    Polynomial operator+=(Polynomial const&o) && {
        add_inplace_impl(o);
        return std::move(*this);
    }
    Polynomial& operator-=(Polynomial const&o) & {
        sub_inplace_impl(o);
        return *this;
    }
    Polynomial operator-=(Polynomial const&o) && {
        sub_inplace_impl(o);
        return std::move(*this);
    }
    Polynomial operator*(Polynomial const&o) const{
        if(empty() || o.empty()) return Polynomial();
        Polynomial ret(size()+o.size()-1);
        poly_mul(ret, *this, o);
        return ret;
    }
    Polynomial& operator*=(Polynomial const&o) & {
        if(empty() || o.empty()) {
            clear();
            return *this;
        }
        poly_mul(*this, *this, o);
        return *this;
    }
    Polynomial operator*=(Polynomial const&o) && {
        if(empty() || o.empty()) {
            clear();
            return std::move(*this);
        }
        poly_mul(*this, *this, o);
        return std::move(*this);
    }

private:
    void left_shift_inplace_impl(size_t const n){
        resize(size()+n, 0);
        std::rotate(begin(), begin() + (size()-n), end());
        normalize();
    }
public:
    Polynomial& operator<<=(size_t const n) & {
        left_shift_inplace_impl(n);
        return *this;
    }
    Polynomial operator<<=(size_t const n) && {
        left_shift_inplace_impl(n);
        return std::move(*this);
    }
    Polynomial operator<<(size_t const n) const {
        Polynomial tmp(*this);
        tmp<<=n;
        return tmp;
    }

private:
    void square_inplace_impl(){
        poly_square(*this, *this);
    }
public:
    Polynomial& square_inplace() & {
        square_inplace_impl();
        return *this;
    }
    Polynomial square_inplace() && {
        square_inplace_impl();
        return std::move(*this);
    }
    Polynomial square() const{
        Polynomial tmp(*this);
        tmp.square_inplace();
        return tmp;
    }

private:
    // should optimize 3 calls to fft(ret) to 1
    template<typename T>
    Polynomial inv_impl(size_t n, T (std::vector<int64_t>::* it_fun) () const) const {
        assert(!empty() && operator[](0) != 0);
        Polynomial ret, tmp_l, tmp_r;
        ret.reserve(2*n);
        tmp_l.reserve(2*n);
        tmp_r.reserve(2*n);
        ret.push_back(NT64<mod>::inv((this->*it_fun)() [0]));
        //cerr << "inv " << *this << "\n";
        for(size_t i=1;i<n;i*=2){
            tmp_l.resize(i);
            tmp_r.resize(i);
            std::copy((this->*it_fun)(), (this->*it_fun)() + std::min(i, size()), tmp_l.begin());
            std::fill(tmp_l.begin() + std::min(i, size()), tmp_l.end(), 0);
            std::copy((this->*it_fun)() + std::min(i, size()), (this->*it_fun)()+std::min(2*i, size()), tmp_r.begin());
            std::fill(tmp_r.begin() + (std::min(2*i, size()) - std::min(i, size())), tmp_r.end(), 0);

            poly_mul(tmp_l, tmp_l, ret);
            poly_mul(tmp_r, tmp_r, ret);
            tmp_r.resize(i);
            for(size_t j = 0;j<i-1;++j){
                NT64<mod>::xadd(tmp_r[j], tmp_l[j+i]);
            }

            poly_mul(tmp_r, tmp_r, ret);
            ret.resize(2*i);
            for(size_t j=0;j<i;++j){
                ret[i+j] = NT64<mod>::sub(0, tmp_r[j]);
            }
        }
        ret.resize(n);
        ret.normalize();
        return ret;
    }

public:
    Polynomial inv(size_t n) const{
        //using cbegin_ptr_t = std::vector<int64_t>::const_iterator (vector<int64_t>::*) () const;
        //cbegin_ptr_t cbegin_fun = &std::vector<int64_t>::cbegin;
        return inv_impl(n, &std::vector<int64_t>::cbegin);
    }
    Polynomial inv_rev(size_t n) const{
        using rbegin_ptr_t = std::vector<int64_t>::const_reverse_iterator (vector<int64_t>::*) () const;
        rbegin_ptr_t rbegin_fun = &std::vector<int64_t>::rbegin;
        return inv_impl(n, rbegin_fun);
    }
    std::pair<Polynomial, Polynomial> divmod(Polynomial const&o, Polynomial const&oinvrev)const{
        if(o.size() > size()) return std::make_pair(Polynomial(1, 0), *this);
        size_t rsize = size() - o.size() + 1;
        assert(oinvrev.size() >= rsize);
        if(oinvrev.size() > 2*rsize){
            return divmod(o, Polynomial(oinvrev.begin(), oinvrev.begin()+rsize));
        }
        Polynomial q = *this;
        std::reverse(q.begin(), q.end());
        q*=oinvrev;
        q.resize(rsize);
        std::reverse(q.begin(), q.end());
        Polynomial r = *this - q * o; // Todo: optimize temporary
        // cerr << *this << " / " << o << " | " << oinvrev << " : " << r << "\n";
        return std::make_pair(std::move(q), std::move(r).normalize());
    }
    std::pair<Polynomial, Polynomial> divmod(Polynomial const&o) const {
        if(o.size() > size()) return std::make_pair(Polynomial(1, 0), *this);
        size_t rsize = size()-o.size()+1;
        Polynomial tmp = o.inv_rev(rsize);
        if(tmp.size() < rsize) tmp.resize(rsize, 0);
        return divmod(o, tmp);
    }
    friend std::ostream& operator<<(std::ostream&o, Polynomial const&p){
        bool first = false;
        for(auto const&e:p){
            if(!first) o << " ";
            first = false;
            o << e;
        }
        return o;
    }

private:
    void derivative_inplace_impl(){
        for(size_t i=1;i<size();++i){
            operator[](i-1) = NT64<mod>::mul(operator[](i), i);
        }
        if(size()>1) pop_back();
        else if(size() == 1) back() = 0;
    }
    void integrate_inplace_impl(int64_t const constant_term){
        emplace_back();
        for(size_t i=size()-1;i>0;--i){
            operator[](i) = NT64<mod>::div(operator[](i-1), i);
        }
        operator[](0) = constant_term;
    }
public:
    Polynomial& derivative_inplace() & {
        derivative_inplace_impl();
        return *this;
    }
    Polynomial derivative_inplace() && {
        derivative_inplace_impl();
        return std::move(*this);
    }
    Polynomial derivative() const{
        return Polynomial(*this).derivative_inplace();
    }
    Polynomial& integrate_inplace(int64_t const constant_term) & {
        integrate_inplace_impl(constant_term);
        return *this;
    }
    Polynomial integrate_inplace(int64_t const constant_term) && {
        integrate_inplace_impl(constant_term);
        return std::move(*this);
    }
    Polynomial integrate(int64_t const constant_term) const{
        return Polynomial(*this).integrate_inplace(constant_term);
    }
    /**
     *  Logarithm of truncated power series
     *  Runtime: O(n log n)
     */
    Polynomial logarithm(size_t const n) const{
        if(n == 0) return Polynomial();
        assert(operator[](0) == 1);

        Polynomial tmp = derivative();
        tmp *= inv(n-1);
        tmp.resize(n-1);
        tmp.integrate_inplace(0);
        return tmp.normalize();
    }
    /**
     *  Exponential of truncated power series
     *  Runtime: O(n log n)
     */
    Polynomial exponential(size_t const n) const{
        assert(!empty() && operator[](0) == 0);
        Polynomial ret, ret_inv;
        ret.reserve(2*n), ret_inv.reserve(2*n);
        ret.push_back(1); ret_inv.push_back(1);
        for(size_t i=1;i<n;i*=2){
            Polynomial tmp = ret_inv.square();
            tmp*= ret;
            tmp.resize(i);
            tmp-= ret_inv;
            ret_inv-= tmp;
            Polynomial q(begin(), begin()+std::min(size(), i));
            q.derivative_inplace();
            q*= ret;
            for(size_t j=i;j<q.size();++j){
                NT64<mod>::xadd(q[j-i], q[j]);
            }
            q.resize(i);
            Polynomial s = ret.derivative();
            s-= q;
            s<<=1;
            for(size_t j=i;j<s.size();++j){
                NT64<mod>::xadd(s[j-i], s[j]);
            }
            s.resize(i);
            s*= ret_inv;
            s.resize(i);
            s<<=i-1;
            s.integrate_inplace(0);
            Polynomial u(2*i, 0);
            std::copy(begin(), begin() + std::min(2*i, size()), u.begin());
            u-= s;
            u.resize(2*i, 0);
            std::copy(u.begin()+i, u.end(), u.begin());
            u.resize(i);
            u*= ret;
            u.resize(i);
            ret.resize(2*i);
            for(size_t j=0;j<i;++j){
                ret[i+j] = u[j];
            }
        }
        ret.resize(n);
        ret.normalize();
        return ret;
    }
};


using fft = Fft<std::complex<double>>;
using fft_c = Fft<Complex<double>>;

using namespace std;
using ll = int64_t;



// n-th term of linear recurrence in O(K log K log N)
signed codechef_rng(){
    constexpr int mod = 104857601;
    using Poly = Polynomial<mod>;
    int k;
    ll n;
    cin >> k >> n;
    --n;
    vector<int> a(k), c(k);
    for(auto &e:a) cin >> e;
    for(auto &e:c) cin >> e;
    Poly p(k+1, 1);
    for(int i=0;i<k;++i){
        p[k-i-1] = (mod-c[i])%mod;;
    }
    Poly b(1, 1), x({0, 1});
    Poly previnv = p.inv_rev(p.size()+2);
    for(ll i=1ll<<61;i;i>>=1){
        if(2*i<=n){
            /// fast
            // b.square_inplace();
            // b = b.divmod(p, previnv).second;
            /// normal
            b = (b*b).divmod(p, previnv).second;
            /// slow
            //b = (b*b).divmod(p).second;
        }
        if(n&i){
            /// faster
            //b<<=1;
            //b = b.divmod(p, previnv).second;
            /// normal
            b = (b*x).divmod(p, previnv).second;
            /// slow
            //b = (b*x).divmod(p).second;
        }
    }
    b.resize(k, 0);
    int64_t res = 0;
    for(int i=0;i<k;++i){
        NT64<mod>::xadd(res, NT64<mod>::mul(b[i], a[i]));
    }
    cout << res << "\n";
    return 0;
}
