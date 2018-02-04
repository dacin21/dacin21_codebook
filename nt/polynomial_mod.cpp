#include <bits/stdc++.h>
/*
 *  FFT implementation with double by dacin21
 *  splits digits into the imaginary part to cope with bigger numbers
*/

using namespace std;
using ll = long long;

namespace fft{
// floored base 2 logarithm
int log2i(unsigned long long a){
    return __builtin_clzll(1) - __builtin_clzll(a);
}
const double PI = 3.1415926535897932384626;
vector<complex<double> > roots;
// pre-calculate complex roots, log(N) calls to sin/cos
void gen_roots(int N){
    if((int)roots.size()!=N){
        roots.clear();
        roots.resize(N);
        for(int i=0;i<N;++i){
            if((i&-i) == i){
                roots[i] = polar(1.0, 2.0*PI*i/N);
            } else {
                roots[i] = roots[i&-i] * roots[i-(i&-i)];
            }
        }
    }
}
void fft(complex<double> const*a, complex<double> *to, int n, bool isInv = false){
    to[0] = a[0];
    for (int i=1, j=0; i<n; ++i) {
        int m = n >> 1;
        for (; j>=m; m >>= 1)
            j -= m;
        j += m;
        to[i] = a[j];
    }
    gen_roots(n);
    for(int iter=1, sh=log2i(n)-1;iter<n;iter*=2, --sh){
        for(int x=0;x<n;x+=2*iter){
            for(int y=0;y<iter;++y){
                complex<double> ome = roots[y<<sh];
                if(isInv) ome = conj(ome);
                complex<double> v = to[x+y], w=to[x+y+iter];
                to[x+y] = v+ome*w;
                to[x+y+iter] = v-ome*w;
            }
        }
    }
}
template<ll mod, typename int_t>
vector<int_t> poly_mul(vector<int_t> const&a, vector<int_t> const&b){
    int logn = log2i(a.size()+b.size()-1)+1;
    int n = 1<<logn;
    vector<complex<double> > x(n), y(n), xx(n), yy(n);
    // split digit into real and imaginary part
    for(int i=0;i<(int)a.size();++i) x[i] = complex<double>(a[i]&((1<<15)-1), a[i]>>15);
    for(int i=0;i<(int)b.size();++i) y[i] = complex<double>(b[i]&((1<<15)-1), b[i]>>15);

    fft(x.data(), xx.data(), n, false);
    fft(y.data(), yy.data(), n, false);
    // use that fft(conj(x)) = reverse(conj(fft(x)))
    // to recover fft(real(x)) and fft(imag(x))
    for(int i=0;i<n;++i){
        int j = (n-i)&(n-1); //reverse index
        complex<double> rx = (xx[i] + conj(xx[j]))*0.5;
        complex<double> ix = (xx[i] - conj(xx[j]))*complex<double>(0, -0.5);
        complex<double> ry = (yy[i] + conj(yy[j]))*0.5;
        complex<double> iy = (yy[i] - conj(yy[j]))*complex<double>(0, -0.5);
        x[i] = (rx*ry + ix*iy*complex<double>(0, 1.0))/(double)n;
        y[i] = (rx*iy + ix*ry)/(double)n;
    }
    fft(x.data(), xx.data(), n, true);
    fft(y.data(), yy.data(), n, true);
    vector<int_t> ret(a.size()+b.size()-1);
    for(int i=0;i<(int)ret.size();++i){
        ll l = llround(xx[i].real()), m = llround(yy[i].real()), r=llround(xx[i].imag());
        ret[i] = (l + (m%mod<<15) + (r%mod<<30))%mod;
    }
    return ret;
}
}
template<ll mod>
struct NT{
    static int add(int const&a, int const&b){
        ll ret = a+b;
        if(ret>=mod) ret-=mod;
        return ret;
    }
    static int& xadd(int& a, int const&b){
        a+=b;
        if(a>=mod) a-=mod;
        return a;
    }
    static int sub(int const&a, int const&b){
        return add(a, mod-b);
    }
    static int& xsub(int& a, int const&b){
        return xadd(a, mod-b);
    }
    static int mul(int const&a, int const&b){
        return a*(ll)b%mod;
    }
    static int& xmul(int &a, int const&b){
        return a=mul(a, b);
    }
    static int inv_rec(int const&a, int const&m){
        assert(a!=0);
        if(a==1) return 1;
        int ret = m+(1-inv_rec(m%a, a)*(ll)m)/a;
        return ret;
    }
	// this is soooo great, can even be used for a sieve
    static int inv_rec_2(int const&a, int const&m){
        assert(a!=0);
        if(a==1) return 1;
        int ret = m-NT<mod>::mul((m/a), inv_rec_2(m%a, m));
        return ret;
    }
    static int inv(int const&a){
        return inv_rec_2(a, mod);
    }
};

template<ll mod>
struct poly : vector<int>{
    poly(size_t a):vector<int>(a){}
    poly(size_t a, int b):vector<int>(a, b){}
    poly(vector<int> const&a):vector<int>(a){}
    poly& normalize(){
        while(size()>1 && back()==0) pop_back();
        return *this;
    }
    poly substr(int l, int r)const{
        if(r>(int)size()) r=size();
        if(l>(int)size()) l=size();
        return poly(vector<int>(begin()+l, begin()+r));
    }
    poly reversed()const{
        return poly(vector<int>(rbegin(), rend()));
    }
    poly operator+(poly const&o)const{
        poly ret(max(size(), o.size()));
        copy(begin(), end(), ret.begin());
        for(int i=0;i<(int)o.size();++i)
            NT<mod>::xadd(ret[i], o[i]);
        return ret.normalize();;
    }
    poly operator-(poly const&o)const{
        poly ret(max(size(), o.size()));
        copy(begin(), end(), ret.begin());
        for(int i=0;i<(int)o.size();++i)
            NT<mod>::xsub(ret[i], o[i]);
        return ret.normalize();
    }
    poly operator*(poly const&o)const{
        poly ret(fft::poly_mul<mod, int>(*this, o));
        return ret.normalize();
    }
    friend ostream& operator<<(ostream&o, poly const&p){
        for(int i=(int)p.size()-1;i>=0;--i){
            o <<(abs(p[i]-mod)<p[i] ? p[i]-mod:p[i]);
            if(i){
                o << "*x";
                if(i>1) o << "^" << i;
                o << " + ";
            }
        }
        return o;
    }
    // inverse mod x^n
    poly inv(int n)const{
        assert(size() && operator[](0));
        if((int)size()>n) return poly(vector<int>(begin(), begin()+n)).inv(n);

        poly ret(1, NT<mod>::inv(operator[](0)));
        ret.reserve(2*n);
        for(int i=1;i<n;i*=2){
            poly l = substr(0, i) * ret; // l[0:i] will be 0
            poly r = substr(i, 2*i) * ret; // r[i:2*i] will be irrelevant
            poly up = (l.substr(i, 2*i) + r.substr(0, i)) * ret;
            ret.resize(2*i);
            for(int j=0;j<i;++j){
                ret[i+j] = NT<mod>::sub(0, up[j]);
            }
        }
        ret.resize(n);
        return ret.normalize();
    }
    pair<poly, poly> div(poly const&o, poly const& oinvrev)const{
        if(o.size()>size()) return {poly(1, 0), *this};
        int rsize = size()-o.size()+1;
        poly q = (reversed()*oinvrev.substr(0, rsize));
        q.resize(rsize);
        reverse(q.begin(), q.end());
        poly r = *this - q*(o);
        return make_pair(q, r.normalize());
    }
    pair<poly, poly> div(poly const&o)const{
        return div(o, o.reversed().inv(size()+2));
    }
};
// n-th term of linear recurrence in O(K log K log N)
signed codechef_rng(){
    const int mod = 104857601;
    int k;
    ll n;
    cin >> k >> n;
    --n;
    vector<int> a(k), c(k);
    for(auto &e:a) cin >> e;
    for(auto &e:c) cin >> e;
    poly<mod> p(k+1, 1);
    for(int i=0;i<k;++i){
        p[k-i-1] = (mod-c[i])%mod;;
    }
    poly<mod> b(1, 1), x(vector<int>({0, 1}));
    poly<mod> previnv = p.reversed().inv(p.size()+2);
    for(ll i=1ll<<61;i;i>>=1){
        if(2*i<=n){
            b = (b*b).div(p, previnv).second;
        }
        if(n&i){
            b = (b*x).div(p, previnv).second; // could be optimized
        }
    }
    b.resize(k, 0);
    int res = 0;
    for(int i=0;i<k;++i){
        NT<mod>::xadd(res, NT<mod>::mul(b[i], a[i]));
    }
    cout << res << "\n";
    return 0;
}



signed main(){
    #ifdef LOCAL_RUN
    freopen("in.txt", "r", stdin);
    #endif
    cin.tie(0); ios_base::sync_with_stdio(false);
    return codechef_rng();
    const int mod = 1e9+7;

    int n, m;
    cin >> n;
    poly<mod> a(n);
    for(auto &e:a) cin >> e;
    cin >> m;
    poly<mod> b(n);
    for(auto &e:b) cin >> e;
    cout << "a = " << a << "\n";
    cout << "b = " << b << "\n";
    cout << "a+b = " << a+b << "\n";
    cout << "a-b = " << a-b << "\n";
    cout << "a*b = " << a*b << "\n";
    cout << "a/b = " << a.div(b).first << "\n";
    cout << "a%b = " << a.div(b).second << "\n";


}
