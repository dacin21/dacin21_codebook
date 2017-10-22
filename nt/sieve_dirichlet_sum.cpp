// compute sum_i={1...n} f(i) where f is multiplicative,
// assuming sum_i={1...n} g(i) and sum_i={1...n} f*g
// can be computed quickly (f*g is dirichlet convolution)
// takes O(n^{2/3}) time in total
// good candidates for f and f*g are:
// f(x) = x^k, f(x) = (x==1), f(x)=1
// FOR BIG NUMERS ADD MOD
struct Dirichlet_Sum{
    using func  = ll(*)(ll);
    ll n, th; //threshold, around n^{2/3}
    func sum_f; // prefix for f (1<=x<=th)
    func sum_g; // prefix for g (0<=x<=n)
    func sum_fg; // prefix for f*g (0<=x<=n)
    unordered_map<ll, ll> cache;

    Dirichlet_Sum(func s_f, func s_g, func s_fg):
        sum_f(s_f), sum_g(s_g), sum_fg(s_fg){}

    ll calc_rec(ll x){
        if(x<=th) return sum_f(x);
        auto it = cache.find(x);
        if(it!=cache.end()) return it->second;
        ll ret = sum_fg(x);
        for(ll i=2, nex;i<=x;i=nex+1){
            nex = x/(x/i);
            ret-= (sum_g(nex) - sum_g(i-1))*calc_rec(x/i);
        }
        ret/=sum_g(1);
        return cache[x] = ret;
    }
    ll get_sum(ll _n, ll _th){
        if(_n<=0) return 0;
        n = _n; th=_th;
        return calc_rec(n);
    }
};
// computes prefix sums of multiplicative function f
// takes O(n log n) time
struct Linear_Sieve{
    using func  = ll(*)(ll, ll, int);
    int n;
    func f; // f(p^k, p, k)
    vector<ll> sum;
    Linear_Sieve(func _f):f(_f){}
    void compute(int _n){
        n=_n; sum.assign(n, 0);
        vector<char> isp(n, 1);
        sum[1] = f(1, 1, 0);
        for(int i=2;i<n;++i){
            if(isp[i]){
                ll pk = i;
                for(int k=1;pk<n;++k, pk*=i){
                    sum[pk] = f(pk, i, k);
                    for(int j=2;j<=(n-1)/pk;++j){
                        isp[j*pk] = 0;
                        sum[j*pk] = sum[j]*sum[pk];
                    }
                }
            }
        }
        partial_sum(sum.begin(), sum.end(), sum.begin());
    }
};
// example for euler totient
ll sg(ll x){return x;}
ll sfg(ll x){return x*(x+1)/2;}
Linear_Sieve ls([](ll pk, ll p, int){return pk==1?1:pk-pk/p;});
ll sf(ll x){return ls.sum[x];}
