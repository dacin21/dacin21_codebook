// given an unknown string s and Q ranges that
// are know to be palindromes, this computes the characters
// that have to be equal in O(Q + n log n)
struct Palindrome{
    Palindrome(){}
    Palindrome(int n_):n(n_), m(3+__lg(n)), qs(m), p(2*n){
        iota(p.begin(), p.end(), 0);
    }

    int f(int i){
        return p[i] == i ? i : p[i] = f(p[i]);
    }
    void u(int a, int b){
        assert(0 <= a && a < 2*n);
        assert(0 <= b && b < 2*n);
        // union with splicing
        // is a bit faster than just path compression
        // also guarantees p[i] <= i
        while(p[a] != p[b]){
            if(p[a] < p[b]) swap(a, b);
            if(p[a] == a){
                p[a] = b;
                return;
            }
            int tmp = p[a];
            p[a] = p[b];
            a = p[tmp];
        }
    }
    int components(){
        int ret = 0;
        for(int i=0;i<(int)p.size();++i){
            if(p[i] == i) ++ret;
        }
        return ret;
    }

    // call this after adding all queries
    void compute(){
        vector<int> p2(2*n);
        for(int l=m-1; l>=0; --l){
            const int s = 1<<l;
            for(int i=0; i<2*n; ++i) p2[i] = f(i);
            for(int i=0; i+s<2*n; ++i){
                const int j = p2[i];
                if(j+s < 2*n) u(i+s, j+s);
            }
            for(auto const&e:qs[l]){
                u(e.first, e.second);
            }
        }
        // link point with mirror-image
        for(int i=0;i<n;++i){
            u(i, 2*n-1 - i);
        }
    }
    // force [l, r] to be a palindrome
    void add_q(int l, int r){
        assert(0 <= l && l <= r && r < n);
        if(l==r) return;
        const int range = r-l+1;
        const int k = __lg(range);
        qs[k].emplace_back(l, 2*n-1 - r);
    }

    int n;
    int m;
    vector<vector<pair<int, int> > > qs;
    vector<int> p;
};
