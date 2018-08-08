#include <bits/stdc++.h>
using ll = long long;
using ull = unsigned long long;
using fl = long double;

template<typename S, typename T>
void xmin(S&a, T const&b){if(b<a) a=b;}
template<typename S, typename T>
void xmax(S&a, T const&b){if(b>a) a=b;}

using namespace std;


void solve(){
    /// SOLVE HERE

}

signed gen(int T){
    mt19937 rng(43151);
    auto get_rand = [&](int64_t l, int64_t r){
        return uniform_int_distribution<int64_t>(l, r)(rng);
    };
    auto get_double = [&](double l, double r){
        return uniform_real_distribution<double>(l, r)(rng);
    };
    ofstream o("gen.txt");
    o << T << "\n";
    for(int cas=0;cas<T;++cas){
        /// GEN HERE

        o << "\n";
    }
    o << endl;
    o.close();
    return 0;
}

int main()
{
    #ifdef LOCAL_RUN
    freopen("inX.txt", "r", stdin);
    //freopen("outX.txt", "w", stdout);
    #endif // LOCAL_RUN
    cin.tie(0); cout.tie(0);
    ios_base::sync_with_stdio(false);
    int TTT; cin >> TTT;
    for(int cas = 1;cas<=TTT;++cas){
        cout << "Case #" << cas << ": ";

        solve();

        //cout << "\n";
        #ifdef LOCAL_RUN
        cout << flush;
        #endif // LOCAL_RUN
    }
    return 0;
}
