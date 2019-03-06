#include <bits/stdc++.h>
using namespace std;

using Point = array<long long, 3>;
int signum(long long x) { return x ? (x>0 ? 1 : -1) : 0; }
long long sq(long long x){ return x*x; }
Point operator-(Point const&a, Point const&b){
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}
Point operator*(Point const&a, Point const&b){
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
long long dot(Point const&a, Point const&b){
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
// side of d relative to plane through a,b,c
int side(Point const&a, Point const&b, Point const&c, Point const&d){
    return signum(dot(a-d, (b-d)* (c-d)));
}
double dist(Point const&a, Point const&b, Point const&c, Point const&d){
    Point n = (b-a)*(c-a);
    return dot(n, d-a)/sqrt(dot(n, n));
}
int slope(Point const&a, Point const&b, Point const&c){
    __int128 x = sq(b[0]-a[0])*(__int128)(sq(c[1]-a[1])+sq(c[2]-a[2])) - sq(c[0]-a[0]) *(__int128)(sq(b[1]-a[1])+sq(b[2]-a[2]));
    return (x>0)-(x<0);
}

// 3d convex hull
vector<array<int, 3> > convexhull3d(vector<Point> const&pts){
    int N = pts.size();
    if(N<3) return vector<array<int, 3> >();
    vector<vector<char> > vis(N, vector<char>(N, 0));
    vector<array<int, 3> > ret;
    // pick first minimal i, pick last maximal k
    int i = min_element(pts.begin(), pts.end())-pts.begin(), j=(i==0);
    for(int k=0;k<N;++k) if(k!=i && make_pair(slope(pts[i], pts[j], pts[k]), pts[k]) >= make_pair(0, pts[j])) j=k;
    if(!any_of(pts.begin(), pts.end(), [&](Point const&a){return (a-pts[i])*(pts[j]-pts[i])!=Point({0, 0 ,0});})) return vector<array<int, 3> >();
    stack<pair<int, int> > s;
    s.emplace(i, j);
    vector<int> cands;
    while(!s.empty()){
        tie(i, j) = s.top();
        s.pop();
        if(vis[i][j]) continue;
        int k = 0;
        while(k<N && ((pts[k]-pts[i])*(pts[j]-pts[i]) == Point({0, 0 ,0}))) ++k;
        assert(k<N);
        for(int it=0;it<2;++it){ // deal with case where we start interior
            for(int l=0;l<N;++l){
                if(l==i || l==j || l==k) continue;
                if(side(pts[i], pts[j], pts[k], pts[l])<0) k=l;
            }
        }
        cands.push_back(i);
        for(int l=0;l<N;++l){
            if(l!=i && side(pts[i], pts[j], pts[k], pts[l])==0) cands.push_back(l);
        }
        cands.push_back(i);
        Point normal = (pts[j]-pts[i])*(pts[k]-pts[j]);
        for(auto &e:normal) e = signum(e);
        auto ccw = [&](Point const&a, Point const&b, Point const&c){
            return signum(dot((b-a)*(c-a), normal));
        };
        auto cmp = [&](int const&aa, int const&bb){
            Point const&a = pts[aa], &b=pts[bb];
            int x = ccw(pts[i], a, b);
            if(x) return x>0;
            return dot(a-pts[i], a-pts[i]) < dot(b-pts[i], b-pts[i]);
        };
        sort(cands.begin()+1, cands.end()-1, cmp);
        int b=0;
        for(int a=1;a<(int)cands.size();++a){
            while(b>0 && ccw(pts[cands[b-1]], pts[cands[b]], pts[cands[a]])<=0)
                --b;
            cands[++b]=cands[a];
        }
        for(int x=0;x<b;++x){
            if(x)ret.push_back({cands[0], cands[x], cands[x+1]});
            vis[cands[x]][cands[x+1]] = 1;
            s.emplace(cands[x+1], cands[x]);
        }
        ret.pop_back();
        cands.clear();
    }
    return ret;
}
vector<Point> strip_collinear(vector<Point> const&pts){
    int n = pts.size();
    vector<Point> ret;
    vector<pair<Point, bool> > tmp;
    for(int i=0;i<n;++i){
        tmp.clear();
        bool f = false;
        for(int j=0;j<n;++j){
            if(j<i && pts[i] == pts[j]) {f=true; break;}
            if(i!=j){
                Point d = pts[j]-pts[i];
                bool neg = (d<Point({0, 0, 0}));
                long long g = __gcd(__gcd(d[0], d[1]), d[2]);
                if(g==0) continue;
                g = (neg?-abs(g):abs(g));
                for(int k:{0, 1, 2}) d[k]/=g;
                tmp.emplace_back(d, neg);
            }
        }
        sort(tmp.begin(), tmp.end());
        for(int i=0;i<(int)tmp.size()-1;++i){
            if(tmp[i].first == tmp[i+1].first && tmp[i].second!=tmp[i+1].second){
                f=true;
            }
        }
        if(!f){
            ret.push_back(pts[i]);
        }
    }
    return ret;
}
// checks if origin is in hull
signed main(){
    freopen("in.txt", "r", stdin);
    freopen("out.txt", "w", stdout);
    cin.tie(0); ios_base::sync_with_stdio(false);
    int T; cin >> T;
    for(int cas=1;cas<=T;++cas){
        cout << "Case #" << cas << ": ";
        int n; cin >> n;
        vector<Point> pts;
        for(int i=0;i<n;++i){
            int a, b, c;
            cin >> a >> b >> c;
            pts.push_back({a, b, c});
        }
        auto faces = convexhull3d(pts);
        bool out = false, in  = false;
        Point orig({0, 0, 0});
        for(auto const&e:faces){
            long long a = dot(orig-pts[e[0]], (pts[e[1]] -pts[e[0]])* (pts[e[2]] - pts[e[0]]))>0;
            if(a>0)
                out = true;
            if(a<0) in = true;
            if(a == 0){
                bool f = false;
                Point normal = (pts[e[1]]-pts[e[0]])*(pts[e[2]]-pts[e[0]]);
                Point normal2 = normal;
                for(int i:{0, 1, 2}) normal2[i] = signum(normal2[i]);
                for(int i:{0, 1, 2}){
                    if(dot(normal2, (orig-pts[e[i]])*(orig-pts[e[(i+1)%3]]))<0) f = true;
                }
                if(!f) in = true;
            }
        }
        if(out || !in) cout << "NO\n";
        else cout << "YES\n";
    }
    return 0;

}
