// convex hull trick in 3D
// requires distinct monotone x-coordinates
// O(log(n)^2) per query or update

#undef _GLIBCXX_DEBUG
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
ll searches = 0;
ll builds = 0;
// persistent treap
struct Perstree{
    static mt19937 RNG;
    struct Persnode{
        array<shared_ptr<Persnode>, 2> ch;
        // my vector, vector of successor
        array<ll, 3> d, d_next;
        ll y;
        Persnode():ch{0, 0}, d{0,0,0}, d_next{0,0,0}, y(RNG()){}
        Persnode(array<ll, 3> const&v):ch{0, 0}, d(v), d_next{0,0,0}, y(RNG()){}
        Persnode(shared_ptr<Persnode> const&n):ch(n->ch), d(n->d), d_next(n->d_next), y(n->y){}
    };
    static pair<shared_ptr<Persnode>, shared_ptr<Persnode> > split(shared_ptr<Persnode> const&t, ll const&x, array<ll, 3> *&link_next, array<ll, 3> *&next_d){
        if(!t) return pair<shared_ptr<Persnode>, shared_ptr<Persnode> > (0, 0);
        shared_ptr<Persnode> ret = make_shared<Persnode>(t);
        if(x<=t->d[0]){
            auto tmp = split(t->ch[0], x, link_next, next_d);
            ret->ch[0] = tmp.second;
            if(tmp.second == nullptr){
                next_d = &(ret->d);
            }
            return {tmp.first, ret};
        } else {
            auto tmp = split(t->ch[1], x, link_next, next_d);
            ret->ch[1] = tmp.first;
            if(tmp.first == nullptr){
                link_next = &(ret->d_next);
            }
            return{ret, tmp.second};
        }
    }
    static shared_ptr<Persnode> merge(shared_ptr<Persnode> const&t, shared_ptr<Persnode> const&o){
        if(!t) return o;
        if(!o) return t;
        if(t->y>=o->y){
            shared_ptr<Persnode> ret = make_shared<Persnode>(t);
            ret->ch[1] = merge(ret->ch[1], o);
            return ret;
        } else {
            shared_ptr<Persnode> ret = make_shared<Persnode>(o);
            ret->ch[0] = merge(t, ret->ch[0]);
            return ret;
        }
    }
    shared_ptr<Persnode> root;
    Perstree():root(0){}
    Perstree(shared_ptr<Persnode> r):root(r){}
    Perstree insert_or_delete(array<ll, 3> const&p){
        array<ll, 3> * link = nullptr, *nextp = nullptr, *trash = nullptr;
        auto tmp = split(root, p[0], link, trash);
        auto tmp2 = split(tmp.second, p[0]+1, trash, nextp);
        if(!tmp2.first){
            if(link) *link= p;
            auto mid = make_shared<Persnode>(p);
            if(nextp) mid->d_next = *nextp;
            return Perstree(merge(merge(tmp.first, mid),tmp2.second));
        }
        if(link && nextp){
            *link = *nextp;
        }
        return Perstree(merge(tmp.first, tmp2.second));
    }
    // extreme point query, done with ternary search
    array<ll, 3> mascal(array<ll, 3> const&dir)const{
        assert(root);
        array<ll, 3> res = root->d;
        shared_ptr<Persnode> cur = root;
        while(cur){
            if((cur->d[0]-res[0])*dir[0] + (cur->d[1]-res[1])*dir[1] + (cur->d[2]-res[2])*dir[2]>=0){
                res = cur->d;
            }
            ++searches;
            shared_ptr<Persnode> r = cur->ch[1];
            if(!r) {
                cur = cur->ch[0];
            } else {
                if((cur->d[0]-cur->d_next[0])*dir[0] + (cur->d[1]-cur->d_next[1])*dir[1] + (cur->d[2]-cur->d_next[2])*dir[2]>=0){
                    cur = cur->ch[0];
                } else {
                    cur = cur->ch[1];
                }
            }
        }
        return res;
    }
    // debug function, prints whole tree
    void print()const{
        stack<pair<shared_ptr<Persnode>, int> > s;
        s.emplace(root, 0);
        while(!s.empty()){
            auto e = s.top();
            s.pop();
            if(e.first==0)continue;
            if(e.second){
                cerr << "  " << e.first->d[0] << "/" << e.first->d[1] << "/" << e.first->d[2] << " ";
                cerr << "|" << e.first->d_next[0] << "/" << e.first->d_next[1] << "/" << e.first->d_next[2] << " ";
            } else {
                s.emplace(e.first->ch[1], 0);
                s.emplace(e.first, 1);
                s.emplace(e.first->ch[0], 0);
            }
        }
        cerr << "\n";
    }
};
mt19937 Perstree::RNG(132132);


struct Chull3D{
    using ll = long long;

    struct Frac{
        ll a, b;
        Frac(ll const&x=0, ll const&y=1):a(x), b(y){if(b<0){a=-a; b=-b;}}
        Frac& reduce(){
            ll g = __gcd(a, b);
            a/=g, b/=g;
            if(a<0){ a=-a, b=-b; }
            return *this;
        }
        Frac& operator+(Frac const&o)const{
            ll g = __gcd(b, o.b);
            return Frac(o.b/g*a + b/g*o.a, b/g*o.b).reduce();
        }
        Frac& operator-(Frac const&o)const{
            ll g = __gcd(b, o.b);
            return Frac(o.b/g*a - b/g*o.a, b/g*o.b).reduce();
        }
        Frac& operator*(Frac const&o)const{
            return Frac(a*o.a, b*o.b).reduce();
        }
        Frac& operator/(Frac const&o)const{
            return Frac(a*o.b, b*o.a).reduce();
        }
        Frac operator-()const{
            return Frac(-a, b);
        }
        bool operator==(Frac const&o)const{
            return ((a<0) == (o.a<0)) && (a*o.b == b*o.a);
        }
        bool operator<(Frac const&o)const{
            //return a/(long double)b < o.a/(long double)o.b;
            return a*(__int128)o.b < o.a*(__int128)b;
        }
        friend ostream& operator<<(ostream&o, Frac const&f){
            return o << (f.a/(long double)f.b);
        }
    };
    struct Point{
        ll x, y, z;
        Point *prev, *next;
        Point():prev(0), next(0){}
        Point(ll _x, ll _y, ll _z, Point*p, Point*n):x(_x), y(_y), z(_z), prev(p), next(n){}
        void act() {
            if (prev->next != this) prev->next = next->prev = this; // insert
            else { prev->next = next; next->prev = prev; } // delete
        }
        tuple<ll, ll, ll> tup(){
            return {x, y, z};
        }
        array<ll,3> arr(){
            return {x, y, z};
        }
        friend ostream& operator<<(ostream&o, Point const&p){
            return o << p.x << "/" << p.y << "/" << p.z;
        }
    };
    Point cross(Point const&a, Point const&b, Point const&c){
        Point x(b.x-a.x, b.y-a.y, b.z-a.z, 0, 0);
        Point y(c.x-a.x, c.y-a.y, c.z-a.z, 0, 0);
        return Point(x.y*y.z-x.z*y.y, x.z*y.x-x.x*y.z, x.x*y.y-x.y*y.x, 0, 0);
    }
    // begin 3D convex hull
    static ll inf;
    Frac INF = Frac(1, 0);
    static Point nil;
    Point *NIL = &nil;

    ll turn0(Point *p, Point *q, Point *r) { // <0 iff cw
        if (p == NIL || q == NIL || r == NIL) return 1.0;
        return (q->x-p->x)*(r->y-p->y) - (r->x-p->x)*(q->y-p->y);
    }
    ll turn(Point *p, Point *q, Point *r) { // <0 iff cw
        ll ans = turn0(p, q, r);
        if(ans == 0){
            ans = ((q->x-p->x)*(r->z-p->z) - (r->x-p->x)*(q->z-p->z));
            //cerr << "tie break " << ans << "\n";
        }
        return ans;
    }
    Frac time(Point *p, Point *q, Point *r) { // when turn changes
        if (p == NIL || q == NIL || r == NIL) return INF;
        return Frac((q->x-p->x)*(r->z-p->z) - (r->x-p->x)*(q->z-p->z) , turn(p,q,r));
    }
    Point *sort(Point P[], int n) { // mergesort
        Point *a, *b, *c, head;
        if (n == 1) { P[0].next = NIL; return P; }
        a = sort(P, n/2);
        b = sort(P+n/2, n-n/2);
        c = &head;
        do
            if (a->tup() < b->tup()) { c = c->next = a; a = a->next; }
            else { c = c->next = b; b = b->next; }
        while (c != NIL);
        return head.next;
    }
    void hull(Point *list, int n, Point **A, Point **B) { // the algorithm
        Point *u, *v, *mid;
        pair<Frac, int> t[6], oldt, newt;
        int i, j, k, l, minl;
        if (n == 1) { A[0] = list->prev = list->next = NIL; return; }

        for (u = list, i = 0; i < n/2-1; u = u->next, i++) ;
        mid = v = u->next;
        hull(list, n/2, B, A); // recurse on left and right sides
        hull(mid, n-n/2, B+n/2*2, A+n/2*2);

        for (;;) // find initial bridge
            if (turn(u, v, v->next) <= 0) v = v->next;
            else if (turn(u->prev, u, v) <= 0) u = u->prev;
            else break;

        //cerr << "::" << n << "\n";
        //cerr << *u << "  " << *v << "\n";
        // merge by tracking bridge uv over time
        int moveu = 0, movev = 0;
        for (i = k = 0, j = n/2*2, oldt = {-INF, -1}; ; oldt = newt) {
            ++builds;
            // second doesn't matter apparently, just varies degenerate faces????
            t[0] = {time(B[i]->prev, B[i], B[i]->next), 0};
            t[1] = {time(B[j]->prev, B[j], B[j]->next), 1};
            t[2] = {time(u->prev, u, v), 2};
            t[3] = {time(u, v, v->next), 2};
            t[4] = {time(u, u->next, v), 3};
            t[5] = {time(u, v->prev, v), 3};
            if(moveu == 1) t[2].first = INF;
            if(movev == 1) t[3].first = INF;
            if(moveu == -1) t[4].first = INF;
            if(movev == -1) t[5].first = INF;
            for(int l=0;l<6;++l){
                if(t[l].first.a == 0 && t[l].first.b == 0)
                    t[l].first = oldt.first;//, cerr << "failsafe\n";
            }
            for (newt = {INF, -1}, l = 0; l < 6; l++)
            if (t[l] >= oldt && t[l] < newt) { minl = l; newt = t[l]; }
            //for(int l=0;l<6;++l){
            //    cerr << setw(14) << t[l].first;
            //}   cerr << " : " << setw(0) << minl << "\n";
            if (newt.first == INF) break;
            if(oldt.first < newt.first){ moveu = movev = 0; }
            if(minl == 2){ moveu = -1;}
            if(minl == 3){ movev = -1;}
            if(minl == 4){ moveu = 1;}
            if(minl == 5){ movev = 1;}
            switch (minl) {
                case 0: if (B[i]->tup() < u->tup()) A[k++] = B[i]; B[i++]->act(); break;
                case 1: if (B[j]->tup() > v->tup()) A[k++] = B[j]; B[j++]->act(); break;
                case 2: A[k++] = u; u = u->prev; break;
                case 3: A[k++] = v; v = v->next; break;
                case 4: A[k++] = u = u->next; break;
                case 5: A[k++] = v = v->prev; break;
            }
        }
        A[k] = NIL;
        u->next = v; v->prev = u; // now go back in time to update pointers
        for (k--; k >= 0; k--) {
            if (A[k]->tup() <= u->tup() || A[k]->tup() >= v->tup()) {
                A[k]->act();
                if (A[k] == u) u = u->prev; else if (A[k] == v) v = v->next;
            } else {
                u->next = A[k];
                A[k]->prev = u;
                v->prev = A[k];
                A[k]->next = v;
                if (A[k]->tup() < mid->tup()) u = A[k];
                else v = A[k];
            }
        }
        //cerr << *u << "  " << *v << "\n";
    }
    // end 3D convex hull
    size_t n;
    vector<array<ll, 3> > points;
    vector<Point> P;
    vector<Point*> A, B;
    vector<pair<Frac, Perstree>> hulls;
    Point* list;
    Chull3D(vector<array<ll, 3> > const&p):n(p.size()), points(p){
        // build linked list
        std::sort(points.begin(), points.end());
        for(auto const&e:points)
            P.emplace_back(e[0], e[1], e[2], (Point*)0, (Point*)0);
        list = sort(P.data(), n);
        A.resize(2*n, 0); B.resize(2*n, 0);
        // construct hull
        hull(list, n, A.data(), B.data());
        // build persistent tree
        build_hull();
    }
    void build_hull(){
        Perstree hul;
        for(Point*l = list;l!=NIL; l = l->next){
            hul = hul.insert_or_delete(l->arr());
        }
        Frac oldt = -INF;
        int i;
        for(i=0; A[i]!=NIL; A[i++]->act()){
            hulls.emplace_back(oldt, hul);
            hul = hul.insert_or_delete(A[i]->arr());
            oldt = time(A[i]->prev, A[i], A[i]->next);
            //cerr << A[i]->arr()[0] << "/" << A[i]->arr()[1] << "/" << A[i]->arr()[2] << "\n";
        }
        for(;i>0;A[--i]->act());
        hulls.emplace_back(oldt, hul);
        //for(auto const&e:hulls){
        //    cerr << e.first << ":  ";
        //    e.second.print();
        //} cerr << "\n";
    }
    // 3D extreme point query
    array<ll, 3> extreme_point(array<ll, 3> const&dir)const{
        // find 2D-hull at correct timepoint
        Frac t(-dir[1], dir[2]);
        auto it = lower_bound(hulls.begin(), hulls.end(), t,[](pair<Frac, Perstree>const&a, Frac const&b){return a.first<b;});
        // query 2D hull
        if(it == hulls.begin()) return it->second.mascal(dir);
        if(it == hulls.end()) return prev(it)->second.mascal(dir);
        array<ll, 3> r1 = it->second.mascal(dir);
        array<ll, 3> r2 = prev(it)->second.mascal(dir);
        return ((r2[0]-r1[0])*dir[0] + (r2[1]-r1[1])*dir[1] + (r2[2]-r1[2])*dir[2])<0 ? r1:r2;
    }
    size_t size(){
        return n;
    }
    // faster than rebuilding, only needs single conquer step
    void merge(Chull3D &o){
        //cerr << "merge: " << size() << "/" << o.size() << "\n";
        n+=o.n;
        points.insert(points.end(), o.points.begin(), o.points.end());
        vector<Point> pres(P.size()+o.P.size());
        copy(P.begin(), P.end(), pres.begin());
        copy(o.P.begin(), o.P.end(), pres.begin()+P.size());
        A.swap(B);
        A.resize(2*pres.size());
        o.A.swap(o.B);
        for(size_t i=0;i<P.size();++i){
            if(pres[i].next!=NIL) pres[i].next = pres[i].next-P.data()+pres.data();
            if(pres[i].prev!=NIL) pres[i].prev = pres[i].prev-P.data()+pres.data();
        }
        for(size_t i=P.size();i<P.size()+o.P.size();++i){
            if(pres[i].next!=NIL) pres[i].next = pres[i].next-o.P.data()+pres.data()+P.size();
            if(pres[i].prev!=NIL) pres[i].prev = pres[i].prev-o.P.data()+pres.data()+P.size();
        }
        for(int i=0;B[i]!=NIL;++i){
            B[i] = B[i]-P.data()+pres.data();
        }
        for(int i=0;o.B[i]!=NIL;++i){
            o.B[i] = o.B[i]-o.P.data()+pres.data()+P.size();
        }
        list = list-P.data()+pres.data();
        o.list = o.list-o.P.data()+pres.data()+P.size();
        Point *u, *v, *mid;
        pair<Frac, int> t[6], oldt, newt;
        int i, j, k, l, minl;

        for (u = list; u->next!=NIL; u = u->next);
        mid = v = o.list;

        for (;;) // find initial bridge
            if (turn(u, v, v->next) <= 0) v = v->next;
            else if (turn(u->prev, u, v) <= 0) u = u->prev;
            else break;

        //cerr << "::" << n << "\n";
        //cerr << *u << "  " << *v << "\n";
        // merge by tracking bridge uv over time
        int moveu = 0, movev = 0;
        for (i = k = 0, j = 0, oldt = {-INF, -1}; ; oldt = newt) {
            ++builds;
            // second doesn't matter apparently, just varies degenerate faces????
            t[0] = {time(B[i]->prev, B[i], B[i]->next), 0};
            t[1] = {time(o.B[j]->prev, o.B[j], o.B[j]->next), 1};
            t[2] = {time(u->prev, u, v), 2};
            t[3] = {time(u, v, v->next), 2};
            t[4] = {time(u, u->next, v), 3};
            t[5] = {time(u, v->prev, v), 3};
            if(moveu == 1) t[2].first = INF;
            if(movev == 1) t[3].first = INF;
            if(moveu == -1) t[4].first = INF;
            if(movev == -1) t[5].first = INF;
            for(int l=0;l<6;++l){
                if(t[l].first.a == 0 && t[l].first.b == 0)
                    t[l].first = oldt.first;//, cerr << "failsafe\n";
            }
            for (newt = {INF, -1}, l = 0; l < 6; l++)
            if (t[l] >= oldt && t[l] < newt) { minl = l; newt = t[l]; }
            //for(int l=0;l<6;++l){
            //    cerr << setw(14) << t[l].first;
            //}   cerr << " : " << setw(0) << minl << "\n";
            if (newt.first == INF) break;
            if(oldt.first < newt.first){ moveu = movev = 0; }
            if(minl == 2){ moveu = -1;}
            if(minl == 3){ movev = -1;}
            if(minl == 4){ moveu = 1;}
            if(minl == 5){ movev = 1;}
            switch (minl) {
                case 0: if (B[i]->tup() < u->tup()) A[k++] = B[i]; B[i++]->act(); break;
                case 1: if (o.B[j]->tup() > v->tup()) A[k++] = o.B[j]; o.B[j++]->act(); break;
                case 2: A[k++] = u; u = u->prev; break;
                case 3: A[k++] = v; v = v->next; break;
                case 4: A[k++] = u = u->next; break;
                case 5: A[k++] = v = v->prev; break;
            }
        }
        A[k] = NIL;
        u->next = v; v->prev = u; // now go back in time to update pointers
        for (k--; k >= 0; k--) {
            if (A[k]->tup() <= u->tup() || A[k]->tup() >= v->tup()) {
                A[k]->act();
                if (A[k] == u) u = u->prev; else if (A[k] == v) v = v->next;
            } else {
                u->next = A[k];
                A[k]->prev = u;
                v->prev = A[k];
                A[k]->next = v;
                if (A[k]->tup() < mid->tup()) u = A[k];
                else v = A[k];
            }
        }
        //cerr << *u << "  " << *v << "\n";

        //cerr << "did"  << "\n";
        hulls.clear();
        // rebuild persistent tree
        build_hull();
        o.B.clear();
        o.A.clear();
        o.points.clear();
        o.P.clear();
        P.swap(pres);
        o.list  = 0;
        o.n = 0;
        //cerr << "done"  << "\n";

    }
};
ll Chull3D::inf = 1e17;
Chull3D::Point Chull3D::nil =  Chull3D::Point(inf, inf, inf, 0, 0);
ll operator*(array<ll, 3> const&a, array<ll, 3> const&b){
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
const ll inf = 1e16;
template<typename T>
void pvec(vector<T> const&v, string const&s){
    string c = (s + string(5, ' ')).substr(0, 5);
    cerr << c;
    for(auto const&e:v) cerr << setw(5) << e << " ";
    cerr << setw(0) << "\n";
}

// DP
vector<ll> DP2(vector<ll> const&v){
    int n = v.size();
    vector<ll> a(n);
    partial_sum(v.begin(), v.end(), a.begin());
    vector<ll> DP(n+1, 0);
    vector<array<ll, 3>> tmp; tmp.reserve(n);
    vector<Chull3D> hulls;
    // to avoid some overhead, buffer this many points before building a hull
    int l = 1+n/150;
    for(int i=0;i<n;++i){
        // merge hulls
        if((int)tmp.size()>=l){
            hulls.emplace_back(tmp);
            tmp.clear();
            while(hulls.size()>1 && hulls.rbegin()[1].size()<=hulls.rbegin()[0].size()){
                hulls.rbegin()[1].merge(hulls.rbegin()[0]);
                hulls.pop_back();
            }
        }
        // query buffered points
        array<ll, 3> pq{-a[i], -i, -1};
        ll val = a[i], tmp2;
        for(auto const&e:tmp){
            if(val < (tmp2 = e*pq)){
                val = tmp2;
            }
        }
        // query hulls
        for(auto const&e:hulls){
            tmp2 = e.extreme_point(pq)*pq;
            if(val < tmp2){
                val = tmp2;
            }
        }
        val+=a[i]*i;
        DP[i+1] = val;
        tmp.push_back({i, a[i], -(DP[i+1]+a[i]*i)});
    }
    //pvec(opt, "opt2:");
    //pvec(DP, "DP2:");
    return DP;

}
char buf[32*1024];
int bpos, bend;
inline char gc(){
    if(bpos>=bend){
        bend = fread(buf, 1, sizeof(buf)-1, stdin);
        buf[bend] = bpos = 0;
    }
    return buf[bpos++];
}
inline void readint(){}
template<typename T, typename ...S>
inline void readint(T&a, S&...b){
    int s=1;
    a=0;
    int c = gc();
    while(c!='-' && (c<'0'|| c>'9')) c = gc();
    if(c=='-') s=-1, c=gc();
    while(c>='0' && c<='9'){
        a = a*10+c-'0';
        c=gc();
    }
    if(s<0) a=-a;
    readint(b...);
}

signed main(){
    //return main_DP();
    #ifdef LOCAL_RUN
    freopen("in_3d.txt", "r", stdin);
    #endif // LOCAL_RUN
    //cin.tie(0); ios_base::sync_with_stdio(false);
    int TT; readint(TT);
    while(TT--){
        int n;readint(n);
        vector<ll> a(n);
        for(auto &e:a) readint(e);

        ll D2 = DP2(a).back();
        cout << D2 << "\n";
        //cerr << "hits: " << tot_hits << " \n";
    }
}
