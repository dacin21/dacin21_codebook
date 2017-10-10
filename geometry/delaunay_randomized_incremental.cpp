// randomized incremental delaunay triangulation implementation
// runs in O(N log N) expected time
// used ~0.7 sec for 10^5 points (random input)
// by Daniel Rutschmann, 2017
// can be used for off-line point location

#include <bits/stdc++.h>
using namespace std;

#if 0
#define asser(x) assert(x)
#else
#define asser(x) do{}while(0)
#endif // 1

struct Delaunay{
    typedef long long geom_t;
    struct Point{
        geom_t x, y;
        Point():x(0), y(0){}
        Point(geom_t _x, geom_t _y):x(_x), y(_y){}
        bool operator<(const Point&b)const{return x==b.x ? y<b.y : x<b.x;}
        bool operator==(const Point&b)const{return x==b.x && y==b.y;}
        friend ostream& operator<<(ostream&f, const Point&p){
            return f << p.x << "/" << p.y;
        }
    };
    struct Face{
        array<Point, 3> corners;
        array<Face*, 3> adj;
        vector<int> bucket;
        Face():adj({0, 0, 0}), bucket(0){}
        Face(Point const&a, Point const&b, Point const&c):corners({a, b, c}), adj({0, 0 ,0}), bucket(0){}
        friend ostream& operator<<(ostream&f, const Face&ff){
            f << "Face: " << &ff << "\n";
            return f<< ff.corners[0] << "\n" << ff.corners[1] << "\n" << ff.corners[2] << "\n";
        }
    };
    static const geom_t INF = numeric_limits<int>::max()/2;
    bool is_infinite(Point const&p){
        return p.x == -INF || p.x == INF || p.y == INF || p.y == -INF;
    }
    bool is_infinite(Face const& f){
        return is_infinite(f.corners[0]) || is_infinite(f.corners[1]) || is_infinite(f.corners[2]);
    }
    Point O, UL, UR;
    vector<Face> faces;
    vector<Point> points;
    vector<Face*> location;

    Delaunay():O(0 , INF), UL(-INF, -INF), UR(INF, -INF){}

    int ccw(Point const&a, Point const&b, Point const&c){
        geom_t val = (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x);
        return (val>0)-(val<0);
    }
    // circumcircle predicate of degree 4 (might overflow/loose precision!)
    bool checkFlip(Face *cur, Point const&d){
        Point a = cur->corners[0], b=cur->corners[1], c=cur->corners[2];
		// predicate is adjusted to avoid overflow with infinite points
        if(is_infinite(a)) return ccw(b, c, d)>0;
        if(is_infinite(b)) return ccw(c, a, d)>0;
        if(is_infinite(c)) return ccw(a, b, d)>0;
        // handle collinear case
        if(ccw(a, b, c) == 0) return ccw(a, b, d)+ccw(b, c, d)+ccw(c, a, d)>0;
		// infinite point is never in any finite circle
        if(is_infinite(d)) return false;

        double a11 = a.x-d.x, a12 = a.y-d.y;
        double a21 = b.x-d.x, a22 = b.y-d.y;
        double a31 = c.x-d.x, a32 = c.y-d.y;
        double a13 = a11*a11+a12*a12, a23 = a21*a21+a22*a22, a33 = a31*a31+a32*a32;
        double det = (a11*(a22*a33-a23*a32) + (a21*(a32*a13-a12*a33)) + (a31*(a12*a23-a22*a13)));
        //cerr << "check: " << *cur << "/" << *other << "->" << det << "\n";
        bool doFlip = det >1e-2;
        return doFlip;
    }
	// self contained block allocator
    Face* getFreeFace(){
        faces.emplace_back();
        return &(faces.back());
    }
    // searches for the face containing a point?
	// adding a point to the bucket without calling insert on it should be faster
    Face* locate(Point const& p, Face* cur){
        asser(cur);
        for(int i=0;i<3;++i){
            if(ccw(cur->corners[i], cur->corners[(i+1)%3], p)<0)return locate(p, cur->adj[i]);
        }
        return cur;
    }

    // get direction of f from other face
    int get_odir(Face*f, int dir, Face* old){
        if(f->adj[dir] == 0) return -1;
        Face*o = f->adj[dir];
        int odir = -1;
        for(int i=0;i<3;++i){
            if(o->adj[i] == old) odir = i;
        }
        asser(odir != -1);
        asser(f->corners[(dir+1)%3] == o->corners[(odir+2)%3]);
        asser(f->corners[(dir+2)%3] == o->corners[(odir+1)%3]);
        return odir;
    }

    // update adj of other face to point to f
    void link_face(Face*f, int dir, Face*old){
        if(f->adj[dir] == 0) return;
        int odir = get_odir(f, dir, old);
        f->adj[dir]->adj[odir] = f;
    }
    // update pointers of points in bucket
    void link_bucket(Face*f){
        if(!f->bucket.empty()) location[f->bucket.front()]=f;
    }
    // check and perform flip in direction dir
    void check_flips(Face*f, int dir){
        if(f->adj[dir]==0) return;
        int odir = get_odir(f, dir, f);
        Face* o = f->adj[dir];
        if(checkFlip(f, o->corners[odir])){
            f->corners[(dir+1)%3] = o->corners[odir];
            o->corners[(odir+1)%3] = f->corners[dir];
            f->adj[dir] = o->adj[(odir+2)%3];
            o->adj[(odir+2)%3] = f;
            o->adj[odir] = f->adj[(dir+2)%3];
            f->adj[(dir+2)%3] = o;
            link_face(f, dir, o);
            link_face(o, odir, f);
            vector<int> tmp(f->bucket.size() + o->bucket.size());
            merge(f->bucket.begin(), f->bucket.end(), o->bucket.begin(), o->bucket.end(), tmp.begin());
            f->bucket.clear();
            o->bucket.clear();
            for(int e:tmp){
                if(ccw(f->corners[dir], f->corners[(dir+1)%3], points[e])>0){
                    asser(ccw(f->corners[0], f->corners[1], points[e])>=0);
                    asser(ccw(f->corners[1], f->corners[2], points[e])>=0);
                    asser(ccw(f->corners[2], f->corners[0], points[e])>=0);
                    f->bucket.push_back(e);
                } else {
                    asser(ccw(o->corners[0], o->corners[1], points[e])>=0);
                    asser(ccw(o->corners[1], o->corners[2], points[e])>=0);
                    asser(ccw(o->corners[2], o->corners[0], points[e])>=0);
                    o->bucket.push_back(e);
                }
            }
            link_bucket(f);
            link_bucket(o);
            check_flips(f, (dir)%3);
            check_flips(f, (dir+1)%3);
            check_flips(o, (odir)%3);
            check_flips(o, (odir+1)%3);
        }
    }
    // split face into 3 triangles, then check flips
    void split(Face*a, int point_index){
        Face*b = new (getFreeFace()) Face(a->corners[0], a->corners[1], points[point_index]);
        Face*c = new (getFreeFace()) Face(a->corners[1], a->corners[2], points[point_index]);
        a->corners[1] = points[point_index];
        b->adj = {c, a, a->adj[2]};
        c->adj = {a, b, a->adj[0]};
        a->adj = {c, a->adj[1], b};

        link_face(b, 2, a);
        link_face(c, 2, a);
        link_face(a, 1, a);
        vector<int> tmpBuck;
        tmpBuck.swap(a->bucket);
        for(int e:tmpBuck){
            if(e==point_index) continue;
            if(ccw(b->corners[1], b->corners[2], points[e])>=0 && ccw(b->corners[2], b->corners[0], points[e])>=0){
                asser(ccw(b->corners[0], b->corners[1], points[e])>=0);
                b->bucket.push_back(e);
            } else if(ccw(c->corners[1], c->corners[2], points[e])>=0 && ccw(c->corners[2], c->corners[0], points[e])>=0){
                asser(ccw(c->corners[0], c->corners[1], points[e])>=0);
                c->bucket.push_back(e);
            } else {
                asser(ccw(a->corners[1], a->corners[2], points[e])>=0 && ccw(a->corners[2], a->corners[0], points[e])>=0);
                asser(ccw(a->corners[0], a->corners[1], points[e])>=0);
                a->bucket.push_back(e);
            }
        }
        link_bucket(a);
        link_bucket(b);
        link_bucket(c);
        check_flips(a, 1);
        check_flips(b, 2);
        check_flips(c, 2);
    }

    Face* locateFace = 0;
    // compute delaunay triangulation
    vector<Face>& triangulate(vector<Point> const&p){
        points = p;
        int N = p.size();
        faces.reserve(3*N);
        random_shuffle(points.begin(), points.end());
		// start with super triangle that contains all points
        locateFace = new (getFreeFace()) Face(O, UL, UR);
        locateFace->bucket.resize(N);
        iota(locateFace->bucket.begin(), locateFace->bucket.end(), 0);
        location.resize(N, locateFace);
        // incremental construction
        for(int i=0;i<N;++i){
            Face* place = location[i];
            split(place, i);
        }
        vector<Face> retFaces;
        retFaces.reserve(faces.size()+1);
        unordered_map<Face*, Face*> decode;
        for(auto &e:faces){
            if(!decode.count(&e)){
                retFaces.push_back(e);
                decode[&e] = &(retFaces.back());
            }
        }
        for(auto &e:retFaces){
            for(int i=0;i<3;++i){
                e.adj[i] = decode[e.adj[i]];
            }
        }
        faces.swap(retFaces);
        locateFace = decode[locateFace];
        return faces;
    }
};

double circum_radius_sq(Delaunay::Face const& t){
    double ab = (t.corners[0].x * t.corners[0].x) + (t.corners[0].y * t.corners[0].y);
    double cd = (t.corners[1].x * t.corners[1].x) + (t.corners[1].y * t.corners[1].y);
    double ef = (t.corners[2].x * t.corners[2].x) + (t.corners[2].y * t.corners[2].y);
    double circum_x = (ab * (t.corners[2].y - t.corners[1].y) + cd * (t.corners[0].y - t.corners[2].y) + ef * (t.corners[1].y - t.corners[0].y)) / (t.corners[0].x * (t.corners[2].y - t.corners[1].y) + t.corners[1].x * (t.corners[0].y - t.corners[2].y) + t.corners[2].x * (t.corners[1].y - t.corners[0].y)) / 2.f;
    double circum_y = (ab * (t.corners[2].x - t.corners[1].x) + cd * (t.corners[0].x - t.corners[2].x) + ef * (t.corners[1].x - t.corners[0].x)) / (t.corners[0].y * (t.corners[2].x - t.corners[1].x) + t.corners[1].y * (t.corners[0].x - t.corners[2].x) + t.corners[2].y * (t.corners[1].x - t.corners[0].x)) / 2.f;
    double circum_radius = ((t.corners[0].x - circum_x) * (t.corners[0].x - circum_x)) + ((t.corners[0].y - circum_y) * (t.corners[0].y - circum_y));
    return circum_radius;
}

long long sqDist(Delaunay::Point const&a, Delaunay::Point const &b){
    if(a.x == -Delaunay::INF || b.x == -Delaunay::INF) return numeric_limits<long long>::max();// auxiliary points
    return (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y);
}

int lander(){
    cin.tie(0);cout.tie(0);ios_base::sync_with_stdio(false);
	int N;cin >> N;
	vector<Delaunay::Point> points(N);
	for(int i = 0; i < N; i++) {
        int x, y;
        cin >> x >> y;
        points[i] = Delaunay::Point(x, y);
	}

	Delaunay triangulation;
	vector<Delaunay::Face> const& triangles = triangulation.triangulate(points);
    double ret=0.0;
    for(auto const &e:triangles){
        if(!triangulation.is_infinite(e))
            ret = max(ret, circum_radius_sq(e));
    }
    cout << fixed << setprecision(10) << sqrt(ret);
	return 0;
}

int main(){
    #ifdef LOCAL_RUN
    freopen("in.txt", "r", stdin);
    #endif // LOCAL_RUN
    return lander();
    int N;cin >> N;
    vector<Delaunay::Point> ps(N);
    for(int i=0;i<N;++i){
        cin >> ps[i].x >> ps[i].y;
    }
    Delaunay dl;
    dl.triangulate(ps);
    for(auto const&e:dl.faces){
        if(!dl.is_infinite(e))
            cout << e << "\n";
    }
}


