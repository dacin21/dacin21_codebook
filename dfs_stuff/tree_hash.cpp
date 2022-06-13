/*
 *  Hash of a tree
 *  Can be used for isomorphism checking
 *  variadic templates are great!!
 *
 *  2022-06: fixed a bug where all levels had the same salt.
 *
 */

typedef uint32_t H_t;
typedef uint64_t HH_t;
template<H_t p>
struct Sethash{
    H_t val, q;
    Sethash():val(1), q(0){}
    Sethash(H_t _q):val(1), q(_q%p){}
    Sethash& operator+=(Sethash const&o){
        val = val*((HH_t)o.val + q)%p;
        return *this;
    }
    // comparing just val should be enough?
    bool operator==(Sethash const&o)const{
        return q==o.q && val==o.val;
    }
    bool operator<(Sethash const&o)const{
        return make_pair(q, val) < make_pair(o.q, o.val);
    }
    friend ostream& operator<<(ostream&o, Sethash const&h){
        return o << h.val << " (" << h.q << ")";
    }
};
class Treehash{
public:
    tuple<Sethash<1000000007u>, Sethash<999999937u>, Sethash<2147483587u>> val;
    Treehash(unsigned int level = -1){salt_tuple(val, level+10);}
private:
    template<size_t I = 0, typename... Tp>
    inline typename enable_if<I == sizeof...(Tp), void>::type
    joinhash(tuple<Tp...>const&){}
    template<size_t I = 0, typename... Tp>
    inline typename enable_if<I < sizeof...(Tp), void>::type
    joinhash(tuple<Tp...>const& t) {
        get<I>(val)+=get<I>(t);
        joinhash<I + 1, Tp...>(t);
    }
    template<size_t I = 0, typename... Tp>
    static inline typename enable_if<I == sizeof...(Tp), void>::type
    print(ostream&, tuple<Tp...>const&){}
    template<size_t I = 0, typename... Tp>
    static inline typename enable_if<I < sizeof...(Tp) && I!=0, void>::type
    print(ostream&o, tuple<Tp...>const& t) {
        o << ", " << get<I>(t);
        print<I + 1, Tp...>(o, t);
    }
    template<size_t I = 0, typename... Tp>
    static inline typename enable_if<I < sizeof...(Tp) && I==0, void>::type
    print(ostream&o, tuple<Tp...>const& t) {
        o << get<I>(t);
        print<I + 1>(o, t);
    }
    template<size_t I = 0, typename... Tp>
    static inline typename enable_if<I == sizeof...(Tp), size_t>::type
    intify(tuple<Tp...>const&){return 0;}
    template<size_t I = 0, typename... Tp>
    static inline typename enable_if<I < sizeof...(Tp), size_t>::type
    intify(tuple<Tp...>const& t){
        return intify<I + 1>(t)*31 + get<I>(t).val;
    }
    template<size_t I = 0, typename... Tp>
    static inline typename enable_if<I == sizeof...(Tp), void>::type
    salt_tuple(tuple<Tp...>&, unsigned int const&){}
    template<size_t I = 0, typename... Tp>
    static inline typename enable_if<I < sizeof...(Tp), void>::type
    salt_tuple(tuple<Tp...>& t, unsigned int const&level){
        get<I>(t) = get_salt((sizeof...(Tp))*level+I);
        salt_tuple<I+1>(t, level);
    }

    static mt19937 rng;
    static vector<H_t> salt_c;
    static H_t get_salt(unsigned int l){
        if(salt_c.size() <= l){
            while(salt_c.size()<l+20){
                salt_c.push_back(uniform_int_distribution<H_t>(0, ~0u>>1)(rng));
            }
        }
        return salt_c[l];
    }
public:
    Treehash& operator+=(Treehash const&o){
        joinhash(o.val);
        return *this;
    }
    bool operator==(Treehash const&o)const{
        return val==o.val;
    }
    bool operator!=(Treehash const&o)const{
        return !operator==(o);
    }
    bool operator<(Treehash const&o) const {
        return val < o.val;
    }
    friend ostream& operator<<(ostream&o, Treehash const&h){
        o << "Treehash: ";
        print(o, h.val);
        return o;
    }
    size_t to_long()const{return intify(val);}
};
mt19937 Treehash::rng;
vector<H_t> Treehash::salt_c;

namespace std{
    template <>
    struct hash<Treehash>{
        size_t operator()(Treehash const&h) const {
            return h.to_long();
        }
    };
}

using Hash = Treehash;

Hash root_hash(vector<vector<int> > const &G, int root){
    int N=G.size(), a, b, d;
    vector<Hash> h(N, 0);
    stack<tuple<int, int, int> > s;
    s.emplace(root, -1, 0);
    while(!s.empty()){
        tie(a, b, d) = s.top();s.pop();
        if(a<0){
            a=-1-a;
            if(~b) h[b]+=h[a];
        } else {
            h[a] = Hash(d);
            s.emplace(-1-a, b, d);
            for(int const&i:G[a]){
                if(i!=b){
                    s.emplace(i, a, d+1);
    }   }   }   }
    return h[root];
}
Hash de_pair(Hash const&a, Hash const&b){
    Hash ret(-2);
    ret += a;
    ret += b;
    return ret;
}
Hash hashify(vector<vector<int> > const&G){
    if(G.empty()) return de_pair(Hash(), Hash());
    queue<int> q;
    vector<int> pre(G.size(), -1);
    q.emplace(0);
    int a;
    while(!q.empty()){
        a=q.front();q.pop();
        for(int const&e:G[a]){
            if(e!=pre[a]){
                q.emplace(e);
                pre[e]=a;
    }   }   }
    q.emplace(a);
    fill(pre.begin(), pre.end(), -1);
    while(!q.empty()){
        a=q.front();q.pop();
        for(int const&e:G[a]){
            if(e!=pre[a]){
                q.emplace(e);
                pre[e]=a;
    }   }   }
    int len=0;
    for(int i=a;i!=-1;i=pre[i])++len;
    for(int i=1;i<(len+1)/2;++i)a=pre[a];
    if(len%2) return de_pair(root_hash(G, a), Hash());
    return de_pair(root_hash(G, a), root_hash(G, pre[a]));
}
Hash hashify(vector<pair<int, int> > const&E){
    int N=E.size()+1;
    vector<vector<int> > G(N);
    for(pair<int, int>const&e:E){
        G[e.first].emplace_back(e.second);
        G[e.second].emplace_back(e.first);
    }
    return hashify(G);
}
