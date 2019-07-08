// Sets on a finite universe [1, U] that can be split or merged in O(log U)
// Also known a sparse segment trees or binary tries.
// This implementation only uses O(n) memory (instead of O(n log U)).

template<typename traits>
struct Mergable_Set{
    using element_t = typename traits::element_t;
    static_assert(is_integral<element_t>::value);
    static constexpr element_t L = traits::min(), R = traits::max();
    static constexpr element_t range = R - L; // use constexpr to test for overflows

    struct Node{
        Node() : l(nullptr), r(nullptr), size(0) {}
        Node(element_t x) : l(nullptr), r(nullptr), a(x), b(x), size(1) {}
        Node(element_t a_, element_t b_) : l(nullptr), r(nullptr), a(a_), b(b_), size(0) {}
        Node(Node*l_, Node*r_, element_t const&a_, element_t const&b_) : l(l_), r(r_), a(a_), b(b_), size(l?l->size():0 + r?r->size:0) {}

        void recalc(){
            assert(a != b);
            size = (l ? l->size : 0) + (r ? r->size : 0);
        }

        Node *l, *r;
        element_t a, b;
        int size;
    };
    struct Set_Handle{
        Node* node = nullptr;
        int size() const {
            return node ? node->size : 0;
        }
    };

    Mergable_Set(int max_total_elements) : pool_size(0), pool(2*max_total_elements) {}

    Set_Handle create_empty_set(){
        return Set_Handle{};
    }
    Set_Handle create_singleton_set(element_t const&x){
        assert(L <= x && x <= R);
        Set_Handle ret{get_free_node(x)};
        return ret;
    }
    Set_Handle merge(Set_Handle u, Set_Handle v){
        Set_Handle ret{merge(u.node, v.node, L, R)};
        return ret;
    }
    pair<Set_Handle, Set_Handle> split(Set_Handle u, element_t x){
        if(!u.node) return {Set_Handle{}, Set_Handle{}};
        auto tmp = split_rec(u.node, x);
        return {Set_Handle{tmp.first}, Set_Handle{tmp.second}};
    }
    pair<Set_Handle, Set_Handle> split_pos(Set_Handle u, int pos){
        assert(0 <= pos && pos <= u.size());
        if(!u.node) return {Set_Handle{}, Set_Handle{}};
        auto tmp = split_pos_rec(u.node, pos);
        return {Set_Handle{tmp.first}, Set_Handle{tmp.second}};
    }
    int lower_count(Set_Handle u, element_t const x){
        int ret = 0;
        for(Node*cur = u;;){
            if(x <= cur->a){
                return ret;
            } else if(x > cur->b){
                return ret + cur->size;
            }
            const int m = cur->a + (cur->b - cur->a)/2;
            if(x <= m){
                cur = cur->l;
            } else {
                ret += cur->l->size;
                cur = cur->r;
            }
        }
    }
    element_t at(Set_Handle u, int pos){
        assert(0 <= pos && pos < u.size());
        for(Node*cur = u.node;;){
            if(cur->a == cur->b) return cur->a;
            const int left_size = cur->l->size;
            if(pos < left_size){
                cur = cur->l;
            } else {
                pos -= left_size;
                cur = cur->r;
            }
        }
    }

private:
    template<typename... Args>
    Node* get_free_node(Args... args){
        const int index = [this](){
            if(free_nodes.empty()){
                assert(pool_size < (int) pool.size());
                return pool_size++;
            }
            int ret = free_nodes.back();
            free_nodes.pop_back();
            return ret;
        }();
        return new (pool.data() + index) Node(std::forward<Args>(args)...);
    }
    void free_node(Node* node){
        free_nodes.push_back(node - pool.data());
    }
    Node* merge(Node*u, Node*v, element_t a, element_t b){
        if(!u) return v;
        if(!v) return u;
        auto m = [&a, &b](){return a + (b-a)/2;};
        for(;;){
            if(a!=b && m() < u->a && m() < v->a){
                a = m()+1;
            } else if (a!=b && u->b <= m() && v->b <= m()){
                b = m();
            } else break;
        }
        if(a==b){
            u->size += v->size;
            free_node(v);
            return u;
        }
        auto decompose = [&](Node*w, Node* &par) -> pair<Node*, Node*>{
            if(w->a == a && w->b == b){
                if(par) free_node(par);
                par = w;
                return {w->l, w->r};
            }
            assert(w->a > m() || w->b <= m());
            if(w->b <= m()) return {w, nullptr};
            else return {nullptr, w};
        };
        Node* ret = get_free_node(a, b);
        auto sub_u = decompose(u, ret), sub_v = decompose(v, ret);
        ret->l = merge(sub_u.first, sub_v.first, a, m());
        ret->r = merge(sub_u.second, sub_v.second, m()+1, b);
        assert(ret->l != nullptr && ret->r != nullptr);
        ret->recalc();
        return ret;
    }
    template<typename T, typename Fun>
    pair<Node*, Node*> split_recurse(Node*u, T const&pivot_l, T const& pivot_r, Fun split_fun, bool go_left){
        if(go_left){
            auto tmp = (this->*split_fun)(u->l, pivot_l);
            if(tmp.second) {
                u->l = tmp.second;
                u->recalc();
                tmp.second = u;
            } else {
                tmp.second = u->r;
                free_node(u);
            }
            return tmp;
        } else {
            auto tmp = (this->*split_fun)(u->r, pivot_r);
            if(tmp.first){
                u->r = tmp.first;
                u->recalc();
                tmp.first = u;
            } else {
                tmp.first = u->l;
                free_node(u);
            }
            return tmp;
        }
    }
    pair<Node*, Node*> split_rec(Node*u, element_t const&x){
        assert(u != nullptr);
        if(x <= u->a) return {nullptr, u};
        if(x > u->b) return {u, nullptr};
        const element_t m = u->a + (u->b - u->a)/2;
        return split_recurse(u, x, x, &Mergable_Set::split_rec, x <= m);
    }

    pair<Node*, Node*> split_pos_rec(Node*u, int const&pos){
        assert(u != nullptr);
        assert(0 <= pos && pos <= u->size);
        if(pos == 0) return {nullptr, u};
        if(pos == u->size) return {u, nullptr};
        if(u->a == u->b){
            Node*v = get_free_node(u->a, u->b);
            v->size = pos;
            u->size -= pos;
            return {v, u};
        }
        const int left_size = u->l->size;
        return split_recurse(u, pos, pos - left_size, &Mergable_Set::split_pos_rec, pos < left_size);
    }

    int pool_size;
    vector<Node> pool;
    vector<int> free_nodes;
};

struct traits{
    using element_t = int;
    static constexpr int min(){ return 0; }
    static constexpr int max(){ return 1e9; }
};


// range sorting queries
signed main(){
    #ifdef LOCAL_RUN
    freopen("in.txt", "r", stdin);
    #endif // LOCAL_RUN
    int n, q;
    cin >> n >> q;

    Mergable_Set<traits> s(n);

    using Sorted_Range = pair<decltype(s)::Set_Handle, bool>;

    // use pair<bool, Set_Handle> for sorted ranges.

    auto split = [&](Sorted_Range range, int pos){
        if(range.second) {
            auto tmp = s.split_pos(range.first, range.first.size() - pos);
            return make_pair(Sorted_Range(tmp.second, true), Sorted_Range(tmp.first, true));
        } else {
            auto tmp = s.split_pos(range.first, pos);
            return make_pair(Sorted_Range(tmp.first, false), Sorted_Range(tmp.second, false));
        }
    };
    auto splitoff = [&](Sorted_Range &range, int pos){
        //cerr << "splitoff: " << range.first.size() << " " << pos << " " << range.second << "\n";
        auto tmp = split(range, pos);
        range = tmp.second;
        //cerr << tmp.first.first.size() << " / " << tmp.second.first.size() << "\n";
        return tmp.first;
    };
    auto at = [&](Sorted_Range &range, int pos){
        if(range.second){
            pos = range.first.size()-pos-1;
            return s.at(range.first, pos);
        } else {
            return s.at(range.first, pos);
        }
    };

    map<int, Sorted_Range> data; // right endpoint, elements
    vector<int> v(n);
    vector<pair<int, int> > sortv;
    copy_n(istream_iterator<int>(cin), n, v.begin());
    for(int i=0;i<n;++i) sortv.emplace_back(v[i], i);
    sort(sortv.begin(), sortv.end());
    for(int i=0;i<n;++i){
        int val = (int)(lower_bound(sortv.begin(), sortv.end(), make_pair(v[i], i))-sortv.begin());
        data.emplace(i, Sorted_Range(s.create_singleton_set(val), false));
    }


    while(q--){
        //runprint();
        // c=1 : increasing, c=2 : decreasing;
        int a, b, c;
        assert(cin >> c >> a >> b);
        assert(0<c&&c<3 && 0<a&&a<=b&&b<=n);
        --a, --b;
        auto it1 = data.lower_bound(a);
        int s1 = it1->second.first.size();
        if(it1->first-a+1!=s1) data.emplace(a-1, splitoff(it1->second, s1-(it1->first-a+1)));
        auto it2 = data.lower_bound(b);
        int s2 = it2->second.first.size();
        if(it2->first!=b){
            data.emplace(b, splitoff(it2->second, s2-(it2->first-b)));
            if(it1==it2) --it1;
        } else ++it2;
        Sorted_Range res(s.create_empty_set(), c==2);
        for(auto it = it1;it!=it2;++it) res.first = s.merge(res.first, it->second.first);
        data.erase(it1, it2);
        data.emplace(b, res);
    }
    for(int i=0;i<n;++i){
        auto it = data.lower_bound(i);
        cout << sortv[at(it->second, it->second.first.size()-1-it->first+i)].first << " ";
    }
    cout << "\n";
}
