// RMQ in <N, 1> for tree rooted at 'root'
// LCA is specialized for cartesian trees (<=2 children per node)
// Algorithm by Baruch Schieber and Uzi Vishkin
struct Lowest_Common_Ancestor{
    uint32_t n, root;
    vector<uint32_t> parent, order;
    struct AI{
        uint32_t ascendant, inlabel;
    };
    vector<AI> ai;
    vector<int> v;
    vector<uint32_t> inlabel_link;

    Lowest_Common_Ancestor(){}
    Lowest_Common_Ancestor(vector<int> const&v_):n(v_.size()), parent(n) {
        v.reserve(n+1);
        v.insert(v.end(), v_.begin(), v_.end());
        // add sentinel
        v.push_back(numeric_limits<int>::min());
        build();
    }
    __attribute__((always_inline))
    inline uint32_t low_bit_mask(uint32_t x){
        return x&-x;
    }
    __attribute__((always_inline))
    inline uint32_t high_bit_mask(uint32_t x){
        return x ? ((~0u>>1)>>__builtin_clz(x)) : 0;
    }
    __attribute__((always_inline))
    inline void label(){
        inlabel_link.resize(n);
        // compute inlabel, i.e. the descendant with most trailing zeros
        for(uint32_t i=n-1;i>0;--i){
            uint32_t v = order[i], p = parent[v];
            if(low_bit_mask(ai[v].inlabel)>low_bit_mask(ai[p].inlabel))
                ai[p].inlabel = ai[v].inlabel;
        }
        // compute ascendant and links of inlabel-paths (link[v] = par[head[inlabel[v]]])
        ai[root].ascendant = low_bit_mask(ai[root].inlabel);
        inlabel_link[ai[root].inlabel-1] = root; // should be unused
        for(uint32_t i=1;i<n;++i){
            uint32_t v = order[i], p = parent[v];
            ai[v].ascendant = ai[p].ascendant | low_bit_mask(ai[v].inlabel);
            if(ai[v].inlabel != ai[p].inlabel) inlabel_link[ai[v].inlabel-1] = p;
        }
    }
    void build(){
        ai.resize(n);
        uint32_t j = n;
        order.resize(n);
        vector<uint32_t> s(1, n);
        s.reserve(n+1);
        for(uint32_t i=n-1;~i;--i){
            uint32_t last = i;
            while(v[i]<v[s.back()]){
                const int u = s.back();
                order[--j] = u;
                ai[u].inlabel = j;
                parent[last] = u;
                last = u;
                s.pop_back();
            }
            parent[last] = i;
            s.push_back(i);
        }
        for(uint32_t i=s.size()-1; i>1; --i){
            const int u = s[i];
            order[--j] = u;
            ai[u].inlabel = j;
            parent[u] = s[i-1];
        }
        root = s[1];
        parent[root] = ~0u;
        order[0] = root;
        ai[root].inlabel =  0;
        
        label();
    }
    uint32_t query(uint32_t a, uint32_t b){
        uint32_t inlabel_a = ai[a].inlabel, ascendant_a = ai[a].ascendant;
        uint32_t inlabel_b = ai[b].inlabel, ascendant_b = ai[b].ascendant;
        uint32_t i_mask = high_bit_mask(inlabel_a^inlabel_b);
        uint32_t j_mask = low_bit_mask(((ascendant_a&ascendant_b)&~i_mask));
        uint32_t inlabel_lca = (inlabel_a==inlabel_b) ? inlabel_a : ((inlabel_a&~((j_mask<<1)-1))|j_mask);
        if(inlabel_a!=inlabel_lca){
            uint32_t k_mask = high_bit_mask(ascendant_a&(j_mask-1));
            uint32_t inlabel_u = (inlabel_a&~k_mask)|(k_mask+1);
            a = inlabel_link[inlabel_u-1];
        }
        if(inlabel_b!=inlabel_lca){
            uint32_t k_mask = high_bit_mask(ascendant_b&(j_mask-1));
            uint32_t inlabel_w = (inlabel_b&~k_mask)|(k_mask+1);
            b = inlabel_link[inlabel_w-1];
        }
        //return v[a]<v[b] ? a:b; // index
        return min(v[a], v[b]);
    }
};
