// LCA in <N, 1> for tree rooted at 'root'
// Algorithm by Baruch Schieber and Uzi Vishkin
struct Lowest_Common_Ancestor{
    uint32_t n, root;
    vector<uint32_t> parent, order, inv_order;
    vector<uint32_t> ascendant, inlabel, inlabel_link;
    Lowest_Common_Ancestor(vector<uint32_t> const&_parent, uint32_t _root):n(_parent.size()), root(_root), parent(_parent){build();}
    Lowest_Common_Ancestor(vector<int>const&_parent, int _root):n(_parent.size()), root(_root), parent(_parent.cbegin(),_parent.cend()){build();}
    uint32_t low_bit_mask(uint32_t x){
        return x&-x;
    }
    uint32_t high_bit_mask(uint32_t x){
        return x ? (~0u>>(__builtin_clz(x))>>1) : 0;
    }
    void preorder(){
        vector<vector<uint32_t> > ch(n);
        for(uint32_t i=0;i<n;++i) if(i!=root) ch[parent[i]].push_back(i);
        order.reserve(n);
        inv_order.resize(n);
        vector<uint32_t> s(1, root);
        s.reserve(n);
        while(!s.empty()){
            uint32_t cur = s.back(); s.pop_back();
            order.push_back(cur);
            inv_order[cur] = order.size();
            for(uint32_t const&c : ch[cur]){
                s.push_back(c);
            }
        }
    }
    void label(){
        ascendant.resize(n);
        inlabel = inv_order;
        inlabel_link.resize(n);
        // compute inlabel, i.e. the descendant with most trailing zeros
        for(uint32_t i=n-1;i>0;--i){
            uint32_t v = order[i], p = parent[v];
            if(low_bit_mask(inlabel[v])>low_bit_mask(inlabel[p]))
                inlabel[p] = inlabel[v];
        }
        // compute ascendant and links of inlabel-paths (link[v] = par[head[inlabel[v]]])
        ascendant[root] = low_bit_mask(inlabel[root]);
        inlabel_link[inlabel[root]-1] = root; // should be unused
        for(uint32_t i=1;i<n;++i){
            uint32_t v = order[i], p = parent[v];
            ascendant[v] = ascendant[p] | low_bit_mask(inlabel[v]);
            if(inlabel[v] != inlabel[p]) inlabel_link[inlabel[v]-1] = p;
        }
    }
    void build(){
        preorder();
        label();
    }
    uint32_t get_lca(uint32_t a, uint32_t b){
        uint32_t inlabel_a = inlabel[a], inlabel_b = inlabel[b];
        uint32_t ascendant_a = ascendant[a], ascendant_b = ascendant[b];
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
        return inv_order[a]<inv_order[b] ? a:b;
    }
	// solves RMQ by LCA
    template<typename T>
    static Lowest_Common_Ancestor build_rmq(vector<T> const&values){
        uint32_t n = values.size();
        vector<uint32_t> p(n), s;
        for(uint32_t i=0;i<values.size();++i){
            uint32_t last = i;
            while(!s.empty() && values[i]<values[s.back()]){
                p[last] = s.back();
                last = s.back();
                s.pop_back();
            }
            p[last] = i;
            s.push_back(i);
        }
        for(uint32_t i=1;i<s.size();++i){
            p[s[i]] = s[i-1];
        }
        p[s[0]] = ~0u;
        return Lowest_Common_Ancestor(p, s[0]);
    }
};