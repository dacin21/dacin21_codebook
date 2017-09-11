//really fast iterative segment-tree implementation
//uses std::function for subtree-merges and leaf-updates
#include <bits/stdc++.h>
using namespace std;

template<typename T, typename S>
struct segTree{
    vector<T> data;
    vector<S> lazy;
    int N, height;
    T idVal;
    S noLazy;
    function<void(T &, S const&)> nodeMod;
    function<T(T const&, T const&)> combine;
    function<S(S const&, S const&)> lazyCombine;
    segTree(vector<T> const&base, T const& identity_value, S const& no_lazy_value, function<void(T &, S const&)> const&node_update, function<T(T const&, T const&)> const&subtree_combine, function<S(S const&, S const&)>combine_lazy_values):N(base.size()), height(__builtin_clz(1)-__builtin_clz(N)), idVal(identity_value),noLazy(no_lazy_value), nodeMod(node_update), combine(subtree_combine), lazyCombine(combine_lazy_values){
        data.resize(2*N);
        lazy.resize(2*N, noLazy);
        copy(base.begin(), base.end(), data.begin()+N);
        for(int i=N-1;i>=0;--i){
            data[i]=combine(data[i<<1], data[i<<1|1]);
        }
    }
    void localUpdate(int pos, S const&val){
        nodeMod(data[pos], val);
        if(pos<N) lazy[pos] = lazyCombine(lazy[pos], val);
    }
    void push(int pos){
        for(int s=height;s>0;--s){
            int i=pos>>s;
            if(lazy[i]!=noLazy){
                localUpdate(i<<1, lazy[i]);
                localUpdate(i<<1|1, lazy[i]);
                lazy[i]=noLazy;
            }
        }
    }
    void rePath(int pos){
        while(pos>>=1) nodeMod(data[pos] = combine(data[pos<<1], data[pos<<1|1]), lazy[pos]);
    }
    void update(int l, int r, S const&val){
        int l2=l+=N, r2=r+=N;
        push(l2); push(r2-1);
        for(;l<r;l>>=1, r>>=1){
            if(l&1) localUpdate(l++, val);
            if(r&1) localUpdate(--r, val);
        }
        rePath(l2);rePath(r2-1);
    }
    T query(int l, int r){
        push(l+N); push(r+N-1);
        T retL=idVal, retR=idVal;
        for(l+=N, r+=N;l<r;l>>=1, r>>=1){
            if(l&1)retL=combine(retL, data[l++]);
            if(r&1)retR=combine(data[--r], retR);
        }
        return combine(retL, retR);
    }
};

