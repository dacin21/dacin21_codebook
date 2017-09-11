//really fast iterative segment-tree implementation
//uses std::function for subtree-merges and leaf-updates
#include <bits/stdc++.h>
using namespace std;

template<typename T>
struct segTree{
    vector<T> data;
    int N;
    T idVal;
    function<T(T const&, T const&)> combine;
    function<void(T &, T const&)> leafMod;
    segTree(vector<T> const&base, T const& identity_value, function<void(T &, T const&)> const&leaf_update, function<T(T const&, T const&)> const&subtree_combine):N(base.size()), idVal(identity_value), combine(subtree_combine), leafMod(leaf_update){
        data.resize(2*N);
        copy(base.begin(), base.end(), data.begin()+N);
        for(int i=N-1;i>=0;--i){
            data[i]=combine(data[i<<1], data[i<<1|1]);
        }
    }
    void update(int pos, T val){
        for(leafMod(data[pos+=N], val);pos>>=1;){
            data[pos]=combine(data[pos<<1], data[pos<<1|1]);
        }
    }
    T query(int l, int r){
        T retL=idVal, retR=idVal;
        for(l+=N, r+=N;l<r;l>>=1, r>>=1){
            if(l&1)retL=combine(retL, data[l++]);
            if(r&1)retR=combine(data[--r], retR);
        }
        return combine(retL, retR);
    }
};
