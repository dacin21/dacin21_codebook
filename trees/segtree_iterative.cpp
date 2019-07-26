//really fast iterative segment-tree implementation
#include <bits/stdc++.h>
using namespace std;

//uses std::function for subtree-merges and leaf-updates
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

// nicer implementation with struct
template<class Segtree_Data>
struct Segment_Tree{
	using T = typename Segtree_Data::node_t;
	int n;
	vector<T>data;
	Segment_Tree(int _n):n(_n), data(2*n, Segtree_Data::node_ne()){
		for(int i=n-1;i>=0;--i) data[i] = Segtree_Data::merge_nodes(data[i<<1], data[i<<1|1]);
	}
	Segment_Tree(vector<T> const&base):n(base.size()), data(2*n, Segtree_Data::node_ne()){
		copy(base.begin(), base.end(), data.begin()+n);
		for(int i=n-1;i>=0;--i) data[i] = Segtree_Data::merge_nodes(data[i<<1], data[i<<1|1]);
	}
	void update(int pos, typename Segtree_Data::update_t const&val){
		for(Segtree_Data::update_node(data[pos+=n], val);pos>>=1;){
			data[pos] = Segtree_Data::merge_nodes(data[pos<<1], data[pos<<1|1]);
		}
	}
	T query(int l, int r)const{
		T retL = Segtree_Data::node_ne(), retR = Segtree_Data::node_ne();
		for(l+=n, r+=n;l<r;l>>=1, r>>=1){
			if(l&1) retL = Segtree_Data::merge_nodes(retL, data[l++]);
			if(r&1) retR = Segtree_Data::merge_nodes(data[--r], retR);
		}
		return Segtree_Data::merge_nodes(retL, retR);
	}
};
struct Segtreedata{
    typedef int node_t;
    typedef int update_t;
    static constexpr node_t node_ne() {return numeric_limits<int>::min();}
    static node_t merge_nodes(node_t const&left, node_t const&right){
        return node_t(max(left,right));
    }
    static void update_node(node_t &node, update_t const&update){
        if(node == node_ne() || update<node) node = update;
    }
};
