//really fast iterative segment-tree implementation
template<class Segtree_Data>
struct Segtree{
	int N, height;
	vector<typename Segtree_Data::node_t> data;
	vector<typename Segtree_Data::update_t> lazy;
	Segtree(vector<typename Segtree_Data::node_t> const&base):N(base.size()), height(__builtin_clz(1)-__builtin_clz(N)), data(2*N, Segtree_Data::node_ne()), lazy(2*N, Segtree_Data::update_ne()){
		copy(base.begin(), base.end(), data.begin()+N);
		for(int i=N-1;i>=0;--i)
			data[i]=Segtree_Data::merge_nodes(data[i<<1], data[i<<1|1]);
	}
	Segtree(int n):N(n), height(__builtin_clz(1)-__builtin_clz(N)), data(2*N, Segtree_Data::node_ne()), lazy(2*N, Segtree_Data::update_ne()){
		for(int i=N-1;i>=0;--i)
			data[i]=Segtree_Data::merge_nodes(data[i<<1], data[i<<1|1]);
	}
	void local_update(int pos, typename Segtree_Data::update_t const&val){
		Segtree_Data::update_node(data[pos], val);
		if(pos<N) lazy[pos] = Segtree_Data::merge_lazy(lazy[pos], val);
	}
	void push(int pos){
		for(int s=height;s>0;--s){
			int i=pos>>s;
			if(lazy[i]!=Segtree_Data::update_ne()){
				local_update(i<<1, lazy[i]);
				local_update(i<<1|1, lazy[i]);
				lazy[i]=Segtree_Data::update_ne();
			}
		}
	}
	void re_path(int pos){
		while(pos>>=1) Segtree_Data::update_node(data[pos] = Segtree_Data::merge_nodes(data[pos<<1], data[pos<<1|1]), lazy[pos]);
	}
	void update(int l, int r, typename Segtree_Data::update_t const&val){
		int l2=l+=N, r2=r+=N;
		push(l2); push(r2-1);
		for(;l<r;l>>=1, r>>=1){
			if(l&1) local_update(l++, val);
			if(r&1) local_update(--r, val);
		}
		re_path(l2);re_path(r2-1);
	}
	typename Segtree_Data::node_t query(int l, int r){
		push(l+N); push(r+N-1);
		typename Segtree_Data::node_t retL=Segtree_Data::node_ne(), retR=Segtree_Data::node_ne();
		for(l+=N, r+=N;l<r;l>>=1, r>>=1){
			if(l&1)retL=Segtree_Data::merge_nodes(retL, data[l++]);
			if(r&1)retR=Segtree_Data::merge_nodes(data[--r], retR);
		}
		return Segtree_Data::merge_nodes(retL, retR);
	}
};
struct Segtreedata{
    typedef int node_t;
    typedef int update_t;
    static constexpr node_t node_ne() {return numeric_limits<int>::min();}
    static constexpr update_t update_ne(){return numeric_limits<int>::max();}
    static node_t merge_nodes(node_t const&left, node_t const&right){
        return node_t(max(left,right));
    }
    static void update_node(node_t &node, update_t const&update){
        if(node == node_ne() || update<node) node = update;
    }
    static update_t merge_lazy(update_t const&old, update_t const&update){
        return min(old, update);
    }
};