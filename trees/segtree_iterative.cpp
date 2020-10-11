//really fast iterative segment-tree implementation

template<class Segtree_Data>
struct Segment_Tree{
	using T = typename Segtree_Data::node_t;
	int n;
	vector<T>data;
	Segment_Tree(int n_) : n(n_), data(2*n, Segtree_Data::node_init()){
		for(int i=n-1;i>=0;--i) data[i] = Segtree_Data::merge_nodes(data[i<<1], data[i<<1|1]);
	}
	Segment_Tree(vector<T> const&base) :  n(base.size()), data(2*n, Segtree_Data::node_init()){
		copy(base.begin(), base.end(), data.begin()+n);
		for(int i=n-1;i>=0;--i) data[i] = Segtree_Data::merge_nodes(data[i<<1], data[i<<1|1]);
	}
	void update(int pos, typename Segtree_Data::update_t const&val){
		for(Segtree_Data::update_node(data[pos+=n], val); pos>>=1;){
			data[pos] = Segtree_Data::merge_nodes(data[pos<<1], data[pos<<1|1]);
		}
	}
    void u(int pos, typename Segtree_Data::update_t const&val){
        return update(pos,val);
    }
	T query(int l, int r) const {
		T retL = Segtree_Data::node_ne(), retR = Segtree_Data::node_ne();
		for(l+=n, r+=n; l<r; l>>=1, r>>=1){
			if(l&1) retL = Segtree_Data::merge_nodes(retL, data[l++]);
			if(r&1) retR = Segtree_Data::merge_nodes(data[--r], retR);
		}
		return Segtree_Data::merge_nodes(retL, retR);
	}
    T q(int l, int r) const {
        return query(l, r);
    }
};
struct Segtreedata{
    typedef int node_t;
    typedef int update_t;
    static constexpr node_t node_init() {return numeric_limits<int>::min();}
    static constexpr node_t node_ne() {return numeric_limits<int>::min();}
    static node_t merge_nodes(node_t const&left, node_t const&right){
        return node_t(max(left,right));
    }
    static void update_node(node_t &node, update_t const&update){
        if(node == node_init() || update<node) node = update;
    }
};
