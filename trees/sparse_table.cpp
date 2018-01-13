// static range query in O(n log n + q)
template<typename op_class>
class Sparse_Table{
    using T = typename op_class::node_t;
    static int lg2(int a){
        return __builtin_clz(1) - __builtin_clz(a);
    }
    int n, lgn;
    vector<T> base, table;

    void build(int l, int r, int lg){
        if(l>r) return;
        if(lg<lgn){
            T*block_start = table.data()+2*n*lg;
            T cur = T(); // prefix table;
            for(int i=l;i<=r;++i){
                cur = op_class::op(cur, base[i]);
                block_start[i] = cur;
            }
            cur = T(); // suffix table;
            for(int i=r;i>=l;--i){
                cur = op_class::op(base[i], cur);
                block_start[n+i] = cur;
            }
        }
        if(l<r){
            int m = min(r+1, l+(1<<(lg-1)));
            build(l, m-1, lg-1);
            build(m, r, lg-1);
        }
    }
    void build_iterative(){
        for(int it=0;it<lgn;++it){
            T*block_start = table.data()+2*n*it;
            for(int l=0, r;l<n;l=r+1){
                r = min(n-1, l+(1<<it)-1);
                T cur = T(); // prefix table;
                for(int i=l;i<=r;++i){
                    cur = op_class::op(cur, base[i]);
                    block_start[i] = cur;
                }
                cur = T(); // suffix table;
                for(int i=r;i>=l;--i){
                    cur = op_class::op(base[i], cur);
                    block_start[n+i] = cur;
                }
            }
        }
    }
public:
    int query(int l, int r){
        if(l==r) return base[l];
        int lg = lg2(r^l);
        T*block_start = table.data()+2*n*lg;
        return op_class::op(block_start[n+l], block_start[r]);
    }
    Sparse_Table(vector<int> const&_base):n(_base.size()),lgn(lg2(2*n-1)), base(_base){
        table.resize(2*n*(lgn));
        build(0, n-1, lgn);
        //build_iterative();
    }
};
struct Op_Class{
    using node_t = int;
    static int op(int a, int b){return a|b;}
};
