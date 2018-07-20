
/*
 *  Recursive algorithm for maximum clique
 *  Takes ~1s in the worst case for n ~ 53 (random with p=0.95)
 *  Faster on sparse graph (p <= 0.85)
 */
template<size_t max_n>
class Clique{
    using bits = bitset<max_n>;
    bits MASK,ZERO,ans;
    const bits *e;
    int N;
    // int64_t calls;
    void bk_init(){
        ans=ZERO;
        MASK=ZERO;
        MASK.flip();
        calls = 0;
    }
    void bk(bits use, bits can_start,bits can_other){
        // ++calls;
        if(can_start.none()&&can_other.none()){
            if(use.count()>ans.count())ans=use;
            return;
        }
        bits r=can_start;
        bool fi=1;
        for(int i=0;i<N;++i){
            if(r[i]){
                if(fi){
                    fi=0;
                    r&=e[i]^MASK;
                }
                use[i]=1;
                bk(use,can_start&e[i],can_other&e[i]);
                use[i]=0;
                can_start[i]=0;
                can_other[i]=1;
            }
        }
    }
    static Clique& get(){
        static Clique c;
        return c;
    }
public:
    static bits find_clique(bits const*g, const int&n){
        Clique& c = get();
        c.e = g;
        c.N = n;
        c.bk_init();
        bits me;
        c.bk(me,c.MASK,c.ZERO);
        // cerr << "Calls: " << c.calls << "\n";
        // c.calls = 0;
        return c.ans;
    }
    static bits find_clique(vector<bits> const&g){
        return find_clique(g.data(), g.size());
    }
    static bits find_clique(array<bits, max_n> const&g, const int&n){
        return find_clique(g.data(), n);
    }
};