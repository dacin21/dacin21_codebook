/*
 * 2-sat in linear time without SCC
 * variables from 0 to n-1, use ~x for negation.
 */
class Two_Sat {
    int n;
    vector<char> val; // assignment: 1 or 0 where x is at val[2x] and ~x at val[2x+1]
    vector<int> q; // queue for bfs
    vector<vector<int> > g; // graph of implications G[x][i] = y means (x -> y)
    int add_variable() {
        g.emplace_back();
        g.emplace_back();
        return n++;
    }

    // x / ~x (0 <= x < n) to indices in 0, ..., 2n-1
    int to_ind(int x) {
        return x<0 ? 2*(~x)+1 : 2*x;
    }
    // do not use! use add_implication instead!
    void add_edge(int a, int b) {
        g[to_ind(a)].push_back(to_ind(b));
    }

    bool bfs(int x) {
        q[0] = x;
        val[x] = 1;
        val[x^1] = 0;
        auto it = q.begin();
        auto it_end = next(it);
        auto fail = [&](){
            for(it = q.begin(); it != it_end; ++it){
                val[*it] = val[*it^1] = -1;
            }
            return false;
        };
        while(it != it_end){
            const int u = *(it++);
            for(int const e : g[u]){
                switch(val[e]){
                    case 0:
                        return fail();
                    case -1:
                        val[e] = 1;
                        val[e^1] = 0;
                        *(it_end++) = e;
                        break;
                    default:
                        break;
                }
            }
        }
        return true;
    }

public:
    Two_Sat(int n_) : n(n_), g(2*n){}
    // Add the or-clause: (a or b)
    void add_or(int a, int b) {
        add_edge(~a,b);
        add_edge(~b,a);
    }
    // Add condition: (not a) or (not b)
    void add_not_both(int a, int b) {
        add_or(~a, ~b);
    }
    // Add the implication: a -> b
    void add_implication(int a, int b) {
        add_or(~a, b);
    }
    // Add condition: x is true
    void add_true(int x) {
        add_or(x,x);
    }
    // Add condition: x is false
    void add_false(int x) {
        add_true(~x);
    }

    // At most one with linear number of clauses
    template<typename T>
    void add_at_most_one(T vars) {
        if(vars.begin() == vars.end()) return;
        auto it = vars.begin();
        int last = -1;
        int cur = *it;
        while(++it != vars.end()){
            if(last != -1){
                int new_cur = add_variable();
                add_implication(cur, new_cur);
                add_implication(last, new_cur);
                cur = new_cur;
            }
            last = cur;
            cur = *it;
            add_not_both(last, cur);
        }
    }

    bool solve() {
        val.assign(2*n, -1);
        q.assign(2*n+1, -1);
        for(int i=0; i<(int)val.size(); i+=2) {
            if(val[i] == -1){
                if(!bfs(i)){
                    if(!bfs(i+1)){
                        return false;
                    }
                }
            }
        }
        return true;
    }
    bool get_val(int x){
        assert(0 <= x && x < n);
        assert(val[2*x] != -1);
        return val[2*x];
    }
};
