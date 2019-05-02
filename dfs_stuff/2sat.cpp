// 2-sat in linear time via backtracking.
class Two_Sat {
    int N; // number of variables
    vector<int> val; // assignment of x is at val[2x] and -x at val[2x+1]
    vector<char> valid; // changes made at time i are kept iff valid[i]
    vector<vector<int> > G; // graph of implications G[x][i] = y means (x -> y)

    Two_Sat(int N_) : N(N_) { // create a formula over N variables (numbered 1 to N)
        G.resize(2*N);
    }

    int add_variable() {
        G.emplace_back();
        G.emplace_back();
        return N++;
    }

private:
    // converts a signed variable index to its position in val[] and G[]
    int to_ind(int x) {
        return 2*(abs(x)-1) + (x<0);
    }

    // Add a directed edge to the graph.
    // You most likely do not want to call this yourself!
    void add_edge(int a, int b) {
        G[to_ind(a)].push_back(to_ind(b));
    }

    int time() {
        return valid.size()-1;
    }

    bool dfs(int x) {
        if(valid[abs(val[x])]) return val[x]>0;
        val[x] = time();
        val[x^1] = -time();
        for(int e:G[x])
            if(!dfs(e))
                return false;
        return true;
    }

public:
    // Add the or-clause: (a or b)
    void add_or(int a, int b) {
        add_edge(-a,b);
        add_edge(-b,a);
    }

    // Add the implication: a -> b
    void add_implication(int a, int b) {
        add_or(-a, b);
    }

    // Add condition: x is true
    void add_true(int x) {
        add_or(x,x);
    }

    // At most one with linear number of clauses
    template<typename T>
    void add_at_most_one(T vars) {
        if(vars.begin() == vars.end()) return;
        int last = *vars.begin();
        int cur = 0;
        for(int const&e:vars){
            if(e == last) continue;
            if(cur == 0) cur = e;
            else {
                add_or(-cur, -e);
                int new_cur = add_variable();
                cur = add_implication(cur, new_cur);
                add_implication(e, new_cur);
                cur = new_cur;
            }
        }
        if(cur != 0){
            add_or(-cur, -last);
        }
    }

    bool solve() {
        val.assign(2*n, 0);
        valid.assign(1, 0);
        for(int i=0; i<val.size(); i+=2) {
            if(!valid[abs(val[i])]) {
                valid.push_back(1);
                if(!dfs(i)) {
                    valid.back()=0;
                    valid.push_back(1);
                    if(!dfs(i+1)) return false;
                }
            }
        }
        return true;
    }
};
