
struct two_sat {
  int N; // number of variables
  vector<int> val; // assignment of x is at val[2x] and -x at val[2x+1]
  vector<char> valid; // changes made at time i are kept iff valid[i]
  vector<vector<int> > G; // graph of implications G[x][i] = y means (x -> y)
  
  two_sat(int N) : N(N) { // create a formula over N variables (numbered 1 to N)
    val.resize(2*N);
    G.resize(2*N);
  }

  int to_ind(int x) { // converts a signed variable index to its position in val[] and G[]
    return 2*(abs(x)-1) + (x<0);
  }

  // Add the implication: a -> b
  void add_implication(int a, int b) {
    G[to_ind(a)].push_back(to_ind(b));
  }

  // Add the or-clause: (a or b)
  void add_or(int a, int b) {
    add_implication(-a,b);
    add_implication(-b,a);
  }

  // Add condition: x is true
  void add_true(int x) {
    add_or(x,x);
  }
  
  int time(){
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

  bool solve() {
	fill(val.begin(), val.end(), 0);
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
