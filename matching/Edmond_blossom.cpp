/*
   Implementation of Edmonds blossom algorithm for the general maximum matching problem
   based on the implementation of the boost_graph_library.
   Runs in O(N M log(N)) worst case and O(N M alpha(N)) average case time.
   Uses O(N+M) memory.
*/

#include <bits/stdc++.h>
using namespace std;

struct Dsu{
  vector<int> p;
  void init(int n){
    p.resize(n);
    iota(p.begin(), p.end(), 0);
  }
  Dsu(int n){
    init(n);
  }
  Dsu(){}
  int find_set(int i){ // iterative path compression
    int j=i;
    int ret=i;
    while(ret!=p[ret]){
      ret=p[ret];
    }
    i=p[j];
    while(i!=j){
      p[j]=ret;
      j=i;
      i=p[i];
    }
    return ret;
  }
  void unite_set(int a, int b){
    p[find_set(a)] = find_set(b);
  }
};

struct Matcher{
  vector<vector<int> > graph; // adjacency list
  vector<int> mate, origin, pred; //vertex vectors
  deque<int> augmenting_path; // vertex list
  vector<pair<int, int> > bridge; // vertex pair vector
  vector<int> vertex_state; // state vectors 0=unreached, 1=odd, 2=even
  vector<int> ancestor_v, ancestor_w; // vertex vectors
  vector<pair<int, int> > even_edges; // edge vector
  Dsu dsu;
  Matcher(int N){
    graph.resize(N);
    mate.resize(N, -1);
    origin.resize(N);
    pred.resize(N);
    bridge.resize(N);
    vertex_state.resize(N);
    ancestor_v.resize(N);
    ancestor_w.resize(N);
  }

  void add_edge(int a, int b){
    assert(a>=0 && a < (int)graph.size());
    assert(b>=0 && b < (int)graph.size());
    graph[a].emplace_back(b);
    graph[b].emplace_back(a);
  }

  void clearMatching(){
    fill(mate.begin(), mate.end(), -1);
  }

  ///parent in DFS tree
  int parent(int vertex){
    if(vertex_state[vertex]==2 && mate[vertex]!=-1)
      return mate[vertex];
    else if(vertex_state[vertex]==1)
      return origin[dsu.find_set(pred[vertex])];
    else
      return vertex;
  }

  void link_set_bridges(int vertex, int LCA, pair<int, int> bridge_edge){ // shrink Blossom to LCA
    for(int v=vertex; v!=LCA;v=parent(v)){
      dsu.unite_set(v, LCA);
      origin[dsu.find_set(LCA)] = LCA;
      if(vertex_state[v] == 1){
        bridge[v] = bridge_edge;
        for(auto const &e:graph[v]){ // repush edges, as this is now a single blossom-vertex
          if(e!=v){
              even_edges.emplace_back(v, e);
          }
        }
      }
    }
  }

  void retrieve_augmenting_path(int v, int w, bool is_reversed){
    if(v==w){
      augmenting_path.push_back(v);
    }
    else if(!is_reversed){
      if(vertex_state[v]==2){ //straight up the blossom
        augmenting_path.push_back(v);
        augmenting_path.push_back(mate[v]);
        retrieve_augmenting_path(pred[mate[v]], w, false);
      } else { // down the blossom and up on the other side
        augmenting_path.push_back(v);
        retrieve_augmenting_path(bridge[v].first, mate[v], true);
        retrieve_augmenting_path(bridge[v].second, w, false);
      }
    } else {
      if(vertex_state[v]==2){ //straight up the blossom
        retrieve_augmenting_path(pred[mate[v]], w, true);
        augmenting_path.push_back(mate[v]);
        augmenting_path.push_back(v);
      } else { // down the blossom and up on the other side
        retrieve_augmenting_path(bridge[v].second, w, true);
        retrieve_augmenting_path(bridge[v].first, mate[v], false);
        augmenting_path.push_back(v);
      }
    }
  }

  /// Searches a single augmenting path
  bool augment(){

    /// INIT
    int timestamp = 0;
    even_edges.clear();
    int N = graph.size();
    dsu.init(N);

    for(int i=0;i<N;++i){
      origin[i]=i;
      pred[i]=i;
      ancestor_v[i] = ancestor_w[i] = 0;
      if(mate[i]==-1){
        vertex_state[i]=2;//Even
        for(auto const &e:graph[i]){
          even_edges.emplace_back(i, e);
        }
      } else {
        vertex_state[i]=0; //unreached
      }
    }

    ///  END OF INIT
    /// AUGMENTING PATH DFS
    int v, w, w_ancestor = -1, v_ancestor = -1;
    bool found_altenating_path=0;

    while(!even_edges.empty() && !found_altenating_path){
      pair<int, int> curEdge = even_edges.back();
      even_edges.pop_back();
      tie(v, w)=curEdge;

      int v2=origin[dsu.find_set(v)];
      int w2=origin[dsu.find_set(w)];

      if(vertex_state[v2]!=2){
        cerr << "swap, should not happen...";
        swap(v, w);
        swap(v2, w2);
      }

      if(vertex_state[w2] == 0){ // Add matched pair to Tree ( v----w2===w2_mate )
          vertex_state[w2]=1;
          int w2_mate = mate[w2];
          vertex_state[w2_mate] = 2;
          for(auto const &e:graph[w2_mate]){
            if(e!=w2_mate){
              even_edges.emplace_back(w2_mate, e);
            }
          }
        pred[w2] = v;
      } else if(vertex_state[w2] == 2 && w2!=v2){ // if w2==v2 this was a inter-blossom edge
        int vUp = v2, wUp = w2; // used to climb up the DFS tree
        int lowest_common_ancestor = -1;
        w_ancestor = -1, v_ancestor=-1;

        // here we either found a blossom, or an augmenting path
        ++timestamp;
        while(lowest_common_ancestor==-1 &&(v_ancestor==-1 || w_ancestor==-1)){
          ancestor_v[vUp] = ancestor_w[wUp] = timestamp;

          if(w_ancestor == -1)
            wUp = parent(wUp);
          if(v_ancestor == -1)
            vUp = parent(vUp);

          if(mate[vUp] == -1)
            v_ancestor = vUp;
          if(mate[wUp] == -1)
            w_ancestor = wUp;

          if(ancestor_w[vUp] == timestamp)
            lowest_common_ancestor = vUp;
          else if(ancestor_v[wUp] == timestamp)
            lowest_common_ancestor = wUp;
          else if(v_ancestor==w_ancestor && v_ancestor != -1)
            lowest_common_ancestor = vUp;
        }

        if (lowest_common_ancestor==-1){ // different DFS tree roots
          found_altenating_path = 1;
        } else {
          // found blossom
          link_set_bridges(w2, lowest_common_ancestor, make_pair(w, v));
          link_set_bridges(v2, lowest_common_ancestor, make_pair(v, w));
        }
      }
    }
    if(!found_altenating_path)
      return false;
    ///  END OF AUGMENTING PATH DFS

    /// AUGMENTATION
    //retrieve augmenting path and use it
    retrieve_augmenting_path(v, v_ancestor, true);
    retrieve_augmenting_path(w, w_ancestor, false);
    //cerr << "augment: "<< v_ancestor << "-" << w_ancestor << " : ";
    while(!augmenting_path.empty()){
      v=augmenting_path.front();augmenting_path.pop_front();
      w=augmenting_path.front();augmenting_path.pop_front();
      mate[v]=w;
      mate[w]=v;
      //cerr << "->" << v << "->" << w;
    }
    //cerr << "\n";
    ///  END OF AUGMENTATION
    return true;
  }

  void get_matching(){
    bool done=false;
    while(!done){
      done = !augment();
    }
  }

  int get_maching_size(){
    get_matching();
    int ret=0;
    for(int i=0;i<(int)graph.size();++i){
      if(mate[i]!=-1 && mate[i]<i)++ret;
    }
    return ret;
  }
};