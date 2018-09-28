// based on 'mwmatching.py' python code from http://jorisvr.nl/article/maximum-matching
/*
 * Based on a Python implementation of the general maximum Matching algorithm
 * uses a primal-dual method.
 * runs in O(V^3)
 * is quite complicated
*/
#include <bits/stdc++.h>
using namespace std;
#define donone(x) do{}while(0)
//#define DBG(x) (cerr << x)
#define DBG(x) donone(x)
#define asser(x) donone(x)
//#define asser(x) assert(x)
struct Matcher{
  using cost_t = long long;
  int N, M; // number of vertices and edges
  cost_t max_edge_weight=0;
  vector<pair<array<int, 2>, cost_t> > edges;//edge list
  vector<vector<int> > graph; // adjacency edge index list
  vector<int> mate; //matching edge vector
  vector<int> vertex_state;//0:unvisited   1:below/-/odd/T   2:above/+/even/S
  vector<int> pred; // predecessor node in DFS tree
  vector<int> blossom_top, blossom_parent, blossom_base; //highest level blossom, parent-level blossom, blossom-representative
  // would use std::optional with c++17 instead of std::unique_ptr
  vector<unique_ptr<vector<int> > > blossom_childs;// list of all nodes contained in this blossom
  vector<unique_ptr<vector<int> > > blossom_endpoints;// = blossom_childs but in edge list ???
  vector<unique_ptr<vector<int> > > best_blossom_edges;
  vector<int> best_edge;
  vector<int> free_blossom_indicies;
  vector<char> is_edge_allowed;
  vector<cost_t> dual;
  queue<int> q;

  int timestamp = 0;

  Matcher(int NN){
    N=NN;
    M=0;
    graph.resize(N);

  }
  void add_edge(int a, int b, cost_t c){
    asser(a>=0 && a<N);
    asser(b>=0 && b<N);
    c*=2;//needed so that blossom shrinking slack is an integer
    max_edge_weight = max(max_edge_weight, c);
    graph[b].emplace_back(2*M);
    graph[a].emplace_back(2*M+1);
    edges.emplace_back(array<int, 2>({a, b}), c);
    ++M;
  }
  int get_endpoint(int index){
    return edges[index/2].first[index%2];
  }
  cost_t slack(int edge_index){
    auto edge = edges[edge_index];
    return dual[edge.first[0]] + dual[edge.first[1]] - edge.second;
  }
  void get_blossom_leaves(int vertex, vector<int> &ret){
    if(vertex < N){
      ret.emplace_back(vertex);
      return;
    } else
      for(auto const &e:*blossom_childs[vertex])
        get_blossom_leaves(e, ret);
  }

  void assign_label(int vertex, int label, int edge_endpoint){
    DBG("assign label:" << vertex << ", " << label << ", " << edge_endpoint << "\n");
    int top = blossom_top[vertex];
    asser(vertex_state[vertex]==0 && vertex_state[top]==0);
    vertex_state[vertex] = vertex_state[top] = label;
    pred[vertex] = pred[top] = edge_endpoint;
    best_edge[vertex] = best_edge[top] = -1;

    if(label==2){
      vector<int> tmp;
      get_blossom_leaves(top, tmp);
      for(auto const &e:tmp){
        q.push(e);
      }
    } else if(label==1){
      int base = blossom_base[top];
      asser(mate[base]!=-1);
      assign_label(get_endpoint(mate[base]), 2, mate[base]^1);
    }
  }
  /// back-tracing to get LCA to find augmenting path or blossom
  int scan_blossom(int v, int w){
    vector<int> path;
    int base=-1;
    while(v!=-1 || w!=-1){
      int v2 = blossom_top[v];
      if(vertex_state[v2]&4){
        base = blossom_base[v2];
        break;
      }
      asser(vertex_state[v2] == 2);
      path.emplace_back(v2);
      vertex_state[v2]|=4;
      asser(pred[v2]==mate[blossom_base[v2]]);
      if(pred[v2]==-1)
        v=-1; // dfs-tree-root found on this path
      else {
        v=get_endpoint(pred[v2]);
        v2 = blossom_top[v];
        asser(vertex_state[v2]==1);
        asser(pred[v2]!=-1);
        v=get_endpoint(pred[v2]);
      }
      if(w!=-1){
        swap(v, w);//next step on other path
      }
    }
    for(auto const &e:path){
      vertex_state[e]&=~4;
      asser(vertex_state[e] == 2);
    }
    return base;
  }

  void add_blossom(int base, int bridge_edge){
    DBG("adding blossom with base: " << base << ", edge:" << edges[bridge_edge].first[0] << "/" << edges[bridge_edge].first[1] << " (" << edges[bridge_edge].second << ")\n");
    int v, w;
    v=edges[bridge_edge].first[0];
    w=edges[bridge_edge].first[1];
    int v2 = blossom_top[v];
    int w2 = blossom_top[w];
    int b2 = blossom_top[base];

    int blossom_index = free_blossom_indicies.back();free_blossom_indicies.pop_back();
    blossom_base[blossom_index] = base;
    blossom_parent[blossom_index] = -1;
    blossom_parent[b2]=blossom_index;
    blossom_childs[blossom_index] = make_unique<vector<int> >();
    blossom_endpoints[blossom_index] = make_unique<vector<int> >();
    /// extract blossom pseudo-nodes
    while(v2!=b2){
        blossom_parent[v2]=blossom_index;
        blossom_childs[blossom_index]->push_back(v2);
        blossom_endpoints[blossom_index]->push_back(pred[v2]);
        asser(vertex_state[v2]==1||(vertex_state[v2]==2 && pred[v2] == mate[blossom_base[v2]]));
        asser(pred[v2]!=-1);
        v = get_endpoint(pred[v2]);
        v2 = blossom_top[v];
    }
    blossom_childs[blossom_index]->push_back(b2);
    reverse(blossom_childs[blossom_index]->begin(), blossom_childs[blossom_index]->end());
    reverse(blossom_endpoints[blossom_index]->begin(), blossom_endpoints[blossom_index]->end());
    blossom_endpoints[blossom_index]->emplace_back(2*bridge_edge);

    while(w2!=b2){
        blossom_parent[w2] = blossom_index;
        blossom_childs[blossom_index]->push_back(w2);
        blossom_endpoints[blossom_index]->push_back(pred[w2]^1);
        asser(vertex_state[w2]==1||(vertex_state[w2]==2 && pred[w2] == mate[blossom_base[w2]]));
        asser(pred[w2]!=-1);
        w = get_endpoint(pred[w2]);
        w2 = blossom_top[w];
    }
    //end of extraction
    //initialize blossom variables
    //state
    asser(vertex_state[b2]==2);
    vertex_state[blossom_index] = 2;
    pred[blossom_index]= pred[b2];
    dual[blossom_index]=0;
    vector<int> tmp;
    //edges
    get_blossom_leaves(blossom_index, tmp);
    // update potential edges
    for(auto const &e:tmp){
      if(vertex_state[blossom_top[e]] == 1){
        q.push(e);
      }
      blossom_top[e] = blossom_index;
    }
    vector<int> best_edge_to(2*N, -1);
    for(auto const &v2:*blossom_childs[blossom_index]){
      vector<int> nblist;
      if(!best_blossom_edges[v2]){
        tmp.clear();
        get_blossom_leaves(v2, tmp);
        for(auto const &vert:tmp)
          for(auto const &p:graph[vert]){
            nblist.push_back(p/2);
          }
      } else {
        nblist = *best_blossom_edges[v2];
      }
      for(auto const &k:nblist){
        int i, j;
        i = edges[k].first[0];
        j = edges[k].first[1];
        if(blossom_top[j]==blossom_index)
          swap(i, j);
        int j2 = blossom_top[j];
        if(j2!=blossom_index && vertex_state[j2]==2 && (best_edge_to[j2]==-1 || slack(k) < slack(best_edge_to[j2]))){
          best_edge_to[j2] = k;
        }
      }
      best_blossom_edges[v2].reset();
      best_edge[v2] = -1;
    }

    best_blossom_edges[blossom_index] = make_unique<vector<int> >();
    for(auto const &e:best_edge_to)
      if(e!=-1)
        best_blossom_edges[blossom_index]->push_back(e);

    best_edge[blossom_index]=-1;
    for(auto const &e:*best_blossom_edges[blossom_index])
      if(best_edge[blossom_index]==-1 || slack(e) < slack(best_edge[blossom_index]))
        best_edge[blossom_index]=e;
  }



  void expand_blossom(int blossom, bool is_endstage){
    DBG("expanding blossom: " << blossom << " (" << is_endstage << ")\n");
    //remove parent relation
    for(auto const &e:*blossom_childs[blossom]){
      blossom_parent[e] = -1;
      if(e<N)
        blossom_top[e]=e;
      else if(is_endstage && dual[e]==0)
        expand_blossom(e, is_endstage);
      else {
        vector<int> tmp;
        get_blossom_leaves(e, tmp);
        for(auto const&f:tmp)
          blossom_top[f]=e;
      }
    }

    if(!is_endstage && vertex_state[blossom]==1){
      //special case of odd blossom
      asser(pred[blossom]!=-1);
      int entry_child = blossom_top[get_endpoint(pred[blossom]^1)];
      int j = find(blossom_childs[blossom]->begin(), blossom_childs[blossom]->end(), entry_child) - blossom_childs[blossom]->begin();
      int jstep, endptrick, jIndex;
      int blossom_child_count = (int)blossom_childs[blossom]->size();
      if(j&1){
        j-=blossom_child_count;
        jstep = 1;
        endptrick = 0;
      } else {
        jstep = -1;
        endptrick = 1;
      }
      int p = pred[blossom];
      while(j!=0){
        vertex_state[get_endpoint(p^1)]=0;
        jIndex = (j - endptrick);
        if(jIndex<0)jIndex+=blossom_child_count;

        vertex_state[get_endpoint((*blossom_endpoints[blossom])[jIndex]^endptrick^1)]=0;
        assign_label(get_endpoint(p^1), 1, p);
        is_edge_allowed[(*blossom_endpoints[blossom])[jIndex]/2] = true;

        j+=jstep;
        jIndex = (j - endptrick);
        if(jIndex<0)jIndex+=blossom_child_count;

        p=(*blossom_endpoints[blossom])[jIndex]^endptrick;
        is_edge_allowed[p/2]=true;
        j+=jstep;
      }


      int bv = (*blossom_childs[blossom])[0];
      vertex_state[get_endpoint(p^1)] = vertex_state[bv] = 1;
      pred[get_endpoint(p^1)] = pred[bv] = p;
      best_edge[bv]=-1;

      j+=jstep;
      jIndex = j;
      if(jIndex<0)jIndex+=blossom_child_count;
      while((*blossom_childs[blossom])[jIndex]!=entry_child){
        bv = (*blossom_childs[blossom])[jIndex];
        if(vertex_state[bv] == 2){
          j+=jstep;
          jIndex = j;
          if(jIndex<0)jIndex+=blossom_child_count;
          continue;
        }
        vector<int> tmp;
        get_blossom_leaves(bv, tmp);
        for(auto const &v:tmp)
          if(vertex_state[v]!=0){
            asser(vertex_state[v]==1);
            asser(blossom_top[v] == bv);
            vertex_state[v]=0;
            vertex_state[get_endpoint(mate[blossom_base[bv]])] = 0;
            assign_label(v, 1, pred[v]);
            break;
          }

        j+=jstep;
        jIndex = j;
        if(jIndex<0)jIndex+=blossom_child_count;
      }
    }
    // free blossom information
    vertex_state[blossom] = pred[blossom] = -1;
    blossom_childs[blossom].reset();
    blossom_endpoints[blossom].reset();
    best_blossom_edges[blossom].reset();
    blossom_base[blossom]=-1;
    best_edge[blossom]=-1;
    free_blossom_indicies.push_back(blossom);
  }

  void augment_blossom(int blossom, int vertex){
    int v2 = vertex;
    while(blossom_parent[v2] != blossom)
      v2 = blossom_parent[v2];

    if(v2 >= N)
      augment_blossom(v2, vertex); // recurse on sub-blossom

    int i, j, jIndex, jstep, endptrick, blossom_child_count = blossom_childs[blossom]->size();
    i=j=find(blossom_childs[blossom]->begin(), blossom_childs[blossom]->end(), v2) - blossom_childs[blossom]->begin();

    if(i&1){
      j-=blossom_child_count;
      jstep=1;
      endptrick=0;
    } else {
      jstep=-1;
      endptrick=1;
    }
    jIndex = j;
    if(jIndex<0)jIndex+=blossom_child_count;
    while(j!=0){
      //recurse on sub-blossoms
      j+=jstep;
      jIndex = j;
      if(jIndex<0)jIndex+=blossom_child_count;
      v2 = (*blossom_childs[blossom])[jIndex];
      jIndex = (j - endptrick);
      if(jIndex<0)jIndex+=blossom_child_count;
      int p = (*blossom_endpoints[blossom])[jIndex] ^ endptrick;
      if(v2>=N)
        augment_blossom(v2, get_endpoint(p));
      j+=jstep;
      jIndex = j;
      if(jIndex<0)jIndex+=blossom_child_count;
      v2 = (*blossom_childs[blossom])[jIndex];
      if(v2>=N)
        augment_blossom(v2, get_endpoint(p^1));

      //match the pseudo-nodes
      mate[get_endpoint(p)]=p^1;
      mate[get_endpoint(p^1)]=p;
    }
    //rotate to get new base at index 0
    rotate(blossom_childs[blossom]->begin(), blossom_childs[blossom]->begin()+i, blossom_childs[blossom]->end());
    rotate(blossom_endpoints[blossom]->begin(), blossom_endpoints[blossom]->begin()+i, blossom_endpoints[blossom]->end());
    blossom_base[blossom] = blossom_base[(*blossom_childs[blossom])[0]];
    asser(blossom_base[blossom]==vertex);
  }

  void augment_matching(int edge_index){
    DBG("augmeting with: " << edges[edge_index].first[0] << "/" << edges[edge_index].first[1] << " (" << edges[edge_index].second << ")\n");
    int v = edges[edge_index].first[0], w = edges[edge_index].first[1];
    for(auto const& e : {make_pair(v, 2*edge_index+1), make_pair(w, 2*edge_index)}){//hack to run for both v and w
      int s, p;
      tie(s, p) = e;
      for(;;){
        int s2 = blossom_top[s];
        asser(vertex_state[s2]==2);
        asser(pred[s2]== mate[blossom_base[s2]]);
        if(s2>=N)
          augment_blossom(s2, s);

        mate[s]=p;
        if(pred[s2]==-1)
          break;

        int t=get_endpoint(pred[s2]);
        int t2 = blossom_top[t];
        asser(vertex_state[t2]==1);
        asser(pred[t2]!=-1);
        s=get_endpoint(pred[t2]);
        int j=get_endpoint(pred[t2]^1);
        asser(blossom_base[t2] == t);
        if(t2 >= N)
          augment_blossom(t2, j);

        mate[j] = pred[t2];

        p = pred[t2]^1;

      }
    }
  }

  void get_matching(bool is_perfect_matching = false){
    mate.resize(N, -1);
    vertex_state.resize(2*N, 0);
    pred.resize(2*N, -1);
    blossom_top.resize(2*N, 0);
    iota(blossom_top.begin(), blossom_top.begin()+N, 0);
    blossom_parent.resize(2*N, -1);
    blossom_base.resize(2*N, -1);
    iota(blossom_base.begin(), blossom_base.begin()+N, 0);
    blossom_childs.resize(2*N);
    blossom_endpoints.resize(2*N);
    best_edge.resize(2*N, -1);
    best_blossom_edges.resize(2*N);
    free_blossom_indicies.resize(N);
    iota(free_blossom_indicies.begin(), free_blossom_indicies.begin()+N, N);
    dual.resize(2*N, 0);
    fill(dual.begin(), dual.begin()+N, max_edge_weight/2);
    is_edge_allowed.clear();
    q = queue<int>();
    is_edge_allowed.resize(M, false);

    // augmenting phases
    for(int t=0;t<N;++t){
      fill(vertex_state.begin(), vertex_state.begin()+2*N, 0);//clear states
      fill(best_edge.begin(), best_edge.begin()+2*N, -1);
      fill(best_blossom_edges.begin()+N, best_blossom_edges.end(), nullptr);
      fill(is_edge_allowed.begin(), is_edge_allowed.begin()+M, false);
      q = queue<int>();

      for(int i=0;i<N;++i)
        if(mate[i]==-1 && vertex_state[blossom_top[i]] == 0)
          assign_label(i, 2, -1);

      bool did_augment = 0;
      for(;;){
        //do primal operations while possible
        while(!did_augment && !q.empty()){
          int v=q.front();q.pop();
          DBG("POP v=" << v << "\n");
          asser(vertex_state[blossom_top[v]]==2);

          for(auto const &e:graph[v]){
            cost_t edge_index_slack;
            int edge_index = e/2;
            int w = get_endpoint(e);
            if(blossom_top[v]==blossom_top[w])
              continue;//internal blossom edge

            if(!is_edge_allowed[edge_index]){
              edge_index_slack = slack(edge_index);
              if(edge_index_slack<=0)// why <= and not ==
                is_edge_allowed[edge_index]=true;
            }
            if(is_edge_allowed[edge_index]){
              if(vertex_state[blossom_top[w]] == 0)
                //Case 1: grow
                assign_label(w, 1, e^1);
              else if(vertex_state[blossom_top[w]] == 2){
                //either blossom or augmenting path
                int base = scan_blossom(v, w);
                if(base!=-1)
                  //found blossom
                  add_blossom(base, edge_index);
                else {
                  //found alternating path
                  augment_matching(edge_index);
                  did_augment = true;
                  break;
                }
              } else if(vertex_state[w]==0){
                //special case: inside T-blossom
                asser(vertex_state[blossom_top[w]] == 1);
                vertex_state[w]=1;
                pred[w] = e^1;
              }
            } else if(vertex_state[blossom_top[w]] == 2){
              //can't rech other blososm yet
              int b = blossom_top[v];
              if(best_edge[b] == -1 || edge_index_slack < slack(best_edge[b]))
                best_edge[b] = edge_index;
            } else if(vertex_state[w] == 0){
              //can't reach w yet (need to modify dual first)
              if(best_edge[w]==-1 || edge_index_slack < slack(best_edge[w]))
                best_edge[w] = edge_index;
            }
          }
        }
        if(did_augment)
          break;
        ///end of primal phase

        DBG("pstate:");
        for(int i=0;i<2*N;++i) DBG( vertex_state[i] << " ");
        DBG("\n");

        int delta_type = -1;
        cost_t delta;
        int delta_edge, delta_blossom;
        delta = delta_edge = delta_blossom = -1; // = NONE???

        //minimum vertex dual
        if(!is_perfect_matching){
          delta_type = 1;
          delta = *min_element(dual.begin(), dual.begin()+N);
        }
        //grow
        //delta2: minimum edge slack between S-vertex and free vertex
        for(int i=0;i<N;++i)
          if(vertex_state[blossom_top[i]] == 0 && best_edge[i]!=-1){
            cost_t d = slack(best_edge[i]);
            if(delta_type==-1 || d < delta){
              delta_type = 2;
              delta = d;
              delta_edge = best_edge[i];
            }
          }
        //blossom shrink or augment
        //delta3: 1/2*minimum edge slack between S-vertices
        for(int i=0;i<2*N;++i)
          if(blossom_parent[i]==-1 && vertex_state[i]==2 && best_edge[i]!=-1){
            cost_t kslack = slack(best_edge[i]);
            asser((kslack/2)*2==kslack);
            cost_t d = kslack/2;
            if(delta_type==-1 || d < delta){
              delta_type = 3;
              delta = d;
              delta_edge = best_edge[i];
            }
          }

        //blossom expand
        //delta4: minimum dual of a blossom
        for(int i=N;i<2*N;++i)
          if(blossom_base[i]>=0 && blossom_parent[i] == -1 && vertex_state[i] == 1)
            if(delta_type==-1 || dual[i]<delta){
              delta_type = 4;
              delta = dual[i];
              delta_blossom = i;
            }

        if(delta_type == -1){
          asser(is_perfect_matching);
          delta_type = 1;
          delta = max((cost_t)0, *min_element(dual.begin(), dual.begin()+N));
        }

        DBG(" -  delta type: " << delta_type << ", delta value: " << delta << "\n");
        //update dual variables
        for(int i=0;i<N;++i){
          if(vertex_state[blossom_top[i]] == 2)
            dual[i]-=delta;
          else if(vertex_state[blossom_top[i]] == 1)
            dual[i]+=delta;
        }
        for(int i=N;i<2*N;++i)
          if(blossom_base[i]!=-1 && blossom_parent[i] == -1){
            if(vertex_state[i] == 2)
              dual[i]+=delta;
            else if(vertex_state[i] == 1)
              dual[i]-=delta;
          }

        if(delta_type == 1) //was improvement, just reduction to 0
          break;
        else if(delta_type == 2){// append new equality edge
          is_edge_allowed[delta_edge] = true;
          int i=edges[delta_edge].first[0], j=edges[delta_edge].first[1];
          if(vertex_state[blossom_top[i]]==0)
            swap(i, j);
          asser(vertex_state[blossom_top[i]] == 2);
          q.push(i);
        } else if(delta_type==3){//append new equality edge
          is_edge_allowed[delta_edge] = true;
          int i = edges[delta_edge].first[0];
          asser(vertex_state[blossom_top[i]] == 2);
          q.push(i);
        } else if(delta_type == 4){//expand blossom
          expand_blossom(delta_blossom, false);
        }

        DBG("duals: ");
        for(int i=0;i<2*N;++i) DBG( dual[i] << " ");
        DBG("\n");
        DBG("state: ");
        for(int i=0;i<2*N;++i) DBG( vertex_state[i] << " ");
        DBG("\n");
        ///end of dual phase
      }
      if(!did_augment)
        break;
      //expand all blossom with dual=0
      for(int i=N;i<2*N;++i)
        if(blossom_parent[i]==-1 && blossom_base[i] != -1 && vertex_state[i] == 2 && dual[i] == 0)
          expand_blossom(i, true);
      //end of stage
    }
    //make mate a vertex vector, not an edge index vector
    for(int i=0;i<N;++i)
      if(mate[i]!=-1)
        mate[i] = get_endpoint(mate[i]);
    //verify that mate is symmetric
    for(int i=0;i<N;++i)
      asser(mate[i]==-1 || mate[mate[i]] == i);
  }
  cost_t get_matching_weight(){
    cost_t ret = 0;
    vector<cost_t> best(graph.size(), numeric_limits<cost_t>::max());
    for(auto &e:edges){
        if(mate[e.first[0]] == e.first[1]){
            best[e.first[0]] = min(best[e.first[0]], e.second);
            best[e.first[1]] = min(best[e.first[1]], e.second);
        }
    }
    for(int i=0;i<N;++i){
        if(mate[i]!=-1 && mate[i]<i) ret+=min(best[i], best[mate[i]]);
    }
    return ret/2;
  }
};
