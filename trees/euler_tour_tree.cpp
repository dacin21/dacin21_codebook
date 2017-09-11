#include <bits/stdc++.h>

using namespace std;

/*
  Euler tour tree,
  by Daniel Rutschmann, Switzerland
  runs in O(N + Q * log(N)) armortized time
  and in O(N) space
*/

struct treeNode{
  treeNode *left, *right, *parent;
  int a, b;
  int subSize;
  long long myAgg, subAgg;
  treeNode(int _a, int _b, treeNode *_par):left(0), right(0), parent(_par), a(_a),b(_b), subSize(1), myAgg(0), subAgg(0){}
  void recalcSubSize(){
    subSize = 1+(left ? left->subSize : 0) + (right ? right->subSize : 0);
    subAgg = myAgg + (left ? left->subAgg :0) + (right ? right->subAgg : 0);
  }
  inline void rotateUp(){
      if(this==parent->right){
          if(parent->parent) (parent->parent->right == parent ? parent->parent->right : parent->parent->left) = this;

          parent->right = left;
          if(left) left->parent = parent;

          left = parent;
          parent = parent->parent;
          left->parent = this;
          left->recalcSubSize();
          recalcSubSize();
      } else {
          if(parent->parent) (parent->parent->left == parent ? parent->parent->left : parent->parent->right) = this;

          parent->left = right;
          if(right) right->parent = parent;

          right = parent;
          parent = parent->parent;
          if(right)right->parent = this;
          right->recalcSubSize();
          recalcSubSize();
      }
    }
};
inline void rotateUp(treeNode *cur){
  cur->rotateUp();
}
inline treeNode *splay(treeNode *cur){
  while(cur->parent){
    if(cur->parent->parent==0){
      rotateUp(cur);
    } else {
      if((cur->parent->left == cur) == (cur->parent->parent->left==cur->parent)){
        rotateUp(cur->parent);
        rotateUp(cur);
      } else {
        rotateUp(cur);
        rotateUp(cur);
      }
    }
  }
  return cur;
}
inline void findLastNodeThenSplayIt(treeNode * &a){
  splay(a);
  while(a->right) a = a->right;
  splay(a);
}
inline void findFirstNodeThenSplayIt(treeNode * &a){
  splay(a);
  while(a->left) a = a->left;
  splay(a);
}

inline pair<treeNode*, treeNode*> split(treeNode *cur){// splits the BST containing cur into the part left of cur, cur itself and the part to the right of cur
  splay(cur);
  if(cur->left){
      cur->left->parent = 0;
      cur->myAgg+=cur->left->subAgg;
      //cur->subAgg+=cur->left->subAgg;
  }
  if(cur->right){
      cur->right->parent = 0;
      findFirstNodeThenSplayIt(cur->right);
      cur->right->myAgg+=cur->myAgg;
      cur->right->recalcSubSize();
  }
  pair<treeNode*, treeNode*> ret = make_pair(cur->left, cur->right);
  cur->left=cur->right=0;
  cur->recalcSubSize();
  return ret;
}

inline treeNode *join(treeNode * a, treeNode * b){
  if(b==0)return a; if(a==0) return b;
  findLastNodeThenSplayIt(a);
 findFirstNodeThenSplayIt(b);
  a->right = b;
  b->parent = a;
  b->myAgg-=a->subAgg;

  b->recalcSubSize();
  a->recalcSubSize();
  return a;
}
inline void join(treeNode *left, treeNode *cur, treeNode *right){//cur has to be a single node
  if(right){
      findFirstNodeThenSplayIt(right);
      right->parent = cur;
      right->myAgg-=cur->myAgg;
      right->recalcSubSize();
  }
  if(left){
      splay(left);
      left->parent = cur;
      cur->myAgg-=left->subAgg;
      //cur->subAgg-=left->subAgg;
  }

  cur->left = left;
  cur->right = right;
  cur->recalcSubSize();
}
inline void evert(treeNode *cur){//wraps the BST so that cur is the first node of the tree
  splay(cur);
  treeNode *l = cur->left;
  cur->myAgg+=l->subAgg;
  cur->left = 0;
  l->parent = 0;
  cur->recalcSubSize();
  join(cur, l);
}

int dist(treeNode *l, treeNode *r){
  if(l==r)return -1;
  splay(r);
  splay(l);
  while(r->parent!=0) r->rotateUp();
  assert(r->left==l);
  return (l->right ? l->right->subSize:0);
}

struct EulerTourTree{
  vector<map<int, treeNode *> > graph;///Pointer to nodes of edges
  void init(int N){
    graph.resize(N);
  }
  EulerTourTree(int N){
    init(N);
    for(int i=0;i<N;++i){
      graph[i][i] = new treeNode(i, i, 0);
    }
  }
  bool connected(int a, int b){
    if(graph[a].empty() || graph[b].empty()) return false; //special case of an empty tree
    treeNode *ptrA = graph[a].begin()->second;
    treeNode *ptrB = graph[b].begin()->second;

    splay(ptrA);
    splay(ptrB);
    while(ptrA->parent) ptrA = ptrA->parent;
    return ptrA == ptrB; //check whether a and b are in the same BST
  }
  void reroot(int a){ //makes a the root of its Euler tour tree
    if(graph[a].empty()) return;//special case of an empty tree
    treeNode* cur = graph[a].begin()->second;
    treeNode* cur2 = cur;
    findLastNodeThenSplayIt(cur2);
    if(cur2->b == a) return; ///a is already the root
    evert(cur);
  }

  void cut(int a, int b){
    treeNode *lEdge = graph[a][b]; //lEdge and rEdge are the nodes of the edges between a and b
    treeNode *rEdge = graph[b][a];
    splay(lEdge), splay(rEdge);

    treeNode *tmp = lEdge;
    while(tmp->parent!=rEdge)tmp = tmp->parent;
    if(tmp==rEdge->right) swap(rEdge, lEdge); //make sure that lEdge is to the left of rEdge (in the BST)
    treeNode *rSub = split(rEdge).second;
    treeNode *lSub = split(lEdge).first;
    join(lSub, rSub);

    delete lEdge;
    delete rEdge;
    graph[a].erase(b);
    graph[b].erase(a);
  }
  void link(int a, int b){
    reroot(b);
    reroot(a);
    treeNode* ab = new treeNode(a, b, 0);
    treeNode* ba = new treeNode(b, a, 0);
    join(graph[a][a], ab, graph[b][b]);
    join(ab, ba);
    graph[a][b] = ab;
    graph[b][a] = ba;
  }
  void updateAgg(int a){
    treeNode* cur = graph[a][a];
    treeNode* cur2 = cur;
    findFirstNodeThenSplayIt(cur2);
    if(cur2!=cur) evert(cur);
    splay(cur);
    int total = (cur->subSize+2)/3;

    cur->myAgg+=total;// update node itself
    cur->recalcSubSize();
    if(cur->right){
      treeNode *next = cur->right;
      while(next->left){
        next=next->left;
      }
      splay(next);
      next->myAgg-=total;
      next->recalcSubSize();
    }
    for(auto &e:graph[a]){
      if(e.first==a) continue;
      treeNode *lEdge = graph[a][e.first];
      treeNode *rEdge = graph[e.first][a];
      int subSize = (dist(lEdge, rEdge)+2)/3;
      splay(lEdge);
      lEdge->myAgg+=total-subSize;
      lEdge->recalcSubSize();
      splay(rEdge);
      rEdge->myAgg-=total-subSize;
      rEdge->recalcSubSize();
    }
    while(graph[a].size() >1){
      int b=graph[a].begin()->first;
      if(b==a) b=graph[a].rbegin()->first;
      cut(a, b);
    }
  }
  long long queryAgg(int a){
    reroot(a);
    treeNode *cur = graph[a][a];
    splay(cur);
    return cur->myAgg + (cur->left ? cur->left->subAgg : 0);
  }

};
inline int dynacon1(){
  cin.tie(0);
  ios_base::sync_with_stdio(false);

  int N, Q, a, b;
  string s;
  cin >> N >> Q;
  EulerTourTree et(N+1);
  for(int i=0;i<Q;++i){
    cin >> s >> a >> b;
    if(s[0]=='c'){
        cout << (et.connected(a, b) ? "YES\n":"NO\n");
    }else if(s[0]=='a'){
      et.link(a, b);
    }else{
      et.cut(a, b);
    }
  }
  return 0;
}

// hackerrank unique-colors, uses agg
int main()
{
  //freopen("in.txt", "r", stdin);
    cin.tie(0);
    vector<vector<int> > G;
    vector<int> c;
    unordered_map<int, vector<int> > colorSets;
    int N, a, b;
    cin >> N;
    c.resize(N);
    G.resize(N);
    EulerTourTree et(N);
    for(int i=0;i<N;++i){
      cin >> c[i];
      colorSets[c[i]].emplace_back(i);
    }
    for(int i=0;i<N-1;++i){
      cin >> a >> b;
      --a, --b;
      et.link(a, b);
      G[a].emplace_back(b);
      G[b].emplace_back(a);
    }
    for(auto &e:colorSets){
      for(auto &f:e.second){
       // cerr << f << ": { ";

        et.updateAgg(f);

        /*for(int i=0;i<N;++i){
          cerr << et.queryAgg(i) << " ";
        }
        cerr << "}\n";*/
      }
      //cerr << " | ";
      for(auto &f:e.second){
        for(auto &g:G[f]){
          if(g<f || c[g]!=c[f]){
            et.link(g, f);
          }
        }
      }

    }
    for(int i=0;i<N;++i){
      cout << et.queryAgg(i) << "\n";
    }

    return 0;
}