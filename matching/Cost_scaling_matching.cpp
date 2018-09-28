// Code by Min_25
// from https://github.com/min-25/scaling-algorithm-for-mwpm/blob/master/scaling_mwpm.cpp

// max weight perfect matching with cost scaling
// runs in O(n m log (n C) log (n))

#include <bits/stdc++.h>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <queue>
#include <functional>

using namespace std;

template <typename T>
class BinaryHeap {
public:
  struct Node {
    bool operator < (const Node& rhs) const { return value < rhs.value; }
    T value; int id;
  };
  BinaryHeap() {}
  BinaryHeap(int N) : size_(0), node(N + 1), index(N, 0) {}
  int size() const { return size_; }
  bool empty() const { return size_ == 0; }
  void clear() { while (size_ > 0) index[node[size_--].id] = 0; }
  T min() const { return node[1].value; }
  int argmin() const { return node[1].id; } // argmin ?
  T get_val(int id) const { return node[index[id]].value; }
  void pop() { if (size_ > 0) pop(1); }
  void erase(int id) { if (index[id]) pop(index[id]); }
  bool has(int id) const { return index[id] != 0; }
  void update(int id, T v) {
    if (!has(id)) return push(id, v);
    bool up = (v < node[index[id]].value);
    node[index[id]].value = v;
    if (up) up_heap(index[id]);
    else down_heap(index[id]);
  }
  void decrease_key(int id, T v) {
    if (!has(id)) return push(id, v);
    if (v < node[index[id]].value) node[index[id]].value = v, up_heap(index[id]);
  }
  void push(int id, T v) {
    // assert(!has(id));
    index[id] = ++size_; node[size_] = {v, id};
    up_heap(size_);
  }
private:
  void pop(int pos) {
    index[node[pos].id] = 0;
    if (pos == size_) { --size_; return; }
    bool up = (node[size_].value < node[pos].value);
    node[pos] = node[size_--]; index[node[pos].id] = pos;
    if (up) up_heap(pos);
    else down_heap(pos);
  }
  void swap_node(int a, int b) {
    swap(node[a], node[b]); index[node[a].id] = a; index[node[b].id] = b;
  }
  void down_heap(int pos) {
    for (int k = pos, nk = k; 2 * k <= size_; k = nk) {
      if (node[2 * k] < node[nk]) nk = 2 * k;
      if (2 * k + 1 <= size_ && node[2 * k + 1] < node[nk]) nk = 2 * k + 1;
      if (nk == k) break;
      swap_node(k, nk);
    }
  }
  void up_heap(int pos) {
    for (int k = pos; k > 1 && node[k] < node[k >> 1]; k >>= 1) swap_node(k, k >> 1);
  }
  int size_;
  vector<Node> node;
  vector<int> index;
};

template <typename Key>
class PairingHeaps {
private:
  struct Node {
    Node() : prev(-1) {} // "prev < 0" means the node is unused.
    Node(Key v) : key(v), child(0), next(0), prev(0) {}
    Key key; int child, next, prev;
  };
public:
  PairingHeaps(int H, int N) : heap(H), node(N) {
    // It consists of `H` Pairing heaps.
    // Each heap-node ID can appear at most 1 time(s) among heaps
    // and should be in [1, N).
  }

  void clear(int h) { if (heap[h]) clear_rec(heap[h]), heap[h] = 0; }
  void clear_all() {
    for (size_t i = 0; i < heap.size(); ++i) heap[i] = 0;
    for (size_t i = 0; i < node.size(); ++i) node[i] = Node();
  }
  bool empty(int h) const { return !heap[h]; }
  bool used(int v) const { return node[v].prev >= 0; }
  Key min(int h) const { return node[heap[h]].key; }
  int argmin(int h) const { return heap[h]; }

  void pop(int h) {
    // assert(!empty(h));
    erase(h, heap[h]);
  }
  void push(int h, int v, Key key) {
    // assert(!used(v));
    node[v] = Node(key);
    heap[h] = merge(heap[h], v);
  }
  void erase(int h, int v) {
    if (!used(v)) return;
    int w = two_pass_pairing(node[v].child);
    if (!node[v].prev) heap[h] = w;
    else {
      cut(v);
      heap[h] = merge(heap[h], w);
    }
    node[v].prev = -1;
  }
  void decrease_key(int h, int v, Key key) {
    if (!used(v)) return push(h, v, key);
    if (!node[v].prev) node[v].key = key;
    else {
      cut(v); node[v].key = key;
      heap[h] = merge(heap[h], v);
    }
  }

private:
  void clear_rec(int v) {
    for (; v; v = node[v].next) {
      if (node[v].child) clear_rec(node[v].child);
      node[v].prev = -1;
    }
  }

  inline void cut(int v) {
    auto& n = node[v]; int pv = n.prev, nv = n.next;
    auto& pn = node[pv];
    if (pn.child == v) pn.child = nv;
    else pn.next = nv;
    node[nv].prev = pv;
    n.next = n.prev = 0;
  }

  int merge(int l, int r) {
    if (!l) return r;
    if (!r) return l;
    if (node[l].key > node[r].key) swap(l, r);
    int lc = node[r].next = node[l].child;
    node[l].child = node[lc].prev = r;
    return node[r].prev = l;
  }

  int two_pass_pairing(int root) {
    if (!root) return 0;
    int a = root; root = 0;
    while (a) {
      int b = node[a].next, na = 0;
      node[a].prev = node[a].next = 0;
      if (b) na = node[b].next, node[b].prev = node[b].next = 0;
      a = merge(a, b);
      node[a].next = root; root = a; a = na;
    }
    int s = node[root].next; node[root].next = 0;
    while (s) {
      int t = node[s].next; node[s].next = 0;
      root = merge(root, s);
      s = t;
    }
    return root;
  }

private:
  vector<int> heap;
  vector<Node> node;
};

template <typename T>
struct PriorityQueue : public priority_queue< T, vector<T>, greater<T> > {
  PriorityQueue() {}
  PriorityQueue(int N) { this->c.reserve(N);}
  T min() const { return this->top(); }
  void clear() { this->c.clear(); }
};

template <typename T>
struct Queue {
  Queue() {}
  Queue(int N) : qh(0), qt(0), data(N) {}
  T operator [] (int i) const { return data[i]; }
  void enqueue(int u) { data[qt++] = u; }
  int dequeue() { return data[qh++]; }
  bool empty() const { return qh == qt; }
  void clear() { qh = qt = 0; }
  int size() const { return qt; }
  int qh, qt;
  vector<T> data;
};

struct DisjointSetUnion {
  DisjointSetUnion() {}
  DisjointSetUnion(int N) : par(N) {
    for (int i = 0; i < N; ++i) par[i] = i;
  }
  int find(int u) { return par[u] == u ? u : (par[u] = find(par[u])); }
  void unite(int u, int v) {
    u = find(u), v = find(v);
    if (u != v) par[v] = u;
  }
  vector<int> par;
};

struct LinkedList {
  LinkedList() {}
  LinkedList(int N, int M) : N(N), next(M) { clear(); }
  void clear() { head.assign(N, -1); }
  void push(int h, int u) { next[u] = head[h], head[h] = u; }
  int N;
  vector<int> head, next;
};
template <typename CostType, typename TotalCostType=int64_t>
class MaximumWeightPerfectMatching {
public:
  using cost_t = CostType;
  using tcost_t = TotalCostType;

private:
  enum Label { kSeparated = -2, kInner = -1, kFree = 0, kOuter = 1 };
  static constexpr cost_t kInf = cost_t(1) << (sizeof(cost_t) * 8 - 2);
  static constexpr int kShift = 2; // 2 or 3 is recommended
  static_assert(kShift > 0, "Shift amount should be positive");

public:
  struct InputEdge { int from, to; cost_t cost; };

private:
  template <typename T> using ModifiableHeap = BinaryHeap<T>;
  template <typename T> using ModifiableHeaps = PairingHeaps<T>;
  template <typename T> using FastHeap = PriorityQueue<T>;

  struct Edge { int to; cost_t cost; };
  struct Link { int from, to; cost_t cost; };
  struct Mate { int v; cost_t cost; };

  struct Node {
    struct NodeLink { int b, v; cost_t cost; };
    Node() {}
    Node(int u) : parent(0), size(1) { link[0] = link[1] = {u, u, 0}; }
    int next_v() const { return link[0].v; }
    int next_b() const { return link[0].b; }
    cost_t next_cost() const { return link[0].cost; }
    int prev_v() const { return link[1].v; }
    int prev_b() const { return link[1].b; }
    cost_t prev_cost() const { return link[1].cost; }
    int parent, size;
    NodeLink link[2];
  };
  struct Event {
    Event() {}
    Event(cost_t time, int id) : time(time), id(id) {}
    bool operator < (const Event& rhs) const { return time < rhs.time; }
    bool operator > (const Event& rhs) const { return time > rhs.time; }
    cost_t time; int id;
  };
  struct EdgeEvent {
    EdgeEvent() {}
    EdgeEvent(cost_t time, int from, int to) : time(time), from(from), to(to) {}
    bool operator > (const EdgeEvent& rhs) const { return time > rhs.time; }
    bool operator < (const EdgeEvent& rhs) const { return time < rhs.time; }
    cost_t time; int from, to;
  };

public:
  MaximumWeightPerfectMatching(int N, const vector<InputEdge>& in)
      : N(N), B((N - 1) / 2), S(N + B + 1), ofs(N + 2), edges(in.size() * 2),
        heap2(S), heap2_2(S), heap2s(S, S), heap3(edges.size()), heap4(S) {

    assert(N % 2 == 0);
    cost_min_ = kInf;
    for (auto& e : in) {
      ofs[e.from + 1]++, ofs[e.to + 1]++;
      cost_min_ = min(cost_min_, e.cost);
    }
    for (int i = 1; i <= N + 1; ++i) ofs[i] += ofs[i - 1];
    const int multiplier = scale_factor();
    for (auto& e : in) {
      edges[ofs[e.from]++] = {e.to, (e.cost - cost_min_) * multiplier};
      edges[ofs[e.to]++] = {e.from, (e.cost - cost_min_) * multiplier};
    }
    for (int i = N + 1; i > 0; --i) ofs[i] = ofs[i - 1];
    ofs[0] = 0;
  }

  tcost_t maximum_weight_perfect_matching() {
    initialize();

    cost_t max_c = 2;
    for (auto& e : edges) max_c = max(max_c, e.cost);

    scale_ = 0;
    while (max_c >> (scale_ + 1)) scale_++;

    while (scale_) {
      const int scale_diff = min(kShift, scale_); scale_ -= scale_diff;
      const int coeff1 = 1 << scale_diff, coeff2 = 2 * coeff1 - 2;

      time_limit_ = (scale_ == 0) ? kInf : sqrt(N);

      for (int u = 1; u <= N; ++u) potential[u] = potential[u] * coeff1 + coeff2;
      for (int b = N + 1; b < S; ++b) potential[b] *= coeff1;

      for (int b = N + 1; b < S; ++b) if (base[b] != b) {
        if (surface[b] == b) distribute_potential(b, 0);
      }
      for (int u = 1; u <= N; ++u) mate[u] = {0, 0}, surface[u] = u;
      for (int u = 1; u <= N; ++u) {
        cost_t p = kInf;
        for (int eid = ofs[u]; eid < ofs[u + 1]; ++eid) {
          auto& e = edges[eid];
          p = min(p, reduced_cost(u, e.to, e));
        }
        potential[u] -= p + 2;
      }
      //fprintf(stderr, "[ Start of Scale %d ]\n", scale_);
      int match = 0;
      for (int c = 1; ; ++c) {
        do_edmonds_search();
        reset_all();
        if (!within_time_limit(time_augment_)) break;
        match += find_maximal();
        reset_flags();
        //fprintf(stderr, "- %2d.%03d: |M|: %d/%d\n", scale_, c, match, N / 2);
        if (2 * match == N) break;
      }
      //fprintf(stderr, "[ End of Scale %d ]\n\n", scale_);
    }
    tcost_t ret = compute_optimal_value();
    return ret;
  }

  vector<int> get_mate(){
    vector<int> ret;
    for(int i=1;i<=N;++i) ret.push_back(mate[i].v-1);
    return ret;
  }

private:
  bool within_time_limit(cost_t t) const {
    return t <= time_limit_;
  }

  int scale_factor() const {
    return 2 * (N / 2 + 1);
  }

  tcost_t compute_optimal_value() const {
    tcost_t ret = 0; const int multiplier = scale_factor();
    for (int u = 1; u <= N; ++u) ret += mate[u].cost / multiplier;
    return (ret >> 1) + tcost_t(N / 2) * cost_min_;
  }

  inline cost_t scaled(cost_t c) const {
    return (c >> scale_) & ~cost_t(1);
  }

  inline cost_t reduced_cost(int u, int v, const Edge& e) const {
    return potential[u] + potential[v] - scaled(e.cost);
  }

  void distribute_potential(int b, cost_t c) {
    if (b <= N) {
      potential[b] += c;
    } else {
      c += potential[b] >> 1;
      for (int beta = base[b], bb = beta; ;) {
        distribute_potential(bb, c);
        if ((bb = node[bb].next_b()) == beta) break;
      }
      free_blossom(b);
    }
  }

  void fix_mate_and_base(int b) {
    if (b <= N) return;
    int bv = base[b], mv = node[bv].link[0].v, bmv = node[bv].link[0].b;
    int d = (node[bmv].link[1].v == mate[mv].v) ? 0 : 1;
    while (1) {
      int mv = node[bv].link[d].v, bmv = node[bv].link[d].b;
      if (node[bmv].link[1 ^ d].v != mate[mv].v) break;
      fix_mate_and_base(bv); fix_mate_and_base(bmv);
      bv = node[bmv].link[d].b;
    }
    fix_mate_and_base(base[b] = bv);
    mate[b] = mate[bv];
  }

  void reset_time() {
    time_current_ = 0;
  }

  void reset_blossom(int b) {
    label[b] = kFree; slack[b] = kInf; lazy[b] = 0;
  }

  void reset_all() {
    for (int v = 1; v <= N; ++v) {
      if (label[v] > 0) potential[v] -= time_current_;
      else {
        int bv = surface[v];
        potential[v] += lazy[bv];
        if (label[bv] == kInner) potential[v] += time_current_ - time_created[bv];
      }
      reset_blossom(v);
    }
    for (int b = N + 1; b < S; ++b) if (base[b] != b) {
      if (surface[b] == b) {
        if (label[b] > 0) potential[b] += (time_current_ - time_created[b]) << 1;
        else if (label[b] == kInner) fix_blossom_potential<kInner>(b);
        else fix_blossom_potential<kFree>(b);
      }
      time_created[b] = 0;
      heap2s.clear(b);
      reset_blossom(b);
    }
    reset_time();
    for (int b = N + 1; b < S; ++b) if (base[b] != b) {
      if (surface[b] == b && potential[b] == 0) expand_free(b, kFree);
    }
    que.clear();
    heap2.clear(); heap2_2.clear(); heap3.clear(); heap4.clear();
  }

  void reset_flags() {
    link[0].from = 0;
    for (int v = 1; v <= N; ++v) label[v] = kFree, link[v].from = 0;
    for (int b = N + 1; b < S; ++b) if (base[b] != b) {
      if (surface[b] == b) fix_mate_and_base(b);
      link[b].from = 0;
    }
  }

  void do_edmonds_search() {
    contract_count_ = 0, time_augment_ = kInf;
    for (int u = 1; u <= N; ++u) if (!mate[u].v) {
      link_blossom(surface[u], {0, 0, 0});
      push_outer_and_fix_potentials(surface[u], u, 0);
    }
    while (1) {
      if (augment()) break;
      if (adjust_dual_variables()) break;
      if (time_current_ == time_augment_) break;
    }
  }

  // heap

  template <Label Lab>
  inline cost_t fix_blossom_potential(int b) {
    // Return the amount.
    // (If v is an atom, the potential[v] will not be changed.)
    cost_t d = lazy[b]; lazy[b] = 0;
    if (Lab == kInner) {
      cost_t dt = time_current_ - time_created[b];
      if (b > N) potential[b] -= dt << 1, assert(potential[b] >= 0);
      d += dt;
    }
    return d;
  }

  template <Label Lab>
  inline void update_heap2(int x, int y, int by, cost_t t) {
    if (t >= slack[y]) return;
    slack[y] = t; best_from[y] = x;
    if (y == by) {
      if (Lab != kInner) heap2.decrease_key(y, EdgeEvent(t + lazy[y], x, y));
    } else {
      int gy = group[y];
      if (gy != y) {
        if (t >= slack[gy]) return;
        slack[gy] = t;
      }
      heap2s.decrease_key(by, gy, EdgeEvent(t, x, y));
      if (Lab == kInner) return;
      EdgeEvent m = heap2s.min(by);
      heap2.decrease_key(by, EdgeEvent(m.time + lazy[by], m.from, m.to));
    }
  }

  void activate_heap2_node(int b) {
    if (b <= N) {
      if (slack[b] < kInf) heap2.push(b, EdgeEvent(slack[b] + lazy[b], best_from[b], b));
    } else {
      if (heap2s.empty(b)) return;
      EdgeEvent m = heap2s.min(b);
      heap2.push(b, EdgeEvent(m.time + lazy[b], m.from, m.to));
    }
  }

  // contract

  void swap_blossom(int a, int b) {
    // Assume that `b` is a maximal blossom.
    swap(base[a], base[b]); if (base[a] == a) base[a] = b;
    swap(heavy[a], heavy[b]); if (heavy[a] == a) heavy[a] = b;
    swap(link[a], link[b]);
    swap(mate[a], mate[b]);
    swap(potential[a], potential[b]); swap(lazy[a], lazy[b]);
    swap(time_created[a], time_created[b]);
    for (int d = 0; d < 2; ++d) node[node[a].link[d].b].link[1 ^ d].b = b;
    swap(node[a], node[b]);
  }

  void set_surface_and_group(int b, int sf, int g) {
    surface[b] = sf, group[b] = g;
    if (b <= N) return;
    for (int bb = base[b]; surface[bb] != sf; bb = node[bb].next_b()) {
      set_surface_and_group(bb, sf, g);
    }
  }

  void merge_smaller_blossoms(int bid) {
    int lb = bid, largest_size = 1;
    for (int beta = base[bid], b = beta; ;) {
      if (node[b].size > largest_size) largest_size = node[b].size, lb = b;
      if ((b = node[b].next_b()) == beta) break;
    }
    for (int beta = base[bid], b = beta; ;) {
      if (b != lb) set_surface_and_group(b, lb, b);
      if ((b = node[b].next_b()) == beta) break;
    }
    group[lb] = lb;
    if (largest_size > 1) {
      surface[bid] = heavy[bid] = lb;
      swap_blossom(lb, bid);
    } else heavy[bid] = 0;
  }

  void contract(int x, int y) {
    int bx = surface[x], by = surface[y]; assert(bx != by);
    const int h = -(++contract_count_);
    link[surface[mate[bx].v]].from = link[surface[mate[by].v]].from = h;

    int lca = -1;
    while (1) {
      if (mate[by].v != 0) swap(bx, by);
      bx = lca = surface[link[bx].from];
      if (link[surface[mate[bx].v]].from == h) break;
      link[surface[mate[bx].v]].from = h;
    }

    const int bid = unused_bid[--unused_last_]; assert(unused_last_ >= 0);
    const cost_t c = potential[x] + potential[y] - (time_current_ << 1) + 2;
    int tree_size = 0;
    for (int d = 0; d < 2; ++d) {
      for (int bv = surface[x]; bv != lca; ) {
        int mv = mate[bv].v, bmv = surface[mv], v = mate[mv].v; cost_t mc = mate[bv].cost;
        int f = link[v].from, t = link[v].to; cost_t lc = link[v].cost;
        tree_size += node[bv].size + node[bmv].size;
        link[mv] = {x, y, c};

        if (bv > N) potential[bv] += (time_current_ - time_created[bv]) << 1;
        if (bmv > N) heap4.erase(bmv);
        push_outer_and_fix_potentials(bmv, label[x], fix_blossom_potential<kInner>(bmv));

        node[bv].link[d] = {bmv, mv, mc};
        node[bmv].link[1 ^ d] = {bv, v, mc}; node[bmv].link[d] = {bv = surface[f], f, lc};
        node[bv].link[1 ^ d] = {bmv, t, lc};
      }
      node[surface[x]].link[1 ^ d] = {surface[y], y, c};
      swap(x, y);
    }
    if (lca > N) potential[lca] += (time_current_ - time_created[lca]) << 1;
    node[bid].size = tree_size + node[lca].size;
    base[bid] = lca; link[bid] = link[lca]; mate[bid] = mate[lca];
    label[bid] = label[lca];
    surface[bid] = bid; time_created[bid] = time_current_;
    potential[bid] = 0; lazy[bid] = 0;

    merge_smaller_blossoms(bid); // O(n log n) time / Edmonds search
  }

  // grow

  void link_blossom(int v, Link l) {
    link[v] = {l.from, l.to, l.cost};
    if (v <= N) return;
    int b = base[v]; link_blossom(b, l);
    int pb = node[b].prev_b();
    l = {node[pb].next_v(), node[b].prev_v(), node[b].prev_cost()};
    for (int bv = b; ; ) {
      int bw = node[bv].next_b();
      if (bw == b) break;
      link_blossom(bw, l);
      Link nl = {node[bw].prev_v(), node[bv].next_v(), node[bv].next_cost()};
      bv = node[bw].next_b();
      link_blossom(bv, nl);
    }
  }

  void push_outer_and_fix_potentials(int v, int l, cost_t d) {
    label[v] = l;
    if (v > N) {
      for (int b = base[v]; label[b] <= 0; b = node[b].next_b()) {
        push_outer_and_fix_potentials(b, l, d);
      }
    } else {
      potential[v] += time_current_ + d;
      que.enqueue(v);
    }
  }

  void grow(int x, int y) {
    int by = surface[y], z = mate[by].v, bz = surface[z];
    if (label[bz] == kInner) {
      label[by] = kInner; time_created[by] = time_current_;
      heap2.erase(by); heap2_2.erase(by);
      heap3.emplace(time_current_ + 1, x, y);
      if (y != by) assert(potential[by] > 0);
    } else {
      int mz = mate[z].v;
      const cost_t p = (potential[mz] + lazy[by])
                     + (potential[z] + lazy[bz]) - mate[z].cost;
      assert(-2 <= p && p <= 0);

      bool visited = (label[by] != kFree);
      if (!visited) link_blossom(by, {0, 0, 0});
      label[by] = kInner; time_created[by] = time_current_; heap2.erase(by);
      if (y != by) heap4.update(by, time_current_ + (potential[by] >> 1));

      const cost_t c = (potential[x] - time_current_) + (potential[y] + lazy[by]) + 2;
      if (!visited) link_blossom(bz, {x, y, c});
      else link[bz] = link[z] = {x, y, c};
      if (p == 0) {
        push_outer_and_fix_potentials(bz, label[x], fix_blossom_potential<kFree>(bz));
        time_created[bz] = time_current_; heap2.erase(bz);
      } else {
        heap2_2.push(bz, time_current_ - p);
      }
    }
  }

  // expand

  void free_blossom(int bid) {
    unused_bid[unused_last_++] = bid;
    base[bid] = bid;
  }

  int recalculate_minimum_slack(int b, int g) {
    // Return the destination of the best edge of blossom `g`.
    if (b <= N) {
      if (slack[b] >= slack[g]) return 0;
      slack[g] = slack[b]; best_from[g] = best_from[b];
      return b;
    }
    int v = 0;
    for (int beta = base[b], bb = beta; ; ) {
      int w = recalculate_minimum_slack(bb, g);
      if (w != 0) v = w;
      if ((bb = node[bb].next_b()) == beta) break;
    }
    return v;
  }

  void construct_smaller_components(int b, int sf, int g) {
    surface[b] = sf, group[b] = g; // `group[b] = g` is unneeded.
    if (b <= N) return;
    for (int bb = base[b]; surface[bb] != sf; bb = node[bb].next_b()) {
      if (bb == heavy[b]) {
        construct_smaller_components(bb, sf, g);
      } else {
        set_surface_and_group(bb, sf, bb);
        int to = 0;
        if (bb > N) slack[bb] = kInf, to = recalculate_minimum_slack(bb, bb);
        else if (slack[bb] < kInf) to = bb;
        if (to > 0) heap2s.push(sf, bb, EdgeEvent(slack[bb], best_from[bb], to));
      }
    }
  }

  void move_to_largest_blossom(int bid) {
    const int h = heavy[bid];
    cost_t d = (time_current_ - time_created[bid]) + lazy[bid]; lazy[bid] = 0;
    assert(d >= 0);
    for (int beta = base[bid], b = beta; ;) {
      time_created[b] = time_current_;
      lazy[b] = d;
      if (b != h) construct_smaller_components(b, b, b), heap2s.erase(bid, b);
      if ((b = node[b].next_b()) == beta) break;
    }
    if (h > 0) swap_blossom(h, bid), bid = h;
    free_blossom(bid);
  }

  void expand(int bid, int l=0) {
    int mv = mate[base[bid]].v;
    move_to_largest_blossom(bid); // O(n log n) time / Edmonds search
    Link old_link = link[mv];
    if (l == 0) l = label[old_link.from], assert(l > 0);
    int old_base = surface[mate[mv].v], root = surface[old_link.to];
    int d = (mate[root].v == node[root].link[0].v) ? 1 : 0;
    for (int b = node[old_base].link[d ^ 1].b; b != root; ) {
      for (int t = 0; t < 2; ++t) {
        label[b] = kSeparated; int nb = node[b].link[d ^ 1].b;
        if (b > N && potential[b] == 0) expand_free(b);
        else activate_heap2_node(b);
        b = nb;
      }
    }
    for (int b = old_base; ; b = node[b].link[d].b) {
      label[b] = kInner;
      int nb = node[b].link[d].b;
      if (b == root) link[mate[b].v] = old_link;
      else link[mate[b].v] = {node[b].link[d].v, node[nb].link[d ^ 1].v, node[nb].link[d ^ 1].cost};
      link[surface[mate[b].v]] = link[mate[b].v]; // fix tree links
      if (b > N) {
        if (potential[b] == 0) expand(b, l);
        else heap4.push(b, time_current_ + (potential[b] >> 1));
      }
      if (b == root) break;
      push_outer_and_fix_potentials(nb, l, fix_blossom_potential<kInner>(b = nb));
    }
  }

  void expand_free(int bid, int l=kSeparated) {
    int old_base = (heavy[bid] == base[bid]) ? bid : base[bid];
    move_to_largest_blossom(bid);
    assert(old_base > 0);
    for (int beta = old_base, b = beta; ; ) {
      label[b] = l; int nb = node[b].next_b();
      if (b > N && potential[b] == 0) expand_free(b, l);
      else activate_heap2_node(b);
      if ((b = nb) == beta) break;
    }
  }

  bool augment() {
    while (!que.empty()) {
      int x = que.dequeue(), bx = surface[x], lx = label[x], mx = mate[x].v;
      cost_t px = potential[x] + 2;
      for (int eid = ofs[x]; eid < ofs[x + 1]; ++eid) {
        auto& e = edges[eid]; int y = e.to, by = surface[y];
        if (bx == by) continue;
        int ly = label[by];
        cost_t time_eligible = px + potential[y] - scaled(e.cost);
        if (ly > 0) {
          time_eligible >>= 1;
          if (time_eligible == time_current_) {
            if (lx != ly) {
              time_augment_ = time_current_; return true;
            } else {
              contract(x, y); bx = surface[x];
            }
          } else if (within_time_limit(time_eligible)) {
            heap3.emplace(time_eligible, x, y);
          }
        } else {
          if (!within_time_limit(time_eligible)) continue;
          if (ly != kInner) {
            if (time_eligible + lazy[by] == time_current_) {
              grow(x, y);
            } else {
              update_heap2<kFree>(x, y, by, time_eligible);
            }
          } else {
            if (mx != y) {
              update_heap2<kInner>(x, y, by, time_eligible);
            }
          }
        }
      }
    }
    return false;
  }

  bool adjust_dual_variables() {
    // delta1 : rematch
    cost_t time1 = kInf;

    // delta2 : grow
    cost_t time2 = kInf;
    if (!heap2.empty()) time2 = heap2.min().time;

    // delta2.2 : grow (ineligible matched edge)
    cost_t time2_2 = kInf;
    if (!heap2_2.empty()) time2_2 = heap2_2.min();

    // delta3 : contract : O(m log n) time / Edmonds search
    cost_t time3 = kInf;
    while (!heap3.empty()) {
      EdgeEvent e = heap3.min();
      int x = e.from, y = e.to;
      if (surface[x] != surface[y]) {
        time3 = e.time;
        break;
      } else heap3.pop();
    }

    // delta4 : expand
    cost_t time4 = kInf;
    if (!heap4.empty()) time4 = heap4.min();

    // -- events --
    cost_t time_next = min(min(time1, min(time2, time2_2)), min(time3, time4));
    if (!within_time_limit(time_next)) return true;

    // printf("%llu: %lld %lld %lld %lld %lld\n",
    //   time_current_, time1, time2, time2_2, time3, time4);
    assert(time_current_ <= time_next && time_next < kInf);
    time_current_ = time_next;

    while (!heap2_2.empty() && heap2_2.min() == time_current_) {
      int bz = heap2_2.argmin(); heap2_2.pop();
      int l = label[link[bz].from];
      push_outer_and_fix_potentials(bz, l, fix_blossom_potential<kFree>(bz));
      time_created[bz] = time_current_;
    }
    while (!heap2.empty() && heap2.min().time == time_current_) {
      int x = heap2.min().from, y = heap2.min().to, by = surface[y];
      if (label[by] == kFree || label[by] == kSeparated) grow(x, y);
      else heap2.pop(), heap3.emplace(time_current_, x, y);
    }
    while (!heap3.empty() && heap3.min().time == time_current_) {
      int x = heap3.min().from, y = heap3.min().to; heap3.pop();
      if (surface[x] == surface[y]) continue;
      int by = surface[y];
      if (label[by] < 0) {
        int l = label[link[by].from];
        push_outer_and_fix_potentials(by, l, fix_blossom_potential<kInner>(by));
        time_created[by] = time_current_;
      }
      if (label[x] == label[y]) contract(x, y);
      else time_augment_ = time_current_;
    }
    while (!heap4.empty() && heap4.min() == time_current_) {
      int b = heap4.argmin(); heap4.pop();
      expand(b);
    }
    return false;
  }

  // maximal augmenting paths

  void rematch(int v, int w, cost_t c) {
    int t = mate[v].v; mate[v] = {w, c};
    if (mate[t].v != v) return;
    if (link[v].to == dsu.find(link[v].to)) {
      mate[t] = {link[v].from, link[v].cost};
      rematch(mate[t].v, t, mate[t].cost);
    } else {
      int x = link[v].from, y = link[v].to; c = link[v].cost;
      rematch(x, y, c); rematch(y, x, c);
    }
  }

  bool dfs_augment(int x, int bx) {
    const cost_t px = potential[x] + 2; int lx = label[bx];
    for (int eid = ofs[x]; eid < ofs[x + 1]; ++eid) {
      auto& e = edges[eid]; int y = e.to; cost_t c = scaled(e.cost);
      if (px + potential[y] != c) continue;
      int by = dsu.find(y), ly = label[by];
      if (ly > 0) { // outer
        if (lx >= ly) continue;
        int stack_beg = stack_last_;
        for (int bv = by; bv != bx; bv = dsu.find(link[bv].from)) {
          int bw = dsu.find(mate[bv].v);
          stack[stack_last_++] = bw; link[bw] = {x, y, c};
          dsu.par[bv] = dsu.par[bw] = bx;
        }
        while (stack_last_ > stack_beg) {
          int bv = stack[--stack_last_];
          for (int v = blossom.head[bv]; v >= 0; v = blossom.next[v]) {
            if (!dfs_augment(v, bx)) continue;
            stack_last_ = stack_beg;
            return true;
          }
        }
      } else if (ly == kFree) {
        label[by] = kInner; int z = mate[by].v, mz = mate[z].v;
        if (z == 0) { rematch(x, y, c); rematch(y, x, c); return true; }
        if (potential[mz] + potential[z] != mate[by].cost) continue;
        int bz = dsu.find(z);
        link[bz] = {x, y, c}; label[bz] = outer_id_++;
        for (int v = blossom.head[bz]; v >= 0; v = blossom.next[v]) {
          if (dfs_augment(v, bz)) return true;
        }
      }
    }
    return false;
  }

  void set_dsu(int b, int p) {
    if (b <= N) dsu.par[b] = p;
    else {
      for (int beta = base[b], bb = beta; ; ) {
        set_dsu(bb, p);
        if ((bb = node[bb].next_b()) == beta) break;
      }
    }
  }

  int find_maximal() {
    for (int b = N + 1; b < S; ++b) if (base[b] != b) {
      if (surface[b] == b) {
        link_blossom(b, {0, 0, 0});
        int beta = base[b]; while (beta > N) beta = base[beta];
        set_dsu(b, beta);
      }
    }
    for (int u = 1; u <= N; ++u) blossom.head[u] = -1;
    for (int u = 1; u <= N; ++u) {
      if (u == surface[u]) dsu.par[u] = u;
      label[u] = kFree;
      blossom.push(dsu.par[u], u);
    }
    int ret = 0; outer_id_ = 1, stack_last_ = 0;
    for (int u = 1; u <= N; ++u) if (!mate[u].v) {
      int bu = dsu.par[u];
      if (label[bu] != kFree) continue;
      label[bu] = outer_id_++;
      for (int v = blossom.head[bu]; v >= 0; v = blossom.next[v]) {
        if (!dfs_augment(v, bu)) continue;
        ret += 1;
        break;
      }
    }
    // assert(ret >= 1);
    return ret;
  }

  // init

  void initialize() {
    que = Queue<int>(N);
    mate.assign(S, {0, 0});
    link.assign(S, {0, 0, 0});
    label.assign(S, kFree);
    base.resize(S); for (int u = 1; u < S; ++u) base[u] = u;
    surface.resize(S); for (int u = 1; u < S; ++u) surface[u] = u;

    potential.resize(S);
    node.resize(S); for (int b = 1; b < S; ++b) node[b] = Node(b);

    stack.resize(N);
    dsu = DisjointSetUnion(N + 1);
    blossom = LinkedList(N + 1, N + 1);

    unused_bid.resize(B); for (int i = 0; i < B; ++i) unused_bid[i] = N + B - i;
    unused_last_ = B;

    // for O(nm log n) implementation
    reset_time();
    time_created.resize(S);
    slack.resize(S); for (int i = 0; i < S; ++i) slack[i] = kInf;
    best_from.assign(S, 0);
    heavy.assign(S, 0);
    lazy.assign(S, 0);
    group.resize(S); for (int i = 0; i < S; ++i) group[i] = i;
  }

private:
  const int N, B, S; // N = |V|, B = (|V| - 1) / 2, S = N + B + 1
  vector<int> ofs;
  vector<Edge> edges;

  Queue<int> que;
  vector<Mate> mate;
  vector<int> surface, base;
  vector<Link> link;
  vector<int> label;
  vector<cost_t> potential;

  vector<int> unused_bid; int unused_last_;
  vector<Node> node;

  cost_t cost_min_;
  int scale_, contract_count_, outer_id_;
  vector<int> stack; int stack_last_;
  DisjointSetUnion dsu;
  LinkedList blossom;

  // for O(nm log n) implementation
  vector<int> heavy, group;
  vector<cost_t> time_created, lazy, slack;
  vector<int> best_from;

  cost_t time_current_, time_augment_, time_limit_;
  ModifiableHeap<EdgeEvent> heap2;
  ModifiableHeap<cost_t> heap2_2;
  ModifiableHeaps<EdgeEvent> heap2s;
  FastHeap<EdgeEvent> heap3;
  ModifiableHeap<cost_t> heap4;
};

template <typename CostType, typename TotalCostType>
constexpr CostType MaximumWeightPerfectMatching<CostType, TotalCostType>::kInf;
template <typename CostType, typename TotalCostType>
constexpr int MaximumWeightPerfectMatching<CostType, TotalCostType>::kShift;
