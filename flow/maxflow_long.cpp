template<typename flow_t = int, typename excess_t = long long>
struct Flow{
    struct Edge{
        int to, next; flow_t f;
        void set(int _to, flow_t _f, int _n){
            tie(to, next, f)=make_tuple(_to, _n, _f);
        }
    };
    int n, ledge, s, t, *fe, *h, *hp, *cnt, htop;
    Edge* ed;
    excess_t *ex, flow;
    bool *isq;
    void add_edge(int a, int b, flow_t cap, flow_t rcap = flow_t(0)){
        assert(a>0&&b>0&&a<n&&b<n);
        if(a==b) return;
        ed[ledge].set(a, rcap, fe[b]);
        fe[b] = ledge++;
        ed[ledge].set(b, cap, fe[a]);
        fe[a] = ledge++;
    }
    void push(int u, int v, int p) {
        flow_t d = ex[u] < ed[p].f ? ex[u] : ed[p].f;
        ed[p].f -= d, ed[p^1].f += d;
        ex[u] -= d, ex[v] += d;
    }
    inline void relabel(int u) {
        if(--cnt[h[u]]){
            int min_h = 2*n;
            for (int p = fe[u]; p; p = ed[p].next)
                if (ed[p].f > 0 && min_h > h[ed[p].to])
                    min_h = h[ed[p].to];
            ++cnt[h[u] = min_h + 1];
        } else {//gap
            for(int i=1;i<n;++i){
                if(i!=s && h[i]>h[u]){
                    --cnt[h[i]];
                    ++cnt[h[i]=h[s]+1];
                }
            }
            ++cnt[h[u]=h[s]+1];
        }
    }
    void hpush(int u) {
        int i;
        for (i = ++htop; h[hp[i >> 1]] < h[u]; i >>= 1)
            hp[i] = hp[i >> 1];
        hp[i] = u;
    }
    void hpop() {
        int last = hp[htop--], i, ch;
        for (i = 1; (i << 1) <= htop; i = ch) {
            ch = i << 1;
            if (ch != htop && h[hp[ch]] < h[hp[ch + 1]]) ++ch;
            if (h[last] >= h[hp[ch]]) break;
            else  hp[i] = hp[ch];
        }
        hp[i] = last;
    }
    inline void init() {
        fill(h, h + n, n+2);
        fill(cnt, cnt + 2*n + 2, 0);
        h[hp[0]=0] = ~0u >> 1;
        int f = 1, b = 1;
        cnt[h[s] = n+1] = cnt[h[t] = 0] = 1;
        hp[b++] = t;
        while (f < b) {
            int u = hp[f++];
            for (int p = fe[u]; p; p = ed[p].next) {
                int v = ed[p].to;
                if (h[v] == n+2 && ed[p^1].f > 0) {
                    ++cnt[h[v] = h[u] + 1];
                    hp[b++] = v;
                }
            }
        }
        flow = htop = 0;
        fill(ex, ex + n, 0);
        fill(isq, isq + n, 0);
        ex[s] = ~0llU >> 1;
        for (int p = fe[s]; p != 0; p = ed[p].next) {
            if (ed[p].f > 0) {
                int v = ed[p].to;
                push(s, v, p);
                if(!isq[v]) {
                    hpush(v);
                    isq[v] = true;
                }
            }
        }
    }
    excess_t max_flow() {
        init();
        while (htop) {
            int u = hp[1];
            if (u == t) {
                hpop();
            } else {
                for (int p = fe[u]; p; p = ed[p].next) {
                    int v = ed[p].to;
                    if (ed[p].f > 0 && h[u] == h[v] + 1) {
                        push(u, v, p);
                        if (isq[v] == false) {
                            if (v != s && v != t)
                                hpush(v);
                            isq[v] = true;
                        }
                        if (ex[u] == 0) {
                            hpop();
                            isq[u] = false;
                            break;
                        }
                    }
                }
                if (ex[u]) relabel(u);
            }
        }
        return ex[t];
    }
    Flow(int N, int M, int S, int T):n(N+1), ledge(2), s(S), t(T), fe(new int[n]), h(new int[n]), hp(new int[n+2]), cnt(new int[2*n+2]), htop(0), ed(new Edge[2*M+2]), ex(new excess_t[n]), isq(new bool[n]){fill(fe, fe+n, 0);}
    Flow(int N, int M):Flow(N+2, M, N+1, N+2){}
    ~Flow(){delete[] fe, delete[]h, delete[]hp, delete[]cnt, delete[]ed, delete[]ex; delete[]isq;}
};