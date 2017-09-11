/*
 * fast maximum cardinality matching in O(\sqrt{V} * E)
 * 0.03 sec on www.spoj.com/MATCHING
 * Implementation by (c) Daniel Peter Rutschmann, 2017. All rights reserved...
 *
 */
#include <cstdio>
#include <iostream>
#include <algorithm>
using namespace std;
const int maxn = 2*100015, maxm = 350005;
int N, M, S, T;
int max_h;

typedef unsigned int flow_t;
typedef unsigned int excess_t;

struct Edge {
	int to, next;
	void set(int _to, int _next){
	  to=_to;
	  next=_next;
	}
};
Edge edges[2*maxm];
int edge_list_begin[maxn], lastEdge = 2;
inline void add_edge(int a, int b) {
    b+=N;
	edges[lastEdge].set(b, edge_list_begin[a]);
	edge_list_begin[a] = lastEdge++;
	edges[lastEdge].set(a, edge_list_begin[b]);
	edge_list_begin[b] = lastEdge++;
}
int h[maxn], e[maxn], match[maxn], q[2*maxn], qbegin, qend, pushcnt;
inline void enqueue(int u) {
    if(qend == 2*maxn) qend=0;
    q[qend++]=u;
}
inline void pop() {
    ++qbegin;
    if(qbegin==2*maxn) qbegin=0;
}
inline void global_relabel(){
    int* qq = q;
    if(qend < maxn) qq = q+maxn;
    int l=0, r=0;
    fill(h, h+N+M, max_h);
    fill(match, match+N, -1);
    for(int i=N;i<N+M;++i){
        if(~match[i]) {
            match[match[i]] = i;
        } else {
            qq[r++]=i;
            h[i] = 0;
        }
    }
    while(l!=r){
        int u = qq[l++];
        for(int p = edge_list_begin[u];p;p = edges[p].next){
            int v = edges[p].to;
            if(h[v] == max_h){
                h[v] = h[u]+1;
                if(~match[v]){
                    h[match[v]] = h[u]+2;
                    qq[r++] = match[v];
    }   }   }   }
}
inline void init_matching() {
	fill(h, h + N, 1);
	fill(h+N, h+N+M, 0);
	fill(e, e+N, 1);
	fill(match, match+N+M, -1);
	qbegin = qend = 0;
	for(int i=0;i<N;++i){
        enqueue(i);
	}
	pushcnt = 0;
	max_h = N+M+5;
}
excess_t max_matching() {
	init_matching();
	while (qbegin!=qend) {
		int u = q[qbegin];
		int minh = max_h;
		for (int p = edge_list_begin[u]; p; p = edges[p].next) {
			int v = edges[p].to;
			minh = min(minh, h[v]+1);
			if (h[u] == h[v] + 1) {
                if(~match[v]){
                    enqueue(match[v]);
                    e[match[v]]=1;
                }
                match[v] = u;
                h[v]+=2;
                e[u]=0;
                pop();
                break;
			}
		}
        ++pushcnt;
		if (e[u]) {
            if(minh >= max_h){
                e[u]=0;
                pop();
            }
			h[u] = minh;
		}
		if(pushcnt == 2*N){
            global_relabel();
            pushcnt = 0;
		}
	}
	excess_t hits=0;
	for(int i=N;i<N+M;++i){
        if(~match[i]){
            ++hits;
            match[match[i]] = i;
        }
	}
	return hits;
}

void reset(int lsize, int rsize){
	n = lsize; m = rsize;
	lastEdge = 2;
	fill(edge_list_begin, edge_list_begin+n+m+4, 0);
}
