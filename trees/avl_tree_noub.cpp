#include <bits/stdc++.h>
using namespace std;

// balanced avl-tree (height may differ by at most 1)
class avlTree{
public:
    struct node{
        node*l,*r;
        int x, height;
        int subsize;
        node(int _x):l(0),r(0),x(_x),height(1),subsize(1){}
        node* recalc(){
            subsize = 1+(l?l->subsize:0)+(r?r->subsize:0);
            height = 1+max(l?l->height:0, r?r->height:0);
            return this;
        }

    }*root;
    avlTree():root(0){}
private:
    /// balanced at most once
    static node* balance(node* a, node* l, node* r){
        if((l?l->height:0) > (r?r->height:0)+1){
            if((l->l?l->l->height:0) >= (l->r?l->r->height:0)){
                a->l = l->r; a->r = r; a->recalc();
                l->r = a; return l->recalc();
            }
            a->l = l->r->r; a->r = r; a->recalc();
            node* ret = l->r;
            l->r = l->r->l; l->recalc();
            ret->r = a; ret->l = l; return ret->recalc();
        } else if((r?r->height:0) > (l?l->height:0)+1){
            if((r->r?r->r->height:0) >= (r->l?r->l->height:0)){
                a->r = r->l; a->l = l; a->recalc();
                r->l = a; return r->recalc();
            }
            a->l = l; a->r = r->l->l; a->recalc();
            node* ret = r->l;
            r->l = r->l->r; r->recalc();
            ret->r = r; ret->l = a; return ret->recalc();
        }
        a->l = l; a->r = r; return a->recalc();
    }
    /// balances at most once
    static node* flat_merge(node*a, node*l, node*r){
        return balance(a, l, r);
    }
    /// inserts singleton-node val into a
    static node* insert(node*a, node* val){
        if(a==0) return val;
        if(val->x < a->x){
            return balance(a, insert(a->l, val), a->r);
        } else {
            return balance(a, a->l, insert(a->r, val));
        }
    }

    static node* push_front(node*a, node*val){
        if(a==0) return val;
        return balance(a, push_front(a->l, val), a->r);
    }

    static node* push_back(node*a, node*val){
        if(a==0) return val;
        return balance(a, a->l, push_back(a->r, val));
    }

    static node* merge(node*a, node*l, node*r){
        if(l==0) return push_front(r, a);
        if(r==0) return push_back(l, a);
        if(l->height > r->height+1){
            return balance(l, l->l, merge(a, l->r, r));
        }
        if(r->height > l->height+1){
            return balance(r, merge(a, l, r->l), r->r);
        }
        a->l = l; a->r = r; return a->recalc();
    }
    static node* front(node*a){
        while(a->l) a=a->l;
        return a;
    }
    static node* back(node*a){
        while(a->r) a=a->r;
        return a;
    }
    static node* pop_front(node*a){
        if(a->l) return balance(a, pop_front(a->l), a->r);
        return a->r;
    }
    static node* pop_back(node*a){
        if(a->r) return balance(a, a->l, pop_back(a->r));
        return a->l;
    }

    static node* merge(node*a, node*b){
        if(a==0) return b;
        if(b==0) return a;
        node*m = front(b);
        return merge(m, a, pop_front(b));
    }

    static node* flat_merge(node*a, node*b){
        if(a==0) return b;
        if(b==0) return a;
        node*m = front(b);
        return balance(m, a, pop_front(b));
    }
    static pair<node*, node*> split(node*a, int pos){
        if(a==0){
            assert(pos==0);
            return make_pair((node*)0, (node*)0);
        }
        if(pos <= (a->l?a->l->subsize:0)){
            pair<node*, node*> L = split(a->l, pos);
            return make_pair(L.first, merge(a, L.second, a->r));
        }
        pair<node*, node*> R = split(a->r, pos -1 -(a->l?a->l->subsize:0));
        return make_pair(merge(a, a->l, R.first), R.second);
    }
    static pair<node*, node*> splitval(node*a, int val){
        if(a==0){
            return make_pair((node*)0, (node*)0);
        }
        if(val < a->x){
            pair<node*, node*> L = splitval(a->l, val);
            return make_pair(L.first, merge(a, L.second, a->r));
        }
        pair<node*, node*> R = splitval(a->r, val);
        return make_pair(merge(a, a->l, R.first), R.second);
    }
    static node* erase(node*a, int pos){
        if(pos == (a->l?a->l->subsize:0)){
            node*ret =  flat_merge(a->l, a->r);
            delete a;
            return ret;
        }
        if(pos < (a->l?a->l->subsize:0)) return balance(a, erase(a->l, pos), a->r);
        return balance(a, a->l, erase(a->r, pos -1 -(a->l?a->l->subsize:0)));
    }
	static node* eraseval(node*a, int val){
        if(val==a->x){
            node*ret =  merge(a->l, a->r);
            delete a;
            return ret;
        }
        if(val<a->x) return balance(a, eraseval(a->l, val), a->r);
        return balance(a, a->l, eraseval(a->r, val));
    }
    static void destroy(node*a){
        if(a==0) return ;
        destroy(a->l);
        destroy(a->r);
        delete a;
    }
public:
    static void check_balance(node * cur){
        if(cur->l) check_balance(cur->l);
        if(cur->r) check_balance(cur->r);
        int lH = (cur->l?cur->l->subsize:0);
        int rH = (cur->r?cur->r->subsize:0);
        assert(lH<=rH+1 && rH<=lH+1);
    }
    static void traverse(ostream&o, node*cur){
        if(cur==0) return;
        o << " ( ";
        if(cur->l) traverse(o, cur->l);
        o << cur->x << "[" << cur->height << "," << cur->subsize << "]";
        if(cur->r) traverse(o, cur->r);
        o << ")";
    }
    void insert(int x){
        root = insert(root, new node(x));
    }
    void erase(int pos){
        root = erase(root, pos);
    }
    void eraseval(int val){
        root = eraseval(root, val);
    }
    int lower_count(int x){
        int ret=0;
        for(node*cur = root;cur;){
            if(x<cur->x) cur=cur->l;
            else {
                ret+=1 + (cur->l?cur->l->subsize:0);
                cur=cur->r;
        }   }
        return ret;
    }
    void destroy(){
        destroy(root);
		root = 0;
    }
};
