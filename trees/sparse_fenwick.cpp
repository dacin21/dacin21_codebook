/*
 *	off-line sparse Fenwick tree
 *	point updates, range queries in O(log(n)^2)
 *	hopefully the compiler optimizes the templates...
 */
template<typename T, size_t dim>
enable_if_t<dim!=0, array<T, dim-1>> pop_front(array<T, dim> const&a){
    array<T, dim-1> ret;
    copy_n(a.begin()+1, dim-1, ret.begin());
    return ret;
}
template<typename T, int dim>
struct Fenwick_Tree{
    vector<Fenwick_Tree<T, dim-1>> trees;
    vector<int> coords;
    vector<array<int, dim> > ups;
    Fenwick_Tree(){}
    unsigned int compress(int const&coord){
        return upper_bound(coords.begin(), coords.end(), coord) - coords.begin();
    }
    void fake_update(array<int, dim> const&pos){
        coords.push_back(pos[0]);
        ups.push_back(pos);
    }
    void compile(){
        sort(coords.begin(), coords.end());
        coords.erase(unique(coords.begin(), coords.end()), coords.end());
        trees.resize(coords.size());
        for(auto &e:ups){
            for(unsigned int x = compress(e[0]);x<=trees.size();x+=x&(-x)){
                trees[x-1].fake_update(pop_front(e));
            }
        }
        ups.clear();
        for(auto &e:trees) e.compile();
    }
    T q(array<int, dim> const&pos){
        T ret = 0; array<int, dim-1> newpos = pop_front(pos);
        for(unsigned int x = compress(pos[0]);x;x-=(x&-x)){
            ret+=trees[x-1].q(newpos);
        }
        return ret;
    }
    T qq(array<int, dim> const&l, array<int, dim> const&r){
        --l[0];
        return q(r)-q(l);
    }
    void up(array<int, dim> const&pos, T const&val){
        array<int, dim-1> newpos =pop_front(pos);
        for(unsigned int x = compress(pos[0]);x<=trees.size();x+=(x&-x)){
            trees[x-1].up(newpos, val);
        }
    }
};
template<typename T>
struct Fenwick_Tree<T, 0>{
    T data;
    Fenwick_Tree():data(){}
    void fake_update(array<int, 0>){}
    void compile(){}
    T q(array<int, 0>){return data;}
    T qq(array<int, 0>, array<int, 0>){return data;}
    void up(array<int, 0>, T const&val){data+=val;}
};


signed main(){
    Fenwick_Tree<int, 2> ft;
    ft.fake_update({1, 4});
    ft.fake_update({2, 2});
    ft.compile();
    cerr << ft.q({3, 3}) << "\n";
    cerr << ft.q({4, 4}) << "\n";
    ft.up({2, 2}, 1);
    cerr << ft.q({2, 1}) << "\n";
    cerr << ft.q({2, 2}) << "\n";
    cerr << ft.q({1, 2}) << "\n";
    cerr << ft.q({4, 4}) << "\n";
    ft.up({1, 4}, 10);
    cerr << ft.q({4, 4}) << "\n";
    cerr << ft.q({3, 4}) << "\n";
    cerr << ft.q({1, 4}) << "\n";
    cerr << ft.q({0, 4}) << "\n";
    cerr << ft.q({4, 3}) << "\n";

    cerr << "done.\n";
}
