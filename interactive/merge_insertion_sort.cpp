// merge-insertion-sort
// Sorting in O(n^2) time with near-optimal number of comparisons
// Useful for interactive problems
// Number of comparisons used is: n lg n - 1.415 n
// The lower bound is: lg (n!) = n lg n - 1.443 n
// Binary search insertion sort would need: n log n - n
struct Merge_Insertion_Sort{
    template<typename F>
    static void apply_permutation(vector<int> const&p, F get){
        const int n = p.size();
        vector<int> q(n);
        for(int i=0;i<n;++i) q[p[i]] = i;
        for(int i=0;i<n;++i){
            while(q[i] != i){
                swap(get(i), get(q[i]));
                swap(q[i], q[q[i]]);
            }
        }
    }
    // ret.first is the sorted vector
    // ret.second[i] is the index of the i-th smallest element in the original vector
    // i.e. ret.second is the permutation that was applied to sort
    template<typename T, typename F>
    static pair<vector<T>, vector<int> > sort(vector<T> v, F comp){
        const int n = v.size();
        if(n <= 1) return {move(v), {{0}}};
        vector<int> preperm(n);
        iota(preperm.begin(), preperm.end(), 0);
        const int M = n-n/2;
        for(int i=0;i<n/2;++i){
            if(comp(v[M+i], v[i])){
                swap(v[M+i], v[i]);
                swap(preperm[M+i], preperm[i]);
            }
        }
        auto ret = sort(vector<T>(v.begin(), v.begin()+n/2), comp);
        apply_permutation(ret.second, [&preperm,M](int const&i)->int&{return preperm[M+i];});
        apply_permutation(ret.second, [&preperm](int const&i)->int&{return preperm[i];});
        apply_permutation(ret.second, [&v,M](int const&i)->T&{return v[M+i];});
        iota(ret.second.begin(), ret.second.end(), 0);
        // insert one element without comparisons
        ret.first.push_back(v.back());
        ret.second.push_back(n-1);
        // now insert the rest in blocks that optimize the binary search sizes
        for(int it=1,r=n-1,s=2; r>n/2; ++it, r-=s, s=(1<<it)-s){
            for(int i=r-s;i<r;++i){
                if(i>=n/2){
                    int a = find(ret.second.begin(), ret.second.end(), i-M) - ret.second.begin();
                    int b = ret.first.size();
                    if(a==b) {
                        assert(i==M-1);
                        a = -1;
                    }
                    while(a+1 < b){
                        const int m = a+(b-a)/2;
                        if(comp(ret.first[m], v[i])){
                            a = m;
                        } else {
                            b = m;
                        }
                    }
                    ret.first.insert(ret.first.begin()+b, v[i]);
                    ret.second.insert(ret.second.begin()+b, i);
                }
            }
        }
        // compose permutations
        apply_permutation(ret.second, [&preperm](int const&i)->int&{return preperm[i];});
        ret.second.swap(preperm);
        return ret;
    }
    template<typename T>
    static pair<vector<T>, vector<int> > sort(vector<T> v){
        return sort(move(v), std::less<T>{});
    }
};
