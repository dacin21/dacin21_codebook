
// Preprocesses in O(n log n) to allow you to query
// the next larger element to the left/right in O(log n)
template<typename Key, typename Value, typename value_comp = less<Value>>
class Cartesian_Ancestors{
public:
    // range maximum query in <n log n, 1>
    template<typename T = Value, typename comp = value_comp>
    class Sparsetable{
        int compind(int i, int j, int l, int sign)const{
            i = RMQ[lgn*i+l], j=RMQ[lgn*(j+sign*(1<<l))+l];
            return comp()(data[i], data[j])? j : i;
        }
        int minind(int l, int r)const{
            assert(l<r);
            return compind(l, r, lg2[r-l],-1);
        }
        void build(){
            lg2.resize(n+1);
            for(int i=2;i<=n;++i) lg2[i]=1+lg2[i/2];
            lgn = lg2.back()+1;
            RMQ.resize(n*lgn);
            for(int i=n-1;i>=0;--i){
                RMQ[lgn*i]=i;
                for(int j=0;j<lg2[n-i];++j) RMQ[lgn*i+j+1] = compind(i,i,j,1);
            }
        }
    public:
        Sparsetable() {}
        Sparsetable(vector<T> v): n(v.size()), data(move(v)) { build(); }
        /// max in [l, r)
        pair<int, T const&> operator()(int const&i, int const&j){
            int k = minind(i, j);
            return pair<int, T const&>(k, data[k]);
        }
        int n, lgn;
        vector<int> RMQ, lg2;
        vector<T> data;
    };
    struct Ret{
        int index;
        Key key;
        Value value;
        bool found;
        friend ostream& operator<<(ostream&o, Ret const&r){
            return o << "(" << r.index << ", " << r.key << ", " << r.value << ", " << r.found << ")";
        }
    };

    Cartesian_Ancestors(vector<Value> values) : n(values.size()), keys(n), data(move(values)), st(data){
        iota(keys.begin(), keys.end(), Key{});
    }
    Cartesian_Ancestors(vector<Key> keys_, vector<Value> values) : n(keys.size()), keys(move(keys_)), data(move(values)), st(data){
        assert(keys.size() == values.size());
        assert(is_sorted(keys.begin(), keys.end()));
    }
    Cartesian_Ancestors(vector<pair<Key, Value> > const&v) : n(v.size()), keys(n), data(n) {
        for(int i=0;i<n;++i){
            tie(keys[i], data[i]) = v[i];
        }
        assert(is_sorted(keys.begin(), keys.end()));
        st = Sparsetable<>(data);
    }

    template<bool allow_equal_key, bool allow_equal_value>
    Ret right_larger(Key const&key, Value const& value){
        auto it = allow_equal_key ? lower_bound(keys.begin(), keys.end(), key) : upper_bound(keys.begin(), keys.end(), key);
        int i = it - keys.begin();
        int l = i-1, r = n;
        while(l+1 < r){
            const int m = l + (r-l)/2;
            auto x = st(i, m+1).second;
            if(comp(x, value)){
                l = m;
            } else if(comp(value, x)){
                r = m;
            } else  {
                if(allow_equal_value){
                    r = m;
                } else {
                    l = m;
                }
            }
        }
        return ret_from_index(r);
    }
    template<bool allow_equal_key, bool allow_equal_value>
    Ret left_larger(Key const&key, Value const& value){
        auto it = allow_equal_key ? upper_bound(keys.begin(), keys.end(), key) : lower_bound(keys.begin(), keys.end(), key);
        int i = it - keys.begin();
        int l = -1, r = i;
        while(l+1 < r){
            const int m = l + (r-l)/2;
            auto x = st(m,i).second;
            if(comp(x, value)){
                r = m;
            } else if(comp(value, x)){
                l = m;
            } else  {
                if(allow_equal_value){
                    l = m;
                } else {
                    r = m;
                }
            }
        }
        return ret_from_index(l);
    }
private:
    Ret ret_from_index(int index){
        if(index == -1 || index == n) return Ret{index, Key{}, Value{}, false};
        return Ret{index, keys[index], data[index], true};
    }
    int n;
    vector<Key> keys;
    vector<Value> data;
    Sparsetable<> st;
    value_comp comp;
};
