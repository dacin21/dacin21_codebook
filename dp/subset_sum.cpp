// Subset sum with integers weights
// Implementation uses std::bitset and compresses equal weights
// For items with total weight S, this runs in O(S^{1.5} / 64)
// (Use Jensen to avoid a factor of log S.)
struct Subset_Sum{
    struct Item{
        int weight;
        int id;
        bool operator<(Item const&o) const {
            return make_pair(weight, id) < make_pair(o.weight, o.id);
        }
    };
    struct Ret{
        bool possible;
        vector<int> ids;
    };
    Subset_Sum(vector<int> w) : Subset_Sum(to_items(move(w))) {}
    Subset_Sum(vector<Item> items_) : n(items_.size()), items(move(items_)) {}

    Ret solve(int goal_weight){
        assert(goal_weight >= 0);
        if(!goal_weight) return Ret{true, {}};
        sort(items.begin(), items.end());
        if(accumulate(items.begin(), items.end(), int64_t{0}, [](int64_t x, Item const&e){return x + e.weight;}) < goal_weight) return ret_impossible();

        return solve_find_impl<>(goal_weight);
    }

private:
    // Used to find a near optimal bitset size
    // Higher values of log_W will exceed the stack limit or hang the compiler.
    template<size_t log_W=29>
    Ret solve_find_impl(int goal_weight){
        if(log_W >= 4 && !(goal_weight >> log_W)){
            return solve_find_impl<log_W - (log_W>=4) >(goal_weight);
        }
        return solve_impl<log_W+1>(goal_weight);
    }
    // implementation in separate function to avoid unused bitsets on the stack frames of solve_find_impl
    template<size_t log_W>
    Ret solve_impl(int goal_weight){
        constexpr int W = 1<<log_W;
        assert(goal_weight < W);
        using bset = bitset<W>;

        bset possible(1);
        vector<pair<int, int> > from(W);

        auto process_step = [&](int weight, int step){
            if(weight*(int64_t)step < W){
                bset added = possible<<(weight*step);
                added &= ~possible;
                possible |= added;
                for(int i=added._Find_first();i<W;i = added._Find_next(i)){
                    from[i] = make_pair(weight, step);
                }
            }
        };
        auto process = [&](int weight, int cnt){
            // compress, keeping at most 2 copies of each power of 2
            for(int step = 1;cnt;){
                cnt-=step;
                process_step(weight, step);
                if(!(cnt&step)) step*=2;
            }
        };
        int run = 0;
        for(int i=0; i<n; ++i){
            ++run;
            // batch equal items
            if(i+1 == n || items[i].weight != items[i+1].weight){
                process(items[i].weight, run);
                run = 0;
            }
        }
        // reconstruct
        if(!possible[goal_weight]) return ret_impossible();
        int cur_weight = goal_weight;
        int i = n-1;
        Ret ret{true, {}};
        while(cur_weight > 0){
            assert(i >= 0);
            auto const&f = from[cur_weight];
            assert(items[i].weight >= f.first);
            if(items[i].weight != f.first) --i;
            else {
                cur_weight -= f.first * f.second;
                assert(i >= f.second-1 && items[i-f.second+1].weight == f.first);
                for(int j=i; j>i-f.second; --j) ret.ids.push_back(items[j].id);
                i-= f.second;
            }
        }
        return ret;
    }
    static vector<Item> to_items(vector<int> w){
        vector<Item> ret(w.size());
        for(int i=0; i<(int)w.size(); ++i){
            ret[i].id = i;
            ret[i].weight = w[i];
        }
        return ret;
    }
    static Ret ret_impossible(){
        return Ret{false, {}};
    }

    int n;
    vector<Item> items;
};
