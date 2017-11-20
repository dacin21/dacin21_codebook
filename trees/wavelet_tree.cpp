// wavelet tree
// supports several orthogonal range queries in O(log N)

template<typename int_t>
struct Wavelet_Tree{
    Wavelet_Tree *l, *r;
    vector<unsigned int> lcnt;
    const int_t low, high;
    const unsigned int size;
    Wavelet_Tree(int_t *begin, int_t *end, int_t const&_low, int_t const&_high, int_t const*sort_begin, int_t const*sort_end):l(0), r(0), low(_low), high(_high), size(end-begin){
        if(size == 0) return;
        if(low == high) return;
        int ia = lower_bound(sort_begin, sort_end, low) - sort_begin;
        int ib = upper_bound(sort_begin, sort_end, high) - sort_begin-1;
        int im = ia+(ib-ia)/2;
        const int_t mid = sort_begin[im]; //low + ((ll)high-low)/2;
        auto is_left = [&](int_t const&a){
            return a<=mid;
        };
        lcnt.resize(size+1);
        lcnt[0] = 0;
        for(size_t i=0;i<size;++i){
            lcnt[i+1] = lcnt[i] + is_left(begin[i]);
        }
        auto midptr = stable_partition(begin, end, is_left);
        l = new Wavelet_Tree(begin, midptr, sort_begin[ia], sort_begin[im], sort_begin, sort_end);
        r = new Wavelet_Tree(midptr, end, sort_begin[im+1], sort_begin[ib], sort_begin, sort_end);
    }
    Wavelet_Tree(vector<int_t> v, vector<int_t> const&v_sorted):Wavelet_Tree(v.data(), v.data()+v.size(), v.empty()?int_t():*min_element(v.begin(), v.end()), v.empty()?int_t():*max_element(v.begin(), v.end()), v_sorted.data(), v_sorted.data()+v_sorted.size()){}
    pair<unsigned int, unsigned int> left_range(unsigned int const&a, unsigned int const&b)const{
        return  make_pair(lcnt[a], lcnt[b+1]-1);
    }
    ~Wavelet_Tree(){delete l; delete r;}
    // k-th smallest element in [a, b]; k is zero based
    int_t kth_element(unsigned int const&a, unsigned int const&b, unsigned int const&k)const{
        assert(k>=0 && k < size);
        if(low == high) return low;
        const unsigned int left_a=left_range(a, b).first, left_b=left_range(a, b).second, left_len = left_b-left_a+1;
        if(left_len > k){
            return l->kth_element(left_a, left_b, k);
        } else {
            return r->kth_element(a-left_a, b-left_b-1, k-left_len);
        }
    }
    // number of occurrences of val in [a, b]
    unsigned int count(unsigned int const&a, unsigned int const&b, int_t const&val)const{
        if(!~b || a>b) return 0u;
        if(val < low) return 0u;
        if(val>high) return 0u;
        if(low == val && high == val)
            return b-a+1;
        const unsigned int left_a=left_range(a, b).first, left_b=left_range(a, b).second;
        return (l?l->count(left_a, left_b, val):0u) + (r?r->count(a-left_a, b-left_b-1, val):0u);
    }
    // number of elements in [a, b] that are <= val
    unsigned int lower_cnt(unsigned int const&a, unsigned int const&b, int_t const&val)const{
        if(!~b || a>b) return 0u;
        if(val < low) return 0u;
        if(val>=high) return b-a+1;
        const unsigned int left_a=left_range(a, b).first, left_b=left_range(a, b).second;
        return (l?l->lower_cnt(left_a, left_b, val):0u) + (r?r->lower_cnt(a-left_a, b-left_b-1, val):0u);
    }
    // number of indices x in [a, b] with y1 <= y[x] <= y2
    unsigned int rectangle_cnt(unsigned int const&x1, unsigned int const&x2, int_t const&y1, int_t const&y2)const{
        if(!~x2 || x1>x2) return 0u;
        if(y1>y2) return 0u;
        if(y1 > high || y2 < low) return 0;
        if(y1<=low && high <=y2) return x2-x1+1;
        const unsigned int left_x1=left_range(x1, x2).first, left_x2=left_range(x1, x2).second;
        return (l?l->rectangle_cnt(left_x1, left_x2, y1, y2):0u) + (r?r->rectangle_cnt(x1-left_x1, x2-left_x2-1, y1, y2):0u);
    }
};
