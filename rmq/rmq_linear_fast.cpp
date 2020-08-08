/*
 * range minimum query with
 *   O(n) build time
 *   O(1) query time
 * assuming log(n)-sized words
 *
 * based on: https://codeforces.com/blog/entry/78931
 */

template<typename T, typename Comp = less<T> >
class RMin : Comp {
public:
    static constexpr int lgB = 6;
    static constexpr int B = 1<<lgB;

    static constexpr int last_bit(uint64_t x){
        return 1+__builtin_ctzll(x);
    }
    static constexpr int first_bit(uint64_t x){
        // return 1+__lg(x);
        return 64 - __builtin_clzll(x);
    }

    RMin(vector<T> v) : data(move(v)){
        build();
    }
    RMin(T const* v, int n_) : data(v, v+n_){
        build();
    }
    /*
     * index of minimum in [i, j)
     * tie break: arbitrary
     */
    int min_index(int i, int j){
        if(j-i <= B){
            return small(j, j-i);
        }
        int ret = smaller_ind(small(i+B), small(j));
        i = (i>>lgB)+1;
        j = j>>lgB;
        if(i != j){
            const int k = first_bit(j-i)-1;
            ret = smaller_ind(ret, sparse[k*nb + i]);
            ret = smaller_ind(ret, sparse[k*nb + j-(1<<k)]);
        }
        return ret;
    }
    /*
     * smallest value in [i, j)
     */

    T const& min(int const i, int const j){
        return data[min_index(i, j)];
    }

    T const& at(int const&i) const {
        return data[i];
    }
private:
    bool comp(int i, int j){
        return Comp::operator()(data[i], data[j]);
    }
    int smaller_ind(int i, int j){
        return comp(j, i) ? j : i;
    }
    // [r-d, r)
    int small(int r, int d=B){
        return r - first_bit(mask[r] & (~0ull>>(B-d)));
    }

    void build(){
        n = data.size();
        nb = n / B;
        mask.resize(n+1);
        //if(n) mask[1] = 1;
        uint64_t m = 0;
        for(int i=1; i<n; ++i){
            m = m<<1 | 1u;
            mask[i] = m;
            if(comp(i, i-1)) {
                m^=1u;
                while(m && comp(i, i-last_bit(m))){
                    m^= m&-m;
                }
            }
        }
        mask[n] = m<<1 | 1u;
        sparse.resize(n);
        for(int i=0; i<nb; ++i){
            sparse[i] = small(B*(i+1), B);
        }
        for(int k=1; (1<<k) < n; ++k){
            for(int i=k*nb, i_end = k*nb + nb-(1<<k); i<=i_end; ++i){
                sparse[i] = smaller_ind(sparse[i-nb], sparse[i-nb+(1<<(k-1))]);
            }
        }
    }

    int n, nb;
    vector<uint64_t> mask;
    vector<T> data;
    vector<int> sparse;
};
