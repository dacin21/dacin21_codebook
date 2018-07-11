/*
 *  Various routines that can be used to generate tests
 *  that break rolling hashes with fixed parameters.
 *
 */

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <iterator>
#include <queue>
#include <random>
#include <tuple>
#include <vector>


namespace anti_hash {
    using ll = int64_t;
    using ull = uint64_t;
    template<typename T>
    using ppair = std::pair<T, T>;
    template<typename T>
    using min_heap = std::priority_queue<T, std::vector<T>, std::greater<T>>;
    using Hash_Base = ppair<ll>;

    template <class T, class S, class C>
    void clear_pq(std::priority_queue<T, S, C>& q) {
        struct Hacked_Queue : private std::priority_queue<T, S, C> {
            static S& container(std::priority_queue<T, S, C>& q) {
                return q.*&Hacked_Queue::c;
            }
        };
        Hacked_Queue::container(q).clear();
    }

    ll pow_mod(ll base, ll exp, ll mod){
        assert(0<=base && base<mod);
        ll ret = 1;
        for(;exp;exp>>=1){
            if(exp&1) ret = ret * (__int128) base % mod;
            base = base * (__int128) base % mod;
        }
        return ret;
    }

    std::mt19937 rng(489321);
    template<typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    inline T get_rand(T const&l, T const&r){
        return std::uniform_int_distribution<T>(l, r)(rng);
    }
    template<typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    inline void fill_random(std::vector<T> &v, T const l, T const r){
        for(T& e : v){
            e = get_rand(l, r);
        }
    }

namespace verify{
    ll calc_hash(Hash_Base b, std::string const&s){
        ll ret = 0;
        for(auto &e:s){
            ret = (ret * (__int128) b.second + e) % b.first;
        }
        return ret;
    }
    ull calc_overflow_hash(Hash_Base b, std::string const&s){
        ull ret = 0;
        for(auto &e:s){
            ret = ret *  b.second + e;
        }
        return ret;
    }
}

namespace overflow_break{
    std::vector<std::string> find_collision(Hash_Base const&base, std::vector<std::string> const&alphabet, int h = 11){
        std::vector<std::string> ret(2);
        for(int i=0;i<(1<<h);++i){
            int x = __builtin_popcount(i)%2;
            ret[0]+=alphabet[x];
            ret[1]+=alphabet[!x];
        }
        return ret;
    }


} // overflow_break

namespace birtday_attack {

    std::vector<std::string> find_collision(Hash_Base const&base, std::vector<std::string> const&alphabet, int count){
        std::cerr << "Birthday-attack for hash collision: " << base.first << "/" << base.second << ", count: " << count;
        std::cerr << ", alphabet: " << alphabet.size() << "\n";
        assert(0 <= base.second && base.second < base.first);
        assert(2 <= count);
        const int alpha_length = alphabet[0].size();
        for(auto const&e:alphabet){
            assert(e.size() == alpha_length);
        }
        const ll base_pow = pow_mod(base.second, alpha_length, base.first);
        std::vector<ll> alphabet_hashes;
        for(auto const&e:alphabet){
            ll ha = 0;
            for(char const&c : e){
                ha = (ha * base.second + c) % base.first;
            }
            alphabet_hashes.push_back(ha);
        }
        // estimate by birthday paradoxon
        const ll needed_words = 2*count*pow(static_cast<double>(base.first), 1.0 - 1.0 / (static_cast<double> (count)));
        long double words_factor = 1;
        int length = 1;
        double total_words = static_cast<double>(alphabet.size());

        for(int it=0;;++it){
            ++length;
            total_words*=alphabet.size();
            if(total_words < needed_words){
                --it;
                continue;
            }
            for(int it_2 = 0; it_2 < 2;++it_2){
                if (total_words > 2 * needed_words){
                    words_factor*=1.2;
                }

                const ll used_words = llround(needed_words * words_factor);

                std::cerr << "  length: " << length << ", words: " << used_words << "\n";

                // store generated words in base <alphabet>
                using alpha_t = int32_t;
                std::vector<alpha_t> chars_storage;
                chars_storage.reserve(used_words * length);

                auto calc_hash = [&base, &base_pow, &alphabet_hashes, &alphabet](std::vector<alpha_t> const&word){
                    ll ret = 0;
                    for(alpha_t const&e : word){
                        ret = (ret * base_pow + alphabet_hashes[e]) % base.first;
                    }
                    return ret;
                };
                std::vector<alpha_t> w(length);
                std::vector<std::pair<uint32_t, int> > hashes;
                auto extract_word_string = [&alphabet, &chars_storage, &length](int i)->std::string{
                    std::string ret;
                    for(int j=0;j<length;++j){
                        ret+=alphabet[chars_storage[i*length + j]];
                    }
                    return ret;
                };
                auto extract_with_equal_hash = [&hashes, &chars_storage, &extract_word_string, &count](int l, int r) -> std::vector<std::string>{
                    std::vector<std::string> ret;
                    for(int i=l;i<r;++i){
                        ret.push_back(extract_word_string(hashes[i].second));
                    }
                    std::sort(ret.begin(), ret.end());
                    ret.erase(std::unique(ret.begin(), ret.end()), ret.end());
                    if((int)ret.size() < count){
                        ret.clear();
                    } else {
                        ret.resize(count);
                        std::cerr << "Collision found, total length: " << ret[0].size() << "\n";
                    }
                    return ret;
                };

                for(ll i=0;i<used_words;++i){
                    fill_random<alpha_t>(w, 0, alphabet.size()-1);
                    chars_storage.insert(chars_storage.end(), w.begin(), w.end());
                    hashes.emplace_back(static_cast<uint32_t>(calc_hash(w)), i);
                }
                std::sort(hashes.begin(), hashes.end());

                int equal_cnt = 1;
                for(int i=1;i<(int)hashes.size();++i){
                    if(hashes[i].first != hashes[i-1].first){
                        if(equal_cnt >= count){
                            auto ret = extract_with_equal_hash(i-equal_cnt, i);
                            if(!ret.empty()){
                                return ret;
                            }
                        }
                        equal_cnt = 0;
                    }
                    ++equal_cnt;
                }
                if(equal_cnt >= count){
                    auto ret = extract_with_equal_hash(hashes.size()-equal_cnt, hashes.size());
                    if(!ret.empty()){
                        return ret;
                    }
                }
            }
        }
    }
} // birtday_attack


namespace tree_attack {
    struct Tree_Node{
        const Tree_Node *l;
        union {
            const Tree_Node *r;
            int id;
        };
        ll val;
        bool operator<(Tree_Node const&o)const{
            return val < o.val;
        }
    };
    void merge_nodes(Tree_Node&ret, Tree_Node const*a, Tree_Node const*b){
        ret.l = a;
        ret.r = b;
        ret.val = a->val - b->val;
        if(ret.val < 0) ret.val = -ret.val;
    }
    std::vector<std::string> find_collision(Hash_Base const&base, std::vector<std::string> const&alphabet, int count){
        std::cerr << "Tree-attack for hash collision: " << base.first << "/" << base.second << ", count: " << count;
        std::cerr << ", alphabet: " << alphabet.size() << "\n";
        assert(count == 2); // only support finding a single collision
        assert(0 <= base.second && base.second < base.first);
        assert(2 <= count);

        const int alpha_length = alphabet[0].size();
        for(auto const&e:alphabet){
            assert(e.size() == alpha_length);
        }
        const ll base_pow = pow_mod(base.second, alpha_length, base.first);
        std::vector<ll> alphabet_hashes;
        for(auto const&e:alphabet){
            ll ha = 0;
            for(char const&c : e){
                ha = (ha * (__int128) base.second + c) % base.first;
            }
            alphabet_hashes.push_back(ha);
        }
        const ll hash_diff = alphabet_hashes[1] - alphabet_hashes[0];
        // TODO: find good estimate for initial depth
        int depth = 1;
        for(;;++depth){
            assert(depth < 20);
            const int n = 1<<depth;
            std::cerr << "  depth: " << depth << ", n: " << n << "\n";
            std::vector<Tree_Node> storage(2*n-1);
            ll base_pow_pow = 1;
            for(int i=0;i<n;++i){
                storage[i] = Tree_Node{nullptr, {.id=i}, static_cast<ll>(base_pow_pow * (__int128) hash_diff % base.first)};
                base_pow_pow = base_pow_pow * (__int128) base_pow % base.first;
            }
            auto cur = storage.begin();
            for(int m=n;m>=1;m/=2){
                std::sort(cur, cur+m);
                if(!cur[0].val){ // found collision
                    std::cerr << "Collision found\n";
                    std::vector<int> out(n, 0);
                    std::function<void(Tree_Node const*, int)> extract = [&out, &extract](Tree_Node const*cur, int sign){
                        if(cur->l == nullptr){
                            int pos = cur->id;
                            out[pos] = sign;
                        } else {
                            if(cur->l->val < cur->r->val){
                                sign = -sign;
                            }
                            extract(cur->l, sign);
                            extract(cur->r, -sign);
                        }
                    };
                    extract(&*cur, 1);

                    std::vector<std::string> ret(2);
                    std::reverse(out.begin(), out.end());
                    for(auto const&e:out){
                        ret[0]+= alphabet[std::max(0, e)];
                        ret[1]+= alphabet[std::max(0, -e)];
                    }
                    return ret;
                }
                std::cerr << m << " : " << cur[m-1].val << "\n";
                for(int i=0;i<m/2;++i){
                    merge_nodes(cur[m+i], &*cur + 2*i, &*cur + 2*i+1);
                    //std::cerr << "merge " << 2*i << " " << 2*i+1 << " -> " << m+i << " | ";
                    //std::cerr << cur[2*i].val << " " << cur[2*i+1].val << " " << cur[m+i].val << "\n";
                }
                cur+= m;
            }
        }
    }
} // tree_attack

namespace multi_tree_attack{
    struct Tree_Node{
        const Tree_Node *l;
        union {
            const Tree_Node *r;
            int id;
        };
        std::vector<ll> val;
        long double avg_val;
        bool operator<(Tree_Node const&o)const{
            return avg_val < o.avg_val;
        }
        void recalc_avg(){
            avg_val = std::accumulate(val.begin(), val.end(), (long double)0) / val.size();
            //avg_val = val.back();
            //avg_val = val.front();
        }
    };
    // should be optimizable to O(n log n)
    void extract_smallest_diffs(std::vector<ll> &out, std::vector<ll> const&a, std::vector<ll> const&b, size_t cnt) {
        assert(std::is_sorted(a.begin(), a.end()));
        assert(std::is_sorted(b.begin(), b.end()));
        out.clear();

        static min_heap<std::tuple<ll, int, bool> > pq;
        clear_pq(pq);
        for(int i=0;i<a.size();++i){
            int j = lower_bound(b.begin(), b.end(), a[i]) - b.begin();
            if(j < (int)b.size()){
                pq.emplace(b[j]-a[i], j+1, 0);
            }
        }
        for(int j=0;j<b.size();++j){
            int i = upper_bound(a.begin(), a.end(), b[j]) - a.begin();
            if(i < (int)a.size()){
                pq.emplace(a[i]-b[j], i+1, 1);
            }
        }
        while(!pq.empty() && out.size() < cnt && (out.size() < 2 || out[1]!=0)){
            auto e = pq.top(); pq.pop();
            if(out.empty() || std::get<0>(e) > out.back() || (std::get<0>(e) == 0 && out.size() == 1) ){
                out.push_back(std::get<0>(e));
            }
            if(std::get<2>(e)){
                int i = std::get<1>(e);
                if(i < a.size()){
                    pq.emplace(std::get<0>(e) - a[i-1] + a[i], i+1, 1);
                }
            } else {
                int j = std::get<1>(e);
                if(j < b.size()){
                    pq.emplace(std::get<0>(e) - b[j-1] + b[j], j+1, 0);
                }
            }
        }
        assert(out[0] >= 0ll);
        assert(std::is_sorted(out.begin(), out.end()));
    }
    void merge_nodes(Tree_Node&ret, Tree_Node const*a, Tree_Node const*b, int cut_off){
        ret.l = a;
        ret.r = b;
        extract_smallest_diffs(ret.val, a->val, b->val, cut_off);
        ret.recalc_avg();
    }
    std::vector<std::string> find_collision(Hash_Base const&base, std::vector<std::string> const&alphabet, int count, int thresh = 5){
        std::cerr << "Multi-Tree-attack for hash collision: " << base.first << "/" << base.second << ", count: " << count;
        std::cerr << ", alphabet: " << alphabet.size() << "\n";
        assert(count == 2); // only support finding a single collision
        assert(0 <= base.second && base.second < base.first);
        assert(2 <= count);

        const int alpha_length = alphabet[0].size();
        for(auto const&e:alphabet){
            assert(e.size() == alpha_length);
        }
        const ll base_pow = pow_mod(base.second, alpha_length, base.first);
        std::vector<ll> alphabet_hashes, hash_diff;
        for(auto const&e:alphabet){
            ll ha = 0;
            for(char const&c : e){
                ha = (ha * (__int128) base.second + c) % base.first;
            }
            alphabet_hashes.push_back(ha);
            hash_diff.push_back((base.first + ha - alphabet_hashes.front()) % base.first);
        }

        // TODO: find good estimate for initial depth
        int depth = 1;
        for(;;++depth){
            assert(depth < 20);
            const int n = 1<<depth;
            std::cerr << "  depth: " << depth << ", n: " << n << "\n";
            std::vector<Tree_Node> storage(2*n-1);
            ll base_pow_pow = 1;
            for(int i=0;i<n;++i){
                storage[i] = Tree_Node{nullptr, {.id=i}, {}};
                for(int j=0;j<std::min((int)hash_diff.size(), thresh+1);++j){
                    storage[i].val.push_back(base_pow_pow * (__int128) hash_diff[j] % base.first);
                }
                std::sort(storage[i].val.begin(), storage[i].val.end());
                storage[i].recalc_avg();
                base_pow_pow = base_pow_pow * (__int128) base_pow % base.first;
            }
            std::cerr << "\n";
            auto cur = storage.begin();
            for(int m=n;m>=1;m/=2){
                std::sort(cur, cur+m);

                for(int i=0;i<m;++i){
                    if(!cur[i].val[1]){ // found collision
                        std::cerr << "Collision found\n";
                        std::vector<int> out(n, 0);
                        std::function<void(Tree_Node const*, ll const&, int)> extract = [&out, &extract, &base_pow, &base, &hash_diff, &thresh](Tree_Node const*cur, ll const& val, int sign){

                            if(cur->l == nullptr){
                                int pos = cur->id;
                                ll base_pow_pow = pow_mod(base_pow, pos, base.first);
                                for(int j=0;j<std::min((int)hash_diff.size(), thresh+1);++j){
                                    if(val == static_cast<ll>(base_pow_pow * (__int128) hash_diff[j] % base.first)){
                                        out[pos] = sign * j;
                                        return;
                                    }
                                }
                                assert(0);
                            } else {
                                //for(ll const& a : cur->l->val) {
                                for(auto it = cur->l->val.rbegin(); it!=cur->l->val.rend();++it){
                                    ll const&a = *it;
                                    if(std::binary_search(cur->r->val.begin(), cur->r->val.end(), a - val) ||
                                        std::binary_search(cur->r->val.begin(), cur->r->val.end(), a + val)){
                                        for(ll const&b : cur->r->val){
                                            ll x = a-b;
                                            if(x == val){
                                                extract(cur->l, a, sign);
                                                extract(cur->r, b, -sign);
                                                return;
                                            } else if(-x == val){
                                                extract(cur->l, a, -sign);
                                                extract(cur->r, b, sign);
                                                return;
                                            }
                                        }
                                        assert(0);
                                    }
                                }
                                assert(0);
                            }
                        };
                        extract(&cur[i], 0ll, 1);

                        std::vector<std::string> ret(2);
                        std::reverse(out.begin(), out.end());
                        for(auto const&e:out){
                            ret[0]+= alphabet[std::max(0, e)];
                            ret[1]+= alphabet[std::max(0, -e)];
                        }
                        return ret;
                    }
                }
                std::cerr << m << " " << cur[0].val.size() << "  : " << cur[m-1].val[0] << " " << cur[m-1].avg_val << "\n";
                for(int i=0;i<m/2;++i){
                    merge_nodes(cur[m+i], &*cur + 2*i, &*cur + 2*i+1, thresh);
                    //std::cerr << "merge " << 2*i << " " << 2*i+1 << " -> " << m+i << " | ";
                    //std::cerr << cur[2*i].val << " " << cur[2*i+1].val << " " << cur[m+i].val << "\n";
                }
                cur+= m;
            }
        }
    }

} // multi_tree_attack

    template<size_t alpha, typename std::enable_if<(2 <= alpha), int>::type = 0>
    ppair<std::string> calc_anti_hash(std::vector<Hash_Base> hashes, std::vector<std::string> alphabet){
        assert(!hashes.empty());
        assert(alphabet.size() >= 2);
        // The final call to find_collision is faster, so we want the largest mod last.
        std::sort(hashes.begin(), hashes.end());

        for(int i=0;i<(int)hashes.size();++i){
            int cnt = (i+1 == (int)hashes.size()) ? 2 : alpha;
            alphabet = birtday_attack::find_collision(hashes[i], alphabet, cnt);
            assert(alphabet.size() == cnt);
        }
        assert(alphabet.size() == 2);
        return std::make_pair(alphabet[0], alphabet[1]);
    }
    template<size_t alpha, typename std::enable_if<(2 <= alpha), int>::type = 0>
    ppair<std::string> calc_anti_hash(std::vector<Hash_Base> hashes, std::string const&alphabet){
        assert(!hashes.empty());
        assert(alphabet.size() >= 2);
        std::vector<std::string> alphabet_2;
        for(auto &e:alphabet){
            alphabet_2.push_back(std::string(1, e));
        }
        return calc_anti_hash<alpha>(std::move(hashes), std::move(alphabet_2));
    }

    const std::string lower = "qwertyuiopasdfghjklzxcvbnm";
    const std::string upper = "QWERTYUIOPASDFGHJKLZXCVBNM";
    const std::string nums = "1234567890";
} // anti_hash

#include <bits/stdc++.h>
using namespace std;


signed rand_test(){
    const int p = 1e9+7;
    for(int i=0;i<(1<<9);++i){
        anti_hash::calc_anti_hash<2>(std::vector<anti_hash::Hash_Base>{make_pair(p, 10000+i)}, anti_hash::lower);
        if(i == ((i>>6)<<6)) cout << i << "\n";
    }
    return 0;
}

signed tree_test(){
    anti_hash::Hash_Base b;
    cin >> b.first >> b.second;
    vector<string> alpha;
    alpha.emplace_back(1, 'A');
    alpha.emplace_back(1, 'B');
    auto res = anti_hash::tree_attack::find_collision(b, alpha, 2);
    assert(res.size() == 2);
    cout << res[0] << " : " << anti_hash::verify::calc_hash(b, res[0]) << "\n";
    cout << res[1] << " : " << anti_hash::verify::calc_hash(b, res[1]) << "\n";
    return 0;
}
signed multi_tree_test(){
    anti_hash::Hash_Base b;
    cin >> b.first >> b.second;
    vector<string> alpha;
    for(int i=0;i<2;++i){
        alpha.emplace_back(1, 'A' + i);
    }
    /*for(int i=0;i<26;++i){
        alpha.emplace_back(1, 'a' + i);
    }*/
    auto res = anti_hash::multi_tree_attack::find_collision(b, alpha, 2, 100000);
    assert(res.size() == 2);
    cout << res[0] << " : " << anti_hash::verify::calc_hash(b, res[0]) << "\n";
    cout << res[1] << " : " << anti_hash::verify::calc_hash(b, res[1]) << "\n";
    return 0;
}

signed multi_tree_rand_test(){
    anti_hash::Hash_Base b;
    cin >> b.first >> b.second;
    vector<string> alpha;
    for(int i=0;i<26;++i){
        alpha.emplace_back(1, 'A' + i);
    }
    for(int i=0;i<26;++i){
        alpha.emplace_back(1, 'a' + i);
    }
    int64_t total = 0;
    for(int i=0;i<(1<<15);++i){
        auto res = anti_hash::multi_tree_attack::find_collision(b, alpha, 2, 10000);
        //auto res = anti_hash::tree_attack::find_collision(b, alpha, 2);
        assert(res.size() == 2);
        assert(anti_hash::verify::calc_hash(b, res[0]) == anti_hash::verify::calc_hash(b, res[1]));
        total += res[0].size() + res[1].size();
        if((i>>8<<8) == i) cout << " " << i << " " << (i * 1.0 / (1<<15)) << " " << total << endl;
        ++b.second;
    }
    cout << " => " << total << endl;
    return 0;
}

signed overflow_test(){
    anti_hash::Hash_Base b;
    cin >> b.second;
    vector<string> alpha;
    alpha.emplace_back(1, 'A');
    //alpha.back().push_back('X');
    alpha.emplace_back(1, 'A'+1);
    //alpha.back().push_back('Y');
    for(int it=0;it<10000000;++it){
        if(it < 6){
            unsigned long long x  = b.second;
            int sum = 0;
            for(int i=0;i<10;++i){
                cout << __builtin_ctzll(x-1) << " ";
                sum+=__builtin_ctzll(x-1);
                x*=x;
            }
            cout << " : " << sum << "\n";
        }
        auto res = anti_hash::overflow_break::find_collision(b, alpha, 10);
        assert(res.size() == 2);
        assert(anti_hash::verify::calc_overflow_hash(b, res[0]) == anti_hash::verify::calc_overflow_hash(b, res[1]));
        //cout << res[0] << " : " << anti_hash::verify::calc_overflow_hash(b, res[0]) << "\n";
        //cout << res[1] << " : " << anti_hash::verify::calc_overflow_hash(b, res[1]) << "\n";
        b.second+=2;
    }
    return 0;
}
