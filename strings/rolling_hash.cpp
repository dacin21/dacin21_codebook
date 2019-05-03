/*
 *  61-bit rolling hashes class with random bases.
 *  Supports multiple hashes by template parameter.
 */

template<size_t len>
class Rolling_Hash_Base : public array<uint64_t, len>{
private:
    static constexpr uint64_t add(uint64_t a, uint64_t b){
        a+=b+1;
        a = (a&mod) + (a>>61);
        return a-1;
    }
    static constexpr uint64_t sub(uint64_t a, uint64_t b){
        return add(a, mod-b);
    }
    static constexpr uint64_t modmul(uint64_t a, uint64_t b){
        uint64_t l1 = (uint32_t)a, h1 = a>>32, l2 = (uint32_t)b, h2 = b>>32;
        uint64_t l = l1*l2, m = l1*h2 + l2*h1, h = h1*h2;
        uint64_t ret = (l&mod) + (l>>61) + (h << 3) + (m >> 29) + (m << 35 >> 3) + 1;
        ret = (ret & mod) + (ret>>61);
        ret = (ret & mod) + (ret>>61);
        return ret-1;
    }
    static const array<uint64_t, len> base;
    static array<uint64_t, len>const& get_base_pow(int exp){
        assert(exp>=0);
        static vector<array<uint64_t, len> > base_pow;
        if((int)base_pow.size() <= exp){
            if(base_pow.empty()){
                base_pow.reserve(1001);
                base_pow.emplace_back();
                for(auto &e:base_pow.back()) e = 1;
            }
            for(int it=base_pow.size(); it<exp+1000;++it){
                base_pow.push_back(base_pow.back());
                auto &e = base_pow.back();
                for(size_t i=0;i<len;++i){
                    e[i] = modmul(e[i], base[i]);
                }
            }
        }
        return base_pow[exp];
    }

public:
    static constexpr uint64_t mod = (1ull<<61) - 1;

    Rolling_Hash_Base() : array<uint64_t, len>{} {};
    Rolling_Hash_Base(Rolling_Hash_Base const&o) = default;
    Rolling_Hash_Base& operator=(Rolling_Hash_Base const&o) = default;

    static array<uint64_t, len> gen_base(){
        seed_seq seed{(uint32_t) chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count(), (uint32_t)random_device()(), (uint32_t)4730921};
        mt19937 rng(seed);
        array<uint64_t, len> ret;
        for(auto &e:ret) e = uniform_int_distribution<uint64_t>(0, mod-1)(rng);
        return ret;
    }
    static vector<Rolling_Hash_Base> hashify(string const&s){
        vector<Rolling_Hash_Base> ret;
        ret.reserve(s.size()+1);
        ret.push_back(Rolling_Hash_Base());
        for(char const&e:s){
            ret.push_back(ret.back() + e);
        }
        return ret;
    }

    template<typename T, typename enable_if<is_integral<T>::value, int>::type = 0>
    Rolling_Hash_Base& operator+=(T const&o){
        for(size_t i=0;i<len;++i){
            (*this)[i] = add(modmul((*this)[i], base[i]), static_cast<uint64_t>(o));
        }
        ++length;
        return *this;
    }
    template<typename T, typename enable_if<is_integral<T>::value, int>::type = 0>
    Rolling_Hash_Base operator+(T const&o)const{
        Rolling_Hash_Base ret(*this);
        ret+=o;
        return ret;
    }
    Rolling_Hash_Base& operator-=(Rolling_Hash_Base const&o){
        assert(length >= o.length);
        auto const&base_pow = get_base_pow(length - o.length);
        //cerr << length << " - " << o.length << "\n";
        for(size_t i=0;i<len;++i){
            (*this)[i] = sub((*this)[i], modmul(o[i], base_pow[i]));
        }
        length-=o.length;
        return *this;
    }
    Rolling_Hash_Base operator-(Rolling_Hash_Base const&o)const{
        Rolling_Hash_Base ret(*this);
        ret-=o;
        return ret;
    }

    size_t length = 0;
};
template<size_t len>
const array<uint64_t, len> Rolling_Hash_Base<len>::base = Rolling_Hash_Base<len>::gen_base();

using Rolling_Hash = Rolling_Hash_Base<2>;
