// Arbitrary precision integers
// (c) Daniel Rutschmann aka dacin21 2018 -

#ifndef BIGNUM_HPP
#define BIGNUM_HPP
#include <vector>
#include <algorithm>
#include <bitset>
#include "utility.hpp"
using namespace std;

/**
 *  signed Bigint in base 10^18, used for Input / Output
 *  don't use for computations, doesn't support most operations
 *
 */
class Bigint_base10{
public:
    static constexpr int64_t BASE = 1e18;
    static constexpr int DIGITS = 18;
private:
    bool is_neg;
    vector<int64_t> data;

public:
    Bigint_base10():is_neg(false), data(1, 0){}
    Bigint_base10(int64_t const&val):is_neg(val<0){
        int64_t abs_val = abs(val);
        if(abs_val < BASE){
            data = {abs_val};
        } else {
            data = {abs_val%BASE, abs_val/BASE};
        }
    }

    Bigint_base10 operator+(Bigint_base10 const&o)const{
        assert_msg(is_neg == o.is_neg, "Addition operands need to have equal sign");
        Bigint_base10 ret;
        ret.is_neg = is_neg;
        ret.data.assign(1+max(data.size(), o.data.size()), 0);
        copy(data.begin(), data.end(), ret.data.begin());
        int64_t carry = 0;
        for(unsigned int i=0;i<o.data.size();++i){
            ret.data[i]+=o.data[i] + carry;
            carry = 0;
            if(ret.data[i] >= BASE){
                carry = 1;
                ret.data[i]-=BASE;
            }
        }
        for(unsigned int i=o.data.size();carry;++i){
            ret.data[i]+=carry;
            carry = 0;
            if(ret.data[i] >= BASE){
                carry = 1;
                ret.data[i]-=BASE;
            }
        }
        return ret.trim();
    }

    Bigint_base10 operator*(int64_t const&o)const{
        if(o == 0){
            return Bigint_base10(0);
        }
        if(o<0){
            return operator*(-o).negate();
        }
        if(o&1){
            return operator+(operator*(o-1));
        }
        return operator+(*this)*(o/2);
    }
    Bigint_base10& operator+=(Bigint_base10 const&o){
        *this = operator+(o);
        return *this;
    }
    Bigint_base10& operator*=(int64_t const&o){
        *this = operator*(o);
        return *this;
    }

    Bigint_base10& trim(){
        while(data.size()>1 && data.back() == 0){
            data.pop_back();
        }
        return *this;
    }

    bool is_zero()const{
        for(auto const&e:data) if(e) return false;
        return true;
    }

    Bigint_base10& negate(){
        is_neg = !is_neg;
        if(is_zero()) is_neg = false;
        return *this;
    }

    friend ostream& operator<<(ostream&o, Bigint_base10 const&b){
        if(b.is_neg) o << '-';
        o << b.data.back();
        o << setfill('0');
        for(auto it = next(b.data.rbegin());it != b.data.rend();++it){
            o << setw(9) << *it;
        }
        o << setw(0) << setfill(' ');
        return o;
    }
    friend istream& operator>>(istream&in, Bigint_base10 &b){
        static string tmp;
        in >> tmp;
        assert_msg(in, "input should be readable as a string");
        if(tmp[0] == '-'){
            b.is_neg = true;
            tmp = tmp.substr(1, -1);
        } else {
            b.is_neg = false;
        }
        assert_msg(all_of(tmp.begin(), tmp.end(), [](char const&c){return '0'<=c && c<='9';}), "Input should consist of digits and possibly a '-'");
        assert_msg(!tmp.empty(), "Input should contain at least one digit");

        b.data.resize((tmp.size()+DIGITS-1)/DIGITS);
        unsigned int i, j;
        for(i=tmp.size()-DIGITS, j=0;i>0;i-=DIGITS, ++j){
            b.data[j] = stoll(tmp.substr(i, DIGITS));
        }
        b.data[j] = stoll(tmp.substr(0, i+DIGITS));
        return in;
    }
};


/**
 *  Biginteger with fixed precision
 *  Has 31*len bits, first bit is sign
 *
 *  Is quite fast
 */
template<size_t len>
struct Bigint_Fixedsize{
    unsigned int data[len];
    static constexpr unsigned int bits = 31;
    Bigint_Fixedsize(){memset(data, 0, sizeof(data));}
    Bigint_Fixedsize(long long const&_a){
        memset(data, 0, sizeof(data));
        unsigned long long a = _a;
        data[0] = a&((1u<<bits)-1);
        data[1] = a>>bits;
        data[1]&=~(1u<<bits);
        if(a>~a){ // negative number, use complement
            for(size_t i=2;i<len;++i){
                data[i] = (1u<<bits)-1;
            }
        }
    }
    //__attribute__((optimize("unroll-loops")))
    Bigint_Fixedsize& operator+=(Bigint_Fixedsize const&o){
        unsigned int carry = 0;
        for(size_t i=0;i<len;++i){
            data[i]+=o.data[i] + carry;
            carry = data[i]>>bits;
            data[i]&=~(1u<<bits);
        }
        return *this;
    }
    //__attribute__((optimize("unroll-loops")))
    Bigint_Fixedsize& operator-=(Bigint_Fixedsize const&o){
        unsigned int carry = 0;
        for(size_t i=0;i<len;++i){
            data[i]-=o.data[i] + carry;
            carry = data[i]>>bits;
            data[i]&=~(1u<<bits);
        }
        return *this;
    }
    Bigint_Fixedsize operator+(Bigint_Fixedsize const&o)const{
        Bigint_Fixedsize ret(*this);
        ret+=o;
        return ret;
    }
    Bigint_Fixedsize operator-(Bigint_Fixedsize const&o)const{
        Bigint_Fixedsize ret(*this);
        ret-=o;
        return ret;
    }
    //__attribute__((optimize("unroll-loops")))
    void multo(Bigint_Fixedsize const&o, Bigint_Fixedsize &ret)const{
        static unsigned int tmp[len+1];
        memset(tmp, 0, sizeof(tmp));
        for(size_t i=0;i<len;++i){
            unsigned long long val = 0;
            for(size_t j=0;j<len-i;++j){
                val+= data[i]*(unsigned long long)o.data[j]+tmp[i+j];
                tmp[i+j] = val&((1u<<bits)-1);
                val>>=bits;
            }
        }
        memcpy(ret.data, tmp, sizeof(ret.data));
    }
    Bigint_Fixedsize& operator*=(Bigint_Fixedsize const&o){
        multo(o, *this);
        return *this;
    }
    Bigint_Fixedsize operator*(Bigint_Fixedsize const&o)const{
        Bigint_Fixedsize ret;
        multo(o, ret);
        return ret;
    }
    Bigint_Fixedsize& negate(){
        unsigned int carry = 0;
        for(size_t i=0;i<len;++i){
            data[i] = ~data[i] + !carry;
            carry = (data[i]>>bits);
            data[i]&=~(1u<<bits);
        }
        return *this;
    }

    Bigint_Fixedsize operator-()const{
        Bigint_Fixedsize ret(*this);
        ret.negate();
        return ret;
    }
    bool operator<(Bigint_Fixedsize const&o)const{
        // treat sign bit
        if(data[len-1]>>(bits-1)!=o.data[len-1]>>(bits-1)){
            return data[len-1]>>(bits-1);
        }
        for(size_t i=len-1;~i;--i){
            if(data[i]!=o.data[i]) return data[i]<o.data[i];
        }
        return false;
    }
    bool operator>(Bigint_Fixedsize const&o)const{
        return o<*this;
    }
    bool operator<=(Bigint_Fixedsize const&o)const{
        return !(operator>(o));
    }
    bool operator>=(Bigint_Fixedsize const&o)const{
        return !(operator<(o));
    }
    bool operator==(Bigint_Fixedsize const&o)const{
        for(size_t i=0;i<len;++i){
            if(data[i] !=o.data[i]) return false;
        }
        return true;
    }
    bool operator!=(Bigint_Fixedsize const&o)const{
        return !operator==(o);
    }
    bool operator!()const{
        for(size_t i=0;i<len;++i){
            if(data[i]) return false;
        }
        return true;
    }
    bool is_negative()const{
        return data[len-1]>>(bits-1);
    }
    void print_binary(ostream&o, Bigint_Fixedsize const&b){
        o << "[";
        for(size_t i=len;i>0;--i){
            o << bitset<bits>(b.data[i-1]);
        }
        o << "]";
    }
    friend ostream& operator<<(ostream&o, Bigint_Fixedsize const&b){
        if(b.is_negative()){
            return o << '-' << -b << "\n";
        }
        Bigint_base10 ret(0);
        int64_t base = 1u<<bits;
        for(int i = len-1;i>=0;--i){
            ret*=base;
            ret+=Bigint_base10(b.data[i]);
        }
        o << ret;

        return o;
    }
    explicit operator long double()const{
        if(is_negative()){
            return (long double)operator-();
        }
        long double ret = 0.0;
        long double base = 1u<<bits;
        for(int i = len-1;i>=0;--i){
            ret = ret * base + data[i];
        }
        return ret;
    }
    /// TODO: implement for larger inputs
    friend istream& operator>>(istream&i, Bigint_Fixedsize &b){
        int64_t tmp;
        i >> tmp;
        b = Bigint_Fixedsize(tmp);
        return i;
    }
};

/**
 *  Biginteger with fixed precision
 *  Has 32*len bits, is signed
 *
 *  Trying out more optimizations.
 */
template<size_t len>
struct Bigint_Fixedsize_Fast{
    unsigned int data[len];
    uint16_t siz;
    bool sign;
    static constexpr unsigned int bits = 32;
    Bigint_Fixedsize_Fast(){
        data[0] = 0;
        siz = 1;
        sign = false;
    }
    Bigint_Fixedsize_Fast(long long a){
        sign = false;
        if(a<0){
            sign = true;
            a=-a;
        }
        siz = 0;
        do{
            long long b = a>>bits;
            data[siz] = a - (b<<bits);
            a = b;
            ++siz;
        } while(a);
    }
    void trim(){
        while(siz>1 && !data[siz-1]) --siz;
        if(siz == 1 && data[0] == 0) sign=false;
    }
    int comp_unsigned(Bigint_Fixedsize_Fast const&o)const{
        uint16_t lim = min(siz, o.siz);
        for(unsigned int i=lim;i<siz;++i){
            if(data[i]){
                return 1;
            }
        }
        for(unsigned int i=lim;i<o.siz;++i){
            if(o.data[i]){
                return -1;
            }
        }
        for(unsigned int i=lim-1;i+1;--i){
            if(data[i]!=o.data[i]){
                return data[i] < o.data[i]?-1:1;
            }
        }
        return 0;
    }
    int comp(Bigint_Fixedsize_Fast const&o)const{
        int sign_mul = 1-2*sign;
        if(sign != o.sign){
            return sign_mul;
        }
        return sign_mul * comp_unsigned(o);
    }
    bool operator<(Bigint_Fixedsize_Fast const&o)const{
        return comp(o)<0;
    }
    bool operator>(Bigint_Fixedsize_Fast const&o)const{
        return comp(o)>0;
    }
    bool operator<=(Bigint_Fixedsize_Fast const&o)const{
        return comp(o)<=0;
    }
    bool operator>=(Bigint_Fixedsize_Fast const&o)const{
        return comp(o)>=0;
    }
    bool operator==(Bigint_Fixedsize_Fast const&o)const{
        return comp(o)==0;
    }
    bool operator!=(Bigint_Fixedsize_Fast const&o)const{
        return comp(o)!=0;
    }
    bool operator!()const{
        return operator==(ZERO);
    }
    Bigint_Fixedsize_Fast operator-()const{
        Bigint_Fixedsize_Fast ret(*this);
        if(!!ret){
            ret.sign ^=1;
        }
        return ret;
    }
    Bigint_Fixedsize_Fast operator*(Bigint_Fixedsize_Fast const&o)const{
        Bigint_Fixedsize_Fast ret;
        ret.siz = min(siz+o.siz, (int)len);
        ret.sign = (sign!=o.sign);
        fill(ret.data, ret.data+ret.siz, 0);
        for(unsigned int i=0;i<siz;++i){
            unsigned long long carry = 0, carry_2;
            for(unsigned int j=0;j<o.siz;++j){
                carry+= data[i]*(unsigned long long)o.data[j] + ret.data[i+j];
                carry_2 = carry >> bits;
                ret.data[i+j] = carry - (carry_2<<bits);
                carry = carry_2;
            }
            for(unsigned int j=i+o.siz;carry;++j){
                carry+= ret.data[j];
                carry_2 = carry >> bits;
                ret.data[j] = carry - (carry_2<<bits);
                carry = carry_2;
            }
        }
        ret.trim();
        return ret;
    }
    Bigint_Fixedsize_Fast& operator*=(Bigint_Fixedsize_Fast const&o){
        *this = operator*(o);
        return *this;
    }
    static void unsigned_add(Bigint_Fixedsize_Fast &ret, Bigint_Fixedsize_Fast const&A, Bigint_Fixedsize_Fast const&B){
        const Bigint_Fixedsize_Fast *a = &A, *b = &B;
        if(a->siz < b->siz) swap(a, b);
        ret.sign = A.sign;
        unsigned long long carry = 0, carry_2;
        unsigned int j;
        for(j=0;j<b->siz;++j){
            carry+=(unsigned long long)a->data[j] + (unsigned long long)b->data[j];
            carry_2 = carry>>bits;
            ret.data[j] = carry - (carry_2<<bits);
            carry = carry_2;
        }
        for(;j<a->siz;++j){
            carry+=a->data[j];
            carry_2 = carry>>bits;
            ret.data[j] = carry - (carry_2<<bits);
            carry = carry_2;
        }
        if(carry){
            ret.data[j++] = carry;
        }
        ret.siz = j;
        ret.trim();
    }
    static void unsigned_subtract(Bigint_Fixedsize_Fast &ret, Bigint_Fixedsize_Fast const&A, Bigint_Fixedsize_Fast const&B){
        int com = A.comp_unsigned(B);
        if(com == 0){
            ret.sign = false;
            ret.siz = 1;
            ret.data[0] = 0;
            return;
        }
        ret.sign = A.sign;
        const Bigint_Fixedsize_Fast *a = &A, *b = &B;
        if(com < 0){
            ret.sign ^= true;
            swap(a, b);
        }
        // deal with case then o is not trimed.
        unsigned int min_siz = min(A.siz, B.siz);
        unsigned long long carry = 0, carry_2;
        unsigned int j;
        for(j=0;j<min_siz;++j){
            carry+=(unsigned long long)a->data[j] - (unsigned long long)b->data[j];
            carry_2 = carry>>bits;
            ret.data[j] = carry - (carry_2<<bits);
            carry = -(carry_2 & 1u);
        }
        for(;carry;++j){
            assert(j < a->siz);
            carry+=a->data[j];
            carry_2 = carry>>bits;
            ret.data[j] = carry - (carry_2<<bits);
            carry = -(carry_2 & 1u);
        }
        copy(a->data+j, a->data+a->siz, ret.data+j);
        ret.siz = a->siz;
        ret.trim();
    }
    static void add(Bigint_Fixedsize_Fast &ret, Bigint_Fixedsize_Fast const&A, Bigint_Fixedsize_Fast const&B){
        if(A.sign == B.sign){
            unsigned_add(ret, A, B);
        } else {
            unsigned_subtract(ret, A, B);
        }
    }
    static void sub(Bigint_Fixedsize_Fast &ret, Bigint_Fixedsize_Fast const&A, Bigint_Fixedsize_Fast const&B){
        if(A.sign != B.sign){
            unsigned_add(ret, A, B);
        } else {
            unsigned_subtract(ret, A, B);
        }
    }
    Bigint_Fixedsize_Fast operator+(Bigint_Fixedsize_Fast const&o)const{
        Bigint_Fixedsize_Fast ret;
        add(ret, *this, o);
        return ret;
    }
    Bigint_Fixedsize_Fast& operator+=(Bigint_Fixedsize_Fast const&o){
        add(*this, *this, o);
        return *this;
    }
    Bigint_Fixedsize_Fast operator-(Bigint_Fixedsize_Fast const&o)const{
        Bigint_Fixedsize_Fast ret;
        sub(ret, *this, o);
        return ret;
    }
    Bigint_Fixedsize_Fast operator-=(Bigint_Fixedsize_Fast const&o){
        sub(*this, *this, o);
        return *this;
    }
    /// TODO: more operators
    void print_binary(ostream&o, Bigint_Fixedsize_Fast const&b){
        o << "[";
        o << sign << "/" << len << "/";
        for(size_t i=siz;i>0;--i){
            o << bitset<bits>(b.data[i-1]);
        }
        o << "]";
    }
    friend ostream& operator<<(ostream&o, Bigint_Fixedsize_Fast const&b){
        if(b.sign){
            return o << '-' << -b;
        }
        Bigint_base10 ret(0);
        int64_t base = 1ll<<bits;
        for(int i = b.siz-1;i>=0;--i){
            ret*=base;
            ret+=Bigint_base10(b.data[i]);
        }
        o << ret;

        return o;
    }
    explicit operator long double()const{
        if(sign){
            return (long double)operator-();
        }
        long double ret = 0.0;
        long double base = 1ll<<bits;
        for(int i = siz-1;i>=0;--i){
            ret = ret * base + data[i];
        }
        return ret;
    }
    /// TODO: implement for larger inputs
    friend istream& operator>>(istream&i, Bigint_Fixedsize_Fast &b){
        int64_t tmp;
        i >> tmp;
        b = Bigint_Fixedsize_Fast(tmp);
        return i;
    }
    static const Bigint_Fixedsize_Fast ZERO;
};
template<size_t len>
const Bigint_Fixedsize_Fast<len> Bigint_Fixedsize_Fast<len>::ZERO(0);


#endif // BIGNUM_HPP
