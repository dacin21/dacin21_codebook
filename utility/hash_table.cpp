/*
 *  Hash table with linear probing
 *  Has compile-time fixed size
 */

template<typename Key, typename Value, size_t log_size, typename Hasher = std::hash<Key>>
class Hash_Table{
    static Hasher hasher;
    static constexpr size_t siz = 1ull<<log_size;
    static constexpr size_t mask = siz-1;
    static constexpr double max_load_factor = 0.9;
    static constexpr size_t max_stored_cnt = siz*max_load_factor;
    const Key dummy;
    vector<Key> keys;
    vector<Value> data;
    size_t stored_cnt;
    size_t hashify(Key const&key){
        size_t val = hasher(key);
        return val * size_t{484391091};
    }
    size_t find(Key const&key){
        size_t pos = hashify(key) & mask;
        while(keys[pos] != key && keys[pos] != dummy){
            pos = (pos + 1) & mask;
        }
        return pos;
    }
public:
    Hash_Table(Key const&dummy_):dummy(dummy_), keys(siz, dummy), data(siz), stored_cnt(0){}
    bool count(Key const&key){
        size_t pos = find(key);
        return keys[pos] == key;
    }
    Value& operator[](Key const&key){
        size_t pos = find(key);
        if(keys[pos] != key){
            keys[pos] = key;
            assert(++stored_cnt <= max_stored_cnt);
        }
        return data[pos];
    }
};
template<typename Key, typename Value, size_t log_size, typename Hasher> Hasher Hash_Table<Key, Value, log_size, Hasher>::hasher;
template<typename Key, typename Value, size_t log_size, typename Hasher> constexpr size_t Hash_Table<Key, Value, log_size, Hasher>::siz;
template<typename Key, typename Value, size_t log_size, typename Hasher> constexpr size_t Hash_Table<Key, Value, log_size, Hasher>::mask;
template<typename Key, typename Value, size_t log_size, typename Hasher> constexpr double Hash_Table<Key, Value, log_size, Hasher>::max_load_factor;
template<typename Key, typename Value, size_t log_size, typename Hasher> constexpr size_t Hash_Table<Key, Value, log_size, Hasher>::max_stored_cnt;
