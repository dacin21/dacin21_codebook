// could maybe be optimized by transposing
template<typename T = int, typename comp = less<T>>
class Sparsetable{
    int compind(int i, int j, int l, int sign)const{
        i = RMQ[lgn*i+l], j=RMQ[lgn*(j+sign*(1<<l))+l];
        return comp()(data[i], data[j])? i : j;
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
    Sparsetable(T const*v, int _N):n(_N), data(v, v+_N){build();}
    Sparsetable(vector<T> const&v):Sparsetable(v.data(), v.size()){};
    // min in [l, r)
    pair<int, T&> operator()(int const&i, int const&j){
        int k = minind(i, j);
        return pair<int, T&>(k, data[k]);
    }
    int n, lgn;
    vector<int> RMQ, lg2;
    vector<T> data;
};
