/*
 *  Short suffix array in O(N log N)
 *  uses MSD-radix sort
 *  0.624 sec for N<500001
 */
class Suffix_Array{
    void make_lcp(int const*s){
        for(int i=0, k=0;i<n;++i)
            if(inv[i]!=0){
                for(int j = sa[inv[i]-1];k<n-max(i, j) && s[i+k]==s[j+k];++k);
                lcp[inv[i]-1]=k;
                if(k)--k;
            } else k=0;
    }
    bool cmp(vector<int> const&rank,int i,int j,int l){
        return i<n-l && j<n-l && rank[i]==rank[j] && rank[i+l]==rank[j+l];
    }
    void build_SA(int const*s){
        int m = max(256, n);
        sa.resize(n);
        vector<int> x(n), y(n), buck(m);
        int i,j,p;
        for(i=0;i<n;++i) buck[x[i]=s[i]]++;
        for(i=1;i<m;++i) buck[i]+=buck[i-1];
        for(i=n-1;i>=0;--i) sa[--buck[x[i]]]=i;
        for(p=0,j=1;p+1<n;m=p+1,j*=2) {
            for(p=0,i=n-j;i<n;++i) y[p++]=i;
            for(i=0;i<n;++i) if (sa[i]>=j) y[p++]=sa[i]-j;
            memset(buck.data(),0,sizeof(int)*(m));
            for(i=0;i<n;++i) ++buck[x[y[i]]];
            for(i=1;i<m;++i) buck[i]+=buck[i-1];
            for(i=n-1;i>=0;--i) sa[--buck[x[y[i]]]]=y[i];
            x.swap(y);
            for(p=0,x[sa[0]]=0,i=1;i<n;i++) x[sa[i]]=cmp(y,sa[i-1],sa[i],j) ? p : ++p;
        }
        x.swap(inv); y.swap(lcp);
        for(int i=0;i<n;++i) inv[sa[i]]=i;
        make_lcp(s);
    }
    public:
    Suffix_Array(int const*s, int _N):n(_N){build_SA(s);}
    Suffix_Array(vector<int> const&s):Suffix_Array(s.data(), s.size()){}
    Suffix_Array(char const*s, int _N):n(_N){build_SA(vector<int> (s, s+n).data());}
    Suffix_Array(string const&s):Suffix_Array(s.data(), s.size()){}
    vector<int> sa, inv, lcp;
    const int n;
};