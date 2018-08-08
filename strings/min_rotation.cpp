// lexicographically smallest rotation
// runs in O(n)
int min_rotation(string s){
    int n = s.size();
    s+=s;
    int ans = 0;
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j){
            if(ans+j == i || s[ans+j] < s[i+j]){
                i+=max(0, j-1); break;
            }
            if(s[ans+j] > s[i+j]){
                ans = i; break;
            }
        }
    }
    return ans;
}
