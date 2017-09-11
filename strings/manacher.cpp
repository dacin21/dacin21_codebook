// UNTESTED
// manacher palyndrome algorithm
vector<int> manacher (string const&s){
	int n = s.size(), m=2*n-1;
	vector<int> ret(m), t(m, -1);
	for(int i=0;i<n;++i) t[2*i] = s[i];
	int x=0;
	for(int i=1;i<m;++i){
		int &y = ret[i] = 0; 
		if(i<x+ret[x]) y = min(ret[2*x-i], x+ret[x]-i); // copy from left side
		while(i-y-1 >=0 && i+y+1 < m && t[i-y-1] == t[i+y+1]) ++y; // scan
		if(i+y>=x+ret[x]) x=i; // move right border
	}
	for(int i=0;i<m;++i) if(i-ret[i]==0 || i+ret[i]== m-1) ++ret[i]; // ???
	for(int&e:ret) e/=2;
	return ret;
}