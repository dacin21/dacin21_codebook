namespace FIO{
    char buf[32*1042|1];
    int bc=0, be=0;
    char gc(){
        if(bc>=be){
            be = fread(buf, 1, sizeof(buf)-1, stdin);
            buf[be] = bc = 0;
        }
        return buf[bc++];
    }
    void readint(){}
    void readuint(){}
    template<typename T, typename ...S>
    void readuint(T &a, S& ...b){
        a=0;
        int x=gc();
        while(x<'0' || x>'9') x=gc();
        while(x>='0' && x<='9'){
            a = a*10+x-'0'; x=gc();
        }
        readuint(b...);
    }
	template<typename T, typename ...S>
    void readint(T &a, S& ...b){
        a=0;
        int x=gc(), s=1;;
        while(x!='-' && (x<'0' || x>'9')) x=gc();
		if(x=='-'){ s=-s; x=gc(); }
        while(x>='0' && x<='9'){
            a = a*10+x-'0'; x=gc();
        }
		if(s<0) a=-a;
        readint(b...);
    }
}
using FIO::readint;
using FIO::readuint;