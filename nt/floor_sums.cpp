// sum K=1...N floor(kp/q)
long long floor_sum(long long N, long long p, long long q){
    long long t = __gcd(p, q);
    p/=t, q/=t;
    long long s=0, z=1;
    while(q>0 && N>0){
        t = p/q;
        s+= N*(N+1)/2*z*t;
        p-= q*t;
        t = N/q;
        s+= z*p*t*(N+1) - z*t*(p*q*t+p+q-1)/2;
        N-= q*t;
        t = N*p/q;
        s+= z*t*N;
        N = t;
        swap(p, q);
        z*= -1;
    }
    return s;
}
// number of integer points the triangle
// ax + by <=c && x, y > 0; where a, b, c > 0
int64_t count_triangle(int64_t a, int64_t b, int64_t c){
    if(b>a) swap(a, b);
    int64_t m = c/a;
    if(a==b) return m*(m-1)/2;
    int64_t k= (a-1)/b, h = (c-a*m)/b;
    return m*(m-1)/2*k + m*h + count_triangle(b, a-b*k, c-b*(k*m+h));
}