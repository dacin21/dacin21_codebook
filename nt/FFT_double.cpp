namespace FFT{
	template <typename dataType>
	struct Complex
	{
		dataType re, im;
		Complex() : re(0), im(0) {}
		Complex (const dataType& a, const dataType& b) : re(a), im(b) {}
		Complex(const complex<dataType>& c) : re(c.real()), im(c.imag()) {}
		Complex& operator = (const complex<dataType>& c){
			this->re = c.real();
			this->im = c.imag();
			return *this;
		}
		dataType real() const{return this->re;}
		dataType imag() const{return this->im;}
		void real(dataType const& r){this->re = r;}
		void imag(dataType const& i){this->im = i;}
		void polar(const dataType& rho, const dataType& theta = 0){
			this->re = rho * cos(theta);
			this->im = rho * sin(theta);
		}
		Complex conj() const{return Complex(this->re, -this->im);}
		dataType norm() const{return sqrt(this->re * this->re + this->im * this->im);}
		dataType normSquared() const{return this->re * this->re + this->im * this->im;}
		Complex inverse() const{return this->conj() / this->normSquared();}
		Complex operator + (const Complex<dataType>& c) const{return Complex(this->re + c.re, this->im + c.im);}
		Complex& operator += (const Complex<dataType>& c){return *this = *this + c;}
		Complex operator - (const Complex<dataType>& c) const{return Complex(this->re - c.re, this->im - c.im);}
		Complex& operator -= (const Complex<dataType>& c){return *this = *this + c;}
		Complex operator * (const Complex<dataType>& c) const{return Complex(this->re * c.real() - this->im * c.imag(), this->re * c.imag() + this->im * c.real());}
		Complex& operator *= (const Complex<dataType>& c){return *this = *this * c;}
		Complex operator * (const dataType& c) const{return Complex(this->re * c, this->im * c);}
		Complex& operator *= (const dataType& c){return *this = *this * c;}
		Complex operator / (const Complex<dataType>& c) const{return *this * c.inverse();}
		Complex& operator /= (const Complex<dataType>& c){return *this = *this / c;}
		Complex operator / (const dataType& c) const{ return Complex(this->re / c, this->im / c);}
		Complex& operator /= (const dataType& c){ return *this = *this / c;}
		friend istream& operator >> (istream& stream, Complex& C){
			return stream >> C.re >> C.im;
		}
		friend ostream& operator << (ostream& stream, const Complex& C){
			return stream << "(" << C.re << "," << C.im << ")";
		}
	};

	using FFT_t = Complex<double>;
	int log2i(unsigned long long a){
		return __builtin_clzll(1) - __builtin_clzll(a);
	}
	int get_size(int n){
		return 1<<(log2i(n)+1);
	}
	bool ispow2(int a){return (a&-a)==a;}
	vector<FFT_t > roots;
	// pre-calculate complex roots, log(N) calls to sin/cos
	void gen_roots(int N){
		if((int)roots.size()!=N){
			roots.clear();
			roots.resize(N);
			for(int i=0;i<N;++i){
				if((i&-i) == i){
					roots[i] = polar(1.0, 2.0*3.1415926535897932384626*i/N);
				} else {
					roots[i] = roots[i&-i] * roots[i-(i&-i)];
				}
			}
		}
	}
	void fft(FFT_t *a, int n, bool isInv = false){
		for (int i=1, j=0; i<n; ++i) {
			int m = n >> 1;
			for (; j>=m; m >>= 1)
				j -= m;
			j += m;
			if (i < j)
				swap(a[i], a[j]);
		}
		gen_roots(n);
		assert((int)roots.size()==n);
		for(int iter=1, sh=log2i(n)-1;iter<n;iter*=2, --sh){
			for(int x=0;x<n;x+=2*iter){
				for(int y=0;y<iter;++y){
					FFT_t ome = roots[y<<sh];
					if(isInv) ome = ome.conj();
					FFT_t v = a[x+y], w=a[x+y+iter];
					a[x+y] = v+ome*w;
					a[x+y+iter] = v-ome*w;
				}
			}
		}
	}
	void fft(vector<FFT_t> &v, bool isInv = false){
		assert(ispow2(v.size()));
		fft(v.data(), v.size(), isInv);
	}
	vector<FFT_t> conv(vector<long long> const&v, int n){
		assert(ispow2(n));
		vector<FFT_t> a(n);
		for(size_t i=0;i<v.size();++i){
			a[i].real(v[i]);
		}
		fft(a, false);
		return a;
	}
	vector<long long> convmul(vector<FFT_t> a, vector<FFT_t> const&b){
		assert(a.size() == b.size());
		for(size_t i=0;i<a.size();++i){
			a[i]*=b[i];
		}
		fft(a, true);
		vector<long long> ret(a.size());
		for(size_t i=0;i<a.size();++i){
			ret[i] = llround(a[i].real()/(int)a.size());
		}
		return ret;
	}
	vector<long long> poly_mul(vector<long long> const&a, vector<long long> const&b){
		int n = get_size(a.size()+b.size()-1);
		vector<FFT_t> x = conv(a, n);
		vector<FFT_t> y = conv(b, n);
		vector<long long> ret = convmul(x, y);
		ret.resize(a.size()+b.size()-1);
		return ret;
	}
	vector<long long> blockify(vector<long long> const&a){
		vector<long long> ret(3*a.size());
		for(size_t i=0;i<a.size();++i){
			ret[3*i] = a[i]%30000;
			ret[3*i+1] = a[i]/30000;
		}
		return ret;
	}
	vector<long long> unblockify(vector<long long> const&a, long long const&mod){
		assert(a.size()%3 == 0);
		vector<long long> ret(a.size()/3);
		for(size_t i=0;i<ret.size();++i){
			ret[i] = (((a[3*i+2])%mod*30000 + a[3*i+1])%mod*30000 + a[3*i])%mod;
		}
		return ret;
	}
	vector<long long> poly_mul_block(vector<long long> const&a, vector<long long> const&b, long long const&mod){
		vector<long long> c = blockify(a), d = blockify(b);
		int n = get_size(c.size()+d.size()-1);
		vector<FFT_t> x = conv(c, n);
		vector<FFT_t> y = conv(d, n);
		vector<long long> ret = convmul(x, y);
		ret.resize(3*(a.size()+b.size()-1));
		vector<long long> r = unblockify(ret, mod);
		return r;
	}
}