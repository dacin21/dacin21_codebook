// different schemes for numerical integration
// approximatively ordered by accuracy
// do NOT use integer types for integration range!

struct Integration_Midpoint{
    template<typename Func, typename S>
    static typename result_of<Func(S)>::type integrate_step(Func f, S l, S r){
        S m = (l+r)/2;
        return f(m) * (r-l);
    }
};
struct Integration_Simpson{
    template<typename Func, typename S>
    static typename result_of<Func(S)>::type integrate_step(Func f, S l, S r){
        S m = (l+r)/2;
        return (f(l) + 4*f(m) + f(r))/6 * (r-l);
    }
};
struct Integration_Gauss_2{
    static constexpr long double A = 1.0l/sqrtl(3)/2, x1=0.5l-A, x2 = 0.5l+A;
    template<typename Func, typename S>
    static typename result_of<Func(S)>::type integrate_step(Func f, S l, S r){
        return (f(l*x1 + r*x2) + f(l*x2+r*x1))/2 * (r-l);
    }
};
struct Integration_NCotes_Open_4{
    template<typename Func, typename S>
    static typename result_of<Func(S)>::type integrate_step(Func f, S l, S r){
        S h = (r-l)/5;
        return (11*f(l+h) + f(l+2*h) + f(r-2*h) + 11*f(r-h))/24 * (r-l);
    }
};
struct Integration_Gauss_3{
    static constexpr long double A = sqrtl(3.0l/5.0l)/2, x1=0.5-A, x2 = 0.5+A;
    template<typename Func, typename S>
    static typename result_of<Func(S)>::type integrate_step(Func f, S l, S r){
        return (5*f(l*x1 + r*x2) + 8*f((l+r)/2) + 5*f(l*x2+r*x1))/18 * (r-l);
    }
};

template<typename Integration_Method>
struct Integrator_Fixedstep{
    template<typename Func, typename S>
    static typename result_of<Func(S)>::type integrate(Func f, S const l, S const r, size_t const steps){
        assert(steps>0);
        typename result_of<Func(S)>::type ret(0);
        S cur_l = l, cur_r;
        for(size_t i=0;i<steps;++i){
            cur_r = (l*(steps-i-1) + r*(i+1))/steps;
            ret+=Integration_Method::integrate_step(f, cur_l, cur_r);
            cur_l = cur_r;
        }
        return ret;
    }
};
template<typename Integration_Method>
class Integrator_Adaptive{
private:
    template<size_t depth_limit, typename Func, typename S>
    static typename result_of<Func(S)>::type integrate(Func f, S const l, S const r, typename result_of<Func(S)>::type const val, typename result_of<Func(S)>::type const eps, const size_t depth){
        if(depth>=depth_limit){
            return val;
        }
        S const m = (l+r)/2;
        typename result_of<Func(S)>::type val_l = Integration_Method::integrate_step(f, l, m);
        typename result_of<Func(S)>::type val_r = Integration_Method::integrate_step(f, m, r);
        typename result_of<Func(S)>::type error = abs(val - val_l - val_r);
        if(error < eps){
            return val_l + val_r;
        }
        return integrate<depth_limit>(f, l, m, val_l, eps/2, depth+1)
            + integrate<depth_limit>(f, m, r, val_r, eps/2, depth+1);
    }
public:
    template<size_t depth_limit, typename Func, typename S>
    static typename result_of<Func(S)>::type integrate(Func f, S const l, S const r, typename result_of<Func(S)>::type const eps){
        return integrate<depth_limit>(f, l, r, Integration_Method::integrate_step(f, l, r), eps, 0);
    }
};
