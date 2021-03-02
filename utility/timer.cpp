template<bool print = false>
struct Timer{
    using hrc = chrono::high_resolution_clock;
    Timer(string name_) : name(move(name_)), start(hrc::now()) {}
    ~Timer(){
        if(print){
            end = hrc::now();
            double time = chrono::duration_cast<chrono::nanoseconds>(end-start).count();
            cerr << name << " : " << time*1e-9 << "\n";
        }
    }
    double get(){
        end = hrc::now();
        double time = chrono::duration_cast<chrono::nanoseconds>(end-start).count();
        return time;
        cerr << name << " : " << time*1e-9 << "\n";
    }

    string name;
    decltype(hrc::now()) start, end;
};
