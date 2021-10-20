#include <toml.hpp>

class InitFunctions {
    private:
        const toml::value table;
    public:
        InitFunctions(toml::value);
        InitFunctions(int);
        double initial_ex(double, double, double);
        double initial_ey(double, double, double);
        double initial_ez(double, double, double);
        double initial_bx(double, double, double);
        double initial_by(double, double, double);
        double initial_bz(double, double, double);
};
