#include <toml.hpp>

class Field_functions {
    private:
        toml::value table;
    public:
        Field_functions(toml::value);
        Field_functions(int);
        double initial_ex(double, double, double, double);
        double initial_ey(double, double, double, double);
        double initial_ez(double, double, double, double);
        double initial_bx(double, double, double, double);
        double initial_by(double, double, double, double);
        double initial_bz(double, double, double, double);
};
