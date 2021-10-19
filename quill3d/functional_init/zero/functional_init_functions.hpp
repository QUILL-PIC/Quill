#include "../../functional_init_interface.hpp"
#include <toml.hpp>

Field_functions::Field_functions(int i) {
}

double Field_functions::initial_ex(double t, double x, double y, double z) {
    return 10;
}
