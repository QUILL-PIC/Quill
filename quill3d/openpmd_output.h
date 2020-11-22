#include <H5Cpp.h>
#include <string>

namespace openpmd {
    void initialize_file(const std::string & filepath);
    void initialize_2d_dataset(const std::string & filepath, const std::string & group, hsize_t n1, hsize_t n2);
    std::string iteration_group_string(int iteration);
    void create_iteration_group(const std::string & filepath, int iteration, double dt, double time, double time_units_SI);
    H5::DataSet open_dataset(const H5::Group & base_group, const std::string & name);
}