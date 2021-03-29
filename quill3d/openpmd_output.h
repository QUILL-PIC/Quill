#ifndef OPENPMD_OUTPUT_H_
#define OPENPMD_OUTPUT_H_

#include <H5Cpp.h>
#include <string>
#include <array>

namespace openpmd {
    enum class Space {
        XY, XZ, YZ, XYZ
    };

    void initialize_file(const std::string & filepath);
    void initialize_2d_dataset(const std::string & filepath, const std::string & name, hsize_t n1, hsize_t n2, Space space, std::array<double, 2> grid_spacing);
    void initialize_3d_dataset(const std::string & filepath, const std::string & name, std::array<hsize_t, 3> size, std::array<double, 3> grid_spacing);
    std::string iteration_meshes_group_string(int iteration);
    void create_iteration_group(const std::string & filepath, int iteration, double dt, double time, double time_units_SI);
    H5::DataSet open_dataset(const H5::Group & base_group, const std::string & name);
}

#endif