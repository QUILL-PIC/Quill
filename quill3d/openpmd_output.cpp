#include "openpmd_output.h"
#include <iostream>
#include <vector>
#include <cstring>
#include <exception>

const std::string MESHES_GROUP = "meshes";

std::string get_parent(const std::string & path) {
    auto last_separator = path.find_last_of('/');
     if (last_separator == std::string::npos) {
         return "/";
     } else {
         return path.substr(0, last_separator);
     }
}

std::string get_name(const std::string & path) {
    auto last_separator = path.find_last_of('/');
     if (last_separator == std::string::npos) {
         return path;
     } else {
         return path.substr(last_separator + 1);
     }
}

H5::Group open_group(const H5::Group & base_group, const std::string & group_name) {
    auto first_separator = group_name.find_first_of('/');

    if (first_separator == std::string::npos) {
        if (base_group.nameExists(group_name)) {
            return base_group.openGroup(group_name);
        } else {
            return base_group.createGroup(group_name);
        }
    } else {
        auto first_group_name = group_name.substr(0, first_separator);
        auto subgroup_name = group_name.substr(first_separator+1);
        auto first_group = open_group(base_group, first_group_name);
        return open_group(first_group, subgroup_name);
    }
}

H5::DataSet create_dataset(const H5::Group & base_group, const std::string & name, const H5::DataSpace & dataspace) {
    auto last_separator = name.find_last_of('/');

    if (last_separator == std::string::npos) {
        return base_group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    } else {
        auto new_base_group = open_group(base_group, name.substr(0, last_separator));
        auto new_name = name.substr(last_separator+1);
        return create_dataset(new_base_group, new_name, dataspace);
    }
}

H5::DataSet openpmd::open_dataset(const H5::Group & base_group, const std::string & name) {
    auto last_separator = name.find_last_of('/');

    if (last_separator == std::string::npos) {
        return base_group.openDataSet(name);
    } else {
        auto new_base_group = open_group(base_group, name.substr(0, last_separator));
        auto new_name = name.substr(last_separator+1);
        return open_dataset(new_base_group, new_name);
    }
}

void create_attribute(H5::H5Object & object, const std::string & name, const std::string & value) {
    H5::DataSpace scalar_dataspace{};
    auto type = H5::StrType(0, value.length());

    auto attribute = object.createAttribute(name, type, scalar_dataspace);

    attribute.write(type, (&value[0]));
}

void create_attribute(H5::H5Object & object, const std::string & name, int value) {
    H5::DataSpace scalar_dataspace{};

    auto attribute = object.createAttribute(name, H5::PredType::NATIVE_UINT32, scalar_dataspace);

    attribute.write(H5::PredType::NATIVE_UINT32, &value);
}

void create_attribute(H5::H5Object & object, const std::string & name, double value) {
    H5::DataSpace scalar_dataspace{};

    auto attribute = object.createAttribute(name, H5::PredType::NATIVE_DOUBLE, scalar_dataspace);

    attribute.write(H5::PredType::NATIVE_DOUBLE, &value);
}

void create_attribute(H5::H5Object & object, const std::string & name, const std::vector<double> & values) {
    auto size = values.size();

    const hsize_t dims[1] { size };
    H5::DataSpace dataspace(1, dims);

    auto attribute = object.createAttribute(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    attribute.write(H5::PredType::NATIVE_DOUBLE, &(values[0]));
}

// An alternative (nicer) implementation using variable-length string arrays.
// Currently not supported by openpmd-validator (see https://github.com/openPMD/openPMD-validator/issues/61)
/*
void create_attribute(H5::H5Object & object, const std::string & name, const std::vector<std::string> & values) {
    auto size = values.size();

    std::vector<const char *> values_cstr(size);
    for (size_t i = 0; i < size; i++) {
        values_cstr[i] = values[i].c_str();
    }

    const hsize_t dims[1] { size };
    H5::DataSpace dataspace(1, dims);

    auto datatype = H5::StrType(H5::PredType::C_S1, H5T_VARIABLE);
    auto attribute = object.createAttribute(name, datatype, dataspace);

    attribute.write(datatype, values_cstr.data());
}
*/

void create_attribute(H5::H5Object & object, const std::string & name, const std::vector<std::string> & values) {
    auto size = values.size();

    size_t max_element_size = 0;
    for (size_t i = 0; i < size; i++) {
        if (values[i].length() > max_element_size) {
            max_element_size = values[i].length() + 1;
        }
    }

    std::vector<char> write_array(size * max_element_size);

    for (size_t i = 0; i < size; i++) {
        strcpy(&(write_array[i * max_element_size]), values[i].c_str());
    }

    const hsize_t dims[1] { size };
    H5::DataSpace dataspace(1, dims);

    auto datatype = H5::StrType(H5::PredType::C_S1, max_element_size);
    auto attribute = object.createAttribute(name, datatype, dataspace);

    for (size_t i = 0; i < size; i++) {
        attribute.write(datatype, write_array.data());
    }
}

std::string iteration_group_string(int iteration) {
    return std::string("data/") + std::to_string(iteration);
}

std::string openpmd::iteration_meshes_group_string(int iteration) {
    return iteration_group_string(iteration) + "/" + MESHES_GROUP;
}

void openpmd::initialize_file(const std::string & filepath) {
    auto file = H5::H5File(filepath, H5F_ACC_TRUNC);

    create_attribute(file, "openPMD", "1.1.0");
    create_attribute(file, "openPMDextension", 0);
    create_attribute(file, "basePath", "/data/%T/");
    create_attribute(file, "software", "Quill");
    create_attribute(file, "iterationEncoding", "groupBased");
    create_attribute(file, "iterationFormat", "/data/%T/");
    create_attribute(file, "meshesPath", MESHES_GROUP + "/");
}

std::vector<std::string> get_axes_labels(openpmd::Space space) {
    if (space == openpmd::Space::XY) {
        return {"x", "y"};
    } else if (space == openpmd::Space::YZ) { 
        return {"y", "z"};
    } else if (space == openpmd::Space::XZ) {
        return {"x", "z"};
    } else if (space == openpmd::Space::XYZ) {
        return {"x", "y", "z"};
    } else {
        throw std::invalid_argument("Space has incorrect value");
    }
}

int get_space_size(openpmd::Space space) {
    if (space == openpmd::Space::XY || space == openpmd::Space::XZ || space == openpmd::Space::YZ) {
        return 2;
    } else if (space == openpmd::Space::XYZ) {
        return 3;
    } else {
        throw std::invalid_argument("Space has incorrect value");
    }
}

void initialize_mesh_attrs(H5::H5Object & obj, openpmd::Space space, const std::vector<double> & grid_spacing) {
    create_attribute(obj, "geometry", "cartesian");
    create_attribute(obj, "dataOrder", "C");
    create_attribute(obj, "gridUnitSI", 1.0);
    create_attribute(obj, "gridSpacing", grid_spacing);
    create_attribute(obj, "gridGlobalOffset", std::vector<double>(get_space_size(space), 0.0));
    create_attribute(obj, "unitDimension", std::vector<double>(7, 0.0));
    create_attribute(obj, "timeOffset", 0.0);
    create_attribute(obj, "axisLabels", get_axes_labels(space));
}

void openpmd::initialize_2d_dataset(const std::string & filepath, const std::string & name, hsize_t n1, hsize_t n2, Space space, std::array<double, 2> grid_spacing) {
    auto file = H5::H5File(filepath, H5F_ACC_RDWR);

    const hsize_t dims[2] {n1, n2};
    H5::DataSpace dataspace(2, dims);

    auto group_name = get_parent(name);
    auto group = open_group(file, group_name);

    auto dataset_name = get_name(name);
    H5::DataSet dataset = create_dataset(group, dataset_name, dataspace);

    std::vector<double> grid_spacing_v{grid_spacing.begin(), grid_spacing.end()};

    if (dataset_name == "x" || dataset_name == "y" || dataset_name == "z") {
        if (!group.attrExists("geometry")) {
            initialize_mesh_attrs(group, space, grid_spacing_v);
        }
    } else {
        initialize_mesh_attrs(dataset, space, grid_spacing_v);
    }

    // TODO proper SI units
    create_attribute(dataset, "unitSI", 1.0);
    // TODO proper position of the grid
    create_attribute(dataset, "position", std::vector<double>{0.0, 0.0});
}

void openpmd::initialize_3d_dataset(const std::string & filepath, const std::string & name, std::array<hsize_t, 3> size, std::array<double, 3> grid_spacing) {
    auto file = H5::H5File(filepath, H5F_ACC_RDWR);

    H5::DataSpace dataspace(3, size.data());

    auto group_name = get_parent(name);
    auto group = open_group(file, group_name);

    auto dataset_name = get_name(name);
    H5::DataSet dataset = create_dataset(group, dataset_name, dataspace);

    std::vector<double> grid_spacing_v{grid_spacing.begin(), grid_spacing.end()};

    if (dataset_name == "x" || dataset_name == "y" || dataset_name == "z") {
        if (!group.attrExists("geometry")) {
            initialize_mesh_attrs(group, openpmd::Space::XYZ, grid_spacing_v);
        }
    } else {
        initialize_mesh_attrs(dataset, openpmd::Space::XYZ, grid_spacing_v);
    }

    // TODO proper SI units
    create_attribute(dataset, "unitSI", 1.0);
    // TODO proper position of the grid
    create_attribute(dataset, "position", std::vector<double>{0.0, 0.0});    
}

void openpmd::create_iteration_group(const std::string & filepath, int iteration, double dt, double time, double time_units_SI) {
    auto file = H5::H5File(filepath, H5F_ACC_RDWR);
    
    auto iteration_group_name = iteration_group_string(iteration);
    auto iteration_group = open_group(file, iteration_group_name);

    create_attribute(iteration_group, "time", time);
    create_attribute(iteration_group, "dt", dt);
    create_attribute(iteration_group, "timeUnitSI", time_units_SI);
}