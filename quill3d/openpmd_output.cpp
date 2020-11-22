#include "openpmd_output.h"
#include <iostream>

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

void create_attribute(H5::Group & group, const std::string & name, const std::string & value) {
    H5::DataSpace scalar_dataspace{};
    auto type = H5::StrType(0, value.length());

    auto attribute = group.createAttribute(name, type, scalar_dataspace);

    attribute.write(type, (&value[0]));
}

void create_attribute(H5::Group & group, const std::string & name, int value) {
    H5::DataSpace scalar_dataspace{};

    auto attribute = group.createAttribute(name, H5::PredType::NATIVE_UINT32, scalar_dataspace);

    attribute.write(H5::PredType::NATIVE_UINT32, &value);
}

void create_attribute(H5::Group & group, const std::string & name, double value) {
    H5::DataSpace scalar_dataspace{};

    auto attribute = group.createAttribute(name, H5::PredType::NATIVE_DOUBLE, scalar_dataspace);

    attribute.write(H5::PredType::NATIVE_DOUBLE, &value);
}

std::string openpmd::iteration_group_string(int iteration) {
    return std::string("data/") + std::to_string(iteration);
}

void openpmd::initialize_file(const std::string & filepath) {
    auto file = H5::H5File(filepath, H5F_ACC_TRUNC);

    create_attribute(file, "openPMD", "1.1.0");
    create_attribute(file, "openPMDextension", 0);
    create_attribute(file, "basePath", "/data/%T/");
    create_attribute(file, "software", "Quill");
    create_attribute(file, "iterationEncoding", "groupBased");
    create_attribute(file, "iterationFormat", "/data/%T/");
}

void openpmd::initialize_2d_dataset(const std::string & filepath, const std::string & name, hsize_t n1, hsize_t n2) {
    auto file = H5::H5File(filepath, H5F_ACC_RDWR);

    const hsize_t dims[2] {n1, n2};
    H5::DataSpace dataspace(2, dims);

    H5::DataSet dataset = create_dataset(file, name, dataspace);
}

void openpmd::create_iteration_group(const std::string & filepath, int iteration, double dt, double time, double time_units_SI) {
    auto file = H5::H5File(filepath, H5F_ACC_RDWR);
    
    auto iteration_group_name = openpmd::iteration_group_string(iteration);
    auto iteration_group = open_group(file, iteration_group_name);

    create_attribute(iteration_group, "time", time);
    create_attribute(iteration_group, "dt", dt);
    create_attribute(iteration_group, "timeUnitSI", time_units_SI);
}