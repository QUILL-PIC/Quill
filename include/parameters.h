#pragma once

#include <vector>
#include <string>
#include <cassert>
#include "physical_constants.h"

class OutputIterations {
public:
    OutputIterations() = default;
    OutputIterations(const std::vector<double> & t_end, const std::vector<double> & output_period, double timestep);
    
    void advance() {
        ++current_iteration_index;
    }

    int get_current_iteration() const {
        return iterations[current_iteration_index];
    }

    double get_current_iteration_time() const {
        return iteration_times[current_iteration_index];
    }

    std::string get_current_iteration_string() const;

    const std::vector<double> & get_iteration_times() const {
        return iteration_times;
    }

    const std::vector<int> & get_iterations() const {
        return iterations;
    }

    int get_last_iteration() const {
        return iterations[size-1];
    }

    double get_last_iteration_time() const {
        return iteration_times[size-1];
    }
private:
    int size;
    std::vector<int> iterations;
    std::vector<double> iteration_times;
    double timestep;
    int current_iteration_index = 0;
};
