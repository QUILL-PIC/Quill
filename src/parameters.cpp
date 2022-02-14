#include "parameters.h"

#include <cmath>

// constructs a series of iterations from a list of itervals ending at t_end each with specified output period at each period
OutputIterations::OutputIterations(const std::vector<double> & t_end, const std::vector<double> & output_period, double timestep) : timestep(timestep) {
    assert(t_end.size() == output_period.size());
    assert(timestep > 0);

    int previous_iteration = -1;

    const int size = t_end.size();
    double t_begin = 0;
    for (int i = 0; i < size; i++) {
        int iteration_count = static_cast<int>((t_end[i] + 0.5 * timestep - t_begin) / output_period[i]) + 1;

        for (int j = 0; j < iteration_count; j++) {
            double time = t_begin + j * output_period[i];
            int iteration = static_cast<int>(round(time / timestep));
            if (iteration > previous_iteration) {
                iteration_times.push_back(time);
                iterations.push_back(iteration);
                previous_iteration = iteration;
            }
        }

        t_begin = t_end[i];
    }

    // the final iteration is used for output
    const double t_final = t_end[size-1];
    const int last_iteration = static_cast<int>(round(t_final / timestep));
    if (last_iteration > previous_iteration) {
        iteration_times.push_back(t_final);
        iterations.push_back(last_iteration);
    }

    this->size = iterations.size();
}

std::string OutputIterations::get_current_iteration_string() const {
    const double ACCURACY = 100.;

    double value = static_cast<int>(round(get_current_iteration_time() / 2 / PI * ACCURACY)) / ACCURACY;

    char str_buffer[100];
    sprintf(str_buffer, "%g", value);
    return std::string(str_buffer);
}