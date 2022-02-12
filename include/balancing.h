#pragma once

#include <vector>

std::vector<int> calculate_optimal_partition(const std::vector<double> & weights, int regions_number, int overlap);
std::vector<double> calculate_partition_weights(const std::vector<double> & weights, const std::vector<int> & partition, int overlap);
double calculate_partition_imbalance(const std::vector<double> & partition_weights);
double calculate_partition_imbalance(const std::vector<double> & weights, const std::vector<int> & partition, int overlap);
std::vector<int> normalize_new_partition(const std::vector<int> & old_partition, const std::vector<int> & new_partition, int overlap);
