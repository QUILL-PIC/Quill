#include "balancing.h"
#include "main.h"
#include <numeric>
#include <algorithm>
#include <cassert>

template <class T>
class array2d {
public:
    array2d(int n1, int n2) : n1(n1), n2(n2) {
        array = std::vector<T>(n1 * n2, 0.0);
    }

    T & operator()(int i, int j) {
        return array[n2 * i + j];
    }
    const T & operator()(int i, int j) const {
        return array[n2 * i + j];
    }
private:
    int n1, n2;
    std::vector<T> array;
};

std::vector<double> spatial_region::calculate_layer_weights(double particle_weight) {
    std::vector<double> weights(nx, ny * nz);

    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            for (int k=0; k<nz; k++) {
                auto current = cp[i][j][k].pl.head;
                while (current!=0) {
                    weights[i] += particle_weight;
                    current = current->next;
                }
            }
        }
    }

    return weights;
}

std::vector<int> retrieve_partition(const array2d<int> & last_border, int k, int i, int overlap) {
    int current_border = last_border(k, i);
    if (k == 0) {
        return std::vector<int>{ current_border };
    } else {
        auto partition = retrieve_partition(last_border, k - 1, current_border, overlap);
        partition.push_back(current_border - overlap + 1);
        return partition;
    }
}

std::vector<int> calculate_optimal_partition(const std::vector<double> & weights, int regions_number, int overlap) {
    int size = weights.size();

    array2d<double> max_weight(regions_number, size);
    array2d<int> last_border(regions_number, size);
    
    int minimum_side_block_size = overlap;
    int minimum_middle_block_size = 2 * overlap;

    max_weight(0, minimum_side_block_size-1) = std::accumulate(weights.begin(), weights.begin() + minimum_side_block_size, 0.0);

    for (int i = minimum_side_block_size; i < size; i++) {
        max_weight(0, i) = max_weight(0, i-1) + weights[i];
        for (int k = 1; k < regions_number - 1; k++) {
            int minimum_i = (k + 1) * overlap - 1;
            if (i < minimum_i) {
                break;
            }
            int initial_left = i - minimum_middle_block_size + 1;
            double current_weight = std::accumulate(weights.begin() + initial_left, weights.begin() + i + 1, 0.0);
            int minimum_left = (k - 1) * overlap;
            max_weight(k, i) = std::max(current_weight, max_weight(k-1, initial_left + overlap - 1));
            last_border(k, i) = initial_left + overlap - 1;
            for (int current_left = initial_left-1; current_left > minimum_left - 1; current_left--) {
                current_weight += weights[current_left];
                double current_max = std::max(current_weight, max_weight(k-1, current_left + overlap - 1));
                if (current_max < max_weight(k, i)) {
                    max_weight(k, i) = current_max;
                    last_border(k, i) = current_left + overlap - 1;
                }
                if (current_weight > max_weight(k, i)) {
                    break;
                }
            }
        }
    }

    int k = regions_number - 1;
    int i = size - 1;
    int initial_left = size - minimum_side_block_size;
    double current_weight = std::accumulate(weights.begin() + initial_left, weights.begin() + size, 0.0);
    int minimum_left = (regions_number - 2) * overlap;
    max_weight(k, i) = std::max(current_weight, max_weight(k-1, initial_left + overlap - 1));
    last_border(k, i) = initial_left + overlap - 1;
    for (int current_left = initial_left-1; current_left > minimum_left - 1; current_left--) {
        current_weight += weights[current_left];
        double current_max = std::max(current_weight, max_weight(k-1, current_left + overlap - 1));
        if (current_max < max_weight(k, i)) {
            max_weight(k, i) = current_max;
            last_border(k, i) = current_left + overlap - 1;
        }
        if (current_weight > max_weight(k, i)) {
            break;
        }
    }

    return retrieve_partition(last_border, k, i, overlap);
}

std::vector<double> calculate_partition_weights(const std::vector<double> & weights, const std::vector<int> & partition, int overlap) {
    int size = partition.size();
    
    std::vector<double> partition_weights(size);
    for (int i = 0; i < size - 1; i++) {
        partition_weights[i] = std::accumulate(weights.begin() + partition[i], weights.begin() + partition[i+1] + overlap, 0.0);
    }
    partition_weights[size-1] = std::accumulate(weights.begin() + partition[size - 1], weights.begin() + weights.size(), 0.0);
    return partition_weights;
}

double calculate_partition_imbalance(const std::vector<double> & partition_weights) {
    return *std::max_element(partition_weights.begin(), partition_weights.end()) / *std::min_element(partition_weights.begin(), partition_weights.end()) - 1.0;
}

double calculate_partition_imbalance(const std::vector<double> & weights, const std::vector<int> & partition, int overlap) {
    return calculate_partition_imbalance(calculate_partition_weights(weights, partition, overlap));
}

std::vector<int> normalize_new_partition(const std::vector<int> & old_partition, const std::vector<int> & new_partition, int overlap) {
    
    // partitions cannot change so much that the new n-th region does not overlap with the old n-th region.

    assert(old_partition[0] == 0);
    assert(new_partition[0] == 0);
    assert(old_partition.size() == new_partition.size());

    size_t size = old_partition.size();

    std::vector<int> normalized_partition(size);
    normalized_partition[0] = new_partition[0];

    for (int i = 1; i < size - 1; i++) {
        normalized_partition[i] = new_partition[i];
        if (new_partition[i] < old_partition[i-1]) {
            normalized_partition[i] = old_partition[i-1];
        }
        if (new_partition[i] > old_partition[i+1]) {
            normalized_partition[i] = old_partition[i+1];
        }
        if (normalized_partition[i] < normalized_partition[i-1] + overlap) {
            normalized_partition[i] = normalized_partition[i-1] + overlap;
        }
    }
    int i = size - 1;
    normalized_partition[i] = new_partition[i];
    if (new_partition[i] < old_partition[i-1]) {
        normalized_partition[i] = old_partition[i-1];
    }
    if (normalized_partition[i] < normalized_partition[i-1] + overlap) {
        normalized_partition[i] = normalized_partition[i-1] + overlap;
    }

    return normalized_partition;
}