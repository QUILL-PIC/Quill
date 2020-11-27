#include "balancing.h"
#include "main.h"
#include <iostream>

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