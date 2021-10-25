#ifndef PUSHER_H_
#define PUSHER_H_

#include "containers.h"

void push_boris(particle & p, vector3d& e_field, vector3d& b_field, double dt);
void push_vay(particle & p, vector3d& e_field, vector3d& b_field, double dt);

#endif /* PUSHER_H_ */
