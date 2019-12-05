#ifndef PUSHER_H_
#define PUSHER_H_

#include "thinparticle.h"

void push_boris(thinparticle &p, double ex, double ey, double ez, double bx, double by, double bz, double dt);
void push_vay(thinparticle &p, double ex, double ey, double ez, double bx, double by, double bz, double dt);

#endif /* PUSHER_H_ */
