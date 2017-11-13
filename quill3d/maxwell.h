#ifndef MAXWELL_H_
#define MAXWELL_H_
#include "containers.h"

class maxwell_solver {
public:
    virtual void fadvance(field3d<celle> & ce, field3d<cellb> & cb, field3d<cellj> & cj, double dt, double dx,
            double dy, double dz) = 0;
    virtual void fzeroing_on_boundaries(field3d<celle> & ce, field3d<cellb> & cb, field3d<cellj> & cj, double dt,
            double dx, double dy, double dz) = 0;
    virtual ~maxwell_solver() {}
};

class ndfx_solver: public maxwell_solver {
    virtual void fadvance(field3d<celle> & ce, field3d<cellb> & cb, field3d<cellj> & cj, double dt, double dx,
            double dy, double dz) override;
    virtual void fzeroing_on_boundaries(field3d<celle> & ce, field3d<cellb> & cb, field3d<cellj> & cj, double dt,
            double dx, double dy, double dz) override;
};

class fdtd_solver: public maxwell_solver {
    virtual void fadvance(field3d<celle> & ce, field3d<cellb> & cb, field3d<cellj> & cj, double dt, double dx,
            double dy, double dz) override;
    virtual void fzeroing_on_boundaries(field3d<celle> & ce, field3d<cellb> & cb, field3d<cellj> & cj, double dt,
            double dx, double dy, double dz) override;
};

#endif /* MAXWELL_H_ */
