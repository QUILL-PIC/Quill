#ifndef MAXWELL_H_
#define MAXWELL_H_

#include <math.h>
#include "containers.h"

class maxwell_solver {
protected:
    field3d<celle> & ce;
    field3d<cellb> & cb;
    field3d<cellj> & cj;
    const double dx;
    const double dy;
    const double dz;
    const double dtdx;
    const double dtdy;
    const double dtdz;
    virtual void advance_b() = 0;
    virtual void advance_e() = 0;
    virtual void advance_b_boundaries() = 0;
    virtual void advance_e_boundaries() = 0;
public:
    maxwell_solver(field3d<celle> & ce0, field3d<cellb> & cb0, field3d<cellj> & cj0, double dt, double dx, double dy,
            double dz);
    virtual void advance();
    virtual void init_boundaries() = 0;
    virtual ~maxwell_solver() {}
};

class fp_solver: public maxwell_solver {
private:
    double ax;
    double ay;
    double az;
    double bx;
    double by;
    double bz;
    virtual void advance_b() override;
    virtual void advance_e() override;
    virtual void advance_b_boundaries() override;
    virtual void advance_e_boundaries() override;
public:
    fp_solver(field3d<celle> & ce0, field3d<cellb> & cb0, field3d<cellj> & cj0, double dt, double dx, double dy,
            double dz);
    virtual void init_boundaries() override;
};

class ndfx_solver: public maxwell_solver {
private:
    const static double kappa;
    const static double bx;
    const static double ax;
    double ay;
    double by;
    double az;
    double bz;
    virtual void advance_b() override;
    virtual void advance_e() override;
    virtual void advance_b_boundaries() override;
    virtual void advance_e_boundaries() override;
public:
    ndfx_solver(field3d<celle> & ce0, field3d<cellb> & cb0, field3d<cellj> & cj0, double dt, double dx, double dy,
            double dz);
    virtual void init_boundaries() override;
};

class fdtd_solver: public maxwell_solver {
private:
    virtual void advance_b() override;
    virtual void advance_e() override;
    virtual void advance_b_boundaries() override;
    virtual void advance_e_boundaries() override;
public:
    fdtd_solver(field3d<celle> & ce0, field3d<cellb> & cb0, field3d<cellj> & cj0, double dt, double dx, double dy,
            double dz);
    virtual void init_boundaries() override;
};

#endif /* MAXWELL_H_ */
