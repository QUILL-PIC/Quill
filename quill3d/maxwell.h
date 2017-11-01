#ifndef MAXWELL_H_
#define MAXWELL_H_

class maxwell_solver {
public:
    virtual void fadvance(spatial_region::celle*** ce, spatial_region::cellb*** cb, spatial_region::cellj*** cj, double dt,
            double dx, double dy, double dz, int nx, int ny, int nz) = 0;
    virtual void fzeroing_on_boundaries(spatial_region::celle*** ce, spatial_region::cellb*** cb, spatial_region::cellj*** cj, double dt,
            double dx, double dy, double dz, int nx, int ny, int nz) = 0;
    virtual ~maxwell_solver() {};
};

class ndfx_solver : public maxwell_solver {
    virtual void fadvance(spatial_region::celle*** ce, spatial_region::cellb*** cb, spatial_region::cellj*** cj,
            double dt, double dx, double dy, double dz, int nx, int ny, int nz) override;
    virtual void fzeroing_on_boundaries(spatial_region::celle*** ce, spatial_region::cellb*** cb,
            spatial_region::cellj*** cj, double dt, double dx, double dy, double dz, int nx, int ny, int nz) override;
};

class fdtd_solver : public maxwell_solver {
    virtual void fadvance(spatial_region::celle*** ce, spatial_region::cellb*** cb, spatial_region::cellj*** cj,
            double dt, double dx, double dy, double dz, int nx, int ny, int nz) override;
    virtual void fzeroing_on_boundaries(spatial_region::celle*** ce, spatial_region::cellb*** cb,
            spatial_region::cellj*** cj, double dt, double dx, double dy, double dz, int nx, int ny, int nz) override;
};

#endif /* MAXWELL_H_ */
