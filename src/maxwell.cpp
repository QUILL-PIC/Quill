#include "maxwell.h"
#include <cmath>

const double PI = 3.141592653589793;

maxwell_solver::maxwell_solver(field3d<celle> & ce0, field3d<cellb> & cb0, field3d<cellj> & cj0, double dt, double dx,
        double dy, double dz) :
        ce(ce0), cb(cb0), cj(cj0), dx(dx), dy(dy), dz(dz), dtdx(dt / dx), dtdy(dt / dy), dtdz(dt / dz) { }

void maxwell_solver::advance() {
    // 1/2 b advance
    advance_b();
    advance_b_boundaries();

    // e advance
    advance_e();
    advance_e_boundaries();

    // 1/2 b advance
    advance_b();
    advance_b_boundaries();
}

five_point_solver::five_point_solver(field3d<celle> & ce0, field3d<cellb> & cb0, field3d<cellj> & cj0, double dt, double dx,
        double dy, double dz) :
        maxwell_solver(ce0, cb0, cj0, dt, dx, dy, dz),
        ax(0.25 * (1 - dx / dt * sin(0.5 * PI * std::min(1.0, dt / dx)))),
        ay(0.25 * (1 - dy / dt * sin(0.5 * PI * std::min(1.0, dt / dy)))),
        az(0.25 * (1 - dz / dt * sin(0.5 * PI * std::min(1.0, dt / dz)))) {
    bx = 1 - 2 * ax;
    by = 1 - 2 * ay;
    bz = 1 - 2 * az;
}

void five_point_solver::advance_b()
{
    const int nx = cb.get_nx();
    const int ny = cb.get_ny();
    const int nz = cb.get_nz();
    // 1/2 b advance
    for(int i=1;i<nx-2;i++) {
        for(int j=1;j<ny-2;j++) {
            for(int k=1;k<nz-2;k++) {
                cb[i][j][k].bx += 0.5*(   -dtdy*( by*(ce[i][j+1][k].ez-ce[i][j][k].ez) + ay*(ce[i][j+2][k].ez-ce[i][j+1][k].ez+ce[i][j][k].ez-ce[i][j-1][k].ez) )    +   dtdz*( bz*(ce[i][j][k+1].ey-ce[i][j][k].ey) + az*(ce[i][j][k+2].ey-ce[i][j][k+1].ey+ce[i][j][k].ey-ce[i][j][k-1].ey) )  );

                cb[i][j][k].by += 0.5*(   dtdx*( bx*(ce[i+1][j][k].ez-ce[i][j][k].ez) + ax*(ce[i+2][j][k].ez-ce[i+1][j][k].ez+ce[i][j][k].ez-ce[i-1][j][k].ez) )    -   dtdz*( bz*(ce[i][j][k+1].ex-ce[i][j][k].ex) + az*(ce[i][j][k+2].ex-ce[i][j][k+1].ex+ce[i][j][k].ex-ce[i][j][k-1].ex) )  );

                cb[i][j][k].bz += 0.5*(   dtdy*( by*(ce[i][j+1][k].ex-ce[i][j][k].ex) + ay*(ce[i][j+2][k].ex-ce[i][j+1][k].ex+ce[i][j][k].ex-ce[i][j-1][k].ex) )    -   dtdx*( bx*(ce[i+1][j][k].ey-ce[i][j][k].ey) + ax*(ce[i+2][j][k].ey-ce[i+1][j][k].ey+ce[i][j][k].ey-ce[i-1][j][k].ey))  );
            }
        }
    }
    {int i=nx-1;
        for(int j=1;j<ny-2;j++) {
            for(int k=1;k<nz-2;k++) {
                cb[i][j][k].bx += 0.5*(   -dtdy*( by*(ce[i][j+1][k].ez-ce[i][j][k].ez) + ay*(ce[i][j+2][k].ez-ce[i][j+1][k].ez+ce[i][j][k].ez-ce[i][j-1][k].ez) )    +   dtdz*( bz*(ce[i][j][k+1].ey-ce[i][j][k].ey) + az*(ce[i][j][k+2].ey-ce[i][j][k+1].ey+ce[i][j][k].ey-ce[i][j][k-1].ey) )  );

                cb[i][j][k].by += 0.5*(   dtdx*( bx*(-ce[i][j][k].ez) + ax*(ce[i][j][k].ez-ce[i-1][j][k].ez) )    -   dtdz*( bz*(ce[i][j][k+1].ex-ce[i][j][k].ex) + az*(ce[i][j][k+2].ex-ce[i][j][k+1].ex+ce[i][j][k].ex-ce[i][j][k-1].ex) )  );

                cb[i][j][k].bz += 0.5*(   dtdy*( by*(ce[i][j+1][k].ex-ce[i][j][k].ex) + ay*(ce[i][j+2][k].ex-ce[i][j+1][k].ex+ce[i][j][k].ex-ce[i][j-1][k].ex) )    -   dtdx*( bx*(-ce[i][j][k].ey) + ax*(ce[i][j][k].ey-ce[i-1][j][k].ey))  );
            }
        }
    }
    for(int i=1;i<nx-2;i++) {
        {int j=ny-1;
            for(int k=1;k<nz-2;k++) {
                cb[i][j][k].bx += 0.5*(   -dtdy*( by*(-ce[i][j][k].ez) + ay*(ce[i][j][k].ez-ce[i][j-1][k].ez) )    +   dtdz*( bz*(ce[i][j][k+1].ey-ce[i][j][k].ey) + az*(ce[i][j][k+2].ey-ce[i][j][k+1].ey+ce[i][j][k].ey-ce[i][j][k-1].ey) )  );

                cb[i][j][k].by += 0.5*(   dtdx*( bx*(ce[i+1][j][k].ez-ce[i][j][k].ez) + ax*(ce[i+2][j][k].ez-ce[i+1][j][k].ez+ce[i][j][k].ez-ce[i-1][j][k].ez) )    -   dtdz*( bz*(ce[i][j][k+1].ex-ce[i][j][k].ex) + az*(ce[i][j][k+2].ex-ce[i][j][k+1].ex+ce[i][j][k].ex-ce[i][j][k-1].ex) )  );

                cb[i][j][k].bz += 0.5*(   dtdy*( by*(-ce[i][j][k].ex) + ay*(ce[i][j][k].ex-ce[i][j-1][k].ex) )    -   dtdx*( bx*(ce[i+1][j][k].ey-ce[i][j][k].ey) + ax*(ce[i+2][j][k].ey-ce[i+1][j][k].ey+ce[i][j][k].ey-ce[i-1][j][k].ey))  );
            }
        }
    }
    for(int i=1;i<nx-2;i++) {
        for(int j=1;j<ny-2;j++) {
            {int k=nz-1;
                cb[i][j][k].bx += 0.5*(   -dtdy*( by*(ce[i][j+1][k].ez-ce[i][j][k].ez) + ay*(ce[i][j+2][k].ez-ce[i][j+1][k].ez+ce[i][j][k].ez-ce[i][j-1][k].ez) )    +   dtdz*( bz*(-ce[i][j][k].ey) + az*(ce[i][j][k].ey-ce[i][j][k-1].ey) )  );

                cb[i][j][k].by += 0.5*(   dtdx*( bx*(ce[i+1][j][k].ez-ce[i][j][k].ez) + ax*(ce[i+2][j][k].ez-ce[i+1][j][k].ez+ce[i][j][k].ez-ce[i-1][j][k].ez) )    -   dtdz*( bz*(-ce[i][j][k].ex) + az*(ce[i][j][k].ex-ce[i][j][k-1].ex) )  );

                cb[i][j][k].bz += 0.5*(   dtdy*( by*(ce[i][j+1][k].ex-ce[i][j][k].ex) + ay*(ce[i][j+2][k].ex-ce[i][j+1][k].ex+ce[i][j][k].ex-ce[i][j-1][k].ex) )    -   dtdx*( bx*(ce[i+1][j][k].ey-ce[i][j][k].ey) + ax*(ce[i+2][j][k].ey-ce[i+1][j][k].ey+ce[i][j][k].ey-ce[i-1][j][k].ey))  );
            }
        }
    }
}

void five_point_solver::advance_e()
{
    const int nx = ce.get_nx();
    const int ny = ce.get_ny();
    const int nz = ce.get_nz();

    for(int i=2;i<nx-2;i++) {
        for(int j=2;j<ny-2;j++) {
            for(int k=2;k<nz-2;k++) {
                ce[i][j][k].ex = ce[i][j][k].ex   +   dtdy*( by*(cb[i][j][k].bz-cb[i][j-1][k].bz) + ay*(cb[i][j+1][k].bz-cb[i][j][k].bz+cb[i][j-1][k].bz-cb[i][j-2][k].bz) )  -   dtdz*( bz*(cb[i][j][k].by-cb[i][j][k-1].by) + az*(cb[i][j][k+1].by-cb[i][j][k].by+cb[i][j][k-1].by-cb[i][j][k-2].by) )  -   cj[i][j][k].jx*dx;

                ce[i][j][k].ey = ce[i][j][k].ey   -   dtdx*( bx*(cb[i][j][k].bz-cb[i-1][j][k].bz) + ax*(cb[i+1][j][k].bz-cb[i][j][k].bz+cb[i-1][j][k].bz-cb[i-2][j][k].bz) )  +   dtdz*( bz*(cb[i][j][k].bx-cb[i][j][k-1].bx) + az*(cb[i][j][k+1].bx-cb[i][j][k].bx+cb[i][j][k-1].bx-cb[i][j][k-2].bx) )  -   cj[i][j][k].jy*dy;

                ce[i][j][k].ez = ce[i][j][k].ez   +   dtdx*( bx*(cb[i][j][k].by-cb[i-1][j][k].by) + ax*(cb[i+1][j][k].by-cb[i][j][k].by+cb[i-1][j][k].by-cb[i-2][j][k].by) )  -   dtdy*( by*(cb[i][j][k].bx-cb[i][j-1][k].bx) + ay*(cb[i][j+1][k].bx-cb[i][j][k].bx+cb[i][j-1][k].bx-cb[i][j-2][k].bx) )  -   cj[i][j][k].jz*dz;
            }
        }
    }
    {int i=1;
        for(int j=2;j<ny-2;j++) {
            for(int k=2;k<nz-2;k++) {
                ce[i][j][k].ex = ce[i][j][k].ex   +   dtdy*( by*(cb[i][j][k].bz-cb[i][j-1][k].bz) + ay*(cb[i][j+1][k].bz-cb[i][j][k].bz+cb[i][j-1][k].bz-cb[i][j-2][k].bz) )  -   dtdz*( bz*(cb[i][j][k].by-cb[i][j][k-1].by) + az*(cb[i][j][k+1].by-cb[i][j][k].by+cb[i][j][k-1].by-cb[i][j][k-2].by) )  -   cj[i][j][k].jx*dx;
            }
        }
    }
    for(int i=2;i<nx-2;i++) {
        {int j=1;
            for(int k=2;k<nz-2;k++) {
                ce[i][j][k].ey = ce[i][j][k].ey   -   dtdx*( bx*(cb[i][j][k].bz-cb[i-1][j][k].bz) + ax*(cb[i+1][j][k].bz-cb[i][j][k].bz+cb[i-1][j][k].bz-cb[i-2][j][k].bz) )  +   dtdz*( bz*(cb[i][j][k].bx-cb[i][j][k-1].bx) + az*(cb[i][j][k+1].bx-cb[i][j][k].bx+cb[i][j][k-1].bx-cb[i][j][k-2].bx) )  -   cj[i][j][k].jy*dy;
            }
        }
    }
    for(int i=2;i<nx-2;i++) {
        for(int j=2;j<ny-2;j++) {
            {int k=1;
                ce[i][j][k].ez = ce[i][j][k].ez   +   dtdx*( bx*(cb[i][j][k].by-cb[i-1][j][k].by) + ax*(cb[i+1][j][k].by-cb[i][j][k].by+cb[i-1][j][k].by-cb[i-2][j][k].by) )  -   dtdy*( by*(cb[i][j][k].bx-cb[i][j-1][k].bx) + ay*(cb[i][j+1][k].bx-cb[i][j][k].bx+cb[i][j-1][k].bx-cb[i][j-2][k].bx) )  -   cj[i][j][k].jz*dz;
            }
        }
    }
}

void five_point_solver::advance_b_boundaries() {
    // заготовка для учета гранусловий
}

void five_point_solver::advance_e_boundaries() {
    // заготовка для учета гранусловий
}

void five_point_solver::init_boundaries()
{
    const int nx = ce.get_nx();
    const int ny = ce.get_ny();
    const int nz = ce.get_nz();
    {int i=0;
        for(int j=0;j<ny-1;j++) {
            for(int k=0;k<nz-1;k++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = 0;
                cb[i][j][k].bz = 0;
            }
        }
    }
    {int j=0;
        for(int i=0;i<nx-1;i++) {
            for(int k=0;k<nz-1;k++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = 0;
                cb[i][j][k].bz = 0;
            }
        }
    }
    {int k=0;
        for(int i=0;i<nx-1;i++) {
            for(int j=0;j<ny-1;j++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = 0;
                cb[i][j][k].bz = 0;
            }
        }
    }
    {int i=nx-1;
        for(int j=0;j<ny-1;j++) {
            for(int k=0;k<nz-1;k++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = 0;
                cb[i][j][k].bz = 0;
            }
        }
    }
    {int j=ny-1;
        for(int i=0;i<nx-1;i++) {
            for(int k=0;k<nz-1;k++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = 0;
                cb[i][j][k].bz = 0;
            }
        }
    }
    {int k=nz-1;
        for(int i=0;i<nx-1;i++) {
            for(int j=0;j<ny-1;j++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = 0;
                cb[i][j][k].bz = 0;
            }
        }
    }
    for(int j=0;j<ny-1;j++) {
        for(int k=0;k<nz-1;k++) {
            cb[nx-2][j][k].by = 0;
            cb[nx-2][j][k].bz = 0;
            ce[1][j][k].ez = 0;
            ce[1][j][k].ey = 0;
        }
    }
    for(int i=0;i<nx-1;i++) {
        for(int k=0;k<nz-1;k++) {
            cb[i][ny-2][k].bx = 0;
            cb[i][ny-2][k].bz = 0;
            ce[i][1][k].ex = 0;
            ce[i][1][k].ez = 0;
        }
    }
    for(int j=0;j<ny-1;j++) {
        for(int i=0;i<nx-1;i++) {
            cb[i][j][nz-2].bx = 0;
            cb[i][j][nz-2].by = 0;
            ce[i][j][1].ex = 0;
            ce[i][j][1].ey = 0;
        }
    }
}

const double ndfx_solver::kappa = 1.25 * (sqrt(3) - 1) * 0.5 / sqrt(3);
const double ndfx_solver::bx = (1 - kappa);
const double ndfx_solver::ax = 0.5 * (1 - bx);

ndfx_solver::ndfx_solver(field3d<celle> & ce0, field3d<cellb> & cb0, field3d<cellj> & cj0, double dt, double dx,
        double dy, double dz) :
        maxwell_solver(ce0, cb0, cj0, dt, dx, dy, dz),
        by(1 - kappa * dx * dx / (dy * dy)),
        bz(1 - kappa * dx * dx / (dz * dz)) {
    ay = 0.5 * (1 - by);
    az = 0.5 * (1 - bz);
}

void ndfx_solver::advance_b()
{
    const int nx = cb.get_nx();
    const int ny = cb.get_ny();
    const int nz = cb.get_nz();
    // 1/2 b advance
    for(int i=2;i<nx-1;i++) {
        for(int j=2;j<ny-1;j++) {
            for(int k=2;k<nz-1;k++) {
                cb[i][j][k].bx += 0.5*(   -dtdy*( bz*(ce[i][j+1][k].ez-ce[i][j][k].ez) + az*(ce[i][j+1][k+1].ez-ce[i][j][k+1].ez+ce[i][j+1][k-1].ez-ce[i][j][k-1].ez) )    +   dtdz*( by*(ce[i][j][k+1].ey-ce[i][j][k].ey) + ay*(ce[i][j+1][k+1].ey-ce[i][j+1][k].ey+ce[i][j-1][k+1].ey-ce[i][j-1][k].ey) )  );
                cb[i][j][k].by += 0.5*(   dtdx*( bz*(ce[i+1][j][k].ez-ce[i][j][k].ez) + az*(ce[i+1][j][k+1].ez-ce[i][j][k+1].ez+ce[i+1][j][k-1].ez-ce[i][j][k-1].ez) )    -   dtdz*( bx*(ce[i][j][k+1].ex-ce[i][j][k].ex) + ax*(ce[i+1][j][k+1].ex-ce[i+1][j][k].ex+ce[i-1][j][k+1].ex-ce[i-1][j][k].ex) )  );
                cb[i][j][k].bz += 0.5*(   dtdy*( bx*(ce[i][j+1][k].ex-ce[i][j][k].ex) + ax*(ce[i+1][j+1][k].ex-ce[i+1][j][k].ex+ce[i-1][j+1][k].ex-ce[i-1][j][k].ex) )    -   dtdx*( by*(ce[i+1][j][k].ey-ce[i][j][k].ey) + ay*(ce[i+1][j+1][k].ey-ce[i][j+1][k].ey+ce[i+1][j-1][k].ey-ce[i][j-1][k].ey) )  );
            }
        }
    }
    {int i=1;
        for(int j=2;j<ny-1;j++) {
            for(int k=2;k<nz-1;k++) {
                cb[i][j][k].by += 0.5*(   dtdx*( bz*(ce[i+1][j][k].ez-ce[i][j][k].ez) + az*(ce[i+1][j][k+1].ez-ce[i][j][k+1].ez+ce[i+1][j][k-1].ez-ce[i][j][k-1].ez) )    -   dtdz*( bx*(ce[i][j][k+1].ex-ce[i][j][k].ex) + ax*(ce[i+1][j][k+1].ex-ce[i+1][j][k].ex+ce[i-1][j][k+1].ex-ce[i-1][j][k].ex) )  );
                cb[i][j][k].bz += 0.5*(   dtdy*( bx*(ce[i][j+1][k].ex-ce[i][j][k].ex) + ax*(ce[i+1][j+1][k].ex-ce[i+1][j][k].ex+ce[i-1][j+1][k].ex-ce[i-1][j][k].ex) )    -   dtdx*( by*(ce[i+1][j][k].ey-ce[i][j][k].ey) + ay*(ce[i+1][j+1][k].ey-ce[i][j+1][k].ey+ce[i+1][j-1][k].ey-ce[i][j-1][k].ey) )  );
            }
        }
    }
    for(int i=2;i<nx-1;i++) {
        {int j=1;
            for(int k=2;k<nz-1;k++) {
                cb[i][j][k].bx += 0.5*(   -dtdy*( bz*(ce[i][j+1][k].ez-ce[i][j][k].ez) + az*(ce[i][j+1][k+1].ez-ce[i][j][k+1].ez+ce[i][j+1][k-1].ez-ce[i][j][k-1].ez) )    +   dtdz*( by*(ce[i][j][k+1].ey-ce[i][j][k].ey) + ay*(ce[i][j+1][k+1].ey-ce[i][j+1][k].ey+ce[i][j-1][k+1].ey-ce[i][j-1][k].ey) )  );
                cb[i][j][k].bz += 0.5*(   dtdy*( bx*(ce[i][j+1][k].ex-ce[i][j][k].ex) + ax*(ce[i+1][j+1][k].ex-ce[i+1][j][k].ex+ce[i-1][j+1][k].ex-ce[i-1][j][k].ex) )    -   dtdx*( by*(ce[i+1][j][k].ey-ce[i][j][k].ey) + ay*(ce[i+1][j+1][k].ey-ce[i][j+1][k].ey+ce[i+1][j-1][k].ey-ce[i][j-1][k].ey) )  );
            }
        }
    }
    for(int i=2;i<nx-1;i++) {
        for(int j=2;j<ny-1;j++) {
            {int k=1;
                cb[i][j][k].bx += 0.5*(   -dtdy*( bz*(ce[i][j+1][k].ez-ce[i][j][k].ez) + az*(ce[i][j+1][k+1].ez-ce[i][j][k+1].ez+ce[i][j+1][k-1].ez-ce[i][j][k-1].ez) )    +   dtdz*( by*(ce[i][j][k+1].ey-ce[i][j][k].ey) + ay*(ce[i][j+1][k+1].ey-ce[i][j+1][k].ey+ce[i][j-1][k+1].ey-ce[i][j-1][k].ey) )  );
                cb[i][j][k].by += 0.5*(   dtdx*( bz*(ce[i+1][j][k].ez-ce[i][j][k].ez) + az*(ce[i+1][j][k+1].ez-ce[i][j][k+1].ez+ce[i+1][j][k-1].ez-ce[i][j][k-1].ez) )    -   dtdz*( bx*(ce[i][j][k+1].ex-ce[i][j][k].ex) + ax*(ce[i+1][j][k+1].ex-ce[i+1][j][k].ex+ce[i-1][j][k+1].ex-ce[i-1][j][k].ex) )  );
            }
        }
    }
}

void ndfx_solver::advance_e()
{
    const int nx = ce.get_nx();
    const int ny = ce.get_ny();
    const int nz = ce.get_nz();

    for(int i=2;i<nx-1;i++) {
        for(int j=2;j<ny-1;j++) {
            for(int k=2;k<nz-1;k++) {
                ce[i][j][k].ex = ce[i][j][k].ex   +   dtdy*( bz*(cb[i][j][k].bz-cb[i][j-1][k].bz) + az*(cb[i][j][k+1].bz-cb[i][j-1][k+1].bz+cb[i][j][k-1].bz-cb[i][j-1][k-1].bz) )  -   dtdz*( by*(cb[i][j][k].by-cb[i][j][k-1].by) + ay*(cb[i][j+1][k].by-cb[i][j+1][k-1].by+cb[i][j-1][k].by-cb[i][j-1][k-1].by) )  -   cj[i][j][k].jx*dx;
                ce[i][j][k].ey = ce[i][j][k].ey   -   dtdx*( bz*(cb[i][j][k].bz-cb[i-1][j][k].bz) + az*(cb[i][j][k+1].bz-cb[i-1][j][k+1].bz+cb[i][j][k-1].bz-cb[i-1][j][k-1].bz) )  +   dtdz*( bx*(cb[i][j][k].bx-cb[i][j][k-1].bx) + ax*(cb[i+1][j][k].bx-cb[i+1][j][k-1].bx+cb[i-1][j][k].bx-cb[i-1][j][k-1].bx) )  -   cj[i][j][k].jy*dy;
                ce[i][j][k].ez = ce[i][j][k].ez   +   dtdx*( by*(cb[i][j][k].by-cb[i-1][j][k].by) + ay*(cb[i][j+1][k].by-cb[i-1][j+1][k].by+cb[i][j-1][k].by-cb[i-1][j-1][k].by) )  -   dtdy*( bx*(cb[i][j][k].bx-cb[i][j-1][k].bx) + ax*(cb[i+1][j][k].bx-cb[i+1][j-1][k].bx+cb[i-1][j][k].bx-cb[i-1][j-1][k].bx) )  -   cj[i][j][k].jz*dz;
            }
        }
    }
    {int i=1;
        for(int j=2;j<ny-1;j++) {
            for(int k=2;k<nz-1;k++) {
                ce[i][j][k].ex = ce[i][j][k].ex   +   dtdy*( bz*(cb[i][j][k].bz-cb[i][j-1][k].bz) + az*(cb[i][j][k+1].bz-cb[i][j-1][k+1].bz+cb[i][j][k-1].bz-cb[i][j-1][k-1].bz) )  -   dtdz*( by*(cb[i][j][k].by-cb[i][j][k-1].by) + ay*(cb[i][j+1][k].by-cb[i][j+1][k-1].by+cb[i][j-1][k].by-cb[i][j-1][k-1].by) )  -   cj[i][j][k].jx*dx;
            }
        }
    }
    for(int i=2;i<nx-1;i++) {
        {int j=1;
            for(int k=2;k<nz-1;k++) {
                ce[i][j][k].ey = ce[i][j][k].ey   -   dtdx*( bz*(cb[i][j][k].bz-cb[i-1][j][k].bz) + az*(cb[i][j][k+1].bz-cb[i-1][j][k+1].bz+cb[i][j][k-1].bz-cb[i-1][j][k-1].bz) )  +   dtdz*( bx*(cb[i][j][k].bx-cb[i][j][k-1].bx) + ax*(cb[i+1][j][k].bx-cb[i+1][j][k-1].bx+cb[i-1][j][k].bx-cb[i-1][j][k-1].bx) )  -   cj[i][j][k].jy*dy;
            }
        }
    }
    for(int i=2;i<nx-1;i++) {
        for(int j=2;j<ny-1;j++) {
            {int k=1;
                ce[i][j][k].ez = ce[i][j][k].ez   +   dtdx*( by*(cb[i][j][k].by-cb[i-1][j][k].by) + ay*(cb[i][j+1][k].by-cb[i-1][j+1][k].by+cb[i][j-1][k].by-cb[i-1][j-1][k].by) )  -   dtdy*( bx*(cb[i][j][k].bx-cb[i][j-1][k].bx) + ax*(cb[i+1][j][k].bx-cb[i+1][j-1][k].bx+cb[i-1][j][k].bx-cb[i-1][j-1][k].bx) )  -   cj[i][j][k].jz*dz;
            }
        }
    }
}

void ndfx_solver::advance_b_boundaries() {
    // заготовка для учета гранусловий
}

void ndfx_solver::advance_e_boundaries() {
    // заготовка для учета гранусловий
}

void ndfx_solver::init_boundaries()
{
    const int nx = ce.get_nx();
    const int ny = ce.get_ny();
    const int nz = ce.get_nz();
    /* обнуление полей на границах (иначе ndfx с металлическими
     * границами не сохраняет энергию) */
    {int i=0;
        for(int j=0;j<ny-1;j++) {
            for(int k=0;k<nz-1;k++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = 0;
                cb[i][j][k].bz = 0;
            }
        }
        i=1;
        for(int j=0;j<ny-1;j++) {
            for(int k=0;k<nz-1;k++) {
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
            }
        }
        i=nx-1;
        for(int j=0;j<ny-1;j++) {
            for(int k=0;k<nz-1;k++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = 0;
                cb[i][j][k].bz = 0;
            }
        }
    }
    {int j=0;
        for(int i=0;i<nx-1;i++) {
            for(int k=0;k<nz-1;k++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = 0;
                cb[i][j][k].bz = 0;
            }
        }
        j=1;
        for(int i=0;i<nx-1;i++) {
            for(int k=0;k<nz-1;k++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].by = 0;
            }
        }
        j=ny-1;
        for(int i=0;i<nx-1;i++) {
            for(int k=0;k<nz-1;k++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = 0;
                cb[i][j][k].bz = 0;
            }
        }
    }
    {int k=0;
        for(int i=0;i<nx-1;i++) {
            for(int j=0;j<ny-1;j++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = 0;
                cb[i][j][k].bz = 0;
            }
        }
        k=1;
        for(int i=0;i<nx-1;i++) {
            for(int j=0;j<ny-1;j++) {
                ce[i][j][k].ey = 0;
                ce[i][j][k].ex = 0;
                cb[i][j][k].bz = 0;
            }
        }
        k=nz-1;
        for(int i=0;i<nx-1;i++) {
            for(int j=0;j<ny-1;j++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = 0;
                cb[i][j][k].bz = 0;
            }
        }
    }
}

fdtd_solver::fdtd_solver(field3d<celle> & ce0, field3d<cellb> & cb0, field3d<cellj> & cj0, double dt, double dx,
        double dy, double dz) :
        maxwell_solver(ce0, cb0, cj0, dt, dx, dy, dz) {
}

void fdtd_solver::advance_b()
{
    const int nx = cb.get_nx();
    const int ny = cb.get_ny();
    const int nz = cb.get_nz();

    for(int i=1;i<nx-1;i++) {
        for(int j=1;j<ny-1;j++) {
            for(int k=1;k<nz-1;k++) {
                cb[i][j][k].bx += 0.5*(-dtdy*(ce[i][j+1][k].ez-ce[i][j][k].ez) + dtdz*(ce[i][j][k+1].ey-ce[i][j][k].ey));
                cb[i][j][k].by += 0.5*( dtdx*(ce[i+1][j][k].ez-ce[i][j][k].ez) - dtdz*(ce[i][j][k+1].ex-ce[i][j][k].ex));
                cb[i][j][k].bz += 0.5*( dtdy*(ce[i][j+1][k].ex-ce[i][j][k].ex) - dtdx*(ce[i+1][j][k].ey-ce[i][j][k].ey));
            }
        }
    }
    {int i=0;
        for(int j=1;j<ny-1;j++) {
            for(int k=1;k<nz-1;k++) {
                cb[i][j][k].by += 0.5*( dtdx*(ce[i+1][j][k].ez-ce[i][j][k].ez) - dtdz*(ce[i][j][k+1].ex-ce[i][j][k].ex));
                cb[i][j][k].bz += 0.5*( dtdy*(ce[i][j+1][k].ex-ce[i][j][k].ex) - dtdx*(ce[i+1][j][k].ey-ce[i][j][k].ey));
            }
        }
    }
    for(int i=1;i<nx-1;i++) {
        {int j=0;
            for(int k=1;k<nz-1;k++) {
                cb[i][j][k].bx += 0.5*(-dtdy*(ce[i][j+1][k].ez-ce[i][j][k].ez) + dtdz*(ce[i][j][k+1].ey-ce[i][j][k].ey));
                cb[i][j][k].bz += 0.5*( dtdy*(ce[i][j+1][k].ex-ce[i][j][k].ex) - dtdx*(ce[i+1][j][k].ey-ce[i][j][k].ey));
            }
        }
    }
    for(int i=1;i<nx-1;i++) {
        for(int j=1;j<ny-1;j++) {
            {int k=0;
                cb[i][j][k].bx += 0.5*(-dtdy*(ce[i][j+1][k].ez-ce[i][j][k].ez) + dtdz*(ce[i][j][k+1].ey-ce[i][j][k].ey));
                cb[i][j][k].by += 0.5*( dtdx*(ce[i+1][j][k].ez-ce[i][j][k].ez) - dtdz*(ce[i][j][k+1].ex-ce[i][j][k].ex));
            }
        }
    }
}

void fdtd_solver::advance_e()
{
    const int nx = ce.get_nx();
    const int ny = ce.get_ny();
    const int nz = ce.get_nz();

    for(int i=1;i<nx-1;i++) {
        for(int j=1;j<ny-1;j++) {
            for(int k=1;k<nz-1;k++) {
                ce[i][j][k].ex +=  dtdy*(cb[i][j][k].bz-cb[i][j-1][k].bz) - dtdz*(cb[i][j][k].by-cb[i][j][k-1].by) - cj[i][j][k].jx*dx;
                ce[i][j][k].ey += -dtdx*(cb[i][j][k].bz-cb[i-1][j][k].bz) + dtdz*(cb[i][j][k].bx-cb[i][j][k-1].bx) - cj[i][j][k].jy*dy;
                ce[i][j][k].ez +=  dtdx*(cb[i][j][k].by-cb[i-1][j][k].by) - dtdy*(cb[i][j][k].bx-cb[i][j-1][k].bx) - cj[i][j][k].jz*dz;
            }
        }
    }
    {int i=0;
        for(int j=1;j<ny-1;j++) {
            for(int k=1;k<nz-1;k++) {
                ce[i][j][k].ex +=  dtdy*(cb[i][j][k].bz-cb[i][j-1][k].bz) - dtdz*(cb[i][j][k].by-cb[i][j][k-1].by) - cj[i][j][k].jx*dx;
            }
        }
    }
    for(int i=1;i<nx-1;i++) {
        {int j=0;
            for(int k=1;k<nz-1;k++) {
                ce[i][j][k].ey += -dtdx*(cb[i][j][k].bz-cb[i-1][j][k].bz) + dtdz*(cb[i][j][k].bx-cb[i][j][k-1].bx) - cj[i][j][k].jy*dy;
            }
        }
    }
    for(int i=1;i<nx-1;i++) {
        for(int j=1;j<ny-1;j++) {
            {int k=0;
                ce[i][j][k].ez +=  dtdx*(cb[i][j][k].by-cb[i-1][j][k].by) - dtdy*(cb[i][j][k].bx-cb[i][j-1][k].bx) - cj[i][j][k].jz*dz;
            }
        }
    }
}

void fdtd_solver::advance_b_boundaries() {
    // заготовка для учета гранусловий
}

void fdtd_solver::advance_e_boundaries() {
    // заготовка для учета гранусловий
}

void fdtd_solver::init_boundaries()
{
    const int nx = ce.get_nx();
    const int ny = ce.get_ny();
    const int nz = ce.get_nz();
    {int i=0;
        for(int j=0;j<ny-1;j++) {
            for(int k=0;k<nz-1;k++) {
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
            }
        }
        i=nx-1;
        for(int j=0;j<ny-1;j++) {
            for(int k=0;k<nz-1;k++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = 0;
                cb[i][j][k].bz = 0;
            }
        }
    }
    {int j=0;
        for(int i=0;i<nx-1;i++) {
            for(int k=0;k<nz-1;k++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].by = 0;
            }
        }
        j=ny-1;
        for(int i=0;i<nx-1;i++) {
            for(int k=0;k<nz-1;k++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = 0;
                cb[i][j][k].bz = 0;
            }
        }
    }
    {int k=0;
        for(int i=0;i<nx-1;i++) {
            for(int j=0;j<ny-1;j++) {
                ce[i][j][k].ey = 0;
                ce[i][j][k].ex = 0;
                cb[i][j][k].bz = 0;
            }
        }
        k=nz-1;
        for(int i=0;i<nx-1;i++) {
            for(int j=0;j<ny-1;j++) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = 0;
                cb[i][j][k].bz = 0;
            }
        }
    }
}

