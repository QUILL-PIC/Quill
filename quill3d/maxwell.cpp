#include "maxwell.h"
#include <cmath>

void ndfx_solver::fadvance(field3d<celle> & ce, field3d<cellb> & cb, field3d<cellj> & cj, double dt, double dx,
        double dy, double dz)
{
    /* «Металлические» границы, наиболее быстрые циклы 'for'; для
     * приблизительного сохранения энергии требуется зануление
     * начальных распределений полей в f_init под металлическими
     * границами */
    constexpr double kappa=1.25*(sqrt(3)-1)*0.5/sqrt(3);
    constexpr double bx = 1 - kappa; // approx 0.74
    const double by = 1 - kappa*dx*dx/(dy*dy);
    const double bz = 1 - kappa*dx*dx/(dz*dz);
    const double ax = 0.5*( 1 - bx );
    const double ay = 0.5*( 1 - by );
    const double az = 0.5*( 1 - bz );
    const double dtdx=dt/dx;
    const double dtdy=dt/dy;
    const double dtdz=dt/dz;
    const int nx = ce.get_nx();
    const int ny = ce.get_ny();
    const int nz = ce.get_nz();
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
    // e advance
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

void ndfx_solver::fzeroing_on_boundaries(field3d<celle> & ce, field3d<cellb> & cb, field3d<cellj> & cj, double dt,
        double dx, double dy, double dz)
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

void fdtd_solver::fadvance(field3d<celle> & ce, field3d<cellb> & cb, field3d<cellj> & cj, double dt, double dx,
        double dy, double dz)
{
    const double dtdx=dt/dx;
    const double dtdy=dt/dy;
    const double dtdz=dt/dz;
    const int nx = ce.get_nx();
    const int ny = ce.get_ny();
    const int nz = ce.get_nz();
    // 1/2 b advance
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
    // e advance
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
    // 1/2 b advance
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

void fdtd_solver::fzeroing_on_boundaries(field3d<celle> & ce, field3d<cellb> & cb, field3d<cellj> & cj, double dt,
        double dx, double dy, double dz)
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

