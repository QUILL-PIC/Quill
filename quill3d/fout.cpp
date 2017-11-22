#include <iostream>
#include <fstream>
#include "main.h"
#include "containers.h"

void spatial_region::compute_rho()
{
    particle* current;
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int k=0;k<nz;k++)
            {
                cj[i][j][k].jx = 0;
                cj[i][j][k].jy = 0;
                cj[i][j][k].jz = 0;
                for (int n=0;n<n_ion_populations;n++)
                    irho[n][i][j][k] = 0;
            }
        }
    }
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int k=0;k<nz;k++)
            {
                current = cp[i][j][k].pl.start;
                while(current!=0)
                {
                    rhodeposition(*current);
                    current = current->next;
                }
            }
        }
    }
}

void spatial_region::compute_N(int n1, int n2, double a)
{
    N_e = 0;
    N_p = 0;
    N_ph = 0;
    N_qp_e = 0;
    N_qp_p = 0;
    N_qp_g = 0;
    for (int i=0;i<n_ion_populations;i++)
        N_qp_i[i] = 0;
    particle* current;
    for (int i=n1; i<nx-n2; i++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int k=0; k<nz; k++)
            {
                current = cp[i][j][k].pl.head;
                while (current!=0)
                {
                    if (current->cmr==-1) {
                        N_e -= current->q;
                        N_qp_e += 1;
                    }
                    else if (current->cmr==1) {
                        N_p += current->q;
                        N_qp_p += 1;
                    }
                    else if (current->cmr==0) {
                        N_ph += current->q;
                        N_qp_g += 1;
                    }
                    for (int ii=0;ii<n_ion_populations;ii++) {
                        if (current->cmr==icmr[ii])
                            N_qp_i[ii] += 1;
                    }
                    current = current->next;
                }
            }
        }
    }
    N_e = N_e*a;
    N_p = N_p*a;
    N_ph = N_ph*a;
}

void spatial_region::compute_energy(int n1, int n2, double a, double b)
{
    energy_f = 0;
    energy_e = 0;
    energy_p = 0;
    energy_ph = 0;
    for (int n=0;n<n_ion_populations;n++)
        ienergy[n] = 0;
    particle* current;
    for (int i=n1; i<nx-n2; i++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int k=0; k<nz; k++)
            {
                energy_f += ce[i][j][k].ex*ce[i][j][k].ex + ce[i][j][k].ey*ce[i][j][k].ey + ce[i][j][k].ez*ce[i][j][k].ez + cb[i][j][k].bx*cb[i][j][k].bx + cb[i][j][k].by*cb[i][j][k].by + cb[i][j][k].bz*cb[i][j][k].bz;
                current = cp[i][j][k].pl.head;
                while (current!=0)
                {
                    if (current->cmr==-1)
                        energy_e -= (current->g - 1)*current->q;
                    else if (current->cmr==1)
                        energy_p += (current->g - 1)*current->q;
                    else if (current->cmr==0)
                        energy_ph += current->g*current->q;
                    else
                    {
                        int n;
                        n = 0;
                        while (n!=n_ion_populations)
                        {
                            if (current->cmr==icmr[n])
                            {
                                ienergy[n] += (current->g-1)*current->q/icmr[n];
                                n = n_ion_populations;
                            }
                            else
                                n++;
                        }
                    }
                    current = current->next;
                }
            }
        }
    }
    energy_f = energy_f*a;
    energy_e = energy_e*b;
    energy_p = energy_p*b;
    energy_ph = energy_ph*b;
    for (int n=0;n<n_ion_populations;n++)
        ienergy[n] = ienergy[n]*b;
}

double getjx(field3d<cellj> & cj, int i, int j, int k) {
    return cj[i][j][k].jx;
}

double getjy(field3d<cellj> & cj, int i, int j, int k) {
    return cj[i][j][k].jy;
}

double getjz(field3d<cellj> & cj, int i, int j, int k) {
    return cj[i][j][k].jz;
}

double getdouble(field3d<double> & a, int i, int j, int k) {
    return a[i][j][k];
}

double getex(field3d<celle> & ce, int i, int j, int k) {
    return ce[i][j][k].ex;
}

double getey(field3d<celle> & ce, int i, int j, int k) {
    return ce[i][j][k].ey;
}

double getez(field3d<celle> & ce, int i, int j, int k) {
    return ce[i][j][k].ez;
}

double getbx(field3d<cellb> & cb, int i, int j, int k) {
    return cb[i][j][k].bx;
}

double getby(field3d<cellb> & cb, int i, int j, int k) {
    return cb[i][j][k].by;
}

double getbz(field3d<cellb> & cb, int i, int j, int k) {
    return cb[i][j][k].bz;
}

template<typename T> void write_binary(ofstream* f, field3d<T> & data, double get(field3d<T>&,
            int, int, int), int n0, int n, int ny, int nz) {
    double* a = new double[(n - n0) * (ny + nz)]; /* for values in xy and xz
                                                     planes */
    for (int i = n0; i < n; ++i) {
        int k = nz / 2;
        for (int j = 0; j < ny; ++j)
            a[(i - n0) * (ny + nz) + j] = get(data, i, j, k);
        int j = ny / 2;
        for (int k = 0; k < nz; ++k)
            a[(i - n0) * (ny + nz) + ny + k] = get(data, i, j, k);
    }
    f->write(reinterpret_cast<char*>(a), sizeof(double) * (n - n0) * (ny +
                nz));
    delete[] a;
}

template<typename T> void write_binary_yzplane(ofstream* f, field3d<T> & data, double
        get(field3d<T>&, int, int, int), int i, int ny, int nz) {
    double* a = new double[ny * nz];
    for(int j=0;j<ny;j++)
        for(int k=0;k<nz;k++)
            a[j * nz + k] = get(data, i, j, k);
    f->write(reinterpret_cast<char*>(a), sizeof(double) * (ny * nz));
    delete[] a;
}

void spatial_region::fout_ex(ofstream* pfout, int n0, int n, ios_base::openmode
        mode)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    if (mode == ios_base::out) {
        for(int i=n0;i<n;i++)
        {
            int k=nz/2;
            for(int j=0;j<ny;j++)
                (*pfout)<<ce[i][j][k].ex<<"\n";
            int j=ny/2;
            for(int k=0;k<nz;k++)
                (*pfout)<<ce[i][j][k].ex<<"\n";
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        write_binary<celle>(pfout, ce, getex, n0, n, ny, nz);
    } else {
        cerr << "fout_ex: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_ex_yzplane(ofstream* pfout, int i, ios_base::openmode
        mode)
{
    if (mode == ios_base::out) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++)
                (*pfout)<<ce[i][j][k].ex<<"\n";
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        write_binary_yzplane<celle>(pfout, ce, getex, i, ny, nz);
    } else {
        cerr << "fout_ex_yzplane: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_ey(ofstream* pfout, int n0, int n, ios_base::openmode
        mode)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    if (mode == ios_base::out) {
        for(int i=n0;i<n;i++)
        {
            int k=nz/2;
            for(int j=0;j<ny;j++)
                (*pfout)<<ce[i][j][k].ey<<"\n";
            int j=ny/2;
            for(int k=0;k<nz;k++)
                (*pfout)<<ce[i][j][k].ey<<"\n";
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        write_binary<celle>(pfout, ce, getey, n0, n, ny, nz);
    } else {
        cerr << "fout_ey: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_ey_yzplane(ofstream* pfout, int i, ios_base::openmode
        mode)
{
    if (mode == ios_base::out) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++)
                (*pfout)<<ce[i][j][k].ey<<"\n";
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        write_binary_yzplane<celle>(pfout, ce, getey, i, ny, nz);
    } else {
        cerr << "fout_ey_yzplane: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_ez(ofstream* pfout, int n0, int n, ios_base::openmode
        mode)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    if (mode == ios_base::out) {
        for(int i=n0;i<n;i++)
        {
            int k=nz/2;
            for(int j=0;j<ny;j++)
                (*pfout)<<ce[i][j][k].ez<<"\n";
            int j=ny/2;
            for(int k=0;k<nz;k++)
                (*pfout)<<ce[i][j][k].ez<<"\n";
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        write_binary<celle>(pfout, ce, getez, n0, n, ny, nz);
    } else {
        cerr << "fout_ez: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_ez_yzplane(ofstream* pfout, int i, ios_base::openmode
        mode)
{
    if (mode == ios_base::out) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++)
                (*pfout)<<ce[i][j][k].ez<<"\n";
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        write_binary_yzplane<celle>(pfout, ce, getez, i, ny, nz);
    } else {
        cerr << "fout_ez_yzplane: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_bx(ofstream* pfout, int n0, int n, ios_base::openmode
        mode)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    if (mode == ios_base::out) {
        for(int i=n0;i<n;i++)
        {
            int k=nz/2;
            for(int j=0;j<ny;j++)
                (*pfout)<<cb[i][j][k].bx<<"\n";
            int j=ny/2;
            for(int k=0;k<nz;k++)
                (*pfout)<<cb[i][j][k].bx<<"\n";
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        write_binary<cellb>(pfout, cb, getbx, n0, n, ny, nz);
    } else {
        cerr << "fout_bx: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_bx_yzplane(ofstream* pfout, int i, ios_base::openmode
        mode)
{
    if (mode == ios_base::out) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++)
                (*pfout)<<cb[i][j][k].bx<<"\n";
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        write_binary_yzplane<cellb>(pfout, cb, getbx, i, ny, nz);
    } else {
        cerr << "fout_bx_yzplane: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_by(ofstream* pfout, int n0, int n, ios_base::openmode
        mode)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    if (mode == ios_base::out) {
        for(int i=n0;i<n;i++)
        {
            int k=nz/2;
            for(int j=0;j<ny;j++)
                (*pfout)<<cb[i][j][k].by<<"\n";
            int j=ny/2;
            for(int k=0;k<nz;k++)
                (*pfout)<<cb[i][j][k].by<<"\n";
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        write_binary<cellb>(pfout, cb, getby, n0, n, ny, nz);
    } else {
        cerr << "fout_by: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_by_yzplane(ofstream* pfout, int i, ios_base::openmode
        mode)
{
    if (mode == ios_base::out) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++)
                (*pfout)<<cb[i][j][k].by<<"\n";
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        write_binary_yzplane<cellb>(pfout, cb, getby, i, ny, nz);
    } else {
        cerr << "fout_by_yzplane: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_bz(ofstream* pfout, int n0, int n, ios_base::openmode
        mode)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    if (mode == ios_base::out) {
        for(int i=n0;i<n;i++)
        {
            int k=nz/2;
            for(int j=0;j<ny;j++)
                (*pfout)<<cb[i][j][k].bz<<"\n";
            int j=ny/2;
            for(int k=0;k<nz;k++)
                (*pfout)<<cb[i][j][k].bz<<"\n";
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        write_binary<cellb>(pfout, cb, getbz, n0, n, ny, nz);
    } else {
        cerr << "fout_bz: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_bz_yzplane(ofstream* pfout, int i, ios_base::openmode
        mode)
{
    if (mode == ios_base::out) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++)
                (*pfout)<<cb[i][j][k].bz<<"\n";
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        write_binary_yzplane<cellb>(pfout, cb, getbz, i, ny, nz);
    } else {
        cerr << "fout_bz_yzplane: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_w(ofstream* pfout, int n0, int n, bool is_last_sr)
{
    double w;
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    double x,y,z;
    vector3d e;
    vector3d b;
    for(int i=n0;i<n;i++)
    {
        int k=nz/2;
        for(int j=0;j<ny;j++)
        {
            x=i;
            y=j;
            z=k;
            if (is_last_sr==1&&x==n-1) x = x - 1;
            if(y==ny-1) y = y-1;
            e = e_to_particle(x,y,z);
            b = b_to_particle(x,y,z);
            w = e.x*e.x + e.y*e.y + e.z*e.z + b.x*b.x + b.y*b.y + b.z*b.z;
            (*pfout)<<w<<"\n";
        }
        int j=ny/2;
        for(int k=0;k<nz;k++)
        {
            x=i;
            y=j;
            z=k;
            if (is_last_sr==1&&x==n-1) x = x - 1;
            if(z==nz-1) z = z - 1;
            e = e_to_particle(x,y,z);
            b = b_to_particle(x,y,z);
            w = e.x*e.x + e.y*e.y + e.z*e.z + b.x*b.x + b.y*b.y + b.z*b.z;
            (*pfout)<<w<<"\n";
        }
    }
}

void spatial_region::fout_w_yzplane(ofstream* pfout, int i)
{
    double w;
    double x,y,z;
    vector3d e;
    vector3d b;
    for(int j=0;j<ny;j++)
    {
        for(int k=0;k<nz;k++)
        {
            x=i;
            y=j;
            z=k;
            if(x==nx-1) x = x-1;
            if(y==ny-1) y = y-1;
            if(z==nz-1) z = z-1;
            e = e_to_particle(x,y,z);
            b = b_to_particle(x,y,z);
            w = e.x*e.x + e.y*e.y + e.z*e.z + b.x*b.x + b.y*b.y + b.z*b.z;
            (*pfout)<<w<<"\n";
        }
    }
}

void spatial_region::fout_inv(ofstream* pfout, int n0, int n, bool is_last_sr)
{
    double inv;
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    double x,y,z;
    vector3d e;
    vector3d b;
    for(int i=n0;i<n;i++)
    {
        int k=nz/2;
        for(int j=0;j<ny;j++)
        {
            x=i;
            y=j;
            z=k;
            if (is_last_sr==1&&x==n-1) x = x - 1;
            if(y==ny-1) y = y - 1;
            e = e_to_particle(x,y,z);
            b = b_to_particle(x,y,z);
            inv = e.x*e.x + e.y*e.y + e.z*e.z - b.x*b.x - b.y*b.y - b.z*b.z;
            (*pfout)<<inv<<"\n";
        }
        int j=ny/2;
        for(int k=0;k<nz;k++)
        {
            x=i;
            y=j;
            z=k;
            if (is_last_sr==1&&x==n-1) x = x - 1;
            if(z==nz-1) z = z - 1;
            e = e_to_particle(x,y,z);
            b = b_to_particle(x,y,z);
            inv = e.x*e.x + e.y*e.y + e.z*e.z - b.x*b.x - b.y*b.y - b.z*b.z;
            (*pfout)<<inv<<"\n";
        }
    }
}

void spatial_region::fout_inv_yzplane(ofstream* pfout, int i)
{
    double inv;
    double x,y,z;
    vector3d e;
    vector3d b;
    for(int j=0;j<ny;j++)
    {
        for(int k=0;k<nz;k++)
        {
            x=i;
            y=j;
            z=k;
            if(x==nx-1) x = x - 1;
            if(y==ny-1) y = y - 1;
            if(z==nz-1) z = z - 1;
            e = e_to_particle(x,y,z);
            b = b_to_particle(x,y,z);
            inv = e.x*e.x + e.y*e.y + e.z*e.z - b.x*b.x - b.y*b.y - b.z*b.z;
            (*pfout)<<inv<<"\n";
        }
    }
}

void spatial_region::fout_rho(ofstream* pfout_x, ofstream* pfout_y,
        ofstream* pfout_z, int n0, int n, ios_base::openmode mode)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    if (mode == ios_base::out) {
        if (pfout_x != 0) {
            for (int i = n0; i < n; i++) {
                int k = nz / 2;
                for (int j = 0; j < ny; j++) {
                    (*pfout_x) << cj[i][j][k].jx << "\n";
                }
                int j = ny / 2;
                for (int k = 0; k < nz; k++) {
                    (*pfout_x) << cj[i][j][k].jx << "\n";
                }
            }
        }
        if (pfout_y!=0)
        {
            for(int i=n0;i<n;i++)
            {
                int k=nz/2;
                for(int j=0;j<ny;j++)
                {
                    (*pfout_y)<<cj[i][j][k].jy<<"\n";
                }
                int j=ny/2;
                for(int k=0;k<nz;k++)
                {
                    (*pfout_y)<<cj[i][j][k].jy<<"\n";
                }
            }
        }
        if (pfout_z!=0)
        {
            for(int i=n0;i<n;i++)
            {
                int k=nz/2;
                for(int j=0;j<ny;j++)
                {
                    (*pfout_z)<<cj[i][j][k].jz<<"\n";
                }
                int j=ny/2;
                for(int k=0;k<nz;k++)
                {
                    (*pfout_z)<<cj[i][j][k].jz<<"\n";
                }
            }
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        if (pfout_x != 0)
            write_binary<cellj>(pfout_x, cj, getjx, n0, n, ny, nz);
        if (pfout_y != 0)
            write_binary<cellj>(pfout_y, cj, getjy, n0, n, ny, nz);
        if (pfout_z != 0)
            write_binary<cellj>(pfout_z, cj, getjz, n0, n, ny, nz);
    } else {
        cerr << "fout_rho: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_rho_yzplane(ofstream* pfout_x, ofstream* pfout_y,
        ofstream* pfout_z, int i, ios_base::openmode mode)
{
    if (mode == ios_base::out) {
        if (pfout_x != 0) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    (*pfout_x) << cj[i][j][k].jx << "\n";
                }
            }
        }
        if (pfout_y != 0) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    (*pfout_y) << cj[i][j][k].jy << "\n";
                }
            }
        }
        if (pfout_z != 0) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    (*pfout_z) << cj[i][j][k].jz << "\n";
                }
            }
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        if (pfout_x != 0)
            write_binary_yzplane<cellj>(pfout_x, cj, getjx, i, ny, nz);
        if (pfout_y != 0)
            write_binary_yzplane<cellj>(pfout_y, cj, getjy, i, ny, nz);
        if (pfout_z != 0)
            write_binary_yzplane<cellj>(pfout_z, cj, getjz, i, ny, nz);
    } else {
        cerr << "fout_rho_yzplane: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_irho(int ion_type, ofstream* pfout, int n0, int n,
        ios_base::openmode mode)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    // ion_type - index for icmr array
    if (mode == ios_base::out) {
        for(int i=n0;i<n;i++)
        {
            int k=nz/2;
            for(int j=0;j<ny;j++)
            {
                (*pfout)<<irho[ion_type][i][j][k]<<"\n";
            }
            int j=ny/2;
            for(int k=0;k<nz;k++)
            {
                (*pfout)<<irho[ion_type][i][j][k]<<"\n";
            }
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        write_binary<double>(pfout, irho[ion_type], getdouble, n0, n, ny, nz);
    } else {
        cerr << "fout_irho: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_irho_yzplane(int ion_type, ofstream* pfout, int i,
        ios_base::openmode mode)
{
    if (mode == ios_base::out) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                (*pfout)<<irho[ion_type][i][j][k]<<"\n";
            }
        }
    } else if (mode == ios_base::out | ios_base::binary) {
        write_binary_yzplane<double>(pfout, irho[ion_type], getdouble, i, ny,
                nz);
    } else {
        cerr << "fout_irho_yzplane: ERROR: wrong mode" << endl;
    }
}

void spatial_region::fout_tracks(double a, int nm) {
    for (int i=0;i<nx;i++) {
        for (int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                particle* current;
                current = cp[i][j][k].pl.head;
                while (current!=0) {
                    if (current->trn!=0 && current->x>=nm && current->x<nx-nm) {
                        std::string file_name;
                        char file_num_pchar[20];
                        file_name = spatial_region::data_folder+"/track_";
                        sprintf(file_num_pchar,"%g",current->cmr);
                        file_name = file_name + file_num_pchar;
                        file_name = file_name + "_";
                        sprintf(file_num_pchar,"%d",current->trn);
                        file_name = file_name + file_num_pchar;
                        ofstream fout(file_name.c_str(),ios::app);
                        fout<<current->q<<'\n'<<current->x*dx/2/PI+a<<'\n'<<current->y*dy/2/PI<<'\n'<<current->z*dz/2/PI<<'\n'<<current->ux<<'\n'<<current->uy<<'\n'<<current->uz<<'\n'<<current->g<<'\n'<<current->chi<<'\n';
                    }
                    current = current->next;
                }
            }
        }
    }
}

void spatial_region::scale_j(const double scale) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                cj[i][j][k].jx *= scale;
                cj[i][j][k].jy *= scale;
                cj[i][j][k].jz *= scale;
            }
        }
    }
}
