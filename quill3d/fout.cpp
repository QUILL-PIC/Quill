#include <fstream>
#include "main.h"
#include "containers.h"
#include "openpmd_output.h"
#include <functional>

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

double w_function(double ex, double ey, double ez, double bx, double by, double bz) {
    return ex*ex + ey*ey + ez*ez + bx*bx + by*by + bz*bz;
};

double inv_function(double ex, double ey, double ez, double bx, double by, double bz) {
    return ex*ex + ey*ey + ez*ez - bx*bx - by*by - bz*bz;
};

void write_4d_array(double * p_first_element, int vector_size, int component, int nx0, int ny0, int nz0, 
    H5::DataSet & xy_dataset, H5::DataSet & xz_dataset, int left, int right, int output_dataset_shift) {
    const hsize_t nx = static_cast<hsize_t>(nx0);
    const hsize_t ny = static_cast<hsize_t>(ny0);
    const hsize_t nz = static_cast<hsize_t>(nz0);

    const hsize_t dims[4] {nx, ny, nz, static_cast<hsize_t>(vector_size)};

    H5::DataSpace memory_dataspace(4, dims);

    const hsize_t x_start = static_cast<hsize_t>(left);
    const hsize_t x_count = static_cast<hsize_t>(right - left);
    
    const hsize_t count_xy[4] {x_count, ny, 1, 1};
    const hsize_t start_xy[4] {x_start, 0, nz/2, static_cast<hsize_t>(component)};
    memory_dataspace.selectHyperslab(H5S_SELECT_SET, count_xy, start_xy);

    const hsize_t count_xy_out[2] {x_count, ny};
    const hsize_t start_xy_out[2] {static_cast<hsize_t>(output_dataset_shift), 0};
    auto write_xy_dataspace = xy_dataset.getSpace();
    write_xy_dataspace.selectHyperslab(H5S_SELECT_SET, count_xy_out, start_xy_out);
    
    xy_dataset.write(p_first_element, H5::PredType::NATIVE_DOUBLE, memory_dataspace, write_xy_dataspace);

    const hsize_t count_xz[4] {x_count, 1, nz, 1};
    const hsize_t start_xz[4] {x_start, ny/2, 0, static_cast<hsize_t>(component)};
    memory_dataspace.selectHyperslab(H5S_SELECT_SET, count_xz, start_xz);

    const hsize_t count_xz_out[2] {x_count, nz};
    const hsize_t start_xz_out[2] {static_cast<hsize_t>(output_dataset_shift), 0};
    auto write_xz_dataspace = xz_dataset.getSpace();
    write_xz_dataspace.selectHyperslab(H5S_SELECT_SET, count_xz_out, start_xz_out);

    xz_dataset.write(p_first_element, H5::PredType::NATIVE_DOUBLE, memory_dataspace, write_xz_dataspace);
}

void write_vector_field(double * p_first_element, int component, int nx0, int ny0, int nz0, 
    H5::DataSet & xy_dataset, H5::DataSet & xz_dataset, int left, int right, int output_dataset_shift) {
    write_4d_array(p_first_element, 3, component, nx0, ny0, nz0, xy_dataset, xz_dataset, left, right, output_dataset_shift);
}

void write_scalar_field(double * p_first_element, int nx0, int ny0, int nz0, H5::DataSet & xy_dataset, 
                        H5::DataSet & xz_dataset, int left, int right, int output_dataset_shift) {
    write_4d_array(p_first_element, 1, 0, nx0, ny0, nz0, xy_dataset, xz_dataset, left, right, output_dataset_shift);
}

void write_4d_array_yz(double * p_first_element, int vector_size, int component, int nx0, int ny0, int nz0, 
                       H5::DataSet & yz_dataset, int position) {
    const hsize_t nx = static_cast<hsize_t>(nx0);
    const hsize_t ny = static_cast<hsize_t>(ny0);
    const hsize_t nz = static_cast<hsize_t>(nz0);

    const hsize_t dims[4] {nx, ny, nz, static_cast<hsize_t>(vector_size)};

    H5::DataSpace memory_dataspace(4, dims);

    const hsize_t count[4] {1, ny, nz, 1};
    const hsize_t start[4] {static_cast<hsize_t>(position), 0, 0, static_cast<hsize_t>(component)};
    memory_dataspace.selectHyperslab(H5S_SELECT_SET, count, start);

    auto write_dataspace = yz_dataset.getSpace();
    
    yz_dataset.write(p_first_element, H5::PredType::NATIVE_DOUBLE, memory_dataspace, write_dataspace);
}

void write_vector_field_yz(double * p_first_element, int component, int nx0, int ny0, int nz0, 
                           H5::DataSet & yz_dataset, int position) {
    write_4d_array_yz(p_first_element, 3, component, nx0, ny0, nz0, yz_dataset, position);
}

void write_scalar_field_yz(double * p_first_element, int nx0, int ny0, int nz0, H5::DataSet & yz_dataset, int position) {
    write_4d_array_yz(p_first_element, 1, 0, nx0, ny0, nz0, yz_dataset, position);
}

void spatial_region::fout_ex_xy_xz(H5::DataSet & xy_dataset, H5::DataSet & xz_dataset, int left, int right, int dataset_shift) {
    write_vector_field(&(ce[0][0][0].ex), 0, nx, ny, nz, xy_dataset, xz_dataset, left, right, dataset_shift);
}

void spatial_region::fout_ey_xy_xz(H5::DataSet & xy_dataset, H5::DataSet & xz_dataset, int left, int right, int dataset_shift) {
    write_vector_field(&(ce[0][0][0].ex), 1, nx, ny, nz, xy_dataset, xz_dataset, left, right, dataset_shift);
}

void spatial_region::fout_ez_xy_xz(H5::DataSet & xy_dataset, H5::DataSet & xz_dataset, int left, int right, int dataset_shift) {
    write_vector_field(&(ce[0][0][0].ex), 2, nx, ny, nz, xy_dataset, xz_dataset, left, right, dataset_shift);
}

void spatial_region::fout_bx_xy_xz(H5::DataSet & xy_dataset, H5::DataSet & xz_dataset, int left, int right, int dataset_shift) {
    write_vector_field(&(cb[0][0][0].bx), 0, nx, ny, nz, xy_dataset, xz_dataset, left, right, dataset_shift);
}

void spatial_region::fout_by_xy_xz(H5::DataSet & xy_dataset, H5::DataSet & xz_dataset, int left, int right, int dataset_shift) {
    write_vector_field(&(cb[0][0][0].bx), 1, nx, ny, nz, xy_dataset, xz_dataset, left, right, dataset_shift);
}

void spatial_region::fout_bz_xy_xz(H5::DataSet & xy_dataset, H5::DataSet & xz_dataset, int left, int right, int dataset_shift) {
    write_vector_field(&(cb[0][0][0].bx), 2, nx, ny, nz, xy_dataset, xz_dataset, left, right, dataset_shift);
}

void spatial_region::fout_jx_xy_xz(H5::DataSet & xy_dataset, H5::DataSet & xz_dataset, int left, int right, int dataset_shift) {
    write_vector_field(&(cj[0][0][0].jx), 0, nx, ny, nz, xy_dataset, xz_dataset, left, right, dataset_shift);
}

void spatial_region::fout_jy_xy_xz(H5::DataSet & xy_dataset, H5::DataSet & xz_dataset, int left, int right, int dataset_shift) {
    write_vector_field(&(cj[0][0][0].jx), 1, nx, ny, nz, xy_dataset, xz_dataset, left, right, dataset_shift);
}

void spatial_region::fout_jz_xy_xz(H5::DataSet & xy_dataset, H5::DataSet & xz_dataset, int left, int right, int dataset_shift) {
    write_vector_field(&(cj[0][0][0].jx), 2, nx, ny, nz, xy_dataset, xz_dataset, left, right, dataset_shift);
}

void spatial_region::fout_irho_xy_xz(int ion_type, H5::DataSet & xy_dataset, H5::DataSet & xz_dataset, int left,
                                     int right, int dataset_shift) {
    write_scalar_field(&(irho[ion_type][0][0][0]), nx, ny, nz, xy_dataset, xz_dataset, left, right, dataset_shift);
}

void spatial_region::fout_field_function(H5::DataSet & xy_dataset, H5::DataSet & xz_dataset, int left, int right, 
                                         int dataset_shift, bool is_last_sr, std::function<double(double, double, double, double, double, double)> func)
{
    const hsize_t x_count = static_cast<hsize_t>(right - left);

    std::vector<double> result(x_count * ny);

    double x,y,z;
    vector3d e;
    vector3d b;
    int index = 0;
    for(int i=left; i<right; i++)
    {
        int k=nz/2;
        for(int j=0;j<ny;j++)
        {
            x=i;
            y=j;
            z=k;
            if (is_last_sr==1&&x==right-1) x = x - 1;
            if(y==ny-1) y = y-1;
            e = e_to_particle(x,y,z);
            b = b_to_particle(x,y,z);
            result[index] = func(e.x, e.y, e.z, b.x, b.y, b.z);
            index++;
        }
    }
    const hsize_t dims_xy[1] {(result.size())};
    H5::DataSpace memory_dataspace_xy(1, dims_xy);

    const hsize_t count_xy_out[2] {x_count, static_cast<hsize_t>(ny)};
    const hsize_t start_xy_out[2] {static_cast<hsize_t>(dataset_shift), 0};
    auto write_dataspace_xy = xy_dataset.getSpace();
    write_dataspace_xy.selectHyperslab(H5S_SELECT_SET, count_xy_out, start_xy_out);

    xy_dataset.write(&(result[0]), H5::PredType::NATIVE_DOUBLE, memory_dataspace_xy, write_dataspace_xy);

    result = std::vector<double>(x_count * nz);
    index = 0;
    for(int i=left; i<right; i++)
    {
        int j=ny/2;
        for(int k=0;k<nz;k++)
        {
            x=i;
            y=j;
            z=k;
            if (is_last_sr==1&&x==right-1) x = x - 1;
            if(z==nz-1) z = z - 1;
            e = e_to_particle(x,y,z);
            b = b_to_particle(x,y,z);
            result[index] = func(e.x, e.y, e.z, b.x, b.y, b.z);
            index++;
        }
    }

    const hsize_t dims_xz[1] {(result.size())};
    H5::DataSpace memory_dataspace_xz(1, dims_xz);

    const hsize_t count_xz_out[2] {x_count, static_cast<hsize_t>(nz)};
    const hsize_t start_xz_out[2] {static_cast<hsize_t>(dataset_shift), 0};
    auto write_dataspace_xz = xz_dataset.getSpace();
    write_dataspace_xz.selectHyperslab(H5S_SELECT_SET, count_xz_out, start_xz_out);

    xz_dataset.write(&(result[0]), H5::PredType::NATIVE_DOUBLE, memory_dataspace_xz, write_dataspace_xz);
}

void spatial_region::fout_w_xy_xz(H5::DataSet & xy_dataset, H5::DataSet & xz_dataset, int left, int right, int dataset_shift, bool is_last_sr) {
    fout_field_function(xy_dataset, xz_dataset, left, right, dataset_shift, is_last_sr, w_function);
}

void spatial_region::fout_inv_xy_xz(H5::DataSet & xy_dataset, H5::DataSet & xz_dataset, int left, int right, int dataset_shift, bool is_last_sr) {
    fout_field_function(xy_dataset, xz_dataset, left, right, dataset_shift, is_last_sr, inv_function);
}

void spatial_region::fout_ex_yz(H5::DataSet & yz_dataset, int position) {
    write_vector_field_yz(&(ce[0][0][0].ex), 0, nx, ny, nz, yz_dataset, position);
}

void spatial_region::fout_ey_yz(H5::DataSet & yz_dataset, int position) {
    write_vector_field_yz(&(ce[0][0][0].ex), 1, nx, ny, nz, yz_dataset, position);
}

void spatial_region::fout_ez_yz(H5::DataSet & yz_dataset, int position) {
    write_vector_field_yz(&(ce[0][0][0].ex), 2, nx, ny, nz, yz_dataset, position);
}

void spatial_region::fout_bx_yz(H5::DataSet & yz_dataset, int position) {
    write_vector_field_yz(&(cb[0][0][0].bx), 0, nx, ny, nz, yz_dataset, position);
}

void spatial_region::fout_by_yz(H5::DataSet & yz_dataset, int position) {
    write_vector_field_yz(&(cb[0][0][0].bx), 1, nx, ny, nz, yz_dataset, position);
}

void spatial_region::fout_bz_yz(H5::DataSet & yz_dataset, int position) {
    write_vector_field_yz(&(cb[0][0][0].bx), 2, nx, ny, nz, yz_dataset, position);
}

void spatial_region::fout_jx_yz(H5::DataSet & yz_dataset, int position) {
    write_vector_field_yz(&(cj[0][0][0].jx), 0, nx, ny, nz, yz_dataset, position);
}

void spatial_region::fout_jy_yz(H5::DataSet & yz_dataset, int position) {
    write_vector_field_yz(&(cj[0][0][0].jx), 1, nx, ny, nz, yz_dataset, position);
}

void spatial_region::fout_jz_yz(H5::DataSet & yz_dataset, int position) {
    write_vector_field_yz(&(cj[0][0][0].jx), 2, nx, ny, nz, yz_dataset, position);
}

void spatial_region::fout_irho_yz(int ion_type, H5::DataSet & yz_dataset, int position) {
    write_scalar_field_yz(&(irho[ion_type][0][0][0]), nx, ny, nz, yz_dataset, position);
}

void spatial_region::fout_field_function_yz(H5::DataSet & yz_dataset, int position,
                                            std::function<double(double, double, double, double, double, double)> func)
{
    std::vector<double> result(ny * nz);
    double x,y,z;
    vector3d e;
    vector3d b;
    int index = 0;
    for(int j=0;j<ny;j++)
    {
        for(int k=0;k<nz;k++)
        {
            x=position;
            y=j;
            z=k;
            if(x==nx-1) x = x-1;
            if(y==ny-1) y = y-1;
            if(z==nz-1) z = z-1;
            e = e_to_particle(x,y,z);
            b = b_to_particle(x,y,z);
            result[index] = func(e.x, e.y, e.z, b.x, b.y, b.z);
            index++;
        }
    }

    const hsize_t dims[1] {(result.size())};
    H5::DataSpace memory_dataspace(1, dims);

    auto write_dataspace = yz_dataset.getSpace();

    yz_dataset.write(&(result[0]), H5::PredType::NATIVE_DOUBLE, memory_dataspace, write_dataspace);
}

void spatial_region::fout_w_yz(H5::DataSet & yz_dataset, int position)
{
    fout_field_function_yz(yz_dataset, position, w_function);
}

void spatial_region::fout_inv_yz(H5::DataSet & yz_dataset, int position)
{
    fout_field_function_yz(yz_dataset, position, inv_function);
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
