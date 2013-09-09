#include <fstream>
#include "main.h"

void spatial_region::compute_rho()
{
    spatial_region::plist::particle* current;
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
    N_qp = 0;
    plist::particle* current;
    for (int i=n1; i<nx-n2; i++)
    {
	for (int j=0; j<ny; j++)
	{
	    for (int k=0; k<nz; k++)
	    {
		current = cp[i][j][k].pl.head;
		while (current!=0)
		{
		    N_qp += 1;
		    if (current->cmr==-1)
			N_e -= current->q;
		    else if (current->cmr==1)
			N_p += current->q;
		    else if (current->cmr==0)
			N_ph += current->q;
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
    plist::particle* current;
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

void spatial_region::fout_ex(ofstream* pfout, int n0, int n)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    for(int i=n0;i<n;i++)
    {
        int k=nz/2;
	for(int j=0;j<ny;j++)
	    (*pfout)<<ce[i][j][k].ex<<"\n";
	int j=ny/2;
	for(int k=0;k<nz;k++)
	    (*pfout)<<ce[i][j][k].ex<<"\n";
    }
}

void spatial_region::fout_ex_yzplane(ofstream* pfout, int i)
{
    for(int j=0;j<ny;j++)
    {
	for(int k=0;k<nz;k++)
	    (*pfout)<<ce[i][j][k].ex<<"\n";
    }
}

void spatial_region::fout_ey(ofstream* pfout, int n0, int n)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    for(int i=n0;i<n;i++)
    {
        int k=nz/2;
	for(int j=0;j<ny;j++)
	    (*pfout)<<ce[i][j][k].ey<<"\n";
	int j=ny/2;
	for(int k=0;k<nz;k++)
	    (*pfout)<<ce[i][j][k].ey<<"\n";
    }
}

void spatial_region::fout_ey_yzplane(ofstream* pfout, int i)
{
    for(int j=0;j<ny;j++)
    {
	for(int k=0;k<nz;k++)
	    (*pfout)<<ce[i][j][k].ey<<"\n";
    }
}

void spatial_region::fout_ez(ofstream* pfout, int n0, int n)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    for(int i=n0;i<n;i++)
    {
        int k=nz/2;
	for(int j=0;j<ny;j++)
	    (*pfout)<<ce[i][j][k].ez<<"\n";
	int j=ny/2;
	for(int k=0;k<nz;k++)
	    (*pfout)<<ce[i][j][k].ez<<"\n";
    }
}

void spatial_region::fout_ez_yzplane(ofstream* pfout, int i)
{
    for(int j=0;j<ny;j++)
    {
	for(int k=0;k<nz;k++)
	    (*pfout)<<ce[i][j][k].ez<<"\n";
    }
}

void spatial_region::fout_bx(ofstream* pfout, int n0, int n)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    for(int i=n0;i<n;i++)
    {
        int k=nz/2;
	for(int j=0;j<ny;j++)
	    (*pfout)<<cb[i][j][k].bx<<"\n";
	int j=ny/2;
	for(int k=0;k<nz;k++)
	    (*pfout)<<cb[i][j][k].bx<<"\n";
    }
}

void spatial_region::fout_bx_yzplane(ofstream* pfout, int i)
{
    for(int j=0;j<ny;j++)
    {
	for(int k=0;k<nz;k++)
	    (*pfout)<<cb[i][j][k].bx<<"\n";
    }
}

void spatial_region::fout_by(ofstream* pfout, int n0, int n)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    for(int i=n0;i<n;i++)
    {
        int k=nz/2;
	for(int j=0;j<ny;j++)
	    (*pfout)<<cb[i][j][k].by<<"\n";
	int j=ny/2;
	for(int k=0;k<nz;k++)
	    (*pfout)<<cb[i][j][k].by<<"\n";
    }
}

void spatial_region::fout_by_yzplane(ofstream* pfout, int i)
{
    for(int j=0;j<ny;j++)
    {
	for(int k=0;k<nz;k++)
	    (*pfout)<<cb[i][j][k].by<<"\n";
    }
}

void spatial_region::fout_bz(ofstream* pfout, int n0, int n)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    for(int i=n0;i<n;i++)
    {
        int k=nz/2;
	for(int j=0;j<ny;j++)
	    (*pfout)<<cb[i][j][k].bz<<"\n";
	int j=ny/2;
	for(int k=0;k<nz;k++)
	    (*pfout)<<cb[i][j][k].bz<<"\n";
    }
}

void spatial_region::fout_bz_yzplane(ofstream* pfout, int i)
{
    for(int j=0;j<ny;j++)
    {
	for(int k=0;k<nz;k++)
	    (*pfout)<<cb[i][j][k].bz<<"\n";
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

void spatial_region::fout_rho(ofstream* pfout, ofstream* pfout_p, ofstream* pfout_ph, int n0, int n)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    for(int i=n0;i<n;i++)
    {
        int k=nz/2;
        for(int j=0;j<ny;j++)
        {
            (*pfout)<<cj[i][j][k].jx<<"\n";
        }
        int j=ny/2;
        for(int k=0;k<nz;k++)
        {
            (*pfout)<<cj[i][j][k].jx<<"\n";
        }
    }
    if (pfout_p!=0)
    {
	for(int i=n0;i<n;i++)
	{
	    int k=nz/2;
	    for(int j=0;j<ny;j++)
	    {
		(*pfout_p)<<cj[i][j][k].jy<<"\n";
	    }
	    int j=ny/2;
	    for(int k=0;k<nz;k++)
	    {
		(*pfout_p)<<cj[i][j][k].jy<<"\n";
	    }
	}
    }
    if (pfout_ph!=0)
    {
	for(int i=n0;i<n;i++)
	{
	    int k=nz/2;
	    for(int j=0;j<ny;j++)
	    {
		(*pfout_ph)<<cj[i][j][k].jz<<"\n";
	    }
	    int j=ny/2;
	    for(int k=0;k<nz;k++)
	    {
		(*pfout_ph)<<cj[i][j][k].jz<<"\n";
	    }
	}
    }
}

void spatial_region::fout_rho_yzplane(ofstream* pfout, ofstream* pfout_p, ofstream* pfout_ph, int i)
{
    for(int j=0;j<ny;j++)
    {
        for(int k=0;k<nz;k++)
        {
            (*pfout)<<cj[i][j][k].jx<<"\n";
        }
    }
    if (pfout_p!=0)
    {
	for(int j=0;j<ny;j++)
	{
	    for(int k=0;k<nz;k++)
	    {
		(*pfout_p)<<cj[i][j][k].jy<<"\n";
	    }
	}
    }
    if (pfout_ph!=0)
    {
	for(int j=0;j<ny;j++)
	{
	    for(int k=0;k<nz;k++)
	    {
		(*pfout_ph)<<cj[i][j][k].jz<<"\n";
	    }
	}
    }
}

void spatial_region::fout_irho(int ion_type, ofstream* pfout, int n0, int n)
{
    /* n-n0 - длина массива по x (для данного spatial_region'а) для
     * вывода сечения плоскостью xy и плоскостью xz */
    // ion_type - index for icmr array
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
}

void spatial_region::fout_irho_yzplane(int ion_type,ofstream* pfout, int i)
{
    for(int j=0;j<ny;j++)
    {
        for(int k=0;k<nz;k++)
        {
            (*pfout)<<irho[ion_type][i][j][k]<<"\n";
        }
    }
}

void spatial_region::fout_tracks(double a, int nm) {
    for (int i=0;i<nx;i++) {
	for (int j=0;j<ny;j++) {
	    for(int k=0;k<nz;k++) {
		plist::particle* current;
		current = cp[i][j][k].pl.head;
		while (current!=0) {
		    if (current->trn!=0 && current->x>=nm && current->x<nx-nm) {
			std::string file_name;
			char file_num_pchar[20];
			file_name = "results/track_";
			sprintf(file_num_pchar,"%g",current->cmr);
			file_name = file_name + file_num_pchar;
			file_name = file_name + "_";
			sprintf(file_num_pchar,"%d",current->trn);
			file_name = file_name + file_num_pchar;
			ofstream fout(file_name,ios::app);
			fout<<current->q<<'\n'<<current->x*dx/2/PI+a<<'\n'<<current->y*dy/2/PI<<'\n'<<current->z*dz/2/PI<<'\n'<<current->ux<<'\n'<<current->uy<<'\n'<<current->uz<<'\n'<<current->g<<'\n'<<current->chi<<'\n';
		    }
		    current = current->next;
		}
	    }
	}
    }
}
