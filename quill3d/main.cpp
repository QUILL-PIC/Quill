#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <sys/times.h>
#include <unistd.h>
#include <memory>
#include <cassert>
#include "mpi.h"
#include "main.h"
#include "maxwell.h"
#include "containers.h"
#include "balancing.h"

using namespace std;

// global variables ------------
double dx,dy,dz,dt;
double xlength,ylength,zlength;
int nx_global, ny_global, nz_global;
double a0y,a0z;
double xsigma,ysigma,zsigma;
double x0;
double xtarget, ytarget, ztarget;
double lambda;
double k0,ne;
double deps,deps_p,deps_ph,deps_i;
int neps,enthp,neps_p,enthp_p,neps_ph,enthp_ph,neps_i,enthp_i;
int xnpic,ynpic,znpic;
double x0fout;
double Nb,epsb,xb,rb,x0b,y0b,phib;
time_t inaccurate_time;
tms tms_struct;
clock_t up_time;
clock_t start_time;
double seconds;
clock_t main_thread_time;
double file_name_accuracy;
bool mwseed;
moving_window mwindow;
double mwtolerance; // the relative level of E^2 + B^2 causing the window to move in the AUTO mode
std::string mwseed_ions;
double mw_mcr;
double crpc;
double* ppd;
double phase,phi,phi_rotate;
double shenergy, shphase; // second harmonic relative energy and phase
ddi* p_last_ddi; // ddi включает t_end, output_period и f - счётчик для вывода данных в файлы
ddi* p_current_ddi;
double t_add_mw;
film* p_last_film;
double x00,y00,z00;
double sigma0,sigma; // sigma0 - радиус пучка в фокальной плоскости
double external_bz;
bool b_sign; // b_sign = 1 соответствует знаку '+', 0 - знаку '-'
int n_ion_populations;
double* icmr;
int n_tracks;
std::string particles_to_track;
bool tr_init;
double tr_start,xtr1,ytr1,ztr1,xtr2,ytr2,ztr2;
double mwspeed,nelflow,vlflow,mcrlflow,Tlflow,nerflow,vrflow,mcrrflow,Trflow;
int i_particle, i_particle_p, i_particle_ph, i_particle_i;  // for "writing down every i-th electron, positron, etc."
std::string e_components_for_output;
std::string b_components_for_output;
std::string j_components_for_output;
std::string f_envelope;
int sscos; // 0 - cos, 1 - super-super cos, 2 - pearl
std::string beam;
std::string beam_particles;
bool write_p;
bool write_ph;
bool write_jx = false, write_jy = false, write_jz = false;
std::string pmerging;
bool pmerging_now;
std::string lp_reflection,f_reflection;
std::string ions;
std::string data_folder;
bool catching_enabled;
bool dump_photons;
bool verbose_logging = true;
bool qed_enabled = true;
ios_base::openmode output_mode;
int init();

std::vector<double> ne_profile_x_coords;
std::vector<double> ne_profile_x_values;
std::vector<double> ne_profile_r_coords;
std::vector<double> ne_profile_r_values;

maxwell_solver_enum solver;
pusher_enum pusher;

bool balancing_enabled;
int balancing_every;
double balancing_threshold;
double balancing_particle_weight;

//------------------------------


int l; // номер шага по времени

unique_ptr<spatial_region> psr; // пространственный регион данного процесса

int n_sr; // число слоёв и mpi-процессов
int mpi_rank; // mpi-ранг текущего процесса
vector<int> nx_sr; // число ячеек (по x) в слое
vector<int> x0_sr; // индексы левых границ слоев
int nx_ich; /* число ячеек для приграничной области слоёв, данные в
               которой замещаются данными соседнего слоя; должно
               быть чётным */
int nm; /* используется при обмене данными между слоями и
           для аккуратного подсчёта спектров, энергии и
           числа частиц */

int nmw = 1;

MPI_Datatype MPI_PARTICLE; // тип для передачи частиц через MPI

// номер слоя, отвечающий индексу i
inline int get_sr_for_index(int i)
{
    //   return lower_bound(x0_sr.begin(), x0_sr.end(), i) - x0_sr.begin();
   for (int j=n_sr-1; j>0; j--) {
       if (i >= x0_sr[j])
           return j;
   }
   return 0;
}

// номер слоя, отвечающий координате x
inline int get_sr_for_x(double x)
{
   return get_sr_for_index(int(x / dx));
}

// индекс, соответствующий координате x в заданном слое 
inline int get_xindex_in_sr(double x, int sr_i) {
    int ix = int(x / dx);
    return ix - x0_sr[sr_i];
}

// индекс, соответствующий координате x в слое, где она находится 
inline int get_xindex_in_sr(double x) {
    return get_xindex_in_sr(x, get_sr_for_x(x));
}

void write_mwcoordinate(ofstream & fout_mwcoordinate) {
    if (mpi_rank == 0) {
        fout_mwcoordinate << ((nmw - 1) * dx / 2 / PI) << endl;
    }
}

void write_N(ofstream& fout_N)
{
    double N_e,N_p,N_ph;
    MPI_Reduce(&(psr->N_e), &N_e, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(psr->N_p), &N_p, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(psr->N_ph), &N_ph, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        fout_N<<N_e<<'\t'<<N_p<<'\t'<<N_ph<<endl;
    }
}

void write_energy(ofstream& fout_energy)
{
    double energy_f, energy_e, energy_p, energy_ph;
    auto ienergy = unique_ptr<double[]>(new double[n_ion_populations]);
    MPI_Reduce(&(psr->energy_f), &energy_f, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(psr->energy_e), &energy_e, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(psr->energy_p), &energy_p, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(psr->energy_ph), &energy_ph, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(psr->ienergy, ienergy.get(), n_ion_populations, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        fout_energy<<energy_f<<'\t'<<energy_e<<'\t'<<energy_p<<'\t'<<energy_ph;
            for (int n=0;n<n_ion_populations;n++)
                fout_energy<<'\t'<<ienergy[n];
            fout_energy<<endl;
    }
}

void write_deleted_particles(bool write_p, bool write_ph)
{
    ios_base::openmode non_binary_mode;
    if (mpi_rank == 0) {
        non_binary_mode = ios_base::out;
    } else {
        non_binary_mode = ios_base::out | ios_base::app;
    }

    std::string file_name;
    char file_num_pchar[20];
    sprintf(file_num_pchar,"%g",
        int([](ddi* a) 
        {
            double b=a->f*a->output_period; 
            if(a->prev!=0) 
                b+=(a->prev)->t_end; 
            return b;
        } (p_current_ddi)/2/PI*file_name_accuracy
        )/file_name_accuracy);
    
    double* spectrum_ph = new double[neps_ph];
    for(int i=0; i<neps_ph; i++) spectrum_ph[i] = 0;
    int i_eps = 0;

    for(int n=0; n<n_sr; n++)
    {
        if (mpi_rank == n) {
            file_name = data_folder + "/deleted" + file_num_pchar;
            ofstream fout_deleted_e(file_name.c_str(), non_binary_mode);

            file_name = data_folder + "/deleted_p" + file_num_pchar;
            ofstream fout_deleted_p(file_name.c_str(), non_binary_mode);

            file_name = data_folder + "/deleted_ph" + file_num_pchar;
            ofstream fout_deleted_ph(file_name.c_str(), non_binary_mode);

            ofstream* fout_deleted_i = new ofstream[n_ion_populations];
            for (int m=0; m<n_ion_populations; ++m)
            {
                char s_cmr[20];
                sprintf(s_cmr,"%g",icmr[m]);
                file_name = data_folder+"/deleted_";
                file_name += s_cmr;
                file_name += "_";
                file_name += file_num_pchar;
                fout_deleted_i[m].open(file_name.c_str(), non_binary_mode);
            }

            vector<spatial_region::deleted_particle>& del_particles = psr->deleted_particles;

            for (auto it = del_particles.begin(); it != del_particles.end(); ++it)
            {
                if ((*it).cmr == -1)
                {
                    i_particle++;
                    if(i_particle > enthp)
                    {
                        i_particle = 0;
                        fout_deleted_e << (*it).q << endl;
                        fout_deleted_e << dx*((*it).x + x0_sr[n])/(2*PI) << endl;
                        fout_deleted_e << dy*((*it).y)/(2*PI) << endl << dz*((*it).z)/(2*PI) << endl;
                        fout_deleted_e << (*it).ux << endl << (*it).uy << endl <<(*it).uz << endl;
                        fout_deleted_e << (*it).g << endl << (*it).chi << endl;
                    }
                }
                else if (write_p && (*it).cmr == 1)
                {
                    i_particle_p++;
                    if(i_particle_p > enthp_p)
                    {
                        i_particle_ph = 0;
                        fout_deleted_p << (*it).q << endl;
                        fout_deleted_p << dx*((*it).x + x0_sr[n])/(2*PI) << endl;
                        fout_deleted_p << dy*((*it).y)/(2*PI) << endl << dz*((*it).z)/(2*PI) << endl;
                        fout_deleted_p << (*it).ux << endl << (*it).uy << endl <<(*it).uz << endl;
                        fout_deleted_p << (*it).g << endl << (*it).chi << endl;
                    }
                }
                else if (write_ph && (*it).cmr == 0)
                {
                    i_eps = (*it).g*0.511/deps_ph;
                    if((i_eps > -1) && (i_eps < neps_ph)) {
                        spectrum_ph[i_eps] = spectrum_ph[i_eps] + (*it).q; // q>0 for photons
                    }
                    i_particle_ph++;
                    if(i_particle_ph > enthp_ph)
                    {
                        i_particle_ph = 0;
                        fout_deleted_ph << (*it).q << endl;
                        fout_deleted_ph << dx*((*it).x + x0_sr[n])/(2*PI) << endl;
                        fout_deleted_ph << dy*((*it).y)/(2*PI) << endl << dz*((*it).z)/(2*PI) << endl;
                        fout_deleted_ph << (*it).ux << endl << (*it).uy << endl <<(*it).uz << endl;
                        fout_deleted_ph << (*it).g << endl << (*it).chi << endl;
                    }
                }
                else
                {
                    for (int j=0; j<n_ion_populations; ++j)
                    {
                        if ((*it).cmr == icmr[j])
                        {
                            i_particle_i++;
                            if(i_particle_i > enthp_i)
                            {
                                i_particle_i = 0;
                                fout_deleted_i[j] << (*it).q << endl;
                                fout_deleted_i[j] << dx*((*it).x + x0_sr[n])/(2*PI) << endl;
                                fout_deleted_i[j] << dy*((*it).y)/(2*PI) << endl << dz*((*it).z)/(2*PI) << endl;
                                fout_deleted_i[j] << (*it).ux << endl << (*it).uy << endl <<(*it).uz << endl;
                                fout_deleted_i[j] << (*it).g << endl << (*it).chi << endl;
                            }
                            break;
                        }
                    }
                }
            }
            del_particles.clear();
            fout_deleted_e.close();
            fout_deleted_p.close();
            fout_deleted_ph.close();
            for (int m=0; m<n_ion_populations; ++m)
            {
                fout_deleted_i[m].close();
            }
            delete[] fout_deleted_i;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (write_ph) {
        //gathering photon spectra on process rank 0
        if (mpi_rank == 0) {
            MPI_Reduce(MPI_IN_PLACE, spectrum_ph, neps_ph, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        } else {
            MPI_Reduce(spectrum_ph, spectrum_ph, neps_ph, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }

        //spectrum output
        if (mpi_rank == 0) {
            file_name = data_folder + "/spectrum_deleted_ph" + file_num_pchar;
            ofstream fout_spectrum_ph(file_name.c_str());

            double spectrum_norm_ph = 1.11485e13 * lambda*dx*dy*dz/(8*PI*PI*PI)/deps_ph;
            for(int i=0; i<neps_ph; i++)
            {
                spectrum_ph[i] = spectrum_ph[i] * spectrum_norm_ph;
                fout_spectrum_ph << spectrum_ph[i] << '\n';
            }

            fout_spectrum_ph.close();
        }
    }

    delete[] spectrum_ph;
}

void write_energy_deleted(ofstream& fout_energy_deleted)
{
    double energy_e_deleted, energy_p_deleted, energy_ph_deleted;
    auto ienergy_deleted = unique_ptr<double[]>(new double[n_ion_populations]);
    MPI_Reduce(&(psr->energy_e_deleted), &energy_e_deleted, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(psr->energy_p_deleted), &energy_p_deleted, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(psr->energy_ph_deleted), &energy_ph_deleted, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(psr->ienergy_deleted, ienergy_deleted.get(), n_ion_populations, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        fout_energy_deleted << energy_e_deleted << '\t' << energy_p_deleted << '\t' << energy_ph_deleted;
        for (int n=0; n<n_ion_populations; n++)
            fout_energy_deleted << '\t' << ienergy_deleted[n];
        fout_energy_deleted << endl;
    }
}

void write_density(bool write_x, bool write_y, bool write_z,
        std::string x_folder, std::string y_folder, std::string z_folder,
        bool write_ions = false, bool scale_j = false)
{
    char file_num_pchar[20];
    ofstream* pof_x = 0;
    ofstream* pof_y = 0;
    ofstream* pof_z = 0;

    sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);

    string file_name_x = data_folder + "/" + x_folder + file_num_pchar;
    string file_name_y = data_folder + "/" + y_folder + file_num_pchar;
    string file_name_z = data_folder + "/" + z_folder + file_num_pchar;

    int onx;
    int onx0;
    for(int i=0;i<n_sr;i++)
    {
        if (mpi_rank == i) {
            if (write_x) {
                pof_x = new ofstream(file_name_x.c_str(), output_mode);
            }

            if (write_y) {
                pof_y = new ofstream(file_name_y.c_str(), output_mode);
            }
            if (write_z) {
                pof_z = new ofstream(file_name_z.c_str(), output_mode);
            }

            if(i==n_sr-1)
                onx = nx_sr[i];
            else
                onx = nx_sr[i]-nx_ich/2;
            if(i==0)
                onx0 = 0;
            else
                onx0 = nx_ich/2;
            if (scale_j) {
                psr->scale_j(dx / dt);
            }
            psr->fout_rho(pof_x,pof_y,pof_z,onx0,onx, output_mode);

            if (pof_x) delete pof_x;
            if (pof_y) delete pof_y;
            if (pof_z) delete pof_z;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    pof_x = 0;
    pof_y = 0;
    pof_z = 0;
    
    int ii = get_sr_for_x(xlength-x0fout);
    if (mpi_rank == ii) {
        ios_base::openmode mode = output_mode | ios_base::app;
        if (write_x) {
            pof_x = new ofstream(file_name_x.c_str(), mode);
        }
        if (write_y) {
            pof_y = new ofstream(file_name_y.c_str(), mode);
        }
        if (write_z) {
            pof_z = new ofstream(file_name_z.c_str(), mode);
        }
        psr->fout_rho_yzplane(pof_x,pof_y,pof_z,get_xindex_in_sr(xlength-x0fout, ii), mode);
        if (pof_x) delete pof_x;
        if (pof_y) delete pof_y;
        if (pof_z) delete pof_z;
    }

    if (write_ions && (ions=="on" || ions=="positrons"))
    {
        char s_cmr[20];
        for (int n=0;n<n_ion_populations;n++)
        {
            sprintf(s_cmr,"%g",icmr[n]);
            string file_name = data_folder+"/irho_";
            file_name += s_cmr;
            file_name += "_";
            file_name += file_num_pchar;
            for(int i=0;i<n_sr;i++)
            {
                if (mpi_rank == i) {
                    ofstream fout_irho(file_name.c_str(), output_mode);
                    if(i==n_sr-1)
                        onx = nx_sr[i];
                    else
                        onx = nx_sr[i]-nx_ich/2;
                    if(i==0)
                        onx0 = 0;
                    else
                        onx0 = nx_ich/2;
                    psr->fout_irho(n,&fout_irho,onx0,onx, output_mode);
                    fout_irho.close();
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
            ii = get_sr_for_x(xlength-x0fout);
            if (mpi_rank == ii) {
                ios_base::openmode mode = output_mode | ios_base::app;
                ofstream fout_irho(file_name.c_str(), mode);
                psr->fout_irho_yzplane(n,&fout_irho,get_xindex_in_sr(xlength-x0fout, ii), mode);
                fout_irho.close();
            }
        }
    }
}

void write_spectrum_phasespace(bool write_p, bool write_ph)
{
    ios_base::openmode non_binary_mode;
    if (mpi_rank == 0) {
        non_binary_mode = ios_base::out;
    } else {
        non_binary_mode = ios_base::out | ios_base::app;
    }

    std::string file_name;
    char file_num_pchar[20];
    
    double* spectrum = new double[neps];
    for(int i=0;i<neps;i++) spectrum[i] = 0;
    double* spectrum_p = new double[neps_p];
    for(int i=0;i<neps_p;i++) spectrum_p[i] = 0;
    double* spectrum_ph = new double[neps_ph];
    for(int i=0;i<neps_ph;i++) spectrum_ph[i] = 0;
    int i_eps;
    double** spectrum_i = new double*[n_ion_populations];
    for (int m=0;m<n_ion_populations;m++) {
        spectrum_i[m] = new double[neps_i];
        for(int i=0;i<neps_i;i++)
            spectrum_i[m][i] = 0;
    }

    // sequential output of phasespace
    particle* current;
    for(int n=0;n<n_sr;n++)
    {
        if (mpi_rank == n) {
            file_name = data_folder+"/phasespace";
            sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
            file_name = file_name + file_num_pchar;
            ofstream fout_phasespace(file_name.c_str(), non_binary_mode);
            file_name = data_folder+"/phasespace_p";
            file_name = file_name + file_num_pchar;
            ofstream fout_phasespace_p(file_name.c_str(), non_binary_mode);
            file_name = data_folder+"/phasespace_ph";
            file_name = file_name + file_num_pchar;
            ofstream fout_phasespace_ph(file_name.c_str(), non_binary_mode);

            ofstream* fout_phasespace_i = new ofstream[n_ion_populations];
            for (int m=0;m<n_ion_populations;m++)
            {
                char s_cmr[20];
                sprintf(s_cmr,"%g",icmr[m]);
                file_name = data_folder+"/phasespace_";
                file_name += s_cmr;
                file_name += "_";
                file_name += file_num_pchar;
                fout_phasespace_i[m].open(file_name.c_str(), non_binary_mode);
            }
            for(int i=nm*(n!=0);i<nx_sr[n]-nm*(n!=n_sr-1);i++)
            {
                for(int j=0;j<ny_global;j++)
                {
                    for(int k=0;k<nz_global;k++)
                    {
                        current = psr->cp[i][j][k].pl.head;
                        while(current!=0)
                        {
                            if (current->cmr==-1)
                            { // electrons
                                i_eps = (current->g-1)*0.511/deps;
                                if((i_eps>-1)&&(i_eps<neps))
                                    spectrum[i_eps] = spectrum[i_eps] - current->q; // q<0 for electrons
                                // phase space
                                i_particle++;
                                if(i_particle>enthp)
                                {
                                    i_particle = 0;
                                    /* координаты выводятся
                                     * нормированными на лазерную
                                     * длину волны */
                                    fout_phasespace<<current->q<<"\n";
                                    fout_phasespace<<dx*(current->x+x0_sr[n])/(2*PI)<<"\n";
                                    fout_phasespace<<dy*(current->y)/(2*PI)<<"\n"<<dz*(current->z)/(2*PI)<<"\n";
                                    fout_phasespace<<current->ux<<"\n"<<current->uy<<"\n"<<current->uz<<"\n";
                                    fout_phasespace<<current->g<<"\n";
                                    fout_phasespace<<current->chi<<"\n";
                                }
                            }
                            else if (current->cmr==1 && write_p)
                            { // positrons
                                i_eps = (current->g-1)*0.511/deps_p;
                                if((i_eps>-1)&&(i_eps<neps_p))
                                    spectrum_p[i_eps] = spectrum_p[i_eps] + current->q; // q>0 for positrons
                                i_particle_p++;
                                if(i_particle_p>enthp_p)
                                {
                                    i_particle_p = 0;
                                    fout_phasespace_p<<current->q<<"\n";
                                    fout_phasespace_p<<dx*(current->x+x0_sr[n])/(2*PI)<<"\n";
                                    fout_phasespace_p<<dy*(current->y)/(2*PI)<<"\n"<<dz*(current->z)/(2*PI)<<"\n";
                                    fout_phasespace_p<<current->ux<<"\n"<<current->uy<<"\n"<<current->uz<<"\n";
                                    fout_phasespace_p<<current->g<<"\n";
                                    fout_phasespace_p<<current->chi<<"\n";
                                }
                            }
                            else if (current->cmr==0 && write_ph)
                            { // photons
                                i_eps = current->g*0.511/deps_ph;
                                if((i_eps>-1)&&(i_eps<neps_ph))
                                    spectrum_ph[i_eps] = spectrum_ph[i_eps] + current->q; // q>0 for photons
                                i_particle_ph++;
                                if(i_particle_ph>enthp_ph)
                                {
                                    i_particle_ph = 0;
                                    fout_phasespace_ph<<current->q<<"\n";
                                    fout_phasespace_ph<<dx*(current->x+x0_sr[n])/(2*PI)<<"\n";
                                    fout_phasespace_ph<<dy*(current->y)/(2*PI)<<"\n"<<dz*(current->z)/(2*PI)<<"\n";
                                    fout_phasespace_ph<<current->ux<<"\n"<<current->uy<<"\n"<<current->uz<<"\n";
                                    fout_phasespace_ph<<current->g<<"\n";
                                    fout_phasespace_ph<<current->chi<<"\n";
                                }
                            }
                            else
                            { // ions
                                int m;
                                m = 0;
                                while (m!=n_ion_populations)
                                {
                                    if (current->cmr==icmr[m])
                                    {
                                        i_eps = (current->g-1)*0.511*proton_mass/deps_i;
                                        if((i_eps>-1)&&(i_eps<neps_i))
                                            spectrum_i[m][i_eps] += current->q;
                                        i_particle_i++;
                                        if(i_particle_i>enthp_i)
                                        {
                                            i_particle_i = 0;
                                            fout_phasespace_i[m]<<current->q<<"\n";
                                            fout_phasespace_i[m]<<dx*(current->x+x0_sr[n])/(2*PI)<<"\n";
                                            fout_phasespace_i[m]<<dy*(current->y)/(2*PI)<<"\n"<<dz*(current->z)/(2*PI)<<"\n";
                                            fout_phasespace_i[m]<<current->ux<<"\n"<<current->uy<<"\n"<<current->uz<<"\n";
                                            fout_phasespace_i[m]<<current->g<<"\n";
                                            fout_phasespace_i[m]<<current->chi<<"\n";
                                        }
                                        m = n_ion_populations;
                                    }
                                    else
                                        m++;
                                }
                            }
                            current = current->next;
                        }
                    }
                }
            }
            fout_phasespace.close();
            fout_phasespace_p.close();
            fout_phasespace_ph.close();
            for (int m=0;m<n_ion_populations;m++)
            {
                fout_phasespace_i[m].close();
            }
            delete [] fout_phasespace_i;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // gathering spectra on rank 0 process
    if (mpi_rank == 0) {
        MPI_Reduce(MPI_IN_PLACE, spectrum, neps, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, spectrum_p, neps_p, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, spectrum_ph, neps_ph, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        for (int i=0; i<n_ion_populations; i++) {
            MPI_Reduce(MPI_IN_PLACE, spectrum_i[i], neps_i, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Reduce(spectrum, spectrum, neps, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(spectrum_p, spectrum_p, neps_p, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(spectrum_ph, spectrum_ph, neps_ph, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        for (int i=0; i<n_ion_populations; i++) {
            MPI_Reduce(spectrum_i[i], spectrum_i[i], neps_i, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
    }

    // output of spectra
    if (mpi_rank == 0) {
        file_name = data_folder+"/spectrum";
        sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
        file_name = file_name + file_num_pchar;
        ofstream fout_spectrum(file_name.c_str());
        file_name = data_folder+"/spectrum_p";
        file_name = file_name + file_num_pchar;
        ofstream fout_spectrum_p(file_name.c_str());
        file_name = data_folder+"/spectrum_ph";
        file_name = file_name + file_num_pchar;
        ofstream fout_spectrum_ph(file_name.c_str());

        ofstream* fout_spectrum_i = new ofstream[n_ion_populations];
        for (int m=0;m<n_ion_populations;m++) {
            char s_cmr[20];
            sprintf(s_cmr,"%g",icmr[m]);
            file_name = data_folder+"/spectrum_";
            file_name += s_cmr;
            file_name += "_";
            file_name += file_num_pchar;
            fout_spectrum_i[m].open(file_name.c_str());
        }

        // spectrum = dN/deps [1/MeV]
        double spectrum_norm = 1.11485e13*lambda*dx*dy*dz/(8*PI*PI*PI)/deps;
        double spectrum_norm_p = 1.11485e13*lambda*dx*dy*dz/(8*PI*PI*PI)/deps_p;
        double spectrum_norm_ph = 1.11485e13*lambda*dx*dy*dz/(8*PI*PI*PI)/deps_ph;
        for(int i=0;i<neps;i++)
        {
            spectrum[i] = spectrum[i]*spectrum_norm;
            fout_spectrum<<spectrum[i]<<"\n";
        }
        if (write_p)
        {
            for(int i=0;i<neps_p;i++)
            {
                spectrum_p[i] = spectrum_p[i]*spectrum_norm_p;
                fout_spectrum_p<<spectrum_p[i]<<"\n";
            }
        }
        if (write_ph)
        {
            for(int i=0;i<neps_ph;i++)
            {
                spectrum_ph[i] = spectrum_ph[i]*spectrum_norm_ph;
                fout_spectrum_ph<<spectrum_ph[i]<<"\n";
            }
        }
        // spectrum_i = dN/deps [1/MeV], but N is the number of nucleons, not ions
        double spectrum_norm_i;
        for (int m=0;m<n_ion_populations;m++)
        {
            spectrum_norm_i =  1.11485e13*lambda*dx*dy*dz/(8*PI*PI*PI)/(icmr[m]*deps_i);
            for(int i=0;i<neps_i;i++)
            {
                spectrum_i[m][i] = spectrum_i[m][i]*spectrum_norm_i;
                fout_spectrum_i[m]<<spectrum_i[m][i]<<"\n";
            }
        }

        fout_spectrum.close();
        fout_spectrum_p.close();
        fout_spectrum_ph.close();
        //
        for (int m=0;m<n_ion_populations;m++)
        {
            fout_spectrum_i[m].close();
        }
        delete[] fout_spectrum_i;
    }

    delete[] spectrum;
    delete[] spectrum_p;
    delete[] spectrum_ph;
    for (int m=0;m<n_ion_populations;m++)
        delete[] spectrum_i[m];
    delete[] spectrum_i;
}

void write_fields()
{
    std::string file_name;
    char file_num_pchar[20];
    int onx;
    int onx0;
    int ii;
    bool is_last_sr; // 0 - not last, 1 - last
    //
    if (e_components_for_output=="x"||e_components_for_output=="xy"||e_components_for_output=="xz"||e_components_for_output=="xyz")
    {
        file_name = data_folder+"/ex";
        sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
        file_name = file_name + file_num_pchar;
        for(int i=0;i<n_sr;i++)
        {
            if (mpi_rank == i) {
                ofstream fout_ex(file_name.c_str(), output_mode);
                if(i==n_sr-1)
                    onx = nx_sr[i];
                else
                    onx = nx_sr[i]-nx_ich/2;
                if(i==0)
                    onx0 = 0;
                else
                    onx0 = nx_ich/2;
                psr->fout_ex(&fout_ex,onx0,onx, output_mode);
                fout_ex.close();
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        ii = get_sr_for_x(xlength-x0fout);
        if (mpi_rank == ii) {
            ios_base::openmode mode = output_mode | ios_base::app;
            ofstream fout_ex(file_name.c_str(), mode);
            psr->fout_ex_yzplane(&fout_ex,get_xindex_in_sr(xlength-x0fout, ii), mode);
            fout_ex.close();
        }

    }
    //
    if (e_components_for_output=="y"||e_components_for_output=="xy"||e_components_for_output=="yz"||e_components_for_output=="xyz")
    {
        file_name = data_folder+"/ey";
        sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
        file_name = file_name + file_num_pchar;
        for(int i=0;i<n_sr;i++)
        {
            if (mpi_rank == i) {
                ofstream fout_ey(file_name.c_str(), output_mode);
                if(i==n_sr-1)
                    onx = nx_sr[i];
                else
                    onx = nx_sr[i]-nx_ich/2;
                if(i==0)
                    onx0 = 0;
                else
                    onx0 = nx_ich/2;
                psr->fout_ey(&fout_ey,onx0,onx, output_mode);
                fout_ey.close();
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        ii = get_sr_for_x(xlength-x0fout);
        if (mpi_rank == ii) {
            ios_base::openmode mode = output_mode | ios_base::app;
            ofstream fout_ey(file_name.c_str(), mode);
            psr->fout_ey_yzplane(&fout_ey,get_xindex_in_sr(xlength-x0fout, ii), mode);
            fout_ey.close();
        }

    }
    //
    if (e_components_for_output=="z"||e_components_for_output=="xz"||e_components_for_output=="yz"||e_components_for_output=="xyz")
    {
        file_name = data_folder+"/ez";
        sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
        file_name = file_name + file_num_pchar;
        for(int i=0;i<n_sr;i++)
        {
            if (mpi_rank == i) {
                ofstream fout_ez(file_name.c_str(), output_mode);
                if(i==n_sr-1)
                    onx = nx_sr[i];
                else
                    onx = nx_sr[i]-nx_ich/2;
                if(i==0)
                    onx0 = 0;
                else
                    onx0 = nx_ich/2;
                psr->fout_ez(&fout_ez,onx0,onx, output_mode);
                fout_ez.close();
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        ii = get_sr_for_x(xlength-x0fout);
        if (mpi_rank == ii) {
            ios_base::openmode mode = output_mode | ios_base::app;
            ofstream fout_ez(file_name.c_str(), mode);
            psr->fout_ez_yzplane(&fout_ez,get_xindex_in_sr(xlength-x0fout, ii), mode);
            fout_ez.close();
        }
    }
    //
    if (b_components_for_output=="x"||b_components_for_output=="xy"||b_components_for_output=="xz"||b_components_for_output=="xyz")
    {
        file_name = data_folder+"/bx";
        sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
        file_name = file_name + file_num_pchar;
        for(int i=0;i<n_sr;i++)
        {
            if (mpi_rank == i) {
                ofstream fout_bx(file_name.c_str(), output_mode);
                if(i==n_sr-1)
                    onx = nx_sr[i];
                else
                    onx = nx_sr[i]-nx_ich/2;
                if(i==0)
                    onx0 = 0;
                else
                    onx0 = nx_ich/2;
                psr->fout_bx(&fout_bx,onx0,onx, output_mode);
                fout_bx.close();
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        ii = get_sr_for_x(xlength-x0fout);
        if (mpi_rank == ii) {
            ios_base::openmode mode = output_mode | ios_base::app;
            ofstream fout_bx(file_name.c_str(), mode);
            psr->fout_bx_yzplane(&fout_bx,get_xindex_in_sr(xlength-x0fout, ii), mode);
            fout_bx.close();
        }
    }
    //
    if (b_components_for_output=="y"||b_components_for_output=="xy"||b_components_for_output=="yz"||b_components_for_output=="xyz")
    {
        file_name = data_folder+"/by";
        sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
        file_name = file_name + file_num_pchar;
        for(int i=0;i<n_sr;i++)
        {
            if (mpi_rank == i) {
                ofstream fout_by(file_name.c_str(), output_mode);
                if(i==n_sr-1)
                    onx = nx_sr[i];
                else
                    onx = nx_sr[i]-nx_ich/2;
                if(i==0)
                    onx0 = 0;
                else
                    onx0 = nx_ich/2;
                psr->fout_by(&fout_by,onx0,onx, output_mode);
                fout_by.close();
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        ii = get_sr_for_x(xlength-x0fout);
        if (mpi_rank == ii) {
            ios_base::openmode mode = output_mode | ios_base::app;
            ofstream fout_by(file_name.c_str(), mode);
            psr->fout_by_yzplane(&fout_by,get_xindex_in_sr(xlength-x0fout, ii), mode);
            fout_by.close();
        }
    }
    //
    if (b_components_for_output=="z"||b_components_for_output=="xz"||b_components_for_output=="yz"||b_components_for_output=="xyz")
    {
        file_name = data_folder+"/bz";
        sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
        file_name = file_name + file_num_pchar;
        for(int i=0;i<n_sr;i++)
        {
            if (mpi_rank == i) {
                ofstream fout_bz(file_name.c_str(), output_mode);
                if(i==n_sr-1)
                    onx = nx_sr[i];
                else
                    onx = nx_sr[i]-nx_ich/2;
                if(i==0)
                    onx0 = 0;
                else
                    onx0 = nx_ich/2;
                psr->fout_bz(&fout_bz,onx0,onx, output_mode);
                fout_bz.close();
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        ii = get_sr_for_x(xlength-x0fout);
        if (mpi_rank == ii) {
            ios_base::openmode mode = output_mode | ios_base::app;
            ofstream fout_bz(file_name.c_str(), mode);
            psr->fout_bz_yzplane(&fout_bz,get_xindex_in_sr(xlength-x0fout, ii), mode);
            fout_bz.close();
        }
    }

    ios_base::openmode non_binary_mode;
    if (mpi_rank == 0) {
        non_binary_mode = ios_base::out;
    } else {
        non_binary_mode = ios_base::out | ios_base::app;
    }
    //
    file_name = data_folder+"/w";
    sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
    file_name = file_name + file_num_pchar;
    for(int i=0;i<n_sr;i++)
    {
        if (mpi_rank == i) {
            ofstream fout_w(file_name.c_str(), non_binary_mode);
            if(i==n_sr-1)
                onx = nx_sr[i];
            else
                onx = nx_sr[i]-nx_ich/2;
            if(i==0)
                onx0 = 0;
            else
                onx0 = nx_ich/2;
            if (i==n_sr-1)
                is_last_sr = 1;
            else
                is_last_sr = 0;
            psr->fout_w(&fout_w,onx0,onx,is_last_sr);
            fout_w.close();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    ii = get_sr_for_x(xlength-x0fout);
    if (mpi_rank == ii) {
        ofstream fout_w(file_name.c_str(), non_binary_mode | ios_base::app);
        psr->fout_w_yzplane(&fout_w,get_xindex_in_sr(xlength-x0fout, ii));
        fout_w.close();
    }

    //
    file_name = data_folder+"/inv";
    sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
    file_name = file_name + file_num_pchar;
    for(int i=0;i<n_sr;i++)
    {
        if (mpi_rank == i) {
            ofstream fout_inv(file_name.c_str(), non_binary_mode);
            if(i==n_sr-1)
                onx = nx_sr[i];
            else
                onx = nx_sr[i]-nx_ich/2;
            if(i==0)
                onx0 = 0;
            else
                onx0 = nx_ich/2;
            if (i==n_sr-1)
                is_last_sr = 1;
            else
                is_last_sr = 0;
            psr->fout_inv(&fout_inv,onx0,onx,is_last_sr);
            fout_inv.close();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    ii = get_sr_for_x(xlength-x0fout);
    if (mpi_rank == ii) {
        ofstream fout_inv(file_name.c_str(), non_binary_mode | ios_base::app);
        psr->fout_inv_yzplane(&fout_inv,get_xindex_in_sr(xlength-x0fout, ii));
        fout_inv.close();
    }
}

vector<double> calculate_global_layer_weights() {
    auto weights = psr->calculate_layer_weights(balancing_particle_weight);

    vector<double> global_weights;

    if (mpi_rank == 0) {
        global_weights = vector<double>(nx_global);

        int left = 0;
        int right = ((n_sr > 0) ? nx_sr[0] - nx_ich / 2 : nx_sr[0]);

        for (int i = left; i < right; i++) {
            global_weights[i] = weights[i];
        }

        for (int n = 1; n < n_sr; n++) {
            left = nx_ich / 2;
            right = (n == n_sr-1 ? nx_sr[n] : nx_sr[n] - nx_ich / 2);
            MPI_Recv(&(global_weights[x0_sr[n] + left]), right-left, MPI_DOUBLE, n, n, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    } else {
        int left = nx_ich / 2;
        int right = (mpi_rank == n_sr-1 ? nx_sr[mpi_rank] : nx_sr[mpi_rank] - nx_ich / 2);
        MPI_Send(&(weights[left]), right - left, MPI_DOUBLE, 0, mpi_rank, MPI_COMM_WORLD);
    }

    return global_weights;
}

void write_layer_weights() {
    auto weights = calculate_global_layer_weights();

    if (mpi_rank == 0) {
        string file_name;
        char file_num_pchar[20];
        
        file_name = data_folder+"/weights";
        sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
        file_name = file_name + file_num_pchar;

        ofstream fout_weights(file_name.c_str(), ios_base::out);

        int length = weights.size();
        for (int i = 0; i < length; i++) {
            fout_weights << weights[i] << "\n";
        }
    }
}

void init_fields()
{
    int i = mpi_rank;
    if (f_envelope=="focused")
    {
        const char* tmpl = lp_reflection.c_str();
        int lp_len = lp_reflection.size();
        const char* tmpf = f_reflection.c_str();
        int f_len = f_reflection.size();
        std::string lp_reflection1 = "";
        std::string lp_reflection2 = "";
        std::string lp_reflection3 = "";
        std::string f_reflection1 = "";
        std::string f_reflection2 = "";
        std::string f_reflection3 = "";
        int jj;
        for( jj=0;jj<lp_len&&tmpl[jj]!='&';jj++ ) {
            lp_reflection1 += tmpl[jj];
        }
        jj++;
        if(jj<lp_len) {
            for( ;jj<lp_len&&tmpl[jj]!='&';jj++ ) {
                lp_reflection2 += tmpl[jj];
            }
        }
        jj++;
        if(jj<lp_len) {
            for( ;jj<lp_len&&tmpl[jj]!='&';jj++ ) {
                lp_reflection3 += tmpl[jj];
            }
        }
        for( jj=0;jj<f_len&&tmpf[jj]!='&';jj++ ) {
            f_reflection1 += tmpf[jj];
        }
        jj++;
        if(jj<f_len) {
            for( ;jj<f_len&&tmpf[jj]!='&';jj++ ) {
                f_reflection2 += tmpf[jj];
            }
        }
        jj++;
        if(jj<f_len) {
            for( ;jj<f_len&&tmpf[jj]!='&';jj++ ) {
                f_reflection3 += tmpf[jj];
            }
        }
        if (shenergy == 0) 
            psr->f_init_focused(a0y,a0z,xsigma,sigma0,xlength/2+x0-dx*x0_sr[i],x0,b_sign,phase,y00,z00,0,0,sscos,1,xtarget,ytarget,ztarget);
        else  // adding second harmonic
        {
            double alpha, beta;
            alpha = sqrt(1 - shenergy);
            beta = sqrt(shenergy);
            psr->f_init_focused(a0y * alpha, a0z * alpha, xsigma,sigma0,xlength/2+x0-dx*x0_sr[i],x0,b_sign,phase,y00,z00,0,0,sscos,1,xtarget,ytarget,ztarget);
            psr->f_init_focused(a0y * beta, a0z * beta, xsigma,sigma0,xlength/2+x0-dx*x0_sr[i],x0,b_sign, shphase,y00,z00, 1, 0,sscos, 2,xtarget,ytarget,ztarget);
        }
        if (phi!=0) {
            psr->f_init_focused(a0y,a0z,xsigma,sigma0,xlength/2-x0-dx*x0_sr[i],-x0,b_sign,phase,-y00,-z00,1,phi,sscos,1,xtarget,ytarget,ztarget);
        }
        else if (lp_reflection1=="xy") {
            psr->f_init_focused((1-2*(f_reflection1=="y"))*a0y,(1-2*(f_reflection1=="z"))*a0z,xsigma,sigma0,xlength/2+x0-dx*x0_sr[i],x0,b_sign,phase,y00,-z00,1,0,sscos,1,xtarget,ytarget,ztarget);
            if (lp_reflection2=="xz") {
                psr->f_init_focused((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2+x0-dx*x0_sr[i],x0,b_sign,phase,-y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                psr->f_init_focused((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2+x0-dx*x0_sr[i],x0,b_sign,phase,-y00,-z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                if (lp_reflection3=="yz") {
                    psr->f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*x0_sr[i],-x0,b_sign,phase,y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                    psr->f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*x0_sr[i],-x0,b_sign,phase,y00,-z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                    psr->f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*x0_sr[i],-x0,b_sign,phase,-y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                    psr->f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*x0_sr[i],-x0,b_sign,phase,-y00,-z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                }
            }
            else if (lp_reflection2=="yz") {
                psr->f_init_focused((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*x0_sr[i],-x0,b_sign,phase,y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                psr->f_init_focused((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*x0_sr[i],-x0,b_sign,phase,y00,-z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                if (lp_reflection3=="xz") {
                    psr->f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*x0_sr[i],-x0,b_sign,phase,-y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                    psr->f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2+x0-dx*x0_sr[i],x0,b_sign,phase,-y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                    psr->f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*x0_sr[i],-x0,b_sign,phase,-y00,-z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                    psr->f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2+x0-dx*x0_sr[i],x0,b_sign,phase,-y00,-z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                }
            }
        }
        else if (lp_reflection1=="xz") {
            psr->f_init_focused((1-2*(f_reflection1=="y"))*a0y,(1-2*(f_reflection1=="z"))*a0z,xsigma,sigma0,xlength/2+x0-dx*x0_sr[i],x0,b_sign,phase,-y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
            if (lp_reflection2=="yz") {
                psr->f_init_focused((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*x0_sr[i],-x0,b_sign,phase,-y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                psr->f_init_focused((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*x0_sr[i],-x0,b_sign,phase,y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
            }
        }
        else if (lp_reflection1=="yz") {
            psr->f_init_focused((1-2*(f_reflection1=="y"))*a0y,(1-2*(f_reflection1=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*x0_sr[i],-x0,b_sign,phase,y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
        }
    } else if (f_envelope == "uniformB") {
        psr->f_init_uniformB(a0y, a0z);
    }
    else // f_envelope == "cos"
    {
        if (f_envelope != "cos") {
            if (mpi_rank == 0) {
                cout << TERM_RED << "Unknown f_envelope value [" << f_envelope << "] during initialization"
                     << TERM_NO_COLOR << endl;
            }
        }
        const char* tmpl = lp_reflection.c_str();
        int lp_len = lp_reflection.size();
        const char* tmpf = f_reflection.c_str();
        int f_len = f_reflection.size();
        std::string lp_reflection1 = "";
        std::string lp_reflection2 = "";
        std::string lp_reflection3 = "";
        std::string f_reflection1 = "";
        std::string f_reflection2 = "";
        std::string f_reflection3 = "";
        int jj;
        for( jj=0;jj<lp_len&&tmpl[jj]!='&';jj++ ) {
            lp_reflection1 += tmpl[jj];
        }
        jj++;
        if(jj<lp_len) {
            for( ;jj<lp_len&&tmpl[jj]!='&';jj++ ) {
                lp_reflection2 += tmpl[jj];
            }
        }
        jj++;
        if(jj<lp_len) {
            for( ;jj<lp_len&&tmpl[jj]!='&';jj++ ) {
                lp_reflection3 += tmpl[jj];
            }
        }
        for( jj=0;jj<f_len&&tmpf[jj]!='&';jj++ ) {
            f_reflection1 += tmpf[jj];
        }
        jj++;
        if(jj<f_len) {
            for( ;jj<f_len&&tmpf[jj]!='&';jj++ ) {
                f_reflection2 += tmpf[jj];
            }
        }
        jj++;
        if(jj<f_len) {
            for( ;jj<f_len&&tmpf[jj]!='&';jj++ ) {
                f_reflection3 += tmpf[jj];
            }
        }

        if (phi_rotate != 0) {
            psr->f_init_cos(a0y,a0z,xsigma,ysigma,zsigma,xlength/2+x0-dx*x0_sr[i],sscos,b_sign,x0,phase,y00,z00,1,phi_rotate,xtarget,ytarget,ztarget);
        } else {
            psr->f_init_cos(a0y,a0z,xsigma,ysigma,zsigma,xlength/2+x0-dx*x0_sr[i],sscos,b_sign,x0,phase,y00,z00,true,0,xtarget,ytarget,ztarget);
            if (phi!=0) {
                psr->f_init_cos(a0y,a0z,xsigma,ysigma,zsigma,xlength/2-x0-dx*x0_sr[i],sscos,b_sign,-x0,phase,-y00,-z00,1,phi,xtarget,ytarget,ztarget);
            }
            else if (lp_reflection1=="xy") {
                psr->f_init_cos((1-2*(f_reflection1=="y"))*a0y,(1-2*(f_reflection1=="z"))*a0z,xsigma,ysigma,zsigma,xlength/2+x0-dx*x0_sr[i],sscos,b_sign,x0,phase,y00,-z00,1,0,xtarget,ytarget,ztarget);
                if (lp_reflection2=="xz") {
                    psr->f_init_cos((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,xlength/2+x0-dx*x0_sr[i],sscos,b_sign,x0,phase,-y00,z00,1,0,xtarget,ytarget,ztarget);
                    psr->f_init_cos((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,xlength/2+x0-dx*x0_sr[i],sscos,b_sign,x0,phase,-y00,-z00,1,0,xtarget,ytarget,ztarget);
                    if (lp_reflection3=="yz") {
                        psr->f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*x0_sr[i],sscos,b_sign,-x0,phase,y00,z00,1,0,xtarget,ytarget,ztarget);
                        psr->f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*x0_sr[i],sscos,b_sign,-x0,phase,y00,-z00,1,0,xtarget,ytarget,ztarget);
                        psr->f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*x0_sr[i],sscos,b_sign,-x0,phase,-y00,z00,1,0,xtarget,ytarget,ztarget);
                        psr->f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*x0_sr[i],sscos,b_sign,-x0,phase,-y00,-z00,1,0,xtarget,ytarget,ztarget);
                    }
                }
                else if (lp_reflection2=="yz") {
                    psr->f_init_cos((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*x0_sr[i],sscos,b_sign,-x0,phase,y00,z00,1,0,xtarget,ytarget,ztarget);
                    psr->f_init_cos((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*x0_sr[i],sscos,b_sign,-x0,phase,y00,-z00,1,0,xtarget,ytarget,ztarget);
                    if (lp_reflection3=="xz") {
                        psr->f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*x0_sr[i],sscos,b_sign,-x0,phase,-y00,z00,1,0,xtarget,ytarget,ztarget);
                        psr->f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,x0+xlength/2-dx*x0_sr[i],sscos,b_sign,x0,phase,-y00,z00,1,0,xtarget,ytarget,ztarget);
                        psr->f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*x0_sr[i],sscos,b_sign,-x0,phase,-y00,-z00,1,0,xtarget,ytarget,ztarget);
                        psr->f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,x0+xlength/2-dx*x0_sr[i],sscos,b_sign,x0,phase,-y00,-z00,1,0,xtarget,ytarget,ztarget);
                    }
                }
            }
            else if (lp_reflection1=="xz") {
                psr->f_init_cos((1-2*(f_reflection1=="y"))*a0y,(1-2*(f_reflection1=="z"))*a0z,xsigma,ysigma,zsigma,xlength/2+x0-dx*x0_sr[i],sscos,b_sign,x0,phase,-y00,z00,1,0,xtarget,ytarget,ztarget);
                if (lp_reflection2=="yz") {
                    psr->f_init_cos((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*x0_sr[i],sscos,b_sign,-x0,phase,-y00,z00,1,0,xtarget,ytarget,ztarget);
                    psr->f_init_cos((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*x0_sr[i],sscos,b_sign,-x0,phase,y00,z00,1,0,xtarget,ytarget,ztarget);
                }
            }
            else if (lp_reflection1=="yz") {
                psr->f_init_cos((1-2*(f_reflection1=="y"))*a0y,(1-2*(f_reflection1=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*x0_sr[i],sscos,b_sign,-x0,phase,y00,z00,1,0,xtarget,ytarget,ztarget);
            }
        }
    }
}

void init_beam()
{
    int i = mpi_rank;
    if (beam_particles=="p")
        psr->add_beam(1,Nb*1.061e-11/(xb*rb*rb*lambda),((epsb>0)-(epsb<0))*sqrt(epsb*epsb/(0.511*0.511)-1),xb,rb,xlength-x0b-dx*x0_sr[i],y0b,phib);
    else if (beam_particles=="ph")
        psr->add_beam(0,Nb*1.061e-11/(xb*rb*rb*lambda),epsb/0.511,xb,rb,xlength-x0b-dx*x0_sr[i],y0b,phib);
    else
        psr->add_beam(-1,Nb*1.061e-11/(xb*rb*rb*lambda),((epsb>0)-(epsb<0))*sqrt(epsb*epsb/(0.511*0.511)-1),xb,rb,xlength-x0b-dx*x0_sr[i],y0b,phib);
}

void init_films()
{
    film* tmp_p_film = p_last_film;
    while (tmp_p_film!=0)
    {
        /* x0film - координата левой границы плёнки, filmwidth - её
         * толщина, gradwidth - толщина части плёнки с линейным ростом
         * плотности */
        if (tmp_p_film->xnpic_film == 0 || tmp_p_film->ynpic_film == 0 || tmp_p_film->znpic_film == 0)
        {
            tmp_p_film->xnpic_film = xnpic;
            tmp_p_film->ynpic_film = ynpic;
            tmp_p_film->znpic_film = znpic;
        }
        if (tmp_p_film->ne != 0)
        {
            tmp_p_film->ne_y0 = tmp_p_film->ne;
            tmp_p_film->ne_y1 = tmp_p_film->ne;
        }

        psr->film(
            tmp_p_film->x0-dx*x0_sr[mpi_rank], 
            tmp_p_film->x0+tmp_p_film->filmwidth-dx*x0_sr[mpi_rank],
            tmp_p_film->ne_y0/(1.11485e+13/lambda/lambda),
            tmp_p_film->ne_y1/(1.11485e+13/lambda/lambda),
            ions == "on",
            1/(proton_mass*tmp_p_film->mcr),
            tmp_p_film->gradwidth,
            tmp_p_film->y0,
            tmp_p_film->y1,
            tmp_p_film->z0,
            tmp_p_film->z1,
            tmp_p_film->T,
            tmp_p_film->vx,
            nelflow != 0 || nerflow != 0,
            tmp_p_film->xnpic_film,
            tmp_p_film->ynpic_film,
            tmp_p_film->znpic_film,
            false,
            tmp_p_film->gradwidth_y
        );
        tmp_p_film = tmp_p_film->prev;
    }
}

void send_field_slice(int left, int width, int destination_rank) {
    int tag = 0;

    const int size = 3 * width * ny_global * nz_global;
    MPI_Send(psr->ce[left][0], size, MPI_DOUBLE, destination_rank, tag++, MPI_COMM_WORLD);
    MPI_Send(psr->cb[left][0], size, MPI_DOUBLE, destination_rank, tag++, MPI_COMM_WORLD);
    MPI_Send(psr->cj[left][0], size, MPI_DOUBLE, destination_rank, tag++, MPI_COMM_WORLD);
    MPI_Send(psr->cbe[left][0], size, MPI_DOUBLE, destination_rank, tag++, MPI_COMM_WORLD);
}

void receive_field_slice(int left, int width, int source_rank) {
    int tag = 0;

    const int size = 3 * width * ny_global * nz_global;
    MPI_Recv(psr->ce[left][0], size, MPI_DOUBLE, source_rank, tag++, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(psr->cb[left][0], size, MPI_DOUBLE, source_rank, tag++, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(psr->cj[left][0], size, MPI_DOUBLE, source_rank, tag++, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(psr->cbe[left][0], size, MPI_DOUBLE, source_rank, tag++, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void exchange_fields(int nm1, int nm2) {
    for (int mod=0; mod<2; mod++) {

        // first even processes send, while odd receive, then vice versa
        if (mpi_rank % 2 == mod) {
            if (mpi_rank > 0) {
                send_field_slice(nm1, nm2, mpi_rank-1);
            }

            if (mpi_rank < n_sr-1) {
                send_field_slice(psr->get_nx() - nm2 - nm1, nm1, mpi_rank+1);
            }
        } else {
            if (mpi_rank < n_sr-1) {
                receive_field_slice(psr->get_nx() - nm2, nm2, mpi_rank+1);
            }

            if (mpi_rank > 0) {
                receive_field_slice(0, nm1, mpi_rank-1);
            }
        }
    }
}

void pack_cell(cellp & p, vector<int> & particle_numbers, size_t & pn_index, vector<particle> & particles, size_t & p_index) {
    particle* current = p.pl.head;
    particle_numbers[pn_index] = 0;
    while (current != 0) {
        if (particles.size() <= p_index) {
            particles.resize(3 * particles.size() / 2 + 1);
        }
        particles[p_index++] = *current;

        current = current->next;
        particle_numbers[pn_index]++;
    }
    pn_index++;
}

void pack_cell(int i, int j, int k, vector<int> & particle_numbers, size_t & pn_index, vector<particle> & particles, size_t & p_index) {
    particle* current = psr->cp[i][j][k].pl.head;
    particle_numbers[pn_index] = 0;
    while (current != 0) {
        if (particles.size() <= p_index) {
            particles.resize(3 * particles.size() / 2 + 1);
        }
        particles[p_index++] = *current;

        current = current->next;
        particle_numbers[pn_index]++;
    }
    pn_index++;
}

void unpack_cell(int i, int j, int k, vector<int> & particle_numbers, size_t & pn_index, vector<particle> & particles, size_t & p_index) {
    plist & pl = psr->cp[i][j][k].pl;
    psr->erase(pl);

    for (int ii=0; ii<particle_numbers[pn_index]; ii++) {
        particle * tmp = psr->new_particle();
        *tmp = particles[p_index++];
        tmp->previous = pl.start;
        pl.start = tmp;
    }
    pl.head = 0;
    while (pl.start != 0) {
        pl.start->next = pl.head;
        pl.head = pl.start;
        pl.start = pl.start->previous;
    }
    pl.start = pl.head;

    pn_index++;
}

int pack_particle_slice(field3d<cellp> & cp, int left, int width, vector<int> & particle_numbers, vector<particle> & particles) {
    size_t pn_index = 0;
    size_t p_index = 0;

    for (int i=left; i<left+width; i++) {
        for (int j=0;j<ny_global;j++) {
            for (int k=0;k<nz_global;k++) {
                pack_cell(cp[i][j][k], particle_numbers, pn_index, particles, p_index);
            }
        }
    }

    return static_cast<int>(p_index);
}

int pack_particle_slice(int left, int width, vector<int> & particle_numbers, vector<particle> & particles) {
    size_t pn_index = 0;
    size_t p_index = 0;

    for (int i=left; i<left+width; i++) {
        for (int j=0;j<ny_global;j++) {
            for (int k=0;k<nz_global;k++) {
                pack_cell(i, j, k, particle_numbers, pn_index, particles, p_index);
            }
        }
    }

    return static_cast<int>(p_index);
}

void unpack_particle_slice(int left, int width, vector<int> & particle_numbers, vector<particle> & particles, int x_diff) {
    size_t pn_index = 0;
    size_t p_index = 0;

    for (int i=left; i<left+width; i++) {
        for (int j=0;j<ny_global;j++) {
            for (int k=0;k<nz_global;k++) {
                unpack_cell(i, j, k, particle_numbers, pn_index, particles, p_index);
                psr->cp[i][j][k].pl.xplus(x_diff);
            }
        }
    }
}

void send_particle_slice(field3d<cellp> & cp, int left, int width, int rank) {
    vector<int> particle_numbers(width * ny_global * nz_global);
    vector<particle> particles(100);
    
    int particles_to_send = pack_particle_slice(cp, left, width, particle_numbers, particles);

    MPI_Send(&particles_to_send, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
    MPI_Send(&(particle_numbers[0]), particle_numbers.size(), MPI_INT, rank, 1, MPI_COMM_WORLD);
    MPI_Send(&(particles[0]), particles_to_send, MPI_PARTICLE, rank, 2, MPI_COMM_WORLD);
}

void receive_particle_slice(int left, int width, int rank, int x_diff) {
    int particles_to_receive;
    MPI_Recv(&particles_to_receive, 1, MPI_INT, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    vector<int> particle_numbers(width * ny_global * nz_global);
    MPI_Recv(&(particle_numbers[0]), particle_numbers.size(), MPI_INT, rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    vector<particle> particles(particles_to_receive);

    MPI_Recv(&(particles[0]), particles_to_receive, MPI_PARTICLE, rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    unpack_particle_slice(left, width, particle_numbers, particles, x_diff);
}

void exchange_particle_slices(int send_left, int send_width, int receive_left, int receive_width, int rank,
        vector<int> & particle_numbers, int pn_size, vector<particle> & particles, int x_diff) {
    int particles_to_send = pack_particle_slice(send_left, send_width, particle_numbers, particles);
    int particles_to_receive = 0;

    MPI_Sendrecv(&particles_to_send, 1, MPI_INT, rank, 0, &particles_to_receive, 1, MPI_INT, rank, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);

    if (particles_to_receive > static_cast<int>(particles.size())) {
        particles.resize(3 * particles_to_receive / 2 + 1);
    }

    MPI_Sendrecv_replace(&(particle_numbers[0]), pn_size, MPI_INT, rank, 1, rank, 1, MPI_COMM_WORLD,
    MPI_STATUSES_IGNORE);

    int particles_size = max(particles_to_receive, particles_to_send);

    MPI_Sendrecv_replace(&(particles[0]), particles_size, MPI_PARTICLE, rank, 2, rank, 2, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);

    unpack_particle_slice(receive_left, receive_width, particle_numbers, particles, x_diff);
}

void exchange_particles(const int nm1, const int nm2) {
    static vector<int> particle_numbers{};
    static vector<particle> particles{};

    const int size1 = ny_global * nz_global * nm1;
    const int size2 = ny_global * nz_global * nm2;
    const int size = max(size1, size2);

    if (static_cast<int>(particle_numbers.size()) < size) {
        particle_numbers.resize(size);
    }

    for (int mod=0; mod<2; mod++) {

        // first even processes send, while odd receive, then vice versa
        if (mpi_rank % 2 == mod) {
            if (mpi_rank > 0) {
                exchange_particle_slices(nm1, nm2, 0, nm1, mpi_rank - 1, particle_numbers, size, particles,
                        -nx_sr[mpi_rank - 1] + nx_ich);
            }
        } else {
            if (mpi_rank < n_sr-1) {
                exchange_particle_slices(psr->get_nx() - nm2 - nm1, nm1, psr->get_nx() - nm2, nm2, mpi_rank + 1,
                        particle_numbers, size, particles, nx_sr[mpi_rank] - nx_ich);
            }
        }
    }

}

void synchronize_regions(bool is_moving_window_iteration) {
    int nm1, nm2;
    if (is_moving_window_iteration) {
        nm1 = nm-1;
        nm2 = nm+1;
    } else {
        nm1 = nm;
        nm2 = nm;
    }

    exchange_fields(nm1, nm2);
    exchange_particles(nm1, nm2);
}

void start_tracking()
{
    if (mpi_rank == 0) {
        cout << "Tracking started for particles: " << particles_to_track << endl;
    }
    if (xtr1 < 0 || xtr2 < 0 || ytr1 < 0 || ytr2 < 0 || ztr1 < 0 || ztr2 < 0 ||
        xtr1 >= xlength || xtr2 >= xlength || ytr1 >= ylength || ytr2 >= ylength || ztr1 >= zlength || ztr2 >= zlength)
    {
        if (mpi_rank == 0) {
            cout << "Error - tracks outside of compulational domain" << endl;
        }
        return;
    }
    long trn = 1; // trn = 0 for untracked particles
    if (tr_init==0) {
        int x1,y1,z1,x2,y2,z2;
        x1 = -1;
        y1 = 0;
        z1 = 0;
        for (int i=0;i<n_tracks;i++) {
            x2 = (xtr1 + (xtr2-xtr1)*i/n_tracks)/dx;
            y2 = (ytr1 + (ytr2-ytr1)*i/n_tracks)/dy;
            z2 = (ztr1 + (ztr2-ztr1)*i/n_tracks)/dz;
            if (x2!=x1 || y2!=y1 || z2!=z1) {
                x1 = x2;
                y1 = y2;
                z1 = z2;
                int n,x;
                n = get_sr_for_index(x2);
                x = x2 - x0_sr[n];
                if (n>0 && x<nm) {
                    n = n - 1;
                    x = x + nx_sr[n] - nx_ich;
                }
                if (mpi_rank == n) {
                    particle* h = psr->cp[x][y2][z2].pl.head;
                    particle* p = h;
                    bool b = 1;
                    if (particles_to_track.find('e') != string::npos)
                    {
                        while (p!=0 && b) {
                            if (p->cmr==-1) {
                                p->trn = trn;
                                trn++;
                                b = 0;
                            }
                            p = p->next;
                        }
                    }
                    if (particles_to_track.find('p') != string::npos)
                    {
                        p = h;
                        b = 1;
                        while (p!=0 && b) {
                            if (p->cmr==1) {
                                p->trn = trn;
                                trn++;
                                b = 0;
                            }
                            p = p->next;
                        }
                    }
                    if (particles_to_track.find('g') != string::npos)
                    {
                        p = h;
                        b = 1;
                        while (p!=0 && b) {
                            if (p->cmr==0) {
                                p->trn = trn;
                                trn++;
                                b = 0;
                            }
                            p = p->next;
                        }
                    }
                    if (particles_to_track.find('i') != string::npos)
                    {
                        for (int m=0;m<n_ion_populations;m++) {
                            p = h;
                            b = 1;
                            while (p!=0 && b) {
                                if (p->cmr==icmr[m]) {
                                    p->trn = trn;
                                    trn++;
                                    b = 0;
                                }
                                p = p->next;
                            }
                        }
                    }
                }
                MPI_Bcast(&trn, 1, MPI_LONG, n, MPI_COMM_WORLD);
            }
        }
    } else {
        for (int i=0;i<n_tracks;i++) {
            int x1 = ( xtr1 + ( xtr2 - xtr1 ) * rand( ) / RAND_MAX ) / dx;
            int y1 = ( ytr1 + ( ytr2 - ytr1 ) * rand( ) / RAND_MAX ) / dy;
            int z1 = ( ztr1 + ( ztr2 - ztr1 ) * rand( ) / RAND_MAX ) / dz;
            MPI_Bcast(&x1, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&y1, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&z1, 1, MPI_INT, 0, MPI_COMM_WORLD);

            int n,x;
            n = get_sr_for_index(x1);
            x = x1 - x0_sr[n];
            if (n!=0 && x<nm) {
                n = n - 1;
                x = x + nx_sr[n] - nx_ich;
            }
            if (mpi_rank == n) {
                particle* h = psr->cp[x][y1][z1].pl.head;
                particle* p = h;
                bool b = 1;
                if (particles_to_track.find('e') != string::npos)
                {
                    while (p!=0 && b) {
                        if (p->cmr==-1) {
                            p->trn = trn;
                            trn++;
                            b = 0;
                        }
                        p = p->next;
                    }
                }
                if (particles_to_track.find('p') != string::npos)
                {
                    p = h;
                    b = 1;
                    while (p!=0 && b) {
                        if (p->cmr==1) {
                            p->trn = trn;
                            trn++;
                            b = 0;
                        }
                        p = p->next;
                    }
                }
                if (particles_to_track.find('g') != string::npos)
                {
                    p = h;
                    b = 1;
                    while (p!=0 && b) {
                        if (p->cmr==0) {
                            p->trn = trn;
                            trn++;
                            b = 0;
                        }
                        p = p->next;
                    }
                }
                if (particles_to_track.find('i') != string::npos)
                {
                    for (int m=0;m<n_ion_populations;m++) {
                        p = h;
                        b = 1;
                        while (p!=0 && b) {
                            if (p->cmr==icmr[m]) {
                                p->trn = trn;
                                trn++;
                                b = 0;
                            }
                            p = p->next;
                        }
                    }
                }
            }
            MPI_Bcast(&trn, 1, MPI_LONG, n, MPI_COMM_WORLD);
        }
    }

    // make the same particles tracked in all regions, so that tracks aren't lost at the first exchange
    synchronize_regions(false);
}

void evaluate_merging_condition()
{
    int N_qp_e, N_qp_p, N_qp_g;
    auto N_qp_i = vector<int>(n_ion_populations);
    int pmerging_flag = 0; // bools are not supported by MPI

    MPI_Reduce(&(psr->N_qp_e), &N_qp_e, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(psr->N_qp_p), &N_qp_p, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(psr->N_qp_g), &N_qp_g, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(psr->N_qp_i, &(N_qp_i[0]), n_ion_populations, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        double crnp = crpc*xlength*ylength*zlength/(dx*dy*dz);
        pmerging_now = false;
        if (pmerging=="ti") {
            int N_qp;
            N_qp = N_qp_e + N_qp_p + N_qp_g;
            for (int i = 0; i < n_ion_populations; i++)
                N_qp += N_qp_i[i];
            if (N_qp>(3+n_ion_populations)*crnp) {
                pmerging_now = true;
                // portion of particles that will be deleted
                ppd[0] = (N_qp - (3+n_ion_populations)*crnp)/N_qp;
                cout<<"\t\033[36m"<<"ppd = "<<ppd[0]<<TERM_NO_COLOR;
            }
        } else if (pmerging=="nl") {
            bool merge = (N_qp_e>crnp)||(N_qp_p>crnp)||(N_qp_g>crnp);
            for (int i = 0; i < n_ion_populations; ++i)
                merge = merge || (N_qp_i[i]>crnp);
            if (merge) {
                pmerging_now = true;
                // portion of particles that will be deleted
                for (int i=0;i<3+n_ion_populations;i++)
                    ppd[i] = 0;
                if (N_qp_e>crnp)
                    ppd[0] = (N_qp_e - crnp)/N_qp_e;
                if (N_qp_p>crnp)
                    ppd[1] = (N_qp_p - crnp)/N_qp_p;
                if (N_qp_g>crnp)
                    ppd[2] = (N_qp_g - crnp)/N_qp_g;
                for (int i=0;i<n_ion_populations;i++) {
                    if (N_qp_i[i]>crnp)
                        ppd[3+i] = (N_qp_i[i] - crnp)/N_qp_i[i];
                }
                cout<<"\t\033[36m"<<"ppd =";
                for (int i=0;i<3+n_ion_populations;i++) {
                    cout<<' '<<ppd[i];
                }
                cout<<TERM_NO_COLOR;
            }
        }
        pmerging_flag = pmerging_now;
    }
    MPI_Bcast(&pmerging_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ppd, 3+n_ion_populations, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    pmerging_now = (pmerging_flag == 1);
}

void add_moving_window_particles()
{
    if (mpi_rank != n_sr-1) {
        return;
    }
    if ((l + 1) * dt < t_add_mw) {
        if (mwseed==1) {
            double n;
            n=1/(k0*k0);
            n *= lin_interpolation((nmw-1)*dx, ne_profile_x_coords, ne_profile_x_values);

            int_vector3d cell_pos;
            cell_pos.i = nx_sr[n_sr-1] - 3;
            int_vector3d v_npic;
            v_npic.i = xnpic;
            v_npic.j = ynpic;
            v_npic.k = znpic;
            for(int j=0;j<int(ylength/dy);j++)
            {
                for(int k=0;k<int(zlength/dz);k++)
                {
                    cell_pos.j=j;
                    cell_pos.k=k;
                    double zcell = k * dz;
                    double ycell = j * dy;
                    double ycenter = ylength / 2.0;
                    double zcenter = zlength / 2.0;
                    double zrel = zcell - zcenter;
                    double yrel = ycell - ycenter;
                    double r = sqrt(zrel * zrel + yrel * yrel);
                    double modifier = lin_interpolation(r, ne_profile_r_coords, ne_profile_r_values);
                    if (modifier != 0.0)
                    {
                        psr->fill_cell_by_particles(-1,cell_pos,v_npic,n*modifier);
                    }
                }
            }
        }
        if (mwseed_ions == "on") {
            double n;
            n=1/(k0*k0);
            n *= lin_interpolation((nmw-1)*dx, ne_profile_x_coords, ne_profile_x_values);

            int_vector3d cell_pos;
            cell_pos.i = nx_sr[n_sr-1] - 3;
            int_vector3d v_npic;
            v_npic.i = xnpic;
            v_npic.j = ynpic;
            v_npic.k = znpic;
            for(int j=0;j<int(ylength/dy);j++)
            {
                for(int k=0;k<int(zlength/dz);k++)
                {
                    cell_pos.j=j;
                    cell_pos.k=k;
                    double zcell = k * dz;
                    double ycell = j * dy;
                    double ycenter = ylength / 2.0;
                    double zcenter = zlength / 2.0;
                    double zrel = zcell - zcenter;
                    double yrel = ycell - ycenter;
                    double r = sqrt(zrel * zrel + yrel * yrel);
                    double modifier = lin_interpolation(r, ne_profile_r_coords, ne_profile_r_values);
                    if (modifier != 0.0)
                    {
                        psr->fill_cell_by_particles( 1 / (proton_mass * mw_mcr)
                                                   , cell_pos, v_npic, n * modifier);
                    }
                }
            }
        }
        
        film* tmp_p_film = p_last_film;
        while (tmp_p_film != 0)
        {
            psr->film(tmp_p_film->x0-dx*x0_sr[n_sr-1]-dx*nmw, tmp_p_film->x0+tmp_p_film->filmwidth-dx*x0_sr[n_sr-1]-dx*nmw,
                tmp_p_film->ne_y0/(1.11485e+13/lambda/lambda), tmp_p_film->ne_y1/(1.11485e+13/lambda/lambda),
                (ions == "on" || ions == "positrons") ? (ions == "on" ? 1 : 2) : 0,
                1/(proton_mass*tmp_p_film->mcr),
                tmp_p_film->gradwidth, tmp_p_film->y0, tmp_p_film->y1, tmp_p_film->z0, tmp_p_film->z1,
                tmp_p_film->T, tmp_p_film->vx, nelflow != 0 || nerflow != 0,
                tmp_p_film->xnpic_film, tmp_p_film->ynpic_film, tmp_p_film->znpic_film, true, tmp_p_film->gradwidth_y);
            tmp_p_film = tmp_p_film->prev;
        }
    }
}

void add_neutral_flow_particles(short direction, double neflow, double vflow, double Tflow, double mcrflow)
{
    static double xlflow = 1;
    static double xrflow = 1;
    
    double& xflow = direction > 0 ? xlflow : xrflow;
    if (xflow >= 1) {
        xflow -= 1;
        double n=neflow/xnpic; // in n_{cr}
        int_vector3d cell_pos, v_npic;
        v_npic.i = 1; // с xnpic приходится обходиться отдельно, см. ниже
        v_npic.j = ynpic;
        v_npic.k = znpic;
        for (int j=0; j<int(ylength/dy); j++) {
            for (int k=0; k<int(zlength/dz); k++) {
                cell_pos.j=j;
                cell_pos.k=k;
                for (int ii=0; ii<xnpic; ii++) {
                    double x0;
                    int index = direction > 0 ? 0 : n_sr - 1;

                    x0 = xrflow - float(ii)/xnpic;
                    if ( x0 >= 0 ) {
                        cell_pos.i = direction > 0 ? 3 : nx_sr[index] - 4;
                    } else {
                        cell_pos.i = direction > 0 ? 2 : nx_sr[index] - 3;
                        x0 += 1;
                    }
                    double tmp = (j * dy - ylength / 2) / (ylength / 2);
                    tmp = cos(0.5 * PI * tmp * tmp * tmp * tmp * tmp);
                    double tr_env = tmp * tmp;
                    tmp =  (k * dz - zlength / 2) / (zlength / 2);
                    tmp = cos(0.5 * PI * tmp * tmp * tmp * tmp * tmp);
                    tr_env *= tmp * tmp;
                    
                    if (mpi_rank == index) {
                        psr->fill_cell_by_particles(-1,cell_pos,v_npic, n * tr_env, direction * vflow/sqrt(1-vflow*vflow), 0, (direction > 0 ? x0 : 1-x0)-0.5,Tflow); // 0.5 - for a compensation in fill_cell... for xnpic = 1
                        if (ions == "on")
                            psr->fill_cell_by_particles(1/(proton_mass*mcrflow),cell_pos,v_npic, n * tr_env, direction * vflow/sqrt(1-vflow*vflow), 0, (direction > 0 ? x0 : 1-x0)-0.5,Tflow / (proton_mass * mcrflow));
                    }
                }
            }
        }
    }

    xflow += dt*vflow/dx;
}

void print_number_of_quasiparticles()
{
    // logging number of pseudoparticles in each thread
    int N_qp_e_total = 0;
    int N_qp_p_total = 0;
    int N_qp_g_total = 0;
    int N_qp_i_total = 0;

    int count = 3 + n_ion_populations;
    int data[count];
    data[0] = psr->N_qp_e;
    data[1] = psr->N_qp_p;
    data[2] = psr->N_qp_g;
    for (int j=0; j<n_ion_populations; j++) {
        data[3+j] = psr->N_qp_i[j];
    }

    if (mpi_rank == 0) {
        for (int i=0; i<n_sr; i++) {
            if (i != 0) {
                MPI_Recv(data, count, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            N_qp_e_total += data[0];
            N_qp_p_total += data[1];
            N_qp_g_total += data[2];
            cout << " Thread #" << i << ":\t" << data[0] << " e, " << data[1] << " p, " << data[2] << " g";
            if (n_ion_populations >= 0)
            {
                int N_qp_i = 0;
                for (int j=0; j<n_ion_populations; j++)
                {
                    N_qp_i += data[3+j];
                }
                N_qp_i_total += N_qp_i;
                cout << ", " << N_qp_i << " i";
            }
            cout << endl;
        }
        cout << "Total:  \t" << N_qp_e_total << " e, " << N_qp_p_total << " p, " << N_qp_g_total << " g";
        if (n_ion_populations >= 0)
        {
            cout << ", " << N_qp_i_total << " i";
        }
        cout << endl;
    } else {
        MPI_Send(data, count, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

void register_mpi_particle() {
    // register MPI_PARTICLE type to be able to send particles
    // struct particle has 10 doubles and 1 int
    const int items = 2;
    int blocklengths[items] = {10, 1};
    MPI_Datatype types[items] = {MPI_DOUBLE, MPI_INT};
    MPI_Aint offsets[items];
    offsets[0] = offsetof(particle, x);
    offsets[1] = offsetof(particle, trn);
    MPI_Datatype tmp_particle_type;
    MPI_Type_create_struct(items, blocklengths, offsets, types, &tmp_particle_type);
    MPI_Type_create_resized(tmp_particle_type, 0, sizeof(particle), &MPI_PARTICLE);
    MPI_Type_commit(&MPI_PARTICLE);
}

bool check_moving_window(int moving_window_iteration, int current_iteration) {
    switch (mwindow) {
    case moving_window::OFF:
        return false;
        break;
    case moving_window::ON:
        return (current_iteration + 1) * dt * mwspeed > moving_window_iteration * dx;
        break;
    case moving_window::AUTO:
        double max_w = psr->get_max_w();
        MPI_Allreduce(MPI_IN_PLACE, &max_w, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        double ratio;
        if (mpi_rank == n_sr-1) {
            double front_slice_max_w = psr->get_max_w(psr->get_nx() - nm - 2, psr->get_nx() - nm);
            if (max_w > 0) {
                ratio = front_slice_max_w / max_w;
            } else {
                ratio = 0.0;
            }
        }
        MPI_Bcast(&ratio, 1, MPI_DOUBLE, n_sr-1, MPI_COMM_WORLD);
        return (ratio > mwtolerance);
        break;
    }
    return false;
}

template <class T>
void resize_field(field3d<T> & field, const vector<int> & partition, const vector<int> & nx_sr_new) {
    field3d<T> previous;
    previous = move(field);
    field = field3d<T>(nx_sr_new[mpi_rank], ny_global, nz_global);

    // copying inside the same 
    int lower = max(x0_sr[mpi_rank], partition[mpi_rank]);
    int upper = min(x0_sr[mpi_rank] + nx_sr[mpi_rank], partition[mpi_rank] + nx_sr_new[mpi_rank]);

    for (int i = lower; i < upper; i++) {
        int new_index = i - partition[mpi_rank];
        int old_index = i - x0_sr[mpi_rank];
        
        for (int j = 0; j < ny_global; j++) {
            for (int k = 0; k < nz_global; k++) {
                field[new_index][j][k] = previous[old_index][j][k];
            }
        }
    }

    // sending to a different process
    for (int mod = 0; mod < 2; mod++) {
        if (mpi_rank % 2 == mod) {
            if (mpi_rank != 0) {
                int lower = max(x0_sr[mpi_rank] + nx_ich, partition[mpi_rank-1]);
                int upper = min(x0_sr[mpi_rank] + nx_sr[mpi_rank], partition[mpi_rank-1] + nx_sr_new[mpi_rank-1]);
                if (lower < upper) {
                    const int left = lower - x0_sr[mpi_rank];
                    const int size = 3 * (upper - lower) * ny_global * nz_global;
                    MPI_Send(previous[left][0], size, MPI_DOUBLE, mpi_rank-1, 0, MPI_COMM_WORLD);
                }
            }
            if (mpi_rank != n_sr-1) {
                int lower = max(x0_sr[mpi_rank], partition[mpi_rank+1]);
                int upper = min(x0_sr[mpi_rank] + nx_sr[mpi_rank] - nx_ich, partition[mpi_rank+1] + nx_sr_new[mpi_rank+1]);
                if (lower < upper) {
                    const int left = lower - x0_sr[mpi_rank];
                    const int size = 3 * (upper - lower) * ny_global * nz_global;
                    MPI_Send(previous[left][0], size, MPI_DOUBLE, mpi_rank+1, 0, MPI_COMM_WORLD);
                }
            }
        } else {
            if (mpi_rank != 0) {
                int lower = max(x0_sr[mpi_rank-1], partition[mpi_rank]);
                int upper = min(x0_sr[mpi_rank-1] + nx_sr[mpi_rank-1] - nx_ich, partition[mpi_rank] + nx_sr_new[mpi_rank]);
                if (lower < upper) {
                    const int left = lower - partition[mpi_rank];
                    const int size = 3 * (upper - lower) * ny_global * nz_global;
                    MPI_Recv(field[left][0], size, MPI_DOUBLE, mpi_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            if (mpi_rank != n_sr-1) {
                int lower = max(x0_sr[mpi_rank+1] + nx_ich, partition[mpi_rank]);
                int upper = min(x0_sr[mpi_rank+1] + nx_sr[mpi_rank+1], partition[mpi_rank] + nx_sr_new[mpi_rank]);
                if (lower < upper) {
                    const int left = lower - partition[mpi_rank];
                    const int size = 3 * (upper - lower) * ny_global * nz_global;
                    MPI_Recv(field[left][0], size, MPI_DOUBLE, mpi_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
    }
}

void resize_particles(const vector<int> & partition, const vector<int> & nx_sr_new) {
    field3d<cellp> particles_previous;
    particles_previous = move(psr->cp);
    psr->cp = field3d<cellp>(nx_sr_new[mpi_rank], ny_global, nz_global);

    // sending to a different process
    for (int mod = 0; mod < 2; mod++) {
        if (mpi_rank % 2 == mod) {
            if (mpi_rank != 0) {
                int lower = max(x0_sr[mpi_rank] + nx_ich, partition[mpi_rank-1]);
                int upper = min(x0_sr[mpi_rank] + nx_sr[mpi_rank], partition[mpi_rank-1] + nx_sr_new[mpi_rank-1]);
                if (lower < upper) {
                    const int left = lower - x0_sr[mpi_rank];
                    send_particle_slice(particles_previous, left, upper - lower, mpi_rank - 1);
                }
            }
            if (mpi_rank != n_sr-1) {
                int lower = max(x0_sr[mpi_rank], partition[mpi_rank+1]);
                int upper = min(x0_sr[mpi_rank] + nx_sr[mpi_rank] - nx_ich, partition[mpi_rank+1] + nx_sr_new[mpi_rank+1]);
                if (lower < upper) {
                    const int left = lower - x0_sr[mpi_rank];
                    send_particle_slice(particles_previous, left, upper - lower, mpi_rank + 1);
                }
            }
        } else {
            if (mpi_rank != 0) {
                int lower = max(x0_sr[mpi_rank-1], partition[mpi_rank]);
                int upper = min(x0_sr[mpi_rank-1] + nx_sr[mpi_rank-1] - nx_ich, partition[mpi_rank] + nx_sr_new[mpi_rank]);
                if (lower < upper) {
                    const int left = lower - partition[mpi_rank];
                    receive_particle_slice(left, upper - lower, mpi_rank - 1, x0_sr[mpi_rank - 1] - partition[mpi_rank]);
                }
            }
            if (mpi_rank != n_sr-1) {
                int lower = max(x0_sr[mpi_rank+1] + nx_ich, partition[mpi_rank]);
                int upper = min(x0_sr[mpi_rank+1] + nx_sr[mpi_rank+1], partition[mpi_rank] + nx_sr_new[mpi_rank]);
                if (lower < upper) {
                    const int left = lower - partition[mpi_rank];
                    receive_particle_slice(left, upper - lower, mpi_rank + 1, x0_sr[mpi_rank + 1] - partition[mpi_rank]);
                }
            }
        }
    }

    // copying inside the same 
    int lower = max(x0_sr[mpi_rank], partition[mpi_rank]);
    int upper = min(x0_sr[mpi_rank] + nx_sr[mpi_rank], partition[mpi_rank] + nx_sr_new[mpi_rank]);

        for (int i = lower; i < upper; i++) {
        int new_index = i - partition[mpi_rank];
        int old_index = i - x0_sr[mpi_rank];
        
        for (int j = 0; j < ny_global; j++) {
            for (int k = 0; k < nz_global; k++) {
                psr->cp[new_index][j][k] = particles_previous[old_index][j][k];
                psr->cp[new_index][j][k].pl.xplus(new_index - old_index);
            }
        }
    }

    // cleaning unused particles

    lower -= x0_sr[mpi_rank];
    upper -= x0_sr[mpi_rank];

    for (int i = 0; i < lower; i++) {
        for (int j = 0; j < ny_global; j++) {
            for (int k = 0; k < nz_global; k++) {
                psr->erase(particles_previous[i][j][k].pl);
            }
        }
    }
    for (int i = upper; i < nx_sr[mpi_rank]; i++) {
        for (int j = 0; j < ny_global; j++) {
            for (int k = 0; k < nz_global; k++) {
                psr->erase(particles_previous[i][j][k].pl);
            }
        }
    }
}

void resize_regions(const vector<int> & partition) {
    assert(partition.size() == x0_sr.size());

    size_t size = partition.size();

    std::vector<int> nx_sr_new(size);

    for (size_t i = 0; i < size - 1; i++) {
        nx_sr_new[i] = partition[i+1] + nx_ich - partition[i];
    }
    nx_sr_new[size-1] = nx_global - partition[size-1];

    if ((nx_sr_new[mpi_rank] != nx_sr[mpi_rank]) || (x0_sr[mpi_rank] != partition[mpi_rank])) {
        resize_field(psr->ce, partition, nx_sr_new);
        resize_field(psr->cb, partition, nx_sr_new);
        resize_field(psr->cbe, partition, nx_sr_new);
        resize_particles(partition, nx_sr_new);
        psr->cj = field3d<cellj>(nx_sr_new[mpi_rank], ny_global, nz_global);
        for (int n = 0; n < n_ion_populations; n++) {
            psr->irho[n] = field3d<double>(nx_sr_new[mpi_rank], ny_global, nz_global);
        }    

        psr->set_nx(nx_sr_new[mpi_rank]);
    }

    x0_sr = partition;
    nx_sr = nx_sr_new;
}

void load_balancing() {
    bool need_to_balance = false;
    
    auto global_weights = calculate_global_layer_weights();
    
    double initial_imbalance;
    if (mpi_rank == 0) {
        initial_imbalance = calculate_partition_imbalance(global_weights, x0_sr, nx_ich);

        need_to_balance = (initial_imbalance > balancing_threshold);
    }

    MPI_Bcast(&need_to_balance, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);

    if (need_to_balance) {
        std::vector<int> new_partition(n_sr);

        if (mpi_rank == 0) {
            auto optimal_partition = calculate_optimal_partition(global_weights, n_sr, nx_ich);
            new_partition = normalize_new_partition(x0_sr, optimal_partition, nx_ich);
            auto normalized_imbalance = calculate_partition_imbalance(global_weights, new_partition, nx_ich);

            cout << "Load balancing: imbalance " << initial_imbalance << ", after balancing " << normalized_imbalance << endl;
        }

        MPI_Bcast(new_partition.data(), new_partition.size(), MPI_INT, 0, MPI_COMM_WORLD);

        resize_regions(new_partition);
    }
}

int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_sr);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if (mpi_rank == 0) {
        cout<<"\n\033[1m"<<"hi!"<<"\033[0m\n"<<endl;
        cout <<"Number of MPI processes: " << n_sr << endl;
    }

    char hostname[MPI_MAX_PROCESSOR_NAME];
    int hostname_len;
    MPI_Get_processor_name(hostname, &hostname_len);
    unsigned long pid = getpid();

    if (mpi_rank == 0) {
        for (int i=0; i<n_sr; i++) {
            if (i != 0) {
                MPI_Recv(&pid, 1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Status status;
                MPI_Probe(i, 1, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_CHAR, &hostname_len);
                MPI_Recv(hostname, hostname_len, MPI_CHAR, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            cout << "Quill process rank = " << i << ", id = " << pid << ", hostname = " << hostname << endl;
        }
    } else {
        MPI_Send(&pid, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
        MPI_Send(hostname, hostname_len + 1, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
    }

    register_mpi_particle();

    MPI_Barrier(MPI_COMM_WORLD);

    up_time = times(&tms_struct);
    start_time = times(&tms_struct);
    inaccurate_time = time(NULL);

    nx_ich = 8;
    nm = nx_ich/2;
    file_name_accuracy = 100;

    if (init()==1) {
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }

    nx_global = (int) (xlength / dx);
    ny_global = (int) (ylength / dy);
    nz_global = (int) (zlength / dz);

    ofstream fout_log;
    if (mpi_rank == 0) {
        fout_log.open(data_folder+"/log",ios::app); // ios:app mode - append to the file
        fout_log<<"start time: "<<ctime(&inaccurate_time);
    }

    ppd = new double[3+n_ion_populations];

    const int nx = nx_global + nx_ich*(n_sr-1); // nx = nx_sr*n_sr - nx_ich*(n_sr-1);
    nx_sr = vector<int>(n_sr, nx / n_sr);
    for (int i=0; i < (nx % n_sr); i++) {
        nx_sr[i]++;
    }

    x0_sr = vector<int>(n_sr);
    for (int i=1; i<n_sr; i++) {
        x0_sr[i] = x0_sr[i-1] + nx_sr[i-1] - nx_ich;
    }

    if(nx_ich*(n_sr-1)>=nx_global) {
        cout<<TERM_RED<<"main: too many slices, aborting..."<<TERM_NO_COLOR<<endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    main_thread_time = times(&tms_struct);
    if (mpi_rank == 0) {
        cout<<"Creating arrays..."<<flush;
    }
    psr = unique_ptr<spatial_region>(new spatial_region());
    psr->init(mpi_rank,dx,dy,dz,dt,lambda/2.4263086e-10,xnpic,ynpic,znpic,n_ion_populations,icmr,data_folder,solver,pusher);
    psr->create_arrays(nx_sr[mpi_rank],ny_global,nz_global,mpi_rank+times(&tms_struct));
    
    init_fields();

    if (beam=="on")
    {
        init_beam();
    }

    init_films();

    psr->f_init_boundaries();
    psr->interpolate_be();

    MPI_Barrier(MPI_COMM_WORLD);

    main_thread_time = times(&tms_struct) - main_thread_time;
    seconds = main_thread_time/100.0;

    if (mpi_rank == 0) {
        cout<<"done!"<<endl;
        fout_log<<"fill arrays: "<<seconds<<"s"<<endl;
    }

    l=0;

    ofstream fout_N;
    ofstream fout_energy;
    ofstream fout_energy_deleted;
    ofstream fout_mwcoordinate;
    if (mpi_rank == 0) {
        fout_N.open(data_folder+"/N");
        fout_energy.open(data_folder+"/energy");

        if (catching_enabled || dump_photons) {
            fout_energy_deleted.open(data_folder+"/energy_deleted");
        }

        if (mwindow != moving_window::OFF) {
            fout_mwcoordinate.open(data_folder+"/mwcoordinate");
        }
    }

    while(l<p_last_ddi->t_end/dt)
    {
        if (pmerging_now)
            psr->pmerging(ppd,pmerging);
        int i = mpi_rank;
        psr->birth_from_vacuum(8*PI*PI/(dx*dy*dz)*2.818e-13/lambda); // 2.818e-13 = e^2/mc^2
        psr->padvance(external_bz);
        psr->compute_N(nm*(i!=0),nm*(i!=n_sr-1),dx*dy*dz*1.11485e13*lambda/(8*PI*PI*PI));
        psr->compute_energy(nm*(i!=0),nm*(i!=n_sr-1),0.5*dx*dy*dz*3.691e4*lambda/1e7,8.2e-14*dx*dy*dz*1.11485e13*lambda/(8*PI*PI*PI)); // энергия в Джоулях
        if (n_tracks > 0) {
            psr->fout_tracks((x0_sr[i]+nmw)*dx/2/PI,nm);
        }
        psr->fadvance();

        bool is_moving_window_iteration = check_moving_window(nmw, l);
        if (is_moving_window_iteration) {
            psr->moving_window();
            nmw++;
        }

        synchronize_regions(is_moving_window_iteration);

        /* вывод плотности, спектра и 'phasespace'-данных для фотонов,
           электронов и позитронов в файлы */

        if(l*dt >= [](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi))
        {
            if (write_jx || write_jy || write_jz) {
                write_density(write_jx, write_jy, write_jz, "jx", "jy", "jz", false, true);
            }
            
            psr->compute_rho();

            write_density(true, write_p, write_ph, "rho", "rho_p", "rho_ph", true);

            write_spectrum_phasespace(write_p, write_ph);            
            
            if (catching_enabled || dump_photons)
            {
                write_deleted_particles(write_p, write_ph);
            }

            write_layer_weights();
        }

        if (balancing_enabled && (l % balancing_every == 0)) {
            load_balancing();
        }


        if (!particles_to_track.empty() && l == int(tr_start/dt)) 
        {
            start_tracking();
        }


        write_N(fout_N);
        write_energy(fout_energy);
        if (mwindow != moving_window::OFF) {
            write_mwcoordinate(fout_mwcoordinate);
        }
        if (catching_enabled || dump_photons)
        {
            write_energy_deleted(fout_energy_deleted);
        }

        if (mpi_rank == 0) {
            cout<<"\033[33m"<<"ct/lambda = "<<l*dt/2/PI<<"\tstep# "<<l<<TERM_NO_COLOR;
        }

        evaluate_merging_condition();

        if (mpi_rank == 0) {
            cout << endl;
        }

        if (is_moving_window_iteration) {
            add_moving_window_particles();
        }


        if (nelflow != 0)
        {
            // neutral flow traveling from the left to the right
            add_neutral_flow_particles(1, nelflow, vlflow, Tlflow, mcrlflow);
        }
        if (nerflow != 0)
        {
            // neutral flow traveling from the right to the left
            add_neutral_flow_particles(-1, nerflow, vrflow, Trflow, mcrrflow);
        }


        // вывод данных в файлы (продолжение)
        if(l*dt>=[](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi))
        {
            if (verbose_logging)
            {
                print_number_of_quasiparticles();
            }
            
            if (mpi_rank == 0) {
                cout<<"output# "<<"\033[1m"<<int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy<<TERM_NO_COLOR<<flush;
            }

            write_fields(); // Ex..Ez, Bx..Bz, w (field energy density), inv (E^2-B^2 - relativistic invariant)
            
            if (mpi_rank == 0) {
                cout<<"\t"<<"up: "<<(times(&tms_struct)-start_time)/100.0<<" s"<<endl;
            }
            p_current_ddi->f++;
        }

        if (l*dt>=p_current_ddi->t_end)
        {
            p_current_ddi = p_current_ddi->next;
        }

        l++;
    }

    delete[] icmr;

    delete[] ppd;

    if (mpi_rank == 0) {
        fout_N.close();
        fout_energy.close();
        if (fout_energy_deleted.is_open())
        {
            fout_energy_deleted.close();
        }
        if (fout_mwcoordinate.is_open()) {
            fout_mwcoordinate.close();
        }
    }

    while(p_last_ddi!=0)
    {
        ddi* tmp = p_last_ddi->prev;
        delete p_last_ddi;
        p_last_ddi = tmp;
    }

    if (mpi_rank == 0) {
        inaccurate_time = time(NULL);
        fout_log<<"stop time: "<<ctime(&inaccurate_time);
        up_time = times(&tms_struct) - up_time;
        seconds = up_time/100.0;
        fout_log<<"uptime: "<<seconds<<"s"<<endl;
        fout_log.close();

        cout<<"\n\033[1mbye!\033[0m\n"<<endl;
    }

    MPI_Finalize();
    return 0;
}

double convert_units(double value, string initial_units, string desired_units) {
    if (desired_units == "1/k") {
        if (initial_units == "um") {
            value *= 2 * PI * 1e-4 / lambda;
        } else if (initial_units == "lambda") {
            value *= 2 * PI;
        } else if (initial_units == "1/k") {
        } else {
            cout << TERM_RED << "Unknown units: " << initial_units << " for conversion to " << desired_units
                    << TERM_NO_COLOR << endl;
        }
    } else if (desired_units == "") {
        if (initial_units != "") {
            cout << TERM_RED << "Cannot convert [" << initial_units << "] to unitless" << TERM_NO_COLOR << endl;
        }
    } else {
        cout << TERM_RED << "Unknown desired units: " << desired_units << TERM_NO_COLOR << endl;
    }
    return value;
}

vector<double> find_array(var * element, string name, string desired_units, string default_units = "",
        vector<double> default_array = { }) {
    var * current = find(name, element);
    string units = current->units == "#" ? default_units : current->units;

    if (!current->input_array.empty()) {
        vector<double> array = current->input_array;
        for (double & v : array) {
            v = convert_units(v, units, desired_units);
        }
        return array;
    } else {
        return default_array;
    }
}

bool find_boolean(var * element, string name, bool default_value) {
    var * current = find(name, element);
    if (current->units == "on" || current->units == "true") {
        return true;
    } else if (current->units == "off" || current->units == "false") {
        return false;
    } else {
        if (current->units != "") {
            cout << TERM_RED << "Value [" << current->units << "] for [" << name << "] is incorrect, using the default value" << TERM_NO_COLOR << endl;
        }
        return default_value;
    }
}

double find_double(var * element, string name, double default_value) {
    var * current = find(name, element);
    if (current->value == 0) {
        return default_value;
    } else {
        return current->value;
    }
}

int find_int(var * element, string name, int default_value) {
    var * current = find(name, element);
    if (current->value == 0) {
        return default_value;
    } else {
        return current->value;
    }
}

int init()
{
    if (mpi_rank == 0)
        cout<<"Initialization"<<'\n'<<flush;
    var* first;
    var* current;
    var* tmp;
    first = new var;
    current = first;
    while (current->read()!=1)
    {
        current->next = new var;
        current = current->next;
    }

    //
    current = find("dx",first);
    dx = current->value*2*PI; // from lambda to c/\omega
    current = find("dy",first);
    dy = current->value*2*PI;
    current = find("dz",first);
    dz = current->value*2*PI;
    current = find("dt",first);
    dt = current->value*2*PI;
    current = find("lambda",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4;
        current->units="cm";
    }
    lambda = current->value;
    current = find("ne",first);
    if (current->units=="ncr")
    {
        current->value = current->value*1.11485e+13/lambda/lambda;
        current->units = "cm^{-3}";
    }
    ne = current->value;
    current = find("lambda",first);
    if (current->units=="lambda_p")
    {
        current->value = current->value*sqrt(1.11485e13/ne);
        current->units="cm";
        lambda = current->value;
    }
    k0 = sqrt(1.11485e+13/lambda/lambda/ne);
    if (dt==0) {
        dt = 0.5*k0*sqrt(4.04);
        dx = 0;
    }
    if (dx==0)
    {
        // dx = dt/( 1 - 1/(k0*k0)*dt*dt/4.04 );
        dx = (1+1e-4)*dt/( 1 - 1/(k0*k0)*dt*dt/4.04 );
        if (dx < 0) {
            cout << TERM_RED << "Cannot calculate stable dx, reduce dt. Aborting..." << TERM_NO_COLOR << endl;
            return 1;
        }
    }

    //
    current = find("xlength",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    if (current->units=="fs")
    {
        current->value = current->value*1e-15*2.99792458e10/lambda;
        current->units="lambda";
    }
    xlength = current->value*2*PI;
    current = find("ylength",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    ylength = current->value*2*PI;
    current = find("zlength",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    zlength = current->value*2*PI;
    var* tmp_last_t_end = first;
    p_last_ddi = 0;
    do {
        double t_end;
        double output_period;
        current = find("t_end",tmp_last_t_end);
        tmp_last_t_end = current->next;
        if (current->units=="um")
        {
            current->value = current->value*1e-4/lambda;
            current->units="lambda";
        }
        if (current->units=="mm")
        {
            current->value = current->value/10/lambda;
            current->units="lambda";
        }
        if (current->units=="cm")
        {
            current->value = current->value/lambda;
            current->units="lambda";
        }
        t_end = current->value*2*PI + dt;
        current = find("output_period",tmp_last_t_end);
        if (current->units=="um")
        {
            current->value = current->value*1e-4/lambda;
            current->units="lambda";
        }
        if (current->units=="mm")
        {
            current->value = current->value/10/lambda;
            current->units="lambda";
        }
        if (current->units=="cm")
        {
            current->value = current->value/lambda;
            current->units="lambda";
        }
        if (current->units=="t_end") {
            current-> value = current->value*(t_end-dt)/2/PI;
        }
        output_period = current->value*2*PI;
        p_current_ddi = p_last_ddi;
        p_last_ddi = new ddi;
        if (p_current_ddi!=0)
            p_current_ddi->next = p_last_ddi;
        p_last_ddi->next = 0;
        p_last_ddi->prev = p_current_ddi;
        p_last_ddi->t_end = t_end;
        p_last_ddi->output_period = output_period;
        p_last_ddi->f = 0;
    } while(find("t_end",tmp_last_t_end)->value!=0);
    //
    current = find("t_add_mw", first);
    if (current->units == "um") {
        current->value = current->value * 1e-4 / lambda;
        current->units = "lambda";
    } else if (current->units == "mm") {
        current->value = current->value / (10 * lambda);
        current->units = "lambda";
    } else if (current->units == "cm") {
        current->value = current->value / lambda;
        current->units="lambda";
    }
    if (current->value != 0) {
        t_add_mw = current->value * 2 * PI;
    } else {
        t_add_mw = p_last_ddi->t_end;
    }
    //
    current = find("xsigma",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    if (current->units=="fs")
    {
        current->value = current->value*1e-15*2.99792458e10/lambda;
        current->units="lambda";
    }
    xsigma = current->value*2*PI;
    current = find("ysigma",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    ysigma = current->value*2*PI;
    current = find("zsigma",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    zsigma = current->value*2*PI;
    current = find("x0",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    if (current->units=="fs")
    {
        current->value = current->value*1e-15*2.99792458e10/lambda;
        current->units="lambda";
    }
    x0 = current->value*2*PI;
    current = find("x0fout",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    if (current->units=="fs")
    {
        current->value = current->value*1e-15*2.99792458e10/lambda;
        current->units="lambda";
    }
    x0fout = current->value*2*PI;
    if (x0fout==0) x0fout = xlength/2; // default value
    a0y = 0;
    a0z = 0;
    current = find("P",first);
    if (current->units=="PW")
    {
        current->value = current->value*1e22;
        current->units="cgs";
    }
    if (current->units=="TW")
    {
        current->value = current->value*1e19;
        current->units="cgs";
    }
    if (current->value!=0)
    {
        a0y = 5.27281e+17*sqrt(2*current->value/(PI*2.42161e+52*ysigma*zsigma));
    }
    current = find("I",first);
    if (current->units=="W/cm^2")
    {
        current->value = current->value*1e7;
        current->units="cgs";
    }
    if (current->value!=0)
    {
        a0y = 5.27281e+17*lambda*sqrt(current->value/(PI*2.42161e+52));
    }
    current = find("W",first);
    if (current->units=="J")
    {
        current->value = current->value*1e7;
        current->units="cgs";
    }
    if (current->value!=0)
    {
        // 3.691e4 = (mc^2)^2/(8 pi^2 e^2)
        a0y = sqrt(current->value/(3.691e4*lambda*PI*sqrt(PI/2)*xsigma*ysigma*zsigma/4));
        current = find("polarization",first);
        if (current->units=="circular")
        {
            a0y = a0y/sqrt(2);
            a0z = a0y;
        }
    }
    current = find("a0",first);
    if (current->value!=0)
    {
        a0y = current->value;
    }
    std::string polarization = "linear"; // default
    current = find("polarization",first);
    if (current->units=="circular")
    {
        a0z = a0y;
        polarization = current->units;
    }
    if (current->units=="elliptic") {
        a0z = current->value*a0y;
        polarization = current->units;
    }
    current = find("a0y",first);
    if (current->value!=0)
    {
        a0y = current->value;
        var* a = find("a0z",first);
        a0z = a->value;
    }
    else
    {
        var* a = find("a0z",first);
        if (a->value!=0)
        {
            a0z = a->value;
            a0y = 0;
        }
    }
    current = find("external_bz",first);
    external_bz = current->value;
    current = find("f_envelope",first);
    f_envelope = current->units;
    if (f_envelope=="") {
        f_envelope = "cos";
        sscos = 0;
    } else if (f_envelope == "sscos") {
        f_envelope = "cos";
        sscos = 1;
    } else if (f_envelope == "pearl") {
        f_envelope = "cos";
        sscos = 2;
    } else if (f_envelope == "tophat") {
        f_envelope = "cos";
        sscos = 3;
    }
    current = find("b_sign",first);
    if (current->value!=-1) b_sign = 1;
    else b_sign = 0;
    current = find("phase",first);
    if (current->units=="pi") current->value = current->value*PI;
    phase = current->value;
    current = find("phi",first);
    if (current->units=="pi") current->value = current->value*PI;
    phi = current->value;
    current = find("phi_rotate",first);
    if (current->units=="pi") current->value = current->value*PI;
    phi_rotate = current->value;
    current = find("y0",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    y00 = current->value*2*PI;
    current = find("z0",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    z00 = current->value*2*PI;

    // positioning the pulse with coordinates (r0, theta, z0) instead of (x0, y0, z0)
    current = find("r0",first);  
    if (current->units != "off" && current->value != 0)  // means that r0 is found (r0 is not zero!)
    {
        if (current->units=="um")
        {
            current->value = current->value*1e-4/lambda;
            current->units="lambda";
        }
        double r0 = current->value*2*PI;
        current = find("theta",first);
        if (current->units=="deg")
        {
            current->value = current->value*PI/180;
            current->units="rad";
        }
        double theta = current->value;
        x0 = -r0*cos(theta) + 0.00001;  // Bug: Quill cannot create laser pulse exactly at xlength/2
        y00 = -r0*sin(theta);
    }
    
    // coordinates (relative to the center of the box) of the point which the laser pulse is proparating to
    current = find("xtarget",first);
    if (current->units != "off" && current->value != 0) 
    {
        if (current->units=="um")
        {
            current->value = current->value*1e-4/lambda;
            current->units="lambda";
        }
        xtarget = current->value * 2*PI;
    }
    current = find("ytarget",first);
    if (current->units != "off" && current->value != 0) 
    {
        if (current->units=="um")
        {
            current->value = current->value*1e-4/lambda;
            current->units="lambda";
        }
        ytarget = current->value * 2*PI;
    }
    current = find("ztarget",first);
    if (current->units != "off" && current->value != 0) 
    {
        if (current->units=="um")
        {
            current->value = current->value*1e-4/lambda;
            current->units="lambda";
        }
        ztarget = current->value * 2*PI;
    }

    if ( f_envelope == "focussedSSC" ) {
        f_envelope = "focused";
        sscos = 1;
    }
    if ( f_envelope == "focused" )
    {
        x00 = sqrt(x0*x0+y00*y00+z00*z00);
        /*sigma = 0.5*(ysigma+zsigma);
          sigma0 = sigma*sigma*sigma*sigma/4 - 4*x00*x00;
          if (sigma0<0)
          {
          cout<<"\n\033[31m"<<"main: improper focused pulse, aborting..."<<TERM_NO_COLOR<<endl;
          return 1;
          }
          sigma0 = sqrt( sigma*sigma/2 - sqrt(sigma0) );*/
        sigma0 = 0.5*(ysigma+zsigma);
        sigma = sqrt( sigma0 * sigma0 + 4 * x00 * x00 / ( sigma0 * sigma0 ) );
        a0y = a0y*sigma0/sigma;
        a0z = a0z*sigma0/sigma;
    }
    current = find("shenergy", first);
    shenergy = current->value;
    current = find("shphase",first);
    if (current->units=="pi")
        current->value = current->value * PI;
    shphase = current->value;
    current = find("mwindow",first);
    mwindow = moving_window::ON;
    if (current->units == "off") {
        mwindow = moving_window::OFF;
    } else if (current->units == "auto") {
        mwindow = moving_window::AUTO;
    }
    current = find("mwtolerance", first);
    mwtolerance = current->value;
    if (mwtolerance == 0.0) {
        mwtolerance = 1e-6;
    }
    current = find("mwspeed",first);
    mwspeed = current->value;
    if (mwindow==moving_window::ON && mwspeed==0)
        mwspeed = 1;
    current = find("mwseed",first);
    mwseed = 1;
    if (current->units=="off") mwseed = 0;
    current = find("mwseed_ions", first);
    mwseed_ions = current->units;
    if (mwseed_ions == "")
        mwseed_ions = "off"; // default
    current = find("mw_mcr", first);
    mw_mcr = current->value;
    if (mw_mcr == 0)
        mw_mcr = 1; // default
    current = find("e_components_for_output", first);
    e_components_for_output = current->units;
    current = find("b_components_for_output", first);
    b_components_for_output = current->units;
    current = find("j_components_for_output", first);
    j_components_for_output = current->units;
    if (e_components_for_output == "")
        e_components_for_output = "none"; // default
    if (b_components_for_output == "")
        b_components_for_output = "none"; // default
    if (j_components_for_output == "")
        j_components_for_output = "none";

    if (j_components_for_output.find_first_of('x') != string::npos) {
        write_jx = true;
    }
    if (j_components_for_output.find_first_of('y') != string::npos) {
        write_jy = true;
    }
    if (j_components_for_output.find_first_of('z') != string::npos) {
        write_jz = true;
    }
    if ((j_components_for_output.find_first_not_of("xyz") != string::npos)
            and (j_components_for_output != "none")) {
        cerr << "WARNING: j_components_for_output=[" << j_components_for_output
                << "] is not a valid option, only x, y, z chars are possible"
                << endl;
    }

    current = find("beam",first);
    beam = current->units;
    if (beam=="") beam = "off";
    current = find("beam_particles",first);
    beam_particles = current->units;
    if (beam_particles=="") beam_particles = "e";
    current = find("Nb",first);
    Nb = current->value;
    current = find("epsb",first);
    if (current->units=="mc^2")
    {
        current->value = current->value*0.511;
        current->units = "MeV";
    }
    epsb = current->value; // в программе epsb - в MeV
    current = find("xb",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    if (current->units=="fs")
    {
        current->value = current->value*1e-15*2.99792458e10/lambda;
        current->units="lambda";
    }
    xb = current->value*2*PI;
    current = find("rb",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    rb = current->value*2*PI;
    current = find("x0b",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    if (current->units=="fs")
    {
        current->value = current->value*1e-15*2.99792458e10/lambda;
        current->units="lambda";
    }
    x0b = current->value*2*PI;
    current = find("y0b",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    if (current->units=="fs")
    {
        current->value = current->value*1e-15*2.99792458e10/lambda;
        current->units="lambda";
    }
    y0b = current->value*2*PI;
    current = find("phib",first);
    if (current->units=="deg")
    {
        current->value = current->value*PI/180;
        current->units="rad";
    }
    phib = current->value;
    current = find("ions",first);
    ions = current->units;
    p_last_film = 0;
    film* tmp_p_film;
    std::string tmp_s;
    tmp = first;
    do
    {
        current = find("film",tmp);
        tmp = current->next;
        tmp_s = current->units;
    } while (tmp_s!="on"&&tmp!=0);
    while (tmp_s=="on")
    {
        tmp_p_film = p_last_film;
        p_last_film = new film;
        p_last_film->prev = tmp_p_film;
        current = find("x0film",tmp);
        if (current->units=="um")
        {
            current->value = current->value*1e-4/lambda;
            current->units="lambda";
        }
        if (current->units=="fs")
        {
            current->value = current->value*1e-15*2.99792458e10/lambda;
            current->units="lambda";
        }
        if (current->units == "dx") {
            current->value *= dx / 2 / PI;
            current->units = "lambda";
        }
        p_last_film->x0 = current->value*2*PI;
        current = find("filmwidth",tmp);
        if (current->units=="um")
        {
            current->value = current->value*1e-4/lambda;
            current->units="lambda";
        }
        if (current->units=="fs")
        {
            current->value = current->value*1e-15*2.99792458e10/lambda;
            current->units="lambda";
        }
        p_last_film->filmwidth = current->value*2*PI;
        current = find("gradwidth",tmp);
        if (current->units=="um")
        {
            current->value = current->value*1e-4/lambda;
            current->units="lambda";
        }
        if (current->units=="fs")
        {
            current->value = current->value*1e-15*2.99792458e10/lambda;
            current->units="lambda";
        }
        p_last_film->gradwidth = current->value*2*PI;
        current = find("gradwidth_y",tmp);
        if (current->units=="um")
        {
            current->value = current->value*1e-4/lambda;
            current->units="lambda";
        }
        if (current->units=="fs")
        {
            current->value = current->value*1e-15*2.99792458e10/lambda;
            current->units="lambda";
        }
        p_last_film->gradwidth_y = current->value*2*PI;
        current = find("y0film",tmp);
        if (current->units=="um")
        {
            current->value = current->value*1e-4/lambda;
            current->units="lambda";
        }
        if (current->units=="fs")
        {
            current->value = current->value*1e-15*2.99792458e10/lambda;
            current->units="lambda";
        }
        p_last_film->y0 = current->value*2*PI;
        current = find("y1film",tmp);
        if (current->units=="um")
        {
            current->value = current->value*1e-4/lambda;
            current->units="lambda";
        }
        if (current->units=="fs")
        {
            current->value = current->value*1e-15*2.99792458e10/lambda;
            current->units="lambda";
        }
        if (current->value == 0)
            p_last_film->y1 = ylength;
        else
            p_last_film->y1 = current->value*2*PI;
        current = find("z0film",tmp);
        if (current->units=="um")
        {
            current->value = current->value*1e-4/lambda;
            current->units="lambda";
        }
        if (current->units=="fs")
        {
            current->value = current->value*1e-15*2.99792458e10/lambda;
            current->units="lambda";
        }
        p_last_film->z0 = current->value*2*PI;
        current = find("z1film",tmp);
        if (current->units=="um")
        {
            current->value = current->value*1e-4/lambda;
            current->units="lambda";
        }
        if (current->units=="fs")
        {
            current->value = current->value*1e-15*2.99792458e10/lambda;
            current->units="lambda";
        }
        if (current->value == 0)
            p_last_film->z1 = zlength;
        else
            p_last_film->z1 = current->value*2*PI;
        current = find("nfilm",tmp);
        if (current->units=="ncr")
        {
            current->value = current->value*1.11485e+13/lambda/lambda;
            current->units = "cm^{-3}";
        }
        else if (current->units=="ne")
        {
            current->value = current->value*ne;
            current->units = "cm^{-3}";
        }
        p_last_film->ne = current->value;
        current = find("nfilm_y0",tmp);
        if (current->units=="ncr")
        {
            current->value = current->value*1.11485e+13/lambda/lambda;
            current->units = "cm^{-3}";
        }
        else if (current->units=="ne")
        {
            current->value = current->value*ne;
            current->units = "cm^{-3}";
        }
        p_last_film->ne_y0 = current->value;
        current = find("nfilm_y1",tmp);
        if (current->units=="ncr")
        {
            current->value = current->value*1.11485e+13/lambda/lambda;
            current->units = "cm^{-3}";
        }
        else if (current->units=="ne")
        {
            current->value = current->value*ne;
            current->units = "cm^{-3}";
        }
        p_last_film->ne_y1 = current->value;
        current = find("mcr",tmp);
        p_last_film->mcr = (current->value == 0 ? 1 : current->value) ;
        current = find("Tfilm",tmp);
        p_last_film->T = current->value;
        current = find("vxfilm",tmp);
        p_last_film->vx = current->value;
        current = find("xnpic_film",tmp);
        p_last_film->xnpic_film = current->value;
        current = find("ynpic_film",tmp);
        p_last_film->ynpic_film = current->value;
        current = find("znpic_film",tmp);
        p_last_film->znpic_film = current->value;
        do
        {
            current = find("film",tmp);
            tmp = current->next;
            tmp_s = current->units;
        } while (tmp_s!="on"&&tmp!=0);
    }
    // neutral flows
    current = find("nelflow",first);
    nelflow = current->value*ne/(1.11485e+13/lambda/lambda);
    current = find("vlflow",first);
    vlflow = current->value;
    current = find("mcrlflow",first);
    mcrlflow = current->value;
    current = find("Tlflow",first);
    Tlflow = current->value;
    current = find("nerflow",first);
    nerflow = current->value*ne/(1.11485e+13/lambda/lambda);
    current = find("vrflow",first);
    vrflow = current->value;
    current = find("mcrrflow",first);
    mcrrflow = current->value;
    current = find("Trflow",first);
    Trflow = current->value;
    if (nelflow!=0 || nerflow!=0)
        mwindow = moving_window::OFF;
    //
    if (ions=="on" || ions=="positrons")
    { // counting of ion populations
        int n;
        n = 0;
        tmp_p_film = p_last_film;
        while (tmp_p_film!=0)
        {
            n++;
            tmp_p_film = tmp_p_film->prev;
        }
        if (nelflow!=0)
            n++;
        if (nerflow!=0)
            n++;
        if (mwindow == moving_window::ON && mwseed_ions == "on") {
            n++;
        }
        double* mcr = new double[n];
        int m;
        bool b;
        m = 0;
        tmp_p_film = p_last_film;
        while (tmp_p_film!=0)
        {
            b = 1;
            for (int i=0;i<m;i++)
            {
                if (tmp_p_film->mcr==mcr[i])
                    b = 0;
            }
            if (b==1)
            {
                mcr[m] = tmp_p_film->mcr;
                m++;
            }
            tmp_p_film = tmp_p_film->prev;
        }
        if (nelflow!=0) {
            b = 1;
            for (int i=0;i<m;i++) {
                if (mcrlflow==mcr[i])
                    b = 0;
            }
            if (b==1) {
                mcr[m] = mcrlflow;
                m++;
            }
        }
        if (nerflow!=0) {
            b = 1;
            for (int i=0;i<m;i++) {
                if (mcrrflow==mcr[i])
                    b = 0;
            }
            if (b==1) {
                mcr[m] = mcrrflow;
                m++;
            }
        }
        if (mwseed_ions == "on") {
            bool b = 1;
            for (int i = 0; i < m; i++) {
                b = b && (mcr[i] != mw_mcr);
            }
            if (b) {
                mcr[m] = mw_mcr;
                m++;
            }
        }
        n_ion_populations = m;
        icmr = new double[m];
        for (int i=0;i<m;i++)
            icmr[i] = 1/(proton_mass*mcr[i]);
        delete [] mcr;
    }
    else
    {
        n_ion_populations = 0;
        icmr = 0;
    }
    current = find("xnpic",first);
    xnpic = current->value;
    current = find("ynpic",first);
    ynpic = current->value;
    current = find("znpic",first);
    znpic = current->value;
    current = find("deps",first);
    if (current->units=="mc^2")
    {
        current->value = current->value*0.511;
        current->units = "MeV";
    }
    deps = current->value; // в программе deps - в MeV
    current = find("deps_p",first);
    if (current->units=="mc^2")
    {
        current->value = current->value*0.511;
        current->units = "MeV";
    }
    deps_p = current->value; // в программе deps - в MeV
    if (deps_p==0) deps_p = deps;
    current = find("deps_ph",first);
    if (current->units=="mc^2")
    {
        current->value = current->value*0.511;
        current->units = "MeV";
    }
    deps_ph = current->value; // в программе deps - в MeV
    if (deps_ph==0) deps_ph = deps;
    current = find("deps_i",first);
    deps_i = current->value; // в программе deps_i - в MeV/нуклон
    if (deps_i==0) deps_i = deps;
    current = find("neps",first);
    neps = current->value;
    current = find("neps_p",first);
    neps_p = current->value;
    if (neps_p==0) neps_p = neps;
    current = find("neps_ph",first);
    neps_ph = current->value;
    if (neps_ph==0) neps_ph = neps;
    current = find("neps_i",first);
    neps_i = current->value;
    if (neps_i==0) neps_i = neps;
    current = find("enthp",first);
    enthp = current->value;
    current = find("enthp_p",first);
    enthp_p = current->value;
    if (enthp_p==0) enthp_p = enthp;
    current = find("enthp_ph",first);
    enthp_ph = current->value;
    if (enthp_ph==0) enthp_ph = enthp;
    current = find("enthp_i",first);
    enthp_i = current->value;
    if (enthp_i==0) enthp_i = enthp;
    
    std::string particles_for_output;
    current = find("particles_for_output",first);
    particles_for_output = current->units;
    if(particles_for_output=="") 
        particles_for_output = "e"; // default; initialize so we could print it to log
    
    if (particles_for_output=="ep" || particles_for_output=="epph")
        write_p = true;
    if (particles_for_output=="eph" || particles_for_output=="epph")
        write_ph = true;
    
    current = find("pmerging",first);
    pmerging = current->units;
    if(pmerging=="") pmerging = "off"; // default
    current = find("crpc",first);
    crpc = current->value;
    current = find("lp_reflection",first);
    lp_reflection = current->units;
    if(lp_reflection=="") lp_reflection = "off";
    current = find("f_reflection",first);
    f_reflection = current->units;
    if(f_reflection=="") f_reflection = "off";
    current = find("n_tracks",first);
    n_tracks = current->value;
    current = find("particles_to_track",first);
    particles_to_track = current->units;
    if (particles_to_track.empty()) particles_to_track = "epgi"; // for backward compatibility
    if (particles_to_track == "none") particles_to_track = "";
    current = find("tr_init",first);
    if (current->units=="volume")
        tr_init = 1;
    else
        tr_init = 0;
    current = find("tr_start",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    if (current->units=="mm")
    {
        current->value = current->value/10/lambda;
        current->units="lambda";
    }
    if (current->units=="cm")
    {
        current->value = current->value/lambda;
        current->units="lambda";
    }
    tr_start = current->value*2*PI;
    current = find("xtr1",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    if (current->units=="fs")
    {
        current->value = current->value*1e-15*2.99792458e10/lambda;
        current->units="lambda";
    }
    xtr1 = current->value*2*PI;
    current = find("ytr1",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    ytr1 = current->value*2*PI;
    current = find("ztr1",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    ztr1 = current->value*2*PI;
    current = find("xtr2",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    if (current->units=="fs")
    {
        current->value = current->value*1e-15*2.99792458e10/lambda;
        current->units="lambda";
    }
    xtr2 = current->value*2*PI;
    current = find("ytr2",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    ytr2 = current->value*2*PI;
    current = find("ztr1",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    ztr1 = current->value*2*PI;
    current = find("ztr2",first);
    if (current->units=="um")
    {
        current->value = current->value*1e-4/lambda;
        current->units="lambda";
    }
    ztr2 = current->value*2*PI;

    ne_profile_x_coords = find_array(first, "ne_profile_x_coords", "1/k", "lambda");
    ne_profile_x_values = find_array(first, "ne_profile_x_values", "");

    if (ne_profile_x_coords.size() != ne_profile_x_values.size()) {
        cout << TERM_RED << "The size of ne_profile_x_coords " << ne_profile_x_coords.size()
                << " is not equal to the size of ne_profile_x_values " << ne_profile_x_values.size() << ". Aborting..."
                << TERM_NO_COLOR << endl;
        return 1;
    }

    if (!is_sorted(ne_profile_x_coords.begin(), ne_profile_x_coords.end())) {
        cout << TERM_RED << "The array ne_profile_x_coords is not sorted. Aborting..." << TERM_NO_COLOR << endl;
        return 1;
    }

    ne_profile_r_coords = find_array(first, "ne_profile_r_coords", "1/k", "lambda");
    ne_profile_r_values = find_array(first, "ne_profile_r_values", "");

    if (ne_profile_r_coords.size() != ne_profile_r_values.size()) {
        cout << TERM_RED << "The size of ne_profile_r_coords " << ne_profile_r_coords.size()
                << " is not equal to the size of ne_profile_r_values " << ne_profile_r_values.size() << ". Aborting..."
                << TERM_NO_COLOR << endl;
        return 1;
    }

    if (!is_sorted(ne_profile_r_coords.begin(), ne_profile_r_coords.end())) {
        cout << TERM_RED << "The array ne_profile_r_coords is not sorted. Aborting..." << TERM_NO_COLOR << endl;
        return 1;
    }

    current = find("data_folder",first);
    data_folder = current->units;
    if (data_folder == "")
        data_folder = "results";

    current = find("catching", first);
    if (current->units == "on")
    {
        catching_enabled = true;
    }

    current = find("qed", first);
    if (current->units == "off") {
        qed_enabled = false;
    }

    current = find("dump_photons", first);
    if (current->units == "on")
    {
        dump_photons = true;
    }

    current = find("output_mode", first);
    if (current->units == "binary")
        output_mode = ios_base::out | ios_base::binary;
    else
        output_mode = ios_base::out;
    //
    if (mpi_rank != 0) {
        output_mode = output_mode | ios_base::app;
    }

    current = find("pusher", first);
    string pusher_str = current->units;

    // the default pusher is Vay
    if (pusher_str.empty()) {
        pusher_str = "vay";
    }

    if (pusher_str == "boris") {
        pusher = pusher_enum::BORIS;
    } else if (pusher_str == "vay") {
        pusher = pusher_enum::VAY;
    } else {
        cout << TERM_RED << "Pusher unknown: " << pusher_str << ". Aborting..." << TERM_NO_COLOR << endl;
        return 1;
    }

    current = find("solver", first);
    string solver_str = current->units;

    //the default is NDFX
    if (solver_str.empty()) {
        solver_str = "ndfx";
    }

    if (solver_str == "ndfx") {
        solver = maxwell_solver_enum::NDFX;
    } else if (solver_str == "fdtd") {
        solver = maxwell_solver_enum::FDTD;
    } else {
        cout << TERM_RED << "Solver unknown: " << current->units << ". Aborting..." << TERM_NO_COLOR << endl;
        return 1;
    }

    balancing_enabled = find_boolean(first, "balancing", false);
    balancing_every = find_int(first, "balancing_every", 20);
    balancing_particle_weight = find_double(first, "balancing_particle_weight", 3.0);
    balancing_threshold = find_double(first, "balancing_threshold", 0.1);

    current = first;
    while (current->next!=0)
    {
        tmp = current->next;
        delete current;
        current = tmp;
    }
    delete current; // удаление последней переменной

    ddi* tmp_ddi = p_last_ddi;
    while (tmp_ddi != 0) {
        p_current_ddi = tmp_ddi; // hence after the loop p_current_ddi points to the first ddi
        tmp_ddi = tmp_ddi->prev;
    }

    if (mpi_rank == 0) {
        ofstream fout_log(data_folder+"/log");
        fout_log<<"dx\n"<<dx/2/PI<<"\n";
        fout_log<<"dy\n"<<dy/2/PI<<"\n";
        fout_log<<"dz\n"<<dz/2/PI<<"\n";
        fout_log<<"dt\n"<<dt/2/PI<<"\n";
        fout_log<<"nx\n"<<int(xlength/dx)<<"\n";
        fout_log<<"ny\n"<<int(ylength/dy)<<"\n";
        fout_log<<"nz\n"<<int(zlength/dz)<<"\n";
        fout_log<<"lambda\n"<<lambda<<"\n";
        fout_log<<"ne\n"<<ne<<"\n";

        ddi* tmp_ddi = p_last_ddi;
        while (tmp_ddi!=0) {
            fout_log<<"t_end\n"<<tmp_ddi->t_end/2/PI<<"\n";
            fout_log<<"output_period\n"<<tmp_ddi->output_period/2/PI<<"\n";
            tmp_ddi = tmp_ddi->prev;
        }
        fout_log << "t_add_mw\n" << t_add_mw / (2 * PI) << '\n';
        fout_log<<"xsigma\n"<<xsigma/2/PI<<"\n";
        fout_log<<"ysigma\n"<<ysigma/2/PI<<"\n";
        fout_log<<"zsigma\n"<<zsigma/2/PI<<"\n";
        fout_log<<"x0fout\n"<<x0fout/2/PI<<"\n";
        fout_log<<"a0y\n"<<a0y<<"\n";
        fout_log<<"a0z\n"<<a0z<<"\n";
        fout_log << "external_bz\n" << external_bz << "\n";
        if (mwindow==moving_window::OFF) {
            fout_log<<"mwindow\n"<<"off"<<"\n";
        } else if (mwindow==moving_window::ON) {
            fout_log<<"mwindow\n"<<"on"<<"\n";
        } else if (mwindow==moving_window::AUTO) {
            fout_log<<"mwindow\n"<<"auto"<<"\n";
        }
        if (mwindow != moving_window::OFF) {
            if (mwseed==0)
                fout_log<<"mwseed\n"<<"off"<<"\n";
            else
                fout_log<<"mwseed\n"<<"on"<<"\n";
            fout_log << "mwseed_ions\n" << mwseed_ions << '\n';
            if (mwseed_ions == "on")
                fout_log << "mw_mcr\n" << mw_mcr << '\n';
        }
        fout_log<<"e_components_for_output\n"<<e_components_for_output<<"\n";
        fout_log<<"b_components_for_output\n"<<b_components_for_output<<"\n";
        if (sscos==0) {
            fout_log<<"f_envelope\n"<<f_envelope<<"\n";
        } else if (sscos == 1) {
            fout_log<<"f_envelope\n"<<"sscos"<<"\n";
        } else if (sscos == 2) {
            fout_log<<"f_envelope\n"<<"pearl"<<"\n";
        } else {
            fout_log<<"f_envelope\n"<<"not known type of envelope!"<<"\n";
        }
        fout_log<<"b_sign\n"<<2.0*b_sign-1<<"\n";
        fout_log<<"x0\n"<<x0/2/PI<<"\n";
        fout_log<<"y0\n"<<y00/2/PI<<"\n";
        fout_log<<"z0\n"<<z00/2/PI<<"\n";
        fout_log<<"xtarget\n"<<xtarget/2/PI<<"\n";
        fout_log<<"ytarget\n"<<ytarget/2/PI<<"\n";
        fout_log<<"ztarget\n"<<ztarget/2/PI<<"\n";
        fout_log<<"phase\n"<<phase<<"\n";
        fout_log<<"lp_reflection\n"<<lp_reflection<<"\n";
        fout_log<<"f_reflection\n"<<f_reflection<<"\n";
        fout_log<<"phi\n"<<phi<<"\n";
        fout_log<<"phi_rotate\n"<<phi_rotate<<"\n";
        fout_log<<"shenergy\n"<<shenergy<<"\n";
        fout_log<<"shphase\n"<<shphase<<"\n";
        fout_log<<"beam\n"<<beam<<"\n";
        fout_log<<"beam_particles\n"<<beam_particles<<"\n";
        fout_log<<"Nb\n"<<Nb<<"\n";
        fout_log<<"epsb\n"<<epsb<<"\n";
        fout_log<<"xb\n"<<xb/2/PI<<"\n";
        fout_log<<"rb\n"<<rb/2/PI<<"\n";
        fout_log<<"x0b\n"<<x0b/2/PI<<"\n";
        fout_log<<"y0b\n"<<y0b/2/PI<<"\n";
        fout_log<<"phib\n"<<phib<<"\n";
        fout_log<<"mwspeed\n"<<mwspeed<<"\n";
        fout_log<<"nelflow\n"<<nelflow<<"\n";
        fout_log<<"vlflow\n"<<vlflow<<"\n";
        fout_log<<"mcrlflow\n"<<mcrlflow<<"\n";
        fout_log<<"Tlflow\n"<<Tlflow<<"\n";
        fout_log<<"nerflow\n"<<nerflow<<"\n";
        fout_log<<"vrflow\n"<<vrflow<<"\n";
        fout_log<<"mcrrflow\n"<<mcrrflow<<"\n";
        fout_log<<"Trflow\n"<<Trflow<<"\n";
        fout_log<<"ions\n"<<ions<<"\n";
        tmp_p_film = p_last_film;
        while (tmp_p_film!=0)
        {
            fout_log<<"film\n"<<"on"<<"\n";
            fout_log<<"x0film\n"<<tmp_p_film->x0/2/PI<<"\n";
            fout_log<<"filmwidth\n"<<tmp_p_film->filmwidth/2/PI<<"\n";
            fout_log<<"gradwidth\n"<<tmp_p_film->gradwidth/2/PI<<"\n";
            fout_log<<"gradwidth_y\n"<<tmp_p_film->gradwidth_y/2/PI<<"\n";
            fout_log<<"y0film\n"<<tmp_p_film->y0/2/PI<<"\n";
            fout_log<<"y1film\n"<<tmp_p_film->y1/2/PI<<"\n";
            fout_log<<"z0film\n"<<tmp_p_film->z0/2/PI<<"\n";
            fout_log<<"z1film\n"<<tmp_p_film->z1/2/PI<<"\n";
            fout_log<<"nfilm\n"<<tmp_p_film->ne<<"\n";
            fout_log<<"mcr\n"<<tmp_p_film->mcr<<"\n";
            fout_log<<"Tfilm\n"<<tmp_p_film->T<<"\n";
            fout_log<<"vxfilm\n"<<tmp_p_film->vx<<"\n";
            fout_log<<"xnpic_film\n"<<tmp_p_film->xnpic_film<<"\n";
            fout_log<<"ynpic_film\n"<<tmp_p_film->ynpic_film<<"\n";
            fout_log<<"znpic_film\n"<<tmp_p_film->znpic_film<<"\n";
            tmp_p_film = tmp_p_film->prev;
        }
        fout_log<<"n_ion_populations\n"<<n_ion_populations<<"\n";
        for (int i=0;i<n_ion_populations;i++)
            fout_log<<"icmr\n"<<icmr[i]<<"\n";
        fout_log<<"xnpic\n"<<xnpic<<"\n";
        fout_log<<"ynpic\n"<<ynpic<<"\n";
        fout_log<<"znpic\n"<<znpic<<"\n";
        fout_log<<"particles_for_output\n"<<particles_for_output<<"\n";
        fout_log<<"deps\n"<<deps<<"\n";
        fout_log<<"neps\n"<<neps<<"\n";
        fout_log<<"enthp\n"<<enthp<<"\n";
        fout_log<<"deps_p\n"<<deps_p<<"\n";
        fout_log<<"neps_p\n"<<neps_p<<"\n";
        fout_log<<"enthp_p\n"<<enthp_p<<"\n";
        fout_log<<"deps_ph\n"<<deps_ph<<"\n";
        fout_log<<"neps_ph\n"<<neps_ph<<"\n";
        fout_log<<"enthp_ph\n"<<enthp_ph<<"\n";
        fout_log<<"deps_i\n"<<deps_i<<"\n";
        fout_log<<"neps_i\n"<<neps_i<<"\n";
        fout_log<<"enthp_i\n"<<enthp_i<<"\n";
        fout_log<<"n_sr\n"<<n_sr<<"\n";
        fout_log<<"n_tracks\n"<<n_tracks<<"\n";
        fout_log<<"particles_to_track\n"<<particles_to_track<<"\n";
        fout_log<<"tr_start\n"<<tr_start/2/PI<<"\n";
        fout_log<<"tr_init\n"<<tr_init<<"\n";
        fout_log<<"xtr1\n"<<xtr1<<"\n";
        fout_log<<"ytr1\n"<<ytr1<<"\n";
        fout_log<<"ztr1\n"<<ztr1<<"\n";
        fout_log<<"xtr2\n"<<xtr2<<"\n";
        fout_log<<"ytr2\n"<<ytr2<<"\n";
        fout_log<<"ztr2\n"<<ztr2<<"\n";
        fout_log<<"pmerging\n"<<pmerging<<"\n";
        fout_log<<"crpc\n"<<crpc<<"\n";
        fout_log<<"$\n";

        fout_log << "catching" << endl << (catching_enabled ? "on" : "off") << endl;
        fout_log << "dump_photons" << endl << (dump_photons ? "on" : "off") << endl;
        fout_log << "qed" << endl;
        fout_log << (qed_enabled ? "on" : "off") << endl;

        fout_log << "balancing" << endl << (balancing_enabled ? "on" : "off") << endl;
        if (balancing_enabled) {
            fout_log << "balancing_every\n" << balancing_every << "\n";
            fout_log << "balancing_threshold\n" << balancing_threshold << "\n";
            fout_log << "balancing_particle_weight\n" << balancing_particle_weight << "\n";
        }

        if (output_mode == (ios_base::out | ios_base::binary))
            fout_log << "output_mode\n" << 1 << '\n';
        else if (output_mode == ios_base::out)
            fout_log << "output_mode\n" << 0 << '\n';
        fout_log<<"#------------------------------\n";
        fout_log<<"polarization = "<<polarization<<"\n";
        fout_log<<"P = "<<(a0y*(a0y>a0z)+a0z*(a0z>=a0y))*(a0y*(a0y>a0z)+a0z*(a0z>=a0y))/8*ysigma*zsigma*8.75e9/1e12<<" TW\n";
        fout_log<<"I = "<<PI*(a0y*(a0y>a0z)+a0z*(a0z>=a0y))*(a0y*(a0y>a0z)+a0z*(a0z>=a0y))*8.75e9/(lambda*lambda)<<" W/cm^2\n";
        fout_log<<"W = "<<(a0y*a0y+a0z*a0z)*xsigma*ysigma*zsigma*PI*sqrt(PI/2)/4*lambda*3.691e4/1e7<<" J\n";
        if (f_envelope=="focused")
        {
            //fout_log<<"a0 in focal plane = "<<(a0y*(a0y>a0z)+a0z*(a0z>=a0y))*sigma/sigma0<<"\n";
            fout_log<<"aperture: F/"<<sigma0/2<<"\n";
        }
        fout_log<<"lambda = "<<lambda*1e4<<" um\n";
        fout_log<<"xsigma = "<<xsigma/2/PI*lambda*1e4<<" um = c*"<<xsigma/2/PI*lambda/2.99792458e10*1e15<<" fs\n";
        fout_log<<"ysigma = "<<ysigma/2/PI*lambda*1e4<<" um\n";
        fout_log<<"zsigma = "<<zsigma/2/PI*lambda*1e4<<" um\n";
        fout_log<<"last_t_end = "<<p_last_ddi->t_end/2/PI*lambda*10<<" mm\n";
        fout_log<<"n_cr = "<<PI/(2.818e-13*lambda*lambda)<<" cm^{-3}\n";
        fout_log<<"ne/n_cr = "<<ne/(PI/(2.818e-13*lambda*lambda))<<"\n";
        fout_log<<"#------------------------------\n";
        fout_log.close();
        cout<<"Initialization done!\n";
    }
    return 0;
}
