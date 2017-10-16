#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <sys/times.h>
#include <numa.h>
#include <unistd.h>
#include "main.h"

using namespace std;

// global variables ------------
double dx,dy,dz,dt;
double xlength,ylength,zlength;
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
bool mwindow,mwseed;
double crpc;
double* ppd;
double phase,phi;
double shenergy, shphase; // second harmonic relative energy and phase
ddi* p_last_ddi; // ddi включает t_end, output_period и f - счётчик для вывода данных в файлы
ddi* p_current_ddi;
double t_add_mw;
film* p_last_film;
double x00,y00,z00;
double sigma0,sigma; // sigma0 - радиус пучка в фокальной плоскости
bool freezing;
bool b_sign; // b_sign = 1 соответствует знаку '+', 0 - знаку '-'
int n_ion_populations;
double* icmr;
int n_tracks;
std::string particles_to_track;
bool tr_init;
double tr_start,xtr1,ytr1,ztr1,xtr2,ytr2,ztr2;
double mwspeed,nelflow,vlflow,mcrlflow,Tlflow,nerflow,vrflow,mcrrflow,Trflow;
double mw_channel_radius;
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
std::string pmerging,pmerging_now;
std::string lp_reflection,f_reflection;
std::string ions;
std::string data_folder;
bool catching_enabled;
bool verbose_logging = true;
bool qed_enabled = true;
ios_base::openmode output_mode;
int init();

std::vector<double> ne_profile_x_coords;
std::vector<double> ne_profile_x_values;

//------------------------------


int l; // номер шага по времени

pthread_mutex_t* sr_mutex_c; // c - create; mutex используется только при создании потоков
pthread_mutex_t* sr_mutex1;
pthread_mutex_t* sr_mutex2;
pthread_mutex_t* sr_mutex_m; // ich - interchange; mutex используется для смешивания на границе слоев
pthread_mutex_t* sr_mutex_m2;
pthread_mutex_t* sr_mutex_rho1;
pthread_mutex_t* sr_mutex_rho2;
pthread_mutex_t* sr_mutex_j1;
pthread_mutex_t* sr_mutex_j2; // mutex для вывода токов

spatial_region* psr; // указатель на массив областей, на которые разбита область вычислений

int n_sr; // число слоёв и потоков
int n_numa_nodes; // number of NUMA nodes
int nx_sr; // число ячеек (по x) в слое
int nx_ich; /* число ячеек для приграничной области слоёв, данные в
               которой замещаются данными соседнего слоя; должно
               быть чётным */
int nm; /* используется при обмене данными между слоями и
           для аккуратного подсчёта спектров, энергии и
           числа частиц */

int nmw = 1;

void* thread_function(void* arg)
{
    int i;
    i = *( (int*) arg );
    pthread_mutex_unlock(&sr_mutex_c[i]);

    numa_run_on_node(int(i*n_numa_nodes/n_sr));
    while(l<p_last_ddi->t_end/dt)
    {
        if (pmerging_now=="on")
            psr[i].pmerging(ppd,pmerging);
        psr[i].birth_from_vacuum(8*PI*PI/(dx*dy*dz)*2.818e-13/lambda); // 2.818e-13 = e^2/mc^2
        psr[i].padvance(freezing);
        psr[i].compute_N(nm*(i!=0),nm*(i!=n_sr-1),dx*dy*dz*1.11485e13*lambda/(8*PI*PI*PI));
        psr[i].compute_energy(nm*(i!=0),nm*(i!=n_sr-1),0.5*dx*dy*dz*3.691e4*lambda/1e7,8.2e-14*dx*dy*dz*1.11485e13*lambda/(8*PI*PI*PI)); // энергия в Джоулях
        psr[i].fout_tracks((i*(nx_sr-nx_ich)+nmw)*dx/2/PI,nm);
        psr[i].fadvance_ndfx();
        if (mwindow==1) psr[i].moving_window(l,nmw,mwspeed);

        // обмен данными на приграничной области слоёв
        if(n_sr>1)
        {
            if(i!=0) pthread_mutex_unlock(&sr_mutex_m[i-1]);
            if(i!=n_sr-1) pthread_mutex_unlock(&sr_mutex_m2[i+1]);
            if(i!=n_sr-1) pthread_mutex_lock(&sr_mutex_m[i]);
            if(i!=0) pthread_mutex_lock(&sr_mutex_m2[i]);
            //
            if (mwindow==1&&(l+1)*dt*mwspeed>nmw*dx)
            {
                /* в этом случае функция moving_window
                 * произвела сдвиг содержимого *psr[i]
                 * на одну ячейку */
                if (i!=n_sr-1) {
                    for (int ii=0;ii<nm+1;ii++) {
                        for (int j=0;j<int(ylength/dy);j++) {
                            for (int k=0;k<int(zlength/dz);k++) {
                                copy(psr[i+1],nm-1+ii,j,k,psr[i],nx_sr-nm-1+ii,j,k);
                                psr[i].copy(psr[i+1].cp[nm-1+ii][j][k].pl,psr[i].cp[nx_sr-nm-1+ii][j][k].pl);
                                psr[i].cp[nx_sr-nm-1+ii][j][k].pl.xplus(nx_sr-nx_ich);
                            }
                        }
                    }
                }
                if (i!=0) {
                    for(int ii=0;ii<nm-1;ii++) {
                        for(int j=0;j<int(ylength/dy);j++) {
                            for(int k=0;k<int(zlength/dz);k++) {
                                copy(psr[i-1],nx_sr-nx_ich+ii,j,k,psr[i],ii,j,k);
                                psr[i].copy(psr[i-1].cp[nx_sr-nx_ich+ii][j][k].pl,psr[i].cp[ii][j][k].pl);
                                psr[i].cp[ii][j][k].pl.xplus(-nx_sr+nx_ich);
                            }
                        }
                    }
                }
            }
            else
            {
                /* в этом случае функция moving_window
                 * не производила сдвига содержимого
                 * *psr[i] на одну ячейку */
                if (i!=n_sr-1) {
                    for (int ii=0;ii<nm;ii++) {
                        for (int j=0;j<int(ylength/dy);j++) {
                            for (int k=0;k<int(zlength/dz);k++) {
                                copy(psr[i+1],nm+ii,j,k,psr[i],nx_sr-nm+ii,j,k);
                                psr[i].copy(psr[i+1].cp[nm+ii][j][k].pl,psr[i].cp[nx_sr-nm+ii][j][k].pl);
                                psr[i].cp[nx_sr-nm+ii][j][k].pl.xplus(nx_sr-nx_ich);
                            }
                        }
                    }
                }
                if (i!=0) {
                    for(int ii=0;ii<nm;ii++) {
                        for(int j=0;j<int(ylength/dy);j++) {
                            for(int k=0;k<int(zlength/dz);k++) {
                                copy(psr[i-1],nx_sr-nx_ich+ii,j,k,psr[i],ii,j,k);
                                psr[i].copy(psr[i-1].cp[nx_sr-nx_ich+ii][j][k].pl,psr[i].cp[ii][j][k].pl);
                                psr[i].cp[ii][j][k].pl.xplus(-nx_sr+nx_ich);
                            }
                        }
                    }
                }
            }
        }
        //
        if(l*dt >= [](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi))
        {
            pthread_mutex_unlock(&sr_mutex_j1[i]);
            pthread_mutex_lock(&sr_mutex_j2[i]);
            // подсчёт плотности частиц
            psr[i].compute_rho();
            pthread_mutex_unlock(&sr_mutex_rho1[i]);
            pthread_mutex_lock(&sr_mutex_rho2[i]);
        }
        //
        pthread_mutex_unlock(&sr_mutex1[i]);
        pthread_mutex_lock(&sr_mutex2[i]);
    }
}

int write_N(ofstream& fout_N)
{
    double N_e,N_p,N_ph;
    N_e = N_p = N_ph = 0;
    
    for (int i=0;i<n_sr;i++) 
        N_e += psr[i].N_e;
    for (int i=0;i<n_sr;i++) 
        N_p += psr[i].N_p;
    for (int i=0;i<n_sr;i++) 
        N_ph += psr[i].N_ph;
        
    fout_N<<N_e<<'\t'<<N_p<<'\t'<<N_ph<<endl;
    
    return N_e + N_p; // for logging in "freezing"; delete when freezing is removed
}

void write_energy(ofstream& fout_energy)
{
    double energy_f,energy_e,energy_p,energy_ph;
    double* ienergy = new double[n_ion_populations];
    energy_f = energy_e = energy_p = energy_ph = 0;
    for (int n=0;n<n_ion_populations;n++)
        ienergy[n] = 0;

    for (int i=0;i<n_sr;i++) 
        energy_f += psr[i].energy_f;
    for (int i=0;i<n_sr;i++) 
        energy_e += psr[i].energy_e;
    for (int i=0;i<n_sr;i++) 
        energy_p += psr[i].energy_p;
    for (int i=0;i<n_sr;i++) 
        energy_ph += psr[i].energy_ph;
    for (int i=0;i<n_sr;i++)
    {
        for (int n=0;n<n_ion_populations;n++)
            ienergy[n] += psr[i].ienergy[n];
    }
    
    fout_energy<<energy_f<<'\t'<<energy_e<<'\t'<<energy_p<<'\t'<<energy_ph;
    for (int n=0;n<n_ion_populations;n++)
        fout_energy<<'\t'<<ienergy[n];
    fout_energy<<endl;
    delete[] ienergy;
}

void write_deleted_particles(ofstream& fout_energy_deleted)
{
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
    
    file_name = data_folder + "/deleted" + file_num_pchar;
    ofstream fout_deleted_e(file_name.c_str());
    
    file_name = data_folder + "/deleted_p" + file_num_pchar;
    ofstream fout_deleted_p(file_name.c_str());
    
    file_name = data_folder + "/deleted_ph" + file_num_pchar;
    ofstream fout_deleted_ph(file_name.c_str());
    
    ofstream* fout_deleted_i = new ofstream[n_ion_populations];
    for (int m=0; m<n_ion_populations; ++m)
    {
        char s_cmr[20];
        sprintf(s_cmr,"%g",icmr[m]);
        file_name = data_folder+"/phasespace_";
        file_name += s_cmr;
        file_name += "_";
        file_name += file_num_pchar;
        fout_deleted_i[m].open(file_name.c_str());
    }
    
    for(int n=0; n<n_sr; n++)
    {
        vector<spatial_region::deleted_particle>& del_particles = psr[n].deleted_particles;
        
        double norm = 8.2e-14*dx*dy*dz*1.11485e13*lambda/(8*PI*PI*PI); // energy in Joules
        for (auto it = del_particles.begin(); it != del_particles.end(); ++it)
        {
            if ((*it).cmr == -1)
            {
                i_particle++;
                if(i_particle > enthp)
                {
                    i_particle = 0;
                    fout_deleted_e << (*it).q << endl;
                    fout_deleted_e << dx*((*it).x + n*(nx_sr-nx_ich))/(2*PI) << endl;
                    fout_deleted_e << dy*((*it).y)/(2*PI) << endl << dz*((*it).z)/(2*PI) << endl;
                    fout_deleted_e << (*it).ux << endl << (*it).uy << endl <<(*it).uz << endl;
                    fout_deleted_e << (*it).g << endl << (*it).chi << endl;
                }
            }
            else if ((*it).cmr == 1)
            {
                i_particle_p++;
                if(i_particle_p > enthp_p)
                {
                    i_particle_ph = 0;
                    fout_deleted_p << (*it).q << endl;
                    fout_deleted_p << dx*((*it).x + n*(nx_sr-nx_ich))/(2*PI) << endl;
                    fout_deleted_p << dy*((*it).y)/(2*PI) << endl << dz*((*it).z)/(2*PI) << endl;
                    fout_deleted_p << (*it).ux << endl << (*it).uy << endl <<(*it).uz << endl;
                    fout_deleted_p << (*it).g << endl << (*it).chi << endl;
                }
            }
            else if ((*it).cmr == 0)
            {
                i_particle_ph++;
                if(i_particle_ph > enthp_ph)
                {
                    i_particle_ph = 0;
                    fout_deleted_ph << (*it).q << endl;
                    fout_deleted_ph << dx*((*it).x + n*(nx_sr-nx_ich))/(2*PI) << endl;
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
                            fout_deleted_i[j] << dx*((*it).x + n*(nx_sr-nx_ich))/(2*PI) << endl;
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
    }

    fout_deleted_e.close();
    fout_deleted_p.close();
    fout_deleted_ph.close();
    for (int m=0; m<n_ion_populations; ++m)
    {
        fout_deleted_i[m].close();
    }
    delete [] fout_deleted_i;
}

void write_energy_deleted(ofstream& fout_energy_deleted)
{
    double energy_e_deleted, energy_p_deleted, energy_ph_deleted;
    double* ienergy_deleted = new double[n_ion_populations];
    energy_e_deleted = energy_p_deleted = energy_ph_deleted = 0;
    for (int n=0; n<n_ion_populations; n++)
        ienergy_deleted[n] = 0;

    for (int i=0; i<n_sr; i++) 
        energy_e_deleted += psr[i].energy_e_deleted;
    for (int i=0; i<n_sr; i++) 
        energy_p_deleted += psr[i].energy_p_deleted;
    for (int i=0; i<n_sr; i++) 
        energy_ph_deleted += psr[i].energy_ph_deleted;
    for (int i=0; i<n_sr; i++)
    {
        for (int n=0;n<n_ion_populations;n++)
            ienergy_deleted[n] += psr[i].ienergy_deleted[n];
    }
    
    fout_energy_deleted << energy_e_deleted << '\t' << energy_p_deleted << '\t' << energy_ph_deleted;
    for (int n=0; n<n_ion_populations; n++)
        fout_energy_deleted << '\t' << ienergy_deleted[n];
    fout_energy_deleted << endl;
    delete[] ienergy_deleted;
}

void write_density(bool write_x, bool write_y, bool write_z,
        std::string x_folder, std::string y_folder, std::string z_folder,
        bool write_ions = false, bool scale_j = false)
{
    std::string file_name;
    char file_num_pchar[20];
    ofstream* pof_x = 0;
    ofstream* pof_y = 0;
    ofstream* pof_z = 0;

    sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);

    if (write_x) {
        file_name = data_folder + "/" + x_folder + file_num_pchar;
        pof_x = new ofstream(file_name.c_str(), output_mode);
    }

    if (write_y) {
        file_name = data_folder + "/" + y_folder + file_num_pchar;
        pof_y = new ofstream(file_name.c_str(), output_mode);
    }
    if (write_z) {
        file_name = data_folder + "/" + z_folder + file_num_pchar;
        pof_z = new ofstream(file_name.c_str(), output_mode);
    }
    
    int onx;
    int onx0;
    for(int i=0;i<n_sr;i++)
    {
        if(i==n_sr-1)
            onx = nx_sr;
        else
            onx = nx_sr-nx_ich/2;
        if(i==0)
            onx0 = 0;
        else
            onx0 = nx_ich/2;
        if (scale_j) {
            psr[i].scale_j(dx / dt);
        }
        psr[i].fout_rho(pof_x,pof_y,pof_z,onx0,onx, output_mode);
    }
    
    int ii = int(((xlength-x0fout)/dx - nx_ich)/(nx_sr - nx_ich));
    if(ii<0) ii = 0;
    psr[ii].fout_rho_yzplane(pof_x,pof_y,pof_z,int((xlength-x0fout)/dx)-ii*(nx_sr-nx_ich), output_mode);

    if (write_ions && (ions=="on"))
    {
        char s_cmr[20];
        for (int n=0;n<n_ion_populations;n++)
        {
            sprintf(s_cmr,"%g",icmr[n]);
            file_name = data_folder+"/irho_";
            file_name += s_cmr;
            file_name += "_";
            file_name += file_num_pchar;
            ofstream fout_irho(file_name.c_str());
            for(int i=0;i<n_sr;i++)
            {
                if(i==n_sr-1)
                    onx = nx_sr;
                else
                    onx = nx_sr-nx_ich/2;
                if(i==0)
                    onx0 = 0;
                else
                    onx0 = nx_ich/2;
                psr[i].fout_irho(n,&fout_irho,onx0,onx, output_mode);
            }
            ii = int(((xlength-x0fout)/dx - nx_ich)/(nx_sr - nx_ich));
            if(ii<0) ii = 0;
            psr[ii].fout_irho_yzplane(n,&fout_irho,int((xlength-x0fout)/dx)-ii*(nx_sr-nx_ich), output_mode);
            fout_irho.close();
        }
    }
    
    if (pof_x) delete pof_x;
    if (pof_y) delete pof_y;
    if (pof_z) delete pof_z;
}

void write_spectrum_phasespace(bool write_p, bool write_ph)
{
    std::string file_name;
    char file_num_pchar[20];
    
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
    file_name = data_folder+"/phasespace";
    sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
    file_name = file_name + file_num_pchar;
    ofstream fout_phasespace(file_name.c_str());
    file_name = data_folder+"/phasespace_p";
    file_name = file_name + file_num_pchar;
    ofstream fout_phasespace_p(file_name.c_str());
    file_name = data_folder+"/phasespace_ph";
    file_name = file_name + file_num_pchar;
    ofstream fout_phasespace_ph(file_name.c_str());
    double* spectrum = new double[neps];
    for(int i=0;i<neps;i++) spectrum[i] = 0;
    double* spectrum_p = new double[neps_p];
    for(int i=0;i<neps_p;i++) spectrum_p[i] = 0;
    double* spectrum_ph = new double[neps_ph];
    for(int i=0;i<neps_ph;i++) spectrum_ph[i] = 0;
    int i_eps;
    ofstream* fout_spectrum_i = new ofstream[n_ion_populations];
    ofstream* fout_phasespace_i = new ofstream[n_ion_populations];
    double** spectrum_i = new double*[n_ion_populations];
    for (int m=0;m<n_ion_populations;m++)
    {
        char s_cmr[20];
        sprintf(s_cmr,"%g",icmr[m]);
        file_name = data_folder+"/spectrum_";
        file_name += s_cmr;
        file_name += "_";
        file_name += file_num_pchar;
        fout_spectrum_i[m].open(file_name.c_str());
        file_name = data_folder+"/phasespace_";
        file_name += s_cmr;
        file_name += "_";
        file_name += file_num_pchar;
        fout_phasespace_i[m].open(file_name.c_str());
        spectrum_i[m] = new double[neps_i];
        for(int i=0;i<neps_i;i++)
            spectrum_i[m][i] = 0;
    }
    spatial_region::plist::particle* current;
    for(int n=0;n<n_sr;n++)
    {
        for(int i=nm*(n!=0);i<nx_sr-nm*(n!=n_sr-1);i++)
        {
            for(int j=0;j<int(ylength/dy);j++)
            {
                for(int k=0;k<int(zlength/dz);k++)
                {
                    current = psr[n].cp[i][j][k].pl.head;
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
                                fout_phasespace<<dx*(current->x+n*(nx_sr-nx_ich))/(2*PI)<<"\n";
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
                                fout_phasespace_p<<dx*(current->x+n*(nx_sr-nx_ich))/(2*PI)<<"\n";
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
                                fout_phasespace_ph<<dx*(current->x+n*(nx_sr-nx_ich))/(2*PI)<<"\n";
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
                                        fout_phasespace_i[m]<<dx*(current->x+n*(nx_sr-nx_ich))/(2*PI)<<"\n";
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
            spectrum_i[m][i] = spectrum_i[m][i]*spectrum_norm_ph;
            fout_spectrum_i[m]<<spectrum_i[m][i]<<"\n";
        }
    }
    delete[] spectrum;
    delete[] spectrum_p;
    delete[] spectrum_ph;
    fout_spectrum.close();
    fout_spectrum_p.close();
    fout_spectrum_ph.close();
    fout_phasespace.close();
    fout_phasespace_p.close();
    fout_phasespace_ph.close();
    //
    for (int m=0;m<n_ion_populations;m++)
    {
        fout_spectrum_i[m].close();
        fout_phasespace_i[m].close();
    }
    delete[] fout_spectrum_i;
    delete[] fout_phasespace_i;
    for (int m=0;m<n_ion_populations;m++)
        delete[] spectrum_i[m];
    delete[] spectrum_i;
}

void write_fields()
{
    std::string file_name;
    char file_num_pchar[20];
    ofstream* pof;
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
        ofstream fout_ex(file_name.c_str());
        pof = &fout_ex;
        for(int i=0;i<n_sr;i++)
        {
            if(i==n_sr-1)
                onx = nx_sr;
            else
                onx = nx_sr-nx_ich/2;
            if(i==0)
                onx0 = 0;
            else
                onx0 = nx_ich/2;
            psr[i].fout_ex(pof,onx0,onx, output_mode);
        }
        ii = int(((xlength-x0fout)/dx - nx_ich)/(nx_sr - nx_ich));
        if(ii<0) ii = 0;
        psr[ii].fout_ex_yzplane(pof,int((xlength-x0fout)/dx)-ii*(nx_sr-nx_ich), output_mode);
        fout_ex.close();
    }
    //
    if (e_components_for_output=="y"||e_components_for_output=="xy"||e_components_for_output=="yz"||e_components_for_output=="xyz")
    {
        file_name = data_folder+"/ey";
        sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
        file_name = file_name + file_num_pchar;
        ofstream fout_ey(file_name.c_str());
        pof = &fout_ey;
        for(int i=0;i<n_sr;i++)
        {
            if(i==n_sr-1)
                onx = nx_sr;
            else
                onx = nx_sr-nx_ich/2;
            if(i==0)
                onx0 = 0;
            else
                onx0 = nx_ich/2;
            psr[i].fout_ey(pof,onx0,onx, output_mode);
        }
        ii = int(((xlength-x0fout)/dx - nx_ich)/(nx_sr - nx_ich));
        if(ii<0) ii = 0;
        psr[ii].fout_ey_yzplane(pof,int((xlength-x0fout)/dx)-ii*(nx_sr-nx_ich), output_mode);
        fout_ey.close();
    }
    //
    if (e_components_for_output=="z"||e_components_for_output=="xz"||e_components_for_output=="yz"||e_components_for_output=="xyz")
    {
        file_name = data_folder+"/ez";
        sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
        file_name = file_name + file_num_pchar;
        ofstream fout_ez(file_name.c_str());
        pof = &fout_ez;
        for(int i=0;i<n_sr;i++)
        {
            if(i==n_sr-1)
                onx = nx_sr;
            else
                onx = nx_sr-nx_ich/2;
            if(i==0)
                onx0 = 0;
            else
                onx0 = nx_ich/2;
            psr[i].fout_ez(pof,onx0,onx, output_mode);
        }
        ii = int(((xlength-x0fout)/dx - nx_ich)/(nx_sr - nx_ich));
        if(ii<0) ii = 0;
        psr[ii].fout_ez_yzplane(pof,int((xlength-x0fout)/dx)-ii*(nx_sr-nx_ich), output_mode);
        fout_ez.close();
    }
    //
    if (b_components_for_output=="x"||b_components_for_output=="xy"||b_components_for_output=="xz"||b_components_for_output=="xyz")
    {
        file_name = data_folder+"/bx";
        sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
        file_name = file_name + file_num_pchar;
        ofstream fout_bx(file_name.c_str(), output_mode);
        pof = &fout_bx;
        for(int i=0;i<n_sr;i++)
        {
            if(i==n_sr-1)
                onx = nx_sr;
            else
                onx = nx_sr-nx_ich/2;
            if(i==0)
                onx0 = 0;
            else
                onx0 = nx_ich/2;
            psr[i].fout_bx(pof,onx0,onx, output_mode);
        }
        ii = int(((xlength-x0fout)/dx - nx_ich)/(nx_sr - nx_ich));
        if(ii<0) ii = 0;
        psr[ii].fout_bx_yzplane(pof,int((xlength-x0fout)/dx)-ii*(nx_sr-nx_ich), output_mode);
        fout_bx.close();
    }
    //
    if (b_components_for_output=="y"||b_components_for_output=="xy"||b_components_for_output=="yz"||b_components_for_output=="xyz")
    {
        file_name = data_folder+"/by";
        sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
        file_name = file_name + file_num_pchar;
        ofstream fout_by(file_name.c_str());
        pof = &fout_by;
        for(int i=0;i<n_sr;i++)
        {
            if(i==n_sr-1)
                onx = nx_sr;
            else
                onx = nx_sr-nx_ich/2;
            if(i==0)
                onx0 = 0;
            else
                onx0 = nx_ich/2;
            psr[i].fout_by(pof,onx0,onx, output_mode);
        }
        ii = int(((xlength-x0fout)/dx - nx_ich)/(nx_sr - nx_ich));
        if(ii<0) ii = 0;
        psr[ii].fout_by_yzplane(pof,int((xlength-x0fout)/dx)-ii*(nx_sr-nx_ich), output_mode);
        fout_by.close();
    }
    //
    if (b_components_for_output=="z"||b_components_for_output=="xz"||b_components_for_output=="yz"||b_components_for_output=="xyz")
    {
        file_name = data_folder+"/bz";
        sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
        file_name = file_name + file_num_pchar;
        ofstream fout_bz(file_name.c_str());
        pof = &fout_bz;
        for(int i=0;i<n_sr;i++)
        {
            if(i==n_sr-1)
                onx = nx_sr;
            else
                onx = nx_sr-nx_ich/2;
            if(i==0)
                onx0 = 0;
            else
                onx0 = nx_ich/2;
            psr[i].fout_bz(pof,onx0,onx, output_mode);
        }
        ii = int(((xlength-x0fout)/dx - nx_ich)/(nx_sr - nx_ich));
        if(ii<0) ii = 0;
        psr[ii].fout_bz_yzplane(pof,int((xlength-x0fout)/dx)-ii*(nx_sr-nx_ich), output_mode);
        fout_bz.close();
    }
    //
    file_name = data_folder+"/w";
    sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
    file_name = file_name + file_num_pchar;
    ofstream fout_w(file_name.c_str());
    pof = &fout_w;
    for(int i=0;i<n_sr;i++)
    {
        if(i==n_sr-1)
            onx = nx_sr;
        else
            onx = nx_sr-nx_ich/2;
        if(i==0)
            onx0 = 0;
        else
            onx0 = nx_ich/2;
        if (i==n_sr-1)
            is_last_sr = 1;
        else
            is_last_sr = 0;
        psr[i].fout_w(pof,onx0,onx,is_last_sr);
    }
    ii = int(((xlength-x0fout)/dx - nx_ich)/(nx_sr - nx_ich));
    if(ii<0) ii = 0;
    psr[ii].fout_w_yzplane(pof,int((xlength-x0fout)/dx)-ii*(nx_sr-nx_ich));
    //
    file_name = data_folder+"/inv";
    sprintf(file_num_pchar,"%g",int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy);
    file_name = file_name + file_num_pchar;
    ofstream fout_inv(file_name.c_str());
    pof = &fout_inv;
    for(int i=0;i<n_sr;i++)
    {
        if(i==n_sr-1)
            onx = nx_sr;
        else
            onx = nx_sr-nx_ich/2;
        if(i==0)
            onx0 = 0;
        else
            onx0 = nx_ich/2;
        if (i==n_sr-1)
            is_last_sr = 1;
        else
            is_last_sr = 0;
        psr[i].fout_inv(pof,onx0,onx,is_last_sr);
    }
    ii = int(((xlength-x0fout)/dx - nx_ich)/(nx_sr - nx_ich));
    if(ii<0) ii = 0;
    psr[ii].fout_inv_yzplane(pof,int((xlength-x0fout)/dx)-ii*(nx_sr-nx_ich));
    fout_w.close();
    fout_inv.close();
}

void init_fields()
{
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
            for (int i=0;i<n_sr;i++) psr[i].f_init_focused(a0y,a0z,xsigma,sigma0,xlength/2+x0-dx*i*(nx_sr - nx_ich),x0,b_sign,phase,y00,z00,0,0,sscos,1,xtarget,ytarget,ztarget);
        else  // adding second harmonic
        {
            double alpha, beta;
            alpha = sqrt(1 - shenergy);
            beta = sqrt(shenergy);
            for (int i=0;i<n_sr;i++)
                psr[i].f_init_focused(a0y * alpha, a0z * alpha, xsigma,sigma0,xlength/2+x0-dx*i*(nx_sr - nx_ich),x0,b_sign,phase,y00,z00,0,0,sscos,1,xtarget,ytarget,ztarget);
            for (int i=0;i<n_sr;i++)
                psr[i].f_init_focused(a0y * beta, a0z * beta, xsigma,sigma0,xlength/2+x0-dx*i*(nx_sr - nx_ich),x0,b_sign, shphase,y00,z00, 1, 0,sscos, 2,xtarget,ytarget,ztarget);
        }
        if (phi!=0) {
            for (int i=0;i<n_sr;i++) psr[i].f_init_focused(a0y,a0z,xsigma,sigma0,xlength/2-x0-dx*i*(nx_sr - nx_ich),-x0,b_sign,phase,-y00,-z00,1,phi,sscos,1,xtarget,ytarget,ztarget);
        }
        else if (lp_reflection1=="xy") {
            for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection1=="y"))*a0y,(1-2*(f_reflection1=="z"))*a0z,xsigma,sigma0,xlength/2+x0-dx*i*(nx_sr - nx_ich),x0,b_sign,phase,y00,-z00,1,0,sscos,1,xtarget,ytarget,ztarget);
            if (lp_reflection2=="xz") {
                for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2+x0-dx*i*(nx_sr - nx_ich),x0,b_sign,phase,-y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2+x0-dx*i*(nx_sr - nx_ich),x0,b_sign,phase,-y00,-z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                if (lp_reflection3=="yz") {
                    for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*i*(nx_sr - nx_ich),-x0,b_sign,phase,y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                    for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*i*(nx_sr - nx_ich),-x0,b_sign,phase,y00,-z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                    for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*i*(nx_sr - nx_ich),-x0,b_sign,phase,-y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                    for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*i*(nx_sr - nx_ich),-x0,b_sign,phase,-y00,-z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                }
            }
            else if (lp_reflection2=="yz") {
                for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*i*(nx_sr - nx_ich),-x0,b_sign,phase,y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*i*(nx_sr - nx_ich),-x0,b_sign,phase,y00,-z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                if (lp_reflection3=="xz") {
                    for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*i*(nx_sr - nx_ich),-x0,b_sign,phase,-y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                    for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2+x0-dx*i*(nx_sr - nx_ich),x0,b_sign,phase,-y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                    for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*i*(nx_sr - nx_ich),-x0,b_sign,phase,-y00,-z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                    for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2+x0-dx*i*(nx_sr - nx_ich),x0,b_sign,phase,-y00,-z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                }
            }
        }
        else if (lp_reflection1=="xz") {
            for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection1=="y"))*a0y,(1-2*(f_reflection1=="z"))*a0z,xsigma,sigma0,xlength/2+x0-dx*i*(nx_sr - nx_ich),x0,b_sign,phase,-y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
            if (lp_reflection2=="yz") {
                for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*i*(nx_sr - nx_ich),-x0,b_sign,phase,-y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
                for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*i*(nx_sr - nx_ich),-x0,b_sign,phase,y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
            }
        }
        else if (lp_reflection1=="yz") {
            for (int i=0;i<n_sr;i++) psr[i].f_init_focused((1-2*(f_reflection1=="y"))*a0y,(1-2*(f_reflection1=="z"))*a0z,xsigma,sigma0,xlength/2-x0-dx*i*(nx_sr - nx_ich),-x0,b_sign,phase,y00,z00,1,0,sscos,1,xtarget,ytarget,ztarget);
        }
    } else if (f_envelope == "uniformB") {
        for (int i = 0; i < n_sr; ++i)
            psr[i].f_init_uniformB(a0y, a0z);
    }
    else // f_envelope == "cos"
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
        for (int i=0;i<n_sr;i++) psr[i].f_init_cos(a0y,a0z,xsigma,ysigma,zsigma,xlength/2+x0-dx*i*(nx_sr - nx_ich),sscos,b_sign,x0,phase,y00,z00,true,0,xtarget,ytarget,ztarget);
        if (phi!=0) {
            for (int i=0;i<n_sr;i++) psr[i].f_init_cos(a0y,a0z,xsigma,ysigma,zsigma,xlength/2-x0-dx*i*(nx_sr - nx_ich),sscos,b_sign,-x0,phase,-y00,-z00,1,phi,xtarget,ytarget,ztarget);
        }
        else if (lp_reflection1=="xy") {
            for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection1=="y"))*a0y,(1-2*(f_reflection1=="z"))*a0z,xsigma,ysigma,zsigma,xlength/2+x0-dx*i*(nx_sr - nx_ich),sscos,b_sign,x0,phase,y00,-z00,1,0,xtarget,ytarget,ztarget);
            if (lp_reflection2=="xz") {
                for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,xlength/2+x0-dx*i*(nx_sr - nx_ich),sscos,b_sign,x0,phase,-y00,z00,1,0,xtarget,ytarget,ztarget);
                for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,xlength/2+x0-dx*i*(nx_sr - nx_ich),sscos,b_sign,x0,phase,-y00,-z00,1,0,xtarget,ytarget,ztarget);
                if (lp_reflection3=="yz") {
                    for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*i*(nx_sr - nx_ich),sscos,b_sign,-x0,phase,y00,z00,1,0,xtarget,ytarget,ztarget);
                    for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*i*(nx_sr - nx_ich),sscos,b_sign,-x0,phase,y00,-z00,1,0,xtarget,ytarget,ztarget);
                    for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*i*(nx_sr - nx_ich),sscos,b_sign,-x0,phase,-y00,z00,1,0,xtarget,ytarget,ztarget);
                    for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*i*(nx_sr - nx_ich),sscos,b_sign,-x0,phase,-y00,-z00,1,0,xtarget,ytarget,ztarget);
                }
            }
            else if (lp_reflection2=="yz") {
                for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*i*(nx_sr - nx_ich),sscos,b_sign,-x0,phase,y00,z00,1,0,xtarget,ytarget,ztarget);
                for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*i*(nx_sr - nx_ich),sscos,b_sign,-x0,phase,y00,-z00,1,0,xtarget,ytarget,ztarget);
                if (lp_reflection3=="xz") {
                    for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*i*(nx_sr - nx_ich),sscos,b_sign,-x0,phase,-y00,z00,1,0,xtarget,ytarget,ztarget);
                    for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,x0+xlength/2-dx*i*(nx_sr - nx_ich),sscos,x0,b_sign,phase,-y00,z00,1,0,xtarget,ytarget,ztarget);
                    for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*i*(nx_sr - nx_ich),sscos,b_sign,-x0,phase,-y00,-z00,1,0,xtarget,ytarget,ztarget);
                    for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection3=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,x0+xlength/2-dx*i*(nx_sr - nx_ich),sscos,b_sign,x0,phase,-y00,-z00,1,0,xtarget,ytarget,ztarget);
                }
            }
        }
        else if (lp_reflection1=="xz") {
            for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection1=="y"))*a0y,(1-2*(f_reflection1=="z"))*a0z,xsigma,ysigma,zsigma,xlength/2+x0-dx*i*(nx_sr - nx_ich),sscos,b_sign,x0,phase,-y00,z00,1,0,xtarget,ytarget,ztarget);
            if (lp_reflection2=="yz") {
                for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*i*(nx_sr - nx_ich),sscos,b_sign,-x0,phase,-y00,z00,1,0,xtarget,ytarget,ztarget);
                for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection2=="y"))*a0y,(1-2*(f_reflection2=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*i*(nx_sr - nx_ich),sscos,b_sign,-x0,phase,y00,z00,1,0,xtarget,ytarget,ztarget);
            }
        }
        else if (lp_reflection1=="yz") {
            for (int i=0;i<n_sr;i++) psr[i].f_init_cos((1-2*(f_reflection1=="y"))*a0y,(1-2*(f_reflection1=="z"))*a0z,xsigma,ysigma,zsigma,-x0+xlength/2-dx*i*(nx_sr - nx_ich),sscos,b_sign,-x0,phase,y00,z00,1,0,xtarget,ytarget,ztarget);
        }
    }
    
    for(int i=0;i<n_sr;i++) psr[i].f_zeroing_on_boundaries();
}

void init_beam()
{
    if (beam_particles=="p")
        for(int i=0;i<n_sr;i++) psr[i].add_beam(1,Nb*1.061e-11/(xb*rb*rb*lambda),((epsb>0)-(epsb<0))*sqrt(epsb*epsb/(0.511*0.511)-1),xb,rb,xlength-x0b-dx*i*(nx_sr - nx_ich),y0b,phib);
    else if (beam_particles=="ph")
        for(int i=0;i<n_sr;i++) psr[i].add_beam(0,Nb*1.061e-11/(xb*rb*rb*lambda),epsb/0.511,xb,rb,xlength-x0b-dx*i*(nx_sr - nx_ich),y0b,phib);
    else
        for(int i=0;i<n_sr;i++) psr[i].add_beam(-1,Nb*1.061e-11/(xb*rb*rb*lambda),((epsb>0)-(epsb<0))*sqrt(epsb*epsb/(0.511*0.511)-1),xb,rb,xlength-x0b-dx*i*(nx_sr - nx_ich),y0b,phib);
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
        for(int i=0;i<n_sr;i++) 
            psr[i].film(tmp_p_film->x0-dx*i*(nx_sr-nx_ich), tmp_p_film->x0+tmp_p_film->filmwidth-dx*i*(nx_sr-nx_ich), 
                tmp_p_film->ne_y0/(1.11485e+13/lambda/lambda), tmp_p_film->ne_y1/(1.11485e+13/lambda/lambda),
                ions=="on", 1/(proton_mass*tmp_p_film->mcr),
                tmp_p_film->gradwidth, tmp_p_film->y0, tmp_p_film->y1, tmp_p_film->z0, tmp_p_film->z1,
                tmp_p_film->T, tmp_p_film->vx, nelflow != 0 || nerflow != 0,
                tmp_p_film->xnpic_film, tmp_p_film->ynpic_film, tmp_p_film->znpic_film, false);
        tmp_p_film = tmp_p_film->prev;
    }
}

void start_tracking()
{
    cout << "Tracking started for particles: " << particles_to_track << endl;
    if (xtr1 < 0 || xtr2 < 0 || ytr1 < 0 || ytr2 < 0 || ztr1 < 0 || ztr2 < 0 ||
        xtr1 >= xlength || xtr2 >= xlength || ytr1 >= ylength || ytr2 >= ylength || ztr1 >= zlength || ztr2 >= zlength)
    {
        cout << "Error - tracks outside of compulational domain" << endl;
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
                n = x2/(nx_sr-nx_ich);
                x = x2%(nx_sr-nx_ich);
                if (n>0 && x<nm) {
                    n = n - 1;
                    x = x + nx_sr - nx_ich;
                }
                spatial_region::plist::particle* h = psr[n].cp[x][y2][z2].pl.head;
                spatial_region::plist::particle* p = h;
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
        }
    } else {
        for (int i=0;i<n_tracks;i++) {
            int x1,y1,z1;
            x1 = ( xtr1 + ( xtr2 - xtr1 ) * rand( ) / RAND_MAX ) / dx;
            y1 = ( ytr1 + ( ytr2 - ytr1 ) * rand( ) / RAND_MAX ) / dy;
            z1 = ( ztr1 + ( ztr2 - ztr1 ) * rand( ) / RAND_MAX ) / dz;
            int n,x;
            n = x1/(nx_sr-nx_ich);
            x = x1%(nx_sr-nx_ich);
            if (i!=0 && x<nm) {
                n = n - 1;
                x = x + nx_sr - nx_ich;
            }
            spatial_region::plist::particle* h = psr[n].cp[x][y1][z1].pl.head;
            spatial_region::plist::particle* p = h;
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
    }
}

void evaluate_merging_condition()
{
    int N_qp_e, N_qp_p, N_qp_g;
    int* N_qp_i;
    N_qp_e = 0;
    N_qp_p = 0;
    N_qp_g = 0;
    N_qp_i = new int[n_ion_populations];
    for (int i=0;i<n_ion_populations;i++)
        N_qp_i[i] = 0;
    for (int i=0;i<n_sr;i++) {
        N_qp_e += psr[i].N_qp_e;
        N_qp_p += psr[i].N_qp_p;
        N_qp_g += psr[i].N_qp_g;
        for (int j=0;j<n_ion_populations;j++)
            N_qp_i[j] += psr[i].N_qp_i[j];
    }
    double crnp = crpc*xlength*ylength*zlength/(dx*dy*dz);
    pmerging_now = "off";
    if (pmerging=="ti") {
        int N_qp;
        N_qp = N_qp_e + N_qp_p + N_qp_g;
        for (int i = 0; i < n_ion_populations; i++)
            N_qp += N_qp_i[i];
        if (N_qp>(3+n_ion_populations)*crnp) {
            pmerging_now = "on";
            // portion of particles that will be deleted
            ppd[0] = (N_qp - (3+n_ion_populations)*crnp)/N_qp;
            cout<<"\t\033[36m"<<"ppd = "<<ppd[0]<<"\033[0m";
        }
    } else if (pmerging=="nl") {
        bool merge = (N_qp_e>crnp)||(N_qp_p>crnp)||(N_qp_g>crnp);
        for (int i = 0; i < n_ion_populations; ++i)
            merge = merge || (N_qp_i[i]>crnp);
        if (merge) {
            pmerging_now = "on";
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
            cout<<"\033[0m";
        }
    }
    delete[] N_qp_i;
}

void add_moving_window_particles()
{
    if (mwindow == 1 && (l + 1) * dt * mwspeed > nmw * dx && (l + 1) * dt < t_add_mw)
    {
        nmw = nmw + 1;
        if (mwseed==1) {
            double n;
            n=1/(k0*k0);
            n *= lin_interpolation((nmw-1)*dx/2.0/PI, ne_profile_x_coords, ne_profile_x_values);

            int_vector3d cell_pos;
            cell_pos.i = nx_sr - 3;
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
                    if (r >= mw_channel_radius)
                    {
                        psr[n_sr-1].fill_cell_by_particles(-1,cell_pos,v_npic,n);
                    }
                    //if (ions=="on")
                    // psr[n_sr-1].fill_cell_by_particles(-1,cell_pos,v_npic,n); // bug??! this adds electrons, not ions; qwe
                }
            }
        }
        
        film* tmp_p_film = p_last_film;
        while (tmp_p_film != 0)
        {
            psr[n_sr-1].film(tmp_p_film->x0-dx*(n_sr-1)*(nx_sr-nx_ich)-dx*nmw, tmp_p_film->x0+tmp_p_film->filmwidth-dx*(n_sr-1)*(nx_sr-nx_ich)-dx*nmw, 
                tmp_p_film->ne_y0/(1.11485e+13/lambda/lambda), tmp_p_film->ne_y1/(1.11485e+13/lambda/lambda),
                ions=="on", 1/(proton_mass*tmp_p_film->mcr),
                tmp_p_film->gradwidth, tmp_p_film->y0, tmp_p_film->y1, tmp_p_film->z0, tmp_p_film->z1,
                tmp_p_film->T, tmp_p_film->vx, nelflow != 0 || nerflow != 0,
                tmp_p_film->xnpic_film, tmp_p_film->ynpic_film, tmp_p_film->znpic_film, true);
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
                    x0 = xrflow - float(ii)/xnpic;
                    if ( x0 >= 0 ) {
                        cell_pos.i = direction > 0 ? 3 : nx_sr - 4;
                    } else {
                        cell_pos.i = direction > 0 ? 2 : nx_sr - 3;
                        x0 += 1;
                    }
                    double tmp = (j * dy - ylength / 2) / (ylength / 2);
                    tmp = cos(0.5 * PI * tmp * tmp * tmp * tmp * tmp);
                    double tr_env = tmp * tmp;
                    tmp =  (k * dz - zlength / 2) / (zlength / 2);
                    tmp = cos(0.5 * PI * tmp * tmp * tmp * tmp * tmp);
                    tr_env *= tmp * tmp;
                    
                    int index = direction > 0 ? 0 : n_sr - 1;
                    psr[index].fill_cell_by_particles(-1,cell_pos,v_npic, n * tr_env, direction * vflow/sqrt(1-vflow*vflow), 0, (direction > 0 ? x0 : 1-x0)-0.5,Tflow); // 0.5 - for a compensation in fill_cell... for xnpic = 1
                    if (ions == "on")
                        psr[index].fill_cell_by_particles(1/(proton_mass*mcrflow),cell_pos,v_npic, n * tr_env, direction * vflow/sqrt(1-vflow*vflow), 0, (direction > 0 ? x0 : 1-x0)-0.5,Tflow / (proton_mass * mcrflow));
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
    for (int i=0; i<n_sr; i++)
    {
        N_qp_e_total += psr[i].N_qp_e;
        N_qp_p_total += psr[i].N_qp_p;
        N_qp_g_total += psr[i].N_qp_g;
        cout << " Thread #" << i << ":\t" << psr[i].N_qp_e << " e, " << psr[i].N_qp_p << " p, " << psr[i].N_qp_g << " g";
        if (n_ion_populations >= 0)
        {
            int N_qp_i = 0;
            for (int j=0; j<n_ion_populations; j++)
            {
                N_qp_i += psr[i].N_qp_i[j];
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
}

int main()
{
    cout<<"\n\033[1m"<<"hi!"<<"\033[0m\n"<<endl;
    
    string hostname;
    ifstream hostname_file("/proc/sys/kernel/hostname");
    if (hostname_file.good())
    {
      getline(hostname_file, hostname);
    }
    hostname_file.close();
    cout << "Quill process id = " << getpid() << ", hostname = " << hostname << endl;
    up_time = times(&tms_struct);
    start_time = times(&tms_struct);
    inaccurate_time = time(NULL);

    nx_ich = 8;
    nm = nx_ich/2;
    file_name_accuracy = 100;

    if (init()==1) return 1;

    ofstream fout_log(data_folder+"/log",ios::app); // ios:app mode - append to the file
    fout_log<<"start time: "<<ctime(&inaccurate_time);

    ppd = new double[3+n_ion_populations];

    nx_sr = ( xlength/dx + nx_ich*(n_sr-1) )/n_sr; // nx = nx_sr*n_sr - nx_ich*(n_sr-1);
    if(nx_ich*(n_sr-1)>=xlength/dx) {cout<<"\033[31m"<<"main: too many slices, aborting..."<<"\033[0m"<<endl; return 1;}

    pthread_t* sr_thread = new pthread_t[n_sr];
    pthread_t* listener_thread = new pthread_t;
    sr_mutex_c = new pthread_mutex_t[n_sr];
    sr_mutex1 = new pthread_mutex_t[n_sr];
    sr_mutex2 = new pthread_mutex_t[n_sr];
    sr_mutex_m = new pthread_mutex_t[n_sr];
    sr_mutex_m2 = new pthread_mutex_t[n_sr];
    sr_mutex_rho1 = new pthread_mutex_t[n_sr];
    sr_mutex_rho2 = new pthread_mutex_t[n_sr];
    sr_mutex_j1 = new pthread_mutex_t[n_sr];
    sr_mutex_j2 = new pthread_mutex_t[n_sr];

    main_thread_time = times(&tms_struct);
    cout<<"Creating arrays..."<<flush;
    psr = new spatial_region[n_sr];
    int node;
    for(int i=0;i<n_sr;i++) 
    {
        node = i*n_numa_nodes/n_sr;
        //
        psr[i].init(i,dx,dy,dz,dt,lambda/2.4263086e-10,xnpic,ynpic,znpic,node,n_ion_populations,icmr,data_folder);
        psr[i].create_arrays(nx_sr,int(ylength/dy),int(zlength/dz),i+times(&tms_struct),node);
        //
    }
    
    init_fields();
    
    if (beam=="on")
    {
        init_beam();
    }

    init_films();
    
    main_thread_time = times(&tms_struct) - main_thread_time;
    seconds = main_thread_time/100.0;
    cout<<"done!"<<endl;
    fout_log<<"fill arrays: "<<seconds<<"s"<<endl;

    for(int i=0;i<n_sr;i++) pthread_mutex_init(&sr_mutex_c[i], 0);
    for(int i=0;i<n_sr;i++) pthread_mutex_lock(&sr_mutex_c[i]);
    for(int i=0;i<n_sr;i++) pthread_mutex_init(&sr_mutex1[i], 0);
    for(int i=0;i<n_sr;i++) pthread_mutex_lock(&sr_mutex1[i]);
    for(int i=0;i<n_sr;i++) pthread_mutex_init(&sr_mutex2[i], 0);
    for(int i=0;i<n_sr;i++) pthread_mutex_lock(&sr_mutex2[i]);
    for(int i=0;i<n_sr;i++) pthread_mutex_init(&sr_mutex_m[i], 0);
    for(int i=0;i<n_sr;i++) pthread_mutex_lock(&sr_mutex_m[i]);
    for(int i=0;i<n_sr;i++) pthread_mutex_init(&sr_mutex_m2[i], 0);
    for(int i=0;i<n_sr;i++) pthread_mutex_lock(&sr_mutex_m2[i]);
    for(int i=0;i<n_sr;i++) pthread_mutex_init(&sr_mutex_rho1[i], 0);
    for(int i=0;i<n_sr;i++) pthread_mutex_lock(&sr_mutex_rho1[i]);
    for(int i=0;i<n_sr;i++) pthread_mutex_init(&sr_mutex_rho2[i], 0);
    for(int i=0;i<n_sr;i++) pthread_mutex_lock(&sr_mutex_rho2[i]);
    for(int i=0;i<n_sr;i++) pthread_mutex_init(&sr_mutex_j1[i], 0);
    for(int i=0;i<n_sr;i++) pthread_mutex_lock(&sr_mutex_j1[i]);
    for(int i=0;i<n_sr;i++) pthread_mutex_init(&sr_mutex_j2[i], 0);
    for(int i=0;i<n_sr;i++) pthread_mutex_lock(&sr_mutex_j2[i]);

    for(int i=0;i<n_sr;i++)
    {
        int* p_i;
        p_i = &i;
        pthread_create(&sr_thread[i],0,thread_function,p_i);
        pthread_mutex_lock(&sr_mutex_c[i]);
    }
    
    pthread_create(listener_thread, 0, listen_for_param_updates, 0);

    l=0;
    ofstream fout_N(data_folder+"/N");
    ofstream fout_energy(data_folder+"/energy");
    ofstream fout_energy_deleted;
    if (catching_enabled)
    {
        fout_energy_deleted.open(data_folder+"/energy_deleted");
    }
    
    while(l<p_last_ddi->t_end/dt)
    {
        /* вывод плотности, спектра и 'phasespace'-данных для фотонов,
           электронов и позитронов в файлы */
        if(l*dt >= [](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi))
        {
            for (int i=0; i<n_sr; i++) pthread_mutex_lock(&sr_mutex_j1[i]);

            if (write_jx || write_jy || write_jz) {
                write_density(write_jx, write_jy, write_jz, "jx", "jy", "jz", false, true);
            }

            for (int i=0; i<n_sr; i++) pthread_mutex_unlock(&sr_mutex_j2[i]);
            for (int i=0; i<n_sr; i++) pthread_mutex_lock(&sr_mutex_rho1[i]);
            
            write_density(true, write_p, write_ph, "rho", "rho_p", "rho_ph", true);
            write_spectrum_phasespace(write_p, write_ph);            
            
            if (catching_enabled)
            {
                write_deleted_particles(fout_energy_deleted);
            }
            
            for(int i=0;i<n_sr;i++) pthread_mutex_unlock(&sr_mutex_rho2[i]);
        }

        for(int i=0;i<n_sr;i++) pthread_mutex_lock(&sr_mutex1[i]);

        if (!particles_to_track.empty() && l == int(tr_start/dt)) 
        {
            start_tracking();
        }

        int N_ep = write_N(fout_N);
        write_energy(fout_energy);
        if (catching_enabled)
        {
            write_energy_deleted(fout_energy_deleted);
        }

        cout<<"\033[33m"<<"ct/lambda = "<<l*dt/2/PI<<"\tstep# "<<l<<"\033[0m";

        evaluate_merging_condition();

        if (freezing==1)
        {
            int N_freezed;
            N_freezed = 0;
            for (int i=0;i<n_sr;i++) N_freezed += psr[i].N_freezed;
            if (N_freezed!=0)
                cout<<"\t\033[36m"<<"N_f/N_c = "<<N_freezed*dx*dy*dz*1.11485e13*lambda/(8*PI*PI*PI)/(N_ep)<<"\033[0m";
        }
        cout<<endl;

        add_moving_window_particles();

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
            
            cout<<"output# "<<"\033[1m"<<int([](ddi* a) {double b=a->f*a->output_period; if(a->prev!=0) b+=(a->prev)->t_end; return b;} (p_current_ddi)/2/PI*file_name_accuracy)/file_name_accuracy<<"\033[0m"<<flush;
            
            write_fields(); // Ex..Ez, Bx..Bz, w (field energy density), inv (E^2-B^2 - relativistic invariant)
            
            cout<<"\t"<<"up: "<<(times(&tms_struct)-start_time)/100.0<<" s"<<endl;
            p_current_ddi->f++;
        }
        
        if (l*dt>=p_current_ddi->t_end)
        {
            p_current_ddi = p_current_ddi->next;
        }

        l++;
        for(int i=0;i<n_sr;i++) pthread_mutex_unlock(&sr_mutex2[i]);
    }

    delete[] sr_thread;
    delete listener_thread;
    delete[] sr_mutex_c;
    delete[] sr_mutex1;
    delete[] sr_mutex2;
    delete[] sr_mutex_m;
    delete[] sr_mutex_m2;
    delete[] sr_mutex_rho1;
    delete[] sr_mutex_rho2;
    delete[] sr_mutex_j1;
    delete[] sr_mutex_j2;
    delete[] icmr;

    delete[] ppd;
    fout_N.close();
    fout_energy.close();
    if (fout_energy_deleted.is_open())
    {
        fout_energy_deleted.close();
    }

    while(p_last_ddi!=0)
    {
        ddi* tmp = p_last_ddi->prev;
        delete p_last_ddi;
        p_last_ddi = tmp;
    }

    inaccurate_time = time(NULL);
    fout_log<<"stop time: "<<ctime(&inaccurate_time);
    up_time = times(&tms_struct) - up_time;
    seconds = up_time/100.0;
    fout_log<<"uptime: "<<seconds<<"s"<<endl;
    fout_log.close();

    cout<<"\n\033[1mbye!\033[0m\n"<<endl;
    return 0;
}

int init()
{
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
        dx = (1+1e-2)*dt/( 1 - 1/(k0*k0)*dt*dt/4.04 );
        if (dx < 0) {
            cout << "\033[31m" << "Cannot calculate stable dx, reduce dt. Aborting..." << "\033[0m" << endl;
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
    current = find("f_envelope",first);
    f_envelope = current->units;
    if (f_envelope=="") {
        f_envelope = "cos";
        sscos = 0;
    } else if (f_envelope == "sscos") {
        f_envelope = "cos";
        sscos = 1;
    } else if (f_envelope == "pearl") {
        f_envelope == "cos";
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
          cout<<"\n\033[31m"<<"main: improper focused pulse, aborting..."<<"\033[0m"<<endl;
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
    mwindow = 1;
    if (current->units=="off") mwindow = 0;
    current = find("mw_channel_radius", first);
    mw_channel_radius = current->value * 2 * PI;
    current = find("mwspeed",first);
    mwspeed = current->value;
    if (mwindow==1 && mwspeed==0)
        mwspeed = 1;
    current = find("mwseed",first);
    mwseed = 1;
    if (current->units=="off") mwseed = 0;
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
    if (ions=="") ions = "off";
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
        if (current->value==0)
            current->value = 1;
        p_last_film->mcr = current->value;
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
        mwindow = 0;
    //
    if (ions=="on")
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
    current = find("n_sr",first);
    n_sr = current->value;
    if (n_sr==0) n_sr = 8;
    current = find("n_numa_nodes",first);
    n_numa_nodes = current->value;
    if (n_numa_nodes==0) n_numa_nodes = 2;
    
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
    current = find("freezing",first);
    if (current->units=="on")
        freezing = 1;
    else
        freezing = 0;
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

    current = find("ne_profile_x_coords",first);
    if (current->units == "um") {
        for (double & v : current->input_array) {
            v *= 1e-4/lambda;
        }
        current->units = "lambda";
    }
    ne_profile_x_coords = current->input_array;
    if (ne_profile_x_coords.empty()) {
        ne_profile_x_coords = {0.0};
    }

    current = find("ne_profile_x_values",first);
    ne_profile_x_values = current->input_array;
    if (ne_profile_x_values.empty()) {
        ne_profile_x_values = {1.0};
    }

    if (ne_profile_x_coords.size() != ne_profile_x_values.size()) {
        cout << "\033[31m" << "The size of ne_profile_x_coords " << ne_profile_x_coords.size()
                << " is not equal to the size of ne_profile_x_values " << ne_profile_x_values.size() << ". Aborting..."
                << "\033[0m" << endl;
        return 1;
    }

    if (!is_sorted(ne_profile_x_coords.begin(), ne_profile_x_coords.end())) {
        cout << "\033[31m" << "The array ne_profile_x_coords is not sorted. Aborting..." "\033[0m" << endl;
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

    current = find("output_mode", first);
    if (current->units == "binary")
        output_mode = ios_base::out | ios_base::binary;
    else
        output_mode = ios_base::out;
    //
    current = first;
    while (current->next!=0)
    {
        tmp = current->next;
        delete current;
        current = tmp;
    }
    delete current; // удаление последней переменной
    //
    int nx_sr1 = ( xlength/dx + nx_ich*(n_sr-1) )/n_sr; // nx = nx_sr*n_sr - nx_ich*(n_sr-1);
    ofstream fout_log(data_folder+"/log");
    fout_log<<"dx\n"<<dx/2/PI<<"\n";
    fout_log<<"dy\n"<<dy/2/PI<<"\n";
    fout_log<<"dz\n"<<dz/2/PI<<"\n";
    fout_log<<"dt\n"<<dt/2/PI<<"\n";
    fout_log<<"nx\n"<<(nx_sr1*n_sr - nx_ich*(n_sr-1))<<"\n";
    fout_log<<"ny\n"<<int(ylength/dy)<<"\n";
    fout_log<<"nz\n"<<int(zlength/dz)<<"\n";
    fout_log<<"lambda\n"<<lambda<<"\n";
    fout_log<<"ne\n"<<ne<<"\n";
    for(int counter=0;counter<(int)ne_profile_x_coords.size();counter++)
        fout_log<<"ne_profile_x_coords "<< counter<< " " << ne_profile_x_coords[counter]<<"\n";
    for(int counter=0;counter<(int)ne_profile_x_values.size();counter++)
        fout_log<<"ne_profile_x_values "<< counter<< " " << ne_profile_x_values[counter]<<"\n";

    ddi* tmp_ddi = p_last_ddi;
    while (tmp_ddi!=0) {
        fout_log<<"t_end\n"<<tmp_ddi->t_end/2/PI<<"\n";
        fout_log<<"output_period\n"<<tmp_ddi->output_period/2/PI<<"\n";
        p_current_ddi = tmp_ddi; // hence after the loop p_current_ddi points to the first ddi
        tmp_ddi = tmp_ddi->prev;
    }
    fout_log << "t_add_mw\n" << t_add_mw / (2 * PI) << '\n';
    fout_log<<"xsigma\n"<<xsigma/2/PI<<"\n";
    fout_log<<"ysigma\n"<<ysigma/2/PI<<"\n";
    fout_log<<"zsigma\n"<<zsigma/2/PI<<"\n";
    fout_log<<"x0fout\n"<<x0fout/2/PI<<"\n";
    fout_log<<"a0y\n"<<a0y<<"\n";
    fout_log<<"a0z\n"<<a0z<<"\n";
    if (mwindow==0)
        fout_log<<"mwindow\n"<<"off"<<"\n";
    else
        fout_log<<"mwindow\n"<<"on"<<"\n";
    if (mwseed==0)
        fout_log<<"mwseed\n"<<"off"<<"\n";
    else
        fout_log<<"mwseed\n"<<"on"<<"\n";
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
    fout_log<<"n_numa_nodes\n"<<n_numa_nodes<<"\n";
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
    fout_log<<"freezing\n";
    if (freezing==1) fout_log<<"on"; else fout_log<<"off";
    fout_log<<'\n';
    
    fout_log << "catching" << endl;
    if (catching_enabled)
        fout_log << "on" << endl;
    else
        fout_log << "off" << endl;
    
    fout_log << "qed" << endl;
    fout_log << (qed_enabled ? "on" : "off") << endl;

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
    return 0;
}
