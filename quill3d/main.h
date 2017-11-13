#include <fstream>
#include <vector>
#include "containers.h"

using namespace std;

const double PI = 3.141592653589793;
const double proton_mass = 1836.1526721; /* 1836... - отношение массы
                                            протона к массе электрона
                                          */

class spatial_region
{
    int nx,ny,nz;
    double dx,dy,dz,dt;
    int xnpic,ynpic,znpic;
    double e_s; // Sauter field normalized to mc\omega/e
    double* random;
    int n_random;
    int node_number;
    int sr_id; // Current spatial region id
    //
    int n_ion_populations;
    double* icmr; // ion charge to mass ratios (=1/proton_mass for proton)
    field3d<double> * irho; // ion densities
    //
    int n_ap; // number of allocated particles
    int n_f; // number of free places
    int npwpa; // number of pages with particles
    int npwpo; // number of pages with pointers
    int page_size;
    int particle_size;
    int pointer_size;
    //
    public:
    double N_e,N_p,N_ph;
    int N_qp_e, N_qp_p, N_qp_g; // N_qp - number of quasiparticles
    int* N_qp_i;
    double energy_f,energy_e,energy_p,energy_ph;
    double* ienergy; // energies of ion populations
    double energy_e_deleted, energy_p_deleted, energy_ph_deleted; // for particles caught at boundaries
    double* ienergy_deleted;
    double N_freezed;
    std::string data_folder;
    //
    class plist
    {
        public:
            class particle
            {
                public:
                    double x,y,z,ux,uy,uz,g,q,cmr,chi; // cmr - charge to mass ratio,
                    int trn; // track name
                    particle* next;
                    particle* previous;
                    particle();
                    vector3d get_displacement(double&);
                    void coordinate_advance(vector3d&);
                    void momentum_advance(vector3d&,vector3d&,double&);
            };
            particle* head;
            particle* start;
            plist();
            void xplus(double);
    };
    //
    class pwpa
    { // page with particles
        public:
            plist::particle* head;
            pwpa* previous;
    };
    class pwpo
    { // page with pointers
        public:
            plist::particle** head;
            pwpo* previous;
    };
    plist::particle* p_lap; // pointer to last allocated particle
    plist::particle** pp_lfp; // pointer to pointer to last free place (particle)
    pwpa* p_lapwpa; // pointer to last page with particles
    pwpo* p_lapwpo; // pointer to last page with pointers
    //
    class cellp
    {
        public:
            plist pl;
            cellp();
    };
    field3d<celle> ce;
    field3d<cellb> cb;
    field3d<cellj> cj;
    field3d<cellbe> cbe;
    cellp*** cp;
    class deleted_particle
    {
        public:
            double cmr, q;
            double x, y, z, ux, uy, uz, g, chi;
            deleted_particle(double cmr, double q, double x, double y, double z, double ux, double uy, double uz, double g, double chi):
                cmr(cmr), q(q), x(x), y(y), z(z), ux(ux), uy(uy), uz(uz), g(g), chi(chi) {};
    };
    vector<deleted_particle> deleted_particles;
    spatial_region();
    void init(int,double,double,double,double,double,int,int,int,int,int,double*,std::string);
    void create_arrays(int,int,int,int,int);
    ~spatial_region();
    void fout_ex(ofstream*,int,int, ios_base::openmode);
    void fout_ex_yzplane(ofstream*,int, ios_base::openmode);
    void fout_ey(ofstream*,int,int, ios_base::openmode);
    void fout_ey_yzplane(ofstream*,int, ios_base::openmode);
    void fout_ez(ofstream*,int,int, ios_base::openmode);
    void fout_ez_yzplane(ofstream*,int, ios_base::openmode);
    void fout_bx(ofstream*,int,int, ios_base::openmode);
    void fout_bx_yzplane(ofstream*,int, ios_base::openmode);
    void fout_by(ofstream*,int,int, ios_base::openmode);
    void fout_by_yzplane(ofstream*,int, ios_base::openmode);
    void fout_bz(ofstream*,int,int, ios_base::openmode);
    void fout_bz_yzplane(ofstream*,int, ios_base::openmode);
    void fout_w(ofstream*,int,int,bool);
    void fout_w_yzplane(ofstream*,int);
    void fout_inv(ofstream*,int,int,bool);
    void fout_inv_yzplane(ofstream*,int);
    void fout_rho(ofstream*,ofstream*,ofstream*,int,int, ios_base::openmode);
    void fout_rho_yzplane(ofstream*,ofstream*,ofstream*,int, ios_base::openmode);
    void fout_irho(int,ofstream*,int,int, ios_base::openmode);
    void fout_irho_yzplane(int,ofstream*,int, ios_base::openmode);
    void fout_tracks(double,int);
    void f_init_cos(double,double,double,double,double,double,int,bool,double,double,double,double,bool,double,double,double,double);
    void f_init_focused(double,double,double,double,double,double,bool,double,double,double,bool,double,int,double,double,double,double);
    void f_init_uniformB(double, double);
    void add_beam(double,double,double,double,double,double,double,double);
    void film(double,double,double,double,bool,double,double,double,double,double,double,double,double,bool,int,int,int,bool);
    void fill_cell_by_particles(double,int_vector3d&,int_vector3d&,double,double=0,double=0,double=0,double=0);
    void fadvance();
    void interpolate_be();
    void f_zeroing_on_boundaries();
    void padvance(bool=0);
    void birth_from_vacuum(double);
    void jdeposition(plist::particle&,vector3d&);
    void rhodeposition(plist::particle&);
    void compute_rho();
    void simple_jdep(plist::particle&,vector3d&,int_vector3d&);
    void place(plist::particle&,int&,int&,int&);
    void place(plist::particle&);
    void p_boundary();
    bool is_inside(int,int,int);
    bool is_inside_global(int, int, int);
    bool is_in_exchange_area(int, int, int);
    vector3d e_to_particle(double&,double&,double&);
    vector3d b_to_particle(double&,double&,double&);
    void moving_window(int,int,double);
    double get_rand();
    plist::particle* bear_particle(double,vector3d&,vector3d&,double,double,double);
    plist::particle* new_particle();
    void delete_particle(plist::particle*,bool=false);
    void erase(plist&);
    void copy(plist&,plist&);
    void compute_N(int,int,double);
    void compute_energy(int,int,double,double);
    double chi(vector3d&,vector3d&,vector3d&,double&);
    double w(const double&,double&,const double&); /* осторожно,
                                                      эта функция может изменить r
                                                      (второй аргумент)! */
    double tilde_w(double&,double&,double&);
    double mathcal_W(vector3d&,vector3d&);
    void pmerging(double*,string);
    double _ppd(double,double); // вспомогательная функция
    void scale_j(double);
    
    private:
    void update_energy_deleted(plist::particle*);
};

class film
{
    public:
        film* prev;
        double x0, filmwidth, gradwidth, y0, y1, z0, z1, ne, ne_y0, ne_y1, mcr, T, vx;
        int xnpic_film, ynpic_film, znpic_film;
        film();
};

class ddi
{ // double, double, int
    public:
        ddi* prev;
        ddi* next;
        double t_end,output_period;
        int f;
};

class var
{
    public:
        std::string name;
        std::vector<double> input_array;
        double value;
        std::string units;
        var* next;
        int read();
        var();
};

vector3d regulate(double&, double&, double&);
var* find(std::string, var*);
void copy(spatial_region&,int,int,int,spatial_region&,int,int,int);
void* listen_for_param_updates(void*);
double lin_interpolation(double, std::vector<double>&, std::vector<double>&);
