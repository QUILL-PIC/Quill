#include <fstream>
#include <vector>
#include <memory>
#include <functional>
#include "containers.h"
#include "maxwell.h"

using namespace std;

const double PI = 3.141592653589793;
const double proton_mass = 1836.1526721; /* 1836... - отношение массы
                                            протона к массе электрона
                                          */

const std::string TERM_RED = "\033[31m";
const std::string TERM_NO_COLOR = "\033[0m";

class spatial_region
{
    int nx,ny,nz;
    double dx,dy,dz,dt;
    int xnpic,ynpic,znpic;
    double e_s; // Sauter field normalized to mc\omega/e
    double* random;
    int n_random;
    int sr_id; // Current spatial region id
    //
    int n_ion_populations;
    double* icmr; // ion charge to mass ratios (=1/proton_mass for proton)
    vector<field3d<double> > irho; // ion densities
    //
    int n_ap; // number of allocated particles
    int n_f; // number of free places
    int npwpa; // number of pages with particles
    int npwpo; // number of pages with pointers
    int page_size;
    int particle_size;
    int pointer_size;
    //
    unique_ptr<maxwell_solver> solver;
    function<void(particle&, vector3d&, vector3d&, double&)> advance_momentum;

    public:
    double N_e,N_p,N_ph;
    int N_qp_e, N_qp_p, N_qp_g; // N_qp - number of quasiparticles
    int* N_qp_i;
    double energy_f,energy_e,energy_p,energy_ph;
    double* ienergy; // energies of ion populations
    double energy_e_deleted, energy_p_deleted, energy_ph_deleted; // for particles caught at boundaries
    double* ienergy_deleted;
    std::string data_folder;
    //
    class pwpa
    { // page with particles
        public:
            particle* head;
            pwpa* previous;
    };
    class pwpo
    { // page with pointers
        public:
            particle** head;
            pwpo* previous;
    };
    particle* p_lap; // pointer to last allocated particle
    particle** pp_lfp; // pointer to pointer to last free place (particle)
    pwpa* p_lapwpa; // pointer to last page with particles
    pwpo* p_lapwpo; // pointer to last page with pointers
    //
    field3d<celle> ce;
    field3d<cellb> cb;
    field3d<cellj> cj;
    field3d<cellbe> cbe;
    field3d<cellp> cp;
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
    void init(int,double,double,double,double,double,int,int,int,int,double*,string,maxwell_solver_enum,pusher_enum);
    void create_arrays(int,int,int,int);
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
    void film(double,double,double,double,bool,double,double,double,double,double,double,double,double,bool,int,int,int,bool,double);
    void fill_cell_by_particles(double,int_vector3d&,int_vector3d&,double,double=0,double=0,double=0,double=0);
    void fadvance();
    void interpolate_be();
    void f_init_boundaries();
    void padvance(double=0);
    void birth_from_vacuum(double);
    void jdeposition(particle&,vector3d&);
    void rhodeposition(particle&);
    void compute_rho();
    void simple_jdep(particle&,vector3d&,int_vector3d&);
    void place(particle&,int&,int&,int&);
    void place(particle&);
    void p_boundary();
    bool is_inside(int,int,int);
    bool is_inside_global(int, int, int);
    bool is_in_exchange_area(int);
    vector3d e_to_particle(double&,double&,double&);
    vector3d b_to_particle(double&,double&,double&);
    void moving_window();
    double get_rand();
    particle* bear_particle(double,vector3d&,vector3d&,double,double,double);
    particle* new_particle();
    void delete_particle(particle*,bool=false);
    void erase(plist&);
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

    int get_nx() {return nx;}
    int get_ny() {return ny;}
    int get_nz() {return nz;}
    
    double get_max_w();
    double get_max_w(double left, double right);

    vector<double> calculate_layer_weights(double particle_weight);

    private:
    void update_energy_deleted(particle*);
};

class film
{
    public:
        film* prev;
        double x0, filmwidth, gradwidth, gradwidth_y, y0, y1, z0, z1, ne, ne_y0, ne_y1, mcr, T, vx;
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
double lin_interpolation(double, std::vector<double>&, std::vector<double>&);
