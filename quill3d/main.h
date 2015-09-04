#include <fstream>

using namespace std;

const double PI = 3.141592653589793;
const double proton_mass = 1836.1526721; /* 1836... - отношение массы
					    протона к массе электрона
					  */

class vector3d
{
public:
    double x,y,z;
    vector3d();
};

class int_vector3d
{
public:
    int i,j,k;
    int_vector3d();
};

class spatial_region
{
    int nx,ny,nz;
    double dx,dy,dz,dt;
    int xnpic,ynpic,znpic;
    double e_s; // Sauter field normalized to mc\omega/e
    double* random;
    int n_random;
    int node_number;
    //
    int n_ion_populations;
    double* icmr; // ion charge to mass ratios (=1/proton_mass for proton)
    double**** irho; // ion densities
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
    class celle
    {
    public:
	double ex,ey,ez;
	celle();
    };
    class cellj
    {
    public:
	double jx,jy,jz;
	cellj();
    };
    class cellb
    {
    public:
	double bx,by,bz;
	cellb();
    };
    class cellbe
    {
    public:
	double bex,bey,bez;
	cellbe();
    };
    class cellp
    {
    public:
	plist pl;
	cellp();
    };
    celle*** ce;
    cellb*** cb;
    cellj*** cj;
    cellbe*** cbe;
    cellp*** cp;
    spatial_region();
    void init(double,double,double,double,double,int,int,int,int,int,double*,std::string);
    void create_arrays(int,int,int,int,int);
    ~spatial_region();
    void fout_ex(ofstream*,int,int);
    void fout_ex_yzplane(ofstream*,int);
    void fout_ey(ofstream*,int,int);
    void fout_ey_yzplane(ofstream*,int);
    void fout_ez(ofstream*,int,int);
    void fout_ez_yzplane(ofstream*,int);
    void fout_bx(ofstream*,int,int);
    void fout_bx_yzplane(ofstream*,int);
    void fout_by(ofstream*,int,int);
    void fout_by_yzplane(ofstream*,int);
    void fout_bz(ofstream*,int,int);
    void fout_bz_yzplane(ofstream*,int);
    void fout_w(ofstream*,int,int,bool);
    void fout_w_yzplane(ofstream*,int);
    void fout_inv(ofstream*,int,int,bool);
    void fout_inv_yzplane(ofstream*,int);
    void fout_rho(ofstream*,ofstream*,ofstream*,int,int);
    void fout_rho_yzplane(ofstream*,ofstream*,ofstream*,int);
    void fout_irho(int,ofstream*,int,int);
    void fout_irho_yzplane(int,ofstream*,int);
    void fout_tracks(double,int);
    void f_init_cos(double,double,double,double,double,double,bool=0,bool=1,double=0,double=0,double=0,double=0,bool=1,double=0);
    void f_init_focused(double,double,double,double,double,double,bool=1,double=0,double=0,double=0,bool=1,double=0,bool=0);
    void f_init_uniformB(double, double);
    void add_beam(double,double,double,double,double,double);
    void film(double,double,double,bool,double,double,double,double,double,double,double);
    void fill_cell_by_particles(double,int_vector3d&,int_vector3d&,double,double=0,double=0,double=0);
    void fadvance_ndfx();
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
    vector3d e_to_particle(double&,double&,double&);
    vector3d b_to_particle(double&,double&,double&);
    void moving_window(int,int,double);
    double get_rand();
    plist::particle* bear_particle(double,vector3d&,vector3d&,double,double,double);
    plist::particle* new_particle();
    void delete_particle(plist::particle*);
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
};

class film
{
public:
    film* prev;
    double x0, filmwidth, gradwidth, y0, y1, z0, z1, ne, mcr, T;
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
    double value;
    std::string units;
    var* next;
    int read();
    var();
};

vector3d regulate(double&, double&, double&);
var* find(std::string, var*);
void copy(spatial_region&,int,int,int,spatial_region&,int,int,int);
