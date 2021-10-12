#include <iostream>
#include <algorithm>
#include <memory>
#include <cstdlib>
#include <unistd.h>
#include "main.h"
#include "maxwell.h"
#include "containers.h"
#include "pusher.h"

extern bool catching_enabled;
extern double lambda;
extern int mpi_rank;

spatial_region::spatial_region()
{
    nx = 0;
    ny = 0;
    nz = 0;
    dx = 1;
    dy = 1;
    dz = 1;
    dt = 0.5;
    xnpic = 1;
    ynpic = 1;
    znpic = 1;
    e_s = 1e10;
    ce = field3d<celle>();
    cj = field3d<cellj>();
    cb = field3d<cellb>();
    cbe = field3d<cellbe>();
    cp = field3d<cellp>();
    n_random = 0;
    energy_e_deleted = 0;
    energy_p_deleted = 0;
    energy_ph_deleted = 0;
    //
    n_ap = 0;
    n_f = 0;
    npwpa = 0;
    npwpo = 0;
    // page_size should not be very small to avoid OS problems with allocation
    const int PS = 64 * 1024 * 1024; // bytes
    if (getpagesize() < PS) {
        page_size = PS;
    } else {
        page_size = getpagesize();
    }
    particle_size = sizeof(particle);
    pointer_size = sizeof(particle*);
    p_lap = 0;
    pp_lfp = 0;
    p_lapwpa = 0;
    p_lapwpo = 0;
    n_ion_populations = 0;
    icmr = 0;
    data_folder = "results";
}

void spatial_region::init(int sr_id0, double dx0, double dy0, double dz0, double dt0, double e_s0, int xnpic0,
        int ynpic0, int znpic0, int n_ion_populations0, double* icmr0, std::string df,
        maxwell_solver_enum solver0, pusher_enum pusher)
{
    sr_id = sr_id0;
    dx = dx0;
    dy = dy0;
    dz = dz0;
    dt = dt0;
    e_s = e_s0;
    xnpic = xnpic0;
    ynpic = ynpic0;
    znpic = znpic0;
    n_ion_populations = n_ion_populations0;
    icmr = icmr0;
    data_folder = df;

    // we should always have at least one page for free particles (to avoid performance issues)
    pwpo* b = new pwpo;
    b->previous = 0;
    b->head = (particle**) malloc(page_size);
    npwpo++;
    p_lapwpo = b;

    switch (solver0) {
    case maxwell_solver_enum::NDFX:
        solver = unique_ptr<maxwell_solver>(new ndfx_solver(ce, cb, cj, dt, dx, dy, dz));
        break;
    case maxwell_solver_enum::FDTD:
        solver = unique_ptr<maxwell_solver>(new fdtd_solver(ce, cb, cj, dt, dx, dy, dz));
        break;
    default:
        solver = nullptr;
    }

    switch (pusher) {
    case pusher_enum::BORIS:
        advance_momentum = push_boris;
        break;
    case pusher_enum::VAY:
        advance_momentum = push_vay;
        break;
    default:
        advance_momentum = nullptr;
    }
}

void spatial_region::create_arrays(int nx0, int ny0, int nz0, int seed)
{
    nx = nx0;
    ny = ny0;
    nz = nz0;

    void* pv;

    ce = field3d<celle>(nx, ny, nz);
    cb = field3d<cellb>(nx, ny, nz);
    cj = field3d<cellj>(nx, ny, nz);
    cbe = field3d<cellbe>(nx, ny, nz);
    cp = field3d<cellp>(nx, ny, nz);

    irho = vector<field3d<double> >(n_ion_populations);
    for (int n=0;n<n_ion_populations;n++)
    {
        irho[n] = field3d<double>(nx, ny, nz);
    }

    pv = malloc(sizeof(double)*55);
    random = (double*) pv;
    srand(seed);
    for(int i=0;i<55;i++)
    {
        random[i] = (double)rand()/(double)RAND_MAX;
    }

    pv = malloc(sizeof(double)*n_ion_populations);
    ienergy = (double*) pv;
    pv = malloc(sizeof(double)*n_ion_populations);
    ienergy_deleted = (double*) pv;
    for (int i=0; i<n_ion_populations; ++i)
    {
        ienergy_deleted[i] = 0.0;
    }

    pv = malloc(sizeof(int)*n_ion_populations);
    N_qp_i = (int*) pv;
}

spatial_region::~spatial_region()
{
    // erasing of particles
    while (p_lapwpo!=0) {
        pwpo* a;
        a = p_lapwpo;
        p_lapwpo = p_lapwpo->previous;
        free((void*) a->head);
        //free(a->head);
        delete a;
    }
    while (p_lapwpa!=0) {
        pwpa* a;
        a = p_lapwpa;
        p_lapwpa = p_lapwpa->previous;
        free((void*) a->head);
        //free(a->head);
        delete a;
    }

    void* pv;

    pv = (void*) random;
    free(pv);

    pv = (void*) ienergy;
    free(pv);
    pv = (void*) ienergy_deleted;
    free(pv);

    pv = (void*) N_qp_i;
    free(pv);
}

particle* spatial_region::new_particle()
{
    if (n_f == 0) { // create particle from scratch
        if (n_ap+1 > npwpa*((int) page_size/particle_size))
        {
            n_ap++;
            npwpa++;
            pwpa* a;
            a = new pwpa;
            a->previous = p_lapwpa;
            a->head = (particle*) malloc(page_size);
            //a->head = (plist::particle*) malloc(page_size);
            p_lapwpa = a;
            p_lap = a->head;
            return p_lap;
        }
        else
        {
            n_ap++;
            p_lap++;
            return p_lap;
        }
    } else { // reuse memory that is left from another particle
        if (n_f-1 == (npwpo-1)*((int) page_size/pointer_size))
        {
            n_f--;
            particle* b;
            b = *pp_lfp;

            if (p_lapwpo == 0)
            {
                cout << "Error - quill didn't create page with free particles on startup" << endl;
            }

            if (npwpo > 1) // don't delete the last page with pointers
            {
                npwpo--;
                free((void*) pp_lfp);
                //free(pp_lfp);
                pwpo* a;
                a = p_lapwpo;
                p_lapwpo = p_lapwpo->previous;
                delete a;
                pp_lfp = p_lapwpo->head + (int) page_size/pointer_size - 1;
            }
            else if (n_f != 0)
            {
                cout << "Error - inconsistent n_f and pp_lfp, program may crash!" << endl;
            }

            return b;
        }
        else
        {
            n_f--;
            pp_lfp--;
            return *(pp_lfp+1);
        }
    }
}

void spatial_region::delete_particle(particle* a, bool force_delete)
{
    if ((catching_enabled && !spatial_region::is_inside_global(a->x, a->y, a->z) &&
        !spatial_region::is_in_exchange_area(a->x)) || force_delete)
    {
        #ifdef QUILL_NOQED
        spatial_region::deleted_particle dparticle(a->cmr, a->q, a->x, a->y, a->z, a->ux, a->uy, a->uz, a->g);
        #else
        spatial_region::deleted_particle dparticle(a->cmr, a->q, a->x, a->y, a->z, a->ux, a->uy, a->uz, a->g, a->chi);
        #endif
        spatial_region::deleted_particles.push_back(dparticle);
        update_energy_deleted(a);
    }
    
    if (n_f == 0)
    {
        n_f = 1;
        pp_lfp = p_lapwpo->head;
        *pp_lfp = a;
    }
    else if (n_f+1 > npwpo*((int) page_size/pointer_size))
    {
        n_f++;
        npwpo++;
        pwpo* b;
        b = new pwpo;
        b->previous = p_lapwpo;
        b->head = (particle**) malloc(page_size);
        //b->head = (plist::particle**) malloc(page_size);
        p_lapwpo = b;
        pp_lfp = b->head;
        *pp_lfp = a;
    }
    else
    {
        n_f++;
        pp_lfp++;
        *pp_lfp = a;
    }
}

void spatial_region::update_energy_deleted(particle* a)
{
    double norm = 8.2e-14*dx*dy*dz*1.11485e13*lambda/(8*PI*PI*PI); // energy in Joules
    if (a->cmr == -1)
    {
        spatial_region::energy_e_deleted -= a->q * (a->g - 1) * norm; // energy should be positive in file
    }
    else if (a->cmr == 1)
    {
        spatial_region::energy_p_deleted += a->q * (a->g - 1) * norm;
    }
    else if (a->cmr == 0)
    {
        spatial_region::energy_ph_deleted += a->q * a->g * norm;
    }
    else
    {
        for (int j=0; j<n_ion_populations; ++j)
        {
            if (a->cmr == icmr[j])
            {
                spatial_region::ienergy_deleted[j] += a->q * (a->g - 1) * norm / icmr[j];
                break;
            }
        }
    }
}

vector3d regulate(double& a, double& b, double& c)
{
    vector3d v;
    if(a>b)
    {
        if(b>c)
        {
            v.x = a;
            v.y = b;
            v.z = c;
        }
        else
        {
            v.z = b;
            if(a>c)
            {
                v.x = a;
                v.y = c;
            }
            else
            {
                v.x = c;
                v.y = a;
            }
        }
    }
    else
    {
        if(c>b)
        {
            v.x = c;
            v.y = b;
            v.z = a;
        }
        else
        {
            v.x = b;
            if(a>c)
            {
                v.y = a;
                v.z = c;
            }
            else
            {
                v.y = c;
                v.z = a;
            }
        }
    }
    return v;
}

var::var()
{
    name="";
    value=0;
    units="";
    next = 0;
}

int var::read()
{
    std::string tmp;
    std::string stmp;
    cin>>name;
    if (name!="$")
    {
        if (mpi_rank == 0) {
            cout<<name;
        }
        cin>>tmp;
        if (tmp[0]=='['){ //detecting the start of an array
            int num=1; //character counter
            int incr=0; //counter for calculating the length of the number
            while (tmp[num]!=']'){ // checking the end of an array
                if (tmp[num]!=';'){ // parsing
                    incr++;
                }
                else {
                    char tmp_1[100]={};
                    for (int i=0;i<incr;i++)
                        tmp_1[i]=tmp[num-incr+i];
                    input_array.push_back(atof(tmp_1)); //adding a new number to a vector
                    incr=0;
                }
                if (tmp[num+1]==']'){
                    char tmp_1[100]={};
                    int i=0;
                    while (i<incr){
                       tmp_1[i]=tmp[num-incr+1+i];
                       i++;
                    }
                    input_array.push_back(atof(tmp_1));
                    incr=0;
                }
                num++;
            }
            if (mpi_rank == 0) {
                for (auto v : input_array)
                    cout << "..." << v;
            }
        }
        else{
            value = (tmp != "#" ? stof(tmp) : 0.0);
            if (mpi_rank == 0) {
                cout<<"..."<<value;
            }
        }
        cin>>units;
        if (mpi_rank == 0) {
            cout<<"..."<<units<<"...";
        }
        cin>>stmp;
        if (stmp=="$")
        {
            if (mpi_rank == 0) {
                cout<<"\033[1m\033[32m"<<"ok"<<"\033[0m\n";
            }
            return 0;
        }
        else
        {
            if (mpi_rank == 0) {
                cout<<"\033[1m\033[31m[!]\033[0m\n";
            }
            return 1;
        }
    }
    else
    {
        if (mpi_rank == 0) {
            cout<<"...reading finished\n";
        }
        return 1;
    }
}

var* find(std::string a, var* b)
{
    var* tmp = b;
    while (tmp->next!=0) /* последняя переменная — пустая, поэтому её
                            проверять не нужно */
    {
        if(tmp->name==a)
            return tmp;
        tmp = tmp->next;
    }
    if (tmp->value != 0.0 || !tmp->units.empty()) {
        cerr << TERM_RED << "Warning! The \"empty\" default config variable was modified to [" << tmp->value <<
                " " << tmp->units << "]. Reverting to [0.0]" << TERM_NO_COLOR << endl;
        tmp->value = 0.0;
        tmp->units = "";
    }
    return tmp; /* «пустая» переменная в силу алгоритма чтения файла
                 */
}

film::film()
{
    prev = 0;
    x0 = 0;
    filmwidth = 0;
    gradwidth = 0;
    y0 = 0;
    y1 = 0;
    z0 = 0;
    z1 = 0;
    ne = 0;
    mcr = 1;
    T = 0;
    vx = 0;
    xnpic_film = 0;
    ynpic_film = 0;
    znpic_film = 0;
}

double lin_interpolation(double coordinate, std::vector<double>& density_coords, std::vector<double>& density_values)
{
    if (density_coords.empty()) {
        return 1.0;
    }

    if (coordinate < density_coords[0]) {
        return density_values[0];
    } else if (coordinate > density_coords[density_coords.size()-1]) {
        return density_values[density_values.size()-1];
    } else {
        int i = lower_bound(density_coords.begin(), density_coords.end(), coordinate) - density_coords.begin();
        const double left = density_values[i-1];
        const double right = density_values[i];
        const double x_rel = (coordinate - density_coords[i-1]);
        const double dx = (density_coords[i] - density_coords[i-1]);
        return left + (right - left) * x_rel / dx;
    }
}

double spatial_region::get_max_w() {
    return get_max_w(0, nx);
}

double spatial_region::get_max_w(double left, double right) {
    double max_w = 0.0;
    for (int i = left; i < right; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                celle & e = ce[i][j][k];
                cellb & b = cb[i][j][k];
                // fields are not interpolated to the same point, but should be fine for estimation
                double w = e.ex * e.ex + e.ey * e.ey + e.ez * e.ez + b.bx * b.bx + b.by * b.by + b.bz * b.bz;
                if (w > max_w) {
                    max_w = w;
                }
            }
        }
    }
    return max_w;
}
