#include <iostream>
#include <algorithm>
#include <numa.h>
#include "main.h"

extern bool catching_enabled;
extern double lambda;

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
    ce = 0;
    cj = 0;
    cb = 0;
    cbe = 0;
    cp = 0;
    n_random = 0;
    energy_e_deleted = 0;
    energy_p_deleted = 0;
    energy_ph_deleted = 0;
    //
    n_ap = 0;
    n_f = 0;
    npwpa = 0;
    npwpo = 0;
    //page_size = numa_pagesize();
    // page_size should not be very small to avoid OS problems with allocation
    const size_t PS = 64 * 1024 * 1024; // bytes
    if (numa_pagesize() < PS) {
        page_size = PS;
    } else {
        page_size = numa_pagesize();
    }
    particle_size = sizeof(plist::particle);
    pointer_size = sizeof(plist::particle*);
    p_lap = 0;
    pp_lfp = 0;
    p_lapwpa = 0;
    p_lapwpo = 0;
    n_ion_populations = 0;
    icmr = 0;
    data_folder = "results";
}

spatial_region::plist::particle::particle()
{
    x=0;
    y=0;
    z=0;
    ux=0;
    uy=0;
    uz=0;
    g=1;
    q=0;
    cmr=-1;
    next=0;
    previous=0;
    chi = 0;
    trn = 0;
}

spatial_region::plist::plist()
{
    head=0;
    start=0;
}

void spatial_region::init(int sr_id0, double dx0, double dy0, double dz0, double dt0, double e_s0, int xnpic0, int ynpic0, int znpic0, int node_number0, int n_ion_populations0, double* icmr0, std::string df)
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
    node_number = node_number0;
    n_ion_populations = n_ion_populations0;
    icmr = icmr0;
    data_folder = df;

    // we should always have at least one page for free particles (to avoid performance issues)
    pwpo* b = new pwpo;
    b->previous = 0;
    b->head = (plist::particle**) numa_alloc_onnode(page_size,node_number);
    npwpo++;
    p_lapwpo = b;
}

void spatial_region::create_arrays(int nx0, int ny0, int nz0, int seed, int node_number)
{
    nx = nx0;
    ny = ny0;
    nz = nz0;

    void* pv;

    pv = numa_alloc_onnode(sizeof(celle**)*nx,node_number);
    ce = (celle***) pv;
    pv = numa_alloc_onnode(sizeof(celle*)*nx*ny,node_number);
    for(int i=0;i<nx;i++)
    {
        ce[i] = ((celle**) pv) + ny*i;
    }
    pv = numa_alloc_onnode(sizeof(celle)*nx*ny*nz,node_number);
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            ce[i][j] = ((celle*) pv) + nz*ny*i + nz*j;
        }
    }

    pv = numa_alloc_onnode(sizeof(cellb**)*nx,node_number);
    cb = (cellb***) pv;
    pv = numa_alloc_onnode(sizeof(cellb*)*nx*ny,node_number);
    for(int i=0;i<nx;i++)
    {
        cb[i] = ((cellb**) pv) + ny*i;
    }
    pv = numa_alloc_onnode(sizeof(cellb)*nx*ny*nz,node_number);
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            cb[i][j] = ((cellb*) pv) + nz*ny*i + nz*j;
        }
    }

    pv = numa_alloc_onnode(sizeof(cellj**)*nx,node_number);
    cj = (cellj***) pv;
    pv = numa_alloc_onnode(sizeof(cellj*)*nx*ny,node_number);
    for(int i=0;i<nx;i++)
    {
        cj[i] = ((cellj**) pv) + ny*i;
    }
    pv = numa_alloc_onnode(sizeof(cellj)*nx*ny*nz,node_number);
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            cj[i][j] = ((cellj*) pv) + nz*ny*i + nz*j;
        }
    }

    pv = numa_alloc_onnode(sizeof(double***)*n_ion_populations,node_number);
    irho = (double****) pv;
    pv = numa_alloc_onnode(sizeof(double**)*n_ion_populations*nx,node_number);
    for (int n=0;n<n_ion_populations;n++)
    {
        irho[n] = ((double***) pv) + n*nx;
    }
    pv = numa_alloc_onnode(sizeof(double*)*n_ion_populations*nx*ny,node_number);
    for (int n=0;n<n_ion_populations;n++)
    {
        for (int i=0;i<nx;i++)
            irho[n][i] = ((double**) pv) + n*nx*ny + i*ny;
    }
    pv = numa_alloc_onnode(sizeof(double)*n_ion_populations*nx*ny*nz,node_number);
    for (int n=0;n<n_ion_populations;n++)
    {
        for (int i=0;i<nx;i++)
        {
            for (int j=0;j<ny;j++)
                irho[n][i][j] = ((double*) pv) + n*nx*ny*nz + i*ny*nz + j*nz;
        }
    }

    pv = numa_alloc_onnode(sizeof(cellbe**)*nx,node_number);
    cbe = (cellbe***) pv;
    pv = numa_alloc_onnode(sizeof(cellbe*)*nx*ny,node_number);
    for(int i=0;i<nx;i++)
    {
        cbe[i] = ((cellbe**) pv) + ny*i;
    }
    pv = numa_alloc_onnode(sizeof(cellbe)*nx*ny*nz,node_number);
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            cbe[i][j] = ((cellbe*) pv) + nz*ny*i + nz*j;
        }
    }

    pv = numa_alloc_onnode(sizeof(cellp**)*nx,node_number);
    cp = (cellp***) pv;
    pv = numa_alloc_onnode(sizeof(cellp*)*nx*ny,node_number);
    for(int i=0;i<nx;i++)
    {
        cp[i] = ((cellp**) pv) + ny*i;
    }
    pv = numa_alloc_onnode(sizeof(cellp)*nx*ny*nz,node_number);
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            cp[i][j] = ((cellp*) pv) + nz*ny*i + nz*j;
        }
    }

    pv = numa_alloc_onnode(sizeof(double)*55,node_number);
    random = (double*) pv;
    srand(seed);
    for(int i=0;i<55;i++)
    {
        random[i] = (double)rand()/(double)RAND_MAX;
    }

    pv = numa_alloc_onnode(sizeof(double)*n_ion_populations,node_number);
    ienergy = (double*) pv;
    pv = numa_alloc_onnode(sizeof(double)*n_ion_populations,node_number);
    ienergy_deleted = (double*) pv;
    for (int i=0; i<n_ion_populations; ++i)
    {
        ienergy_deleted[i] = 0.0;
    }

    pv = numa_alloc_onnode(sizeof(int)*n_ion_populations,node_number);
    N_qp_i = (int*) pv;
}

spatial_region::~spatial_region()
{
    // erasing of particles
    while (p_lapwpo!=0) {
        pwpo* a;
        a = p_lapwpo;
        p_lapwpo = p_lapwpo->previous;
        numa_free((void*) a->head,page_size);
        //free(a->head);
        delete a;
    }
    while (p_lapwpa!=0) {
        pwpa* a;
        a = p_lapwpa;
        p_lapwpa = p_lapwpa->previous;
        numa_free((void*) a->head,page_size);
        //free(a->head);
        delete a;
    }

    void* pv;

    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            pv = (void*) ce[i][j];
            numa_free(pv, sizeof(celle)*nz);
        }
        pv = (void*) ce[i];
        numa_free(pv, sizeof(celle*)*ny);
    }
    pv = (void*) ce;
    numa_free(pv, sizeof(celle**)*nx);

    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            pv = (void*) cb[i][j];
            numa_free(pv, sizeof(cellb)*nz);
        }
        pv = (void*) cb[i];
        numa_free(pv, sizeof(cellb*)*ny);
    }
    pv = (void*) cb;
    numa_free(pv, sizeof(cellb**)*nx);

    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            pv = (void*) cj[i][j];
            numa_free(pv, sizeof(cellj)*nz);
        }
        pv = (void*) cj[i];
        numa_free(pv, sizeof(cellj*)*ny);
    }
    pv = (void*) cj;
    numa_free(pv, sizeof(cellj**)*nx);

    for (int n=0;n<n_ion_populations;n++)
    {
        for (int i=0;i<nx;i++)
        {
            for (int j=0;j<ny;j++)
            {
                pv = (void*) irho[n][i][j];
                numa_free(pv, sizeof(double)*nz);
            }
        }
    }
    for (int n=0;n<n_ion_populations;n++)
    {
        for (int i=0;i<nx;i++)
        {
            pv = (void*) irho[n][i];
            numa_free(pv, sizeof(double*)*ny);
        }
    }
    for (int n=0;n<n_ion_populations;n++)
    {
        pv = (void*) irho[n];
        numa_free(pv, sizeof(double**)*nx);
    }
    pv = (void*) irho;
    numa_free(pv, sizeof(double***)*n_ion_populations);

    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            pv = (void*) cbe[i][j];
            numa_free(pv, sizeof(cellbe)*nz);
        }
        pv = (void*) cbe[i];
        numa_free(pv, sizeof(cellbe*)*ny);
    }
    pv = (void*) cbe;
    numa_free(pv, sizeof(cellbe**)*nx);

    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            pv = (void*) cp[i][j];
            numa_free(pv, sizeof(cellp)*nz);
        }
        pv = (void*) cp[i];
        numa_free(pv, sizeof(cellp*)*ny);
    }
    pv = (void*) cp;
    numa_free(pv, sizeof(cellp**)*nx);

    pv = (void*) random;
    numa_free(pv, sizeof(double)*55);

    pv = (void*) ienergy;
    numa_free(pv, sizeof(double*)*n_ion_populations);
    pv = (void*) ienergy_deleted;
    numa_free(pv, sizeof(double*)*n_ion_populations);

    pv = (void*) N_qp_i;
    numa_free(pv, sizeof(int*)*n_ion_populations);
}

spatial_region::plist::particle* spatial_region::new_particle()
{
    if (n_f == 0) { // create particle from scratch
        if (n_ap+1 > npwpa*((int) page_size/particle_size))
        {
            n_ap++;
            npwpa++;
            pwpa* a;
            a = new pwpa;
            a->previous = p_lapwpa;
            a->head = (plist::particle*) numa_alloc_onnode(page_size,node_number);
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
            plist::particle* b;
            b = *pp_lfp;

            if (p_lapwpo == 0)
            {
                cout << "Error - quill didn't create page with free particles on startup" << endl;
            }

            if (npwpo > 1) // don't delete the last page with pointers
            {
                npwpo--;
                numa_free((void*) pp_lfp,page_size);
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

void spatial_region::delete_particle(plist::particle* a)
{
    if (catching_enabled)
    {
        if (!spatial_region::is_inside_global(a->x, a->y, a->z) && !spatial_region::is_in_exchange_area(a->x, a->y, a->z))
        {
            spatial_region::deleted_particle dparticle(a->cmr, a->q, a->x, a->y, a->z, a->ux, a->uy, a->uz, a->g, a->chi);
            spatial_region::deleted_particles.push_back(dparticle);
            update_energy_deleted(a);
        }
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
        b->head = (plist::particle**) numa_alloc_onnode(page_size,node_number);
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

void spatial_region::update_energy_deleted(plist::particle* a)
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

vector3d::vector3d()
{
    x=0;
    y=0;
    z=0;
}

int_vector3d::int_vector3d()
{
    i=0;
    j=0;
    k=0;
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
        cout<<name;
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
            for (int i=0;i<input_array.size();i++)
                    cout<<"..."<<input_array.at(i);
        }
        else{
            value = (tmp != "#" ? stof(tmp) : 0.0);
            cout<<"..."<<value;
        }
        cin>>units;
        cout<<"..."<<units<<"...";
        cin>>stmp;
        if (stmp=="$")
        {
            cout<<"\033[1m\033[32m"<<"ok"<<"\033[0m\n";
            return 0;
        }
        else
        {
            cout<<"\033[1m\033[31m[!]\033[0m\n";
            return 1;
        }
    }
    else
    {
        cout<<"...reading finished\n";
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
