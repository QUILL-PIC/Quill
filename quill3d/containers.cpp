#include <cstdlib>
#include "containers.h"

template <typename T>
field3d<T>::field3d(int nx0, int ny0, int nz0) : nx(nx0), ny(ny0), nz(nz0) {

    void* pv;

    pv = malloc(sizeof(T**)*nx);
    p = (T***) pv;

    pv = malloc(sizeof(T*)*nx*ny);
    for(int i=0; i<nx; i++)
    {
        p[i] = ((T**) pv) + ny*i;
    }

    pv = malloc(sizeof(T)*nx*ny*nz);
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            p[i][j] = ((T*) pv) + nz*ny*i + nz*j;
        }
    }
}

template <typename T>
field3d<T>::~field3d() {
    free_memory();
}

template <typename T>
field3d<T> & field3d<T>::operator =(field3d<T> && other) {
    free_memory();
    nx = other.nx;
    ny = other.ny;
    nz = other.nz;
    p = other.p;
    other.p = nullptr;
    other.nx = 0;
    other.ny = 0;
    other.nz = 0;

    return *this;
}

template <typename T>
void field3d<T>::free_memory() {
    if (!p) {
        return;
    }

    void* pv;

    pv = (void*) p[0][0];
    free(pv);

    pv = (void*) p[0];
    free(pv);

    pv = (void*) p;
    free(pv);
}

particle::particle()
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

vector3d particle::get_displacement(double dt)
{
    vector3d displacement;
    double a;
    a = dt/g;
    displacement.x = ux*a;
    displacement.y = uy*a;
    displacement.z = uz*a;
    return displacement;
}

void particle::coordinate_advance(vector3d& a)
{
    x = x + a.x;
    y = y + a.y;
    z = z + a.z;
}

plist::plist() : head(nullptr), start(nullptr) { }

void plist::xplus(double x_adjunct)
{
    particle* current;
    current = head;
    while(current!=0)
    {
        current->x = current->x + x_adjunct;
        current = current->next;
    }
}

template class field3d<double>;
template class field3d<celle>;
template class field3d<cellb>;
template class field3d<cellj>;
template class field3d<cellbe>;
template class field3d<cellp>;
