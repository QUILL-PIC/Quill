#include "containers.h"
#include "compilation_defines.h"

template <typename T>
field3d<T>::field3d(int nx0, int ny0, int nz0) : nx(nx0), ny(ny0), nz(nz0) {

    p = std::unique_ptr<T**[]>(new T**[nx]());

    p2 = std::unique_ptr<T*[]>(new T*[nx * ny]());
    for(int i=0; i<nx; i++)
    {
        p[i] = &(p2[ny*i]);
    }

    p3 = std::unique_ptr<T[]>(new T[nx * ny * nz]());
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            p[i][j] = &(p3[nz*ny*i + nz*j]);
        }
    }
}

template <typename T>
field3d<T> & field3d<T>::operator =(field3d<T> && other) {
    nx = other.nx;
    ny = other.ny;
    nz = other.nz;
    p = std::move(other.p);
    p2 = std::move(other.p2);
    p3 = std::move(other.p3);
    other.nx = 0;
    other.ny = 0;
    other.nz = 0;

    return *this;
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
    #ifndef QUILL_NOQED
    chi = 0;
    #endif
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
