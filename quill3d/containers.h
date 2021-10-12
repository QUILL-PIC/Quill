#ifndef CONTAINERS_H_
#define CONTAINERS_H_

#include <memory>
#include "compilation_defines.h"

struct vector3d
{
    double x,y,z;
    vector3d() : x(0), y(0), z(0) {}
};

struct int_vector3d
{
    int i,j,k;
    int_vector3d() : i(0), j(0), k(0) {}
};

struct celle
{
    double ex,ey,ez;
};
struct cellj
{
    double jx,jy,jz;
};
struct cellb
{
    double bx,by,bz;
};
struct cellbe
{
    double bex,bey,bez;
};

template <typename T>
class field3d
{
public:
    field3d() : nx(0), ny(0), nz(0), p(nullptr), p2(nullptr), p3(nullptr) {}
    field3d(field3d<T>&) = delete;
    field3d(int nx, int ny, int nz);
    inline T** const & operator[](int i) const { return p[i]; }
    field3d<T> & operator=(field3d<T> &) = delete;
    field3d<T> & operator=(field3d<T> &&);
    int get_nx() { return nx; }
    int get_ny() { return ny; }
    int get_nz() { return nz; }
private:
    int nx, ny, nz;
    std::unique_ptr<T**[]> p;
    std::unique_ptr<T*[]> p2;
    std::unique_ptr<T[]> p3;
};

struct particle
{
    double x,y,z,ux,uy,uz,g,q,cmr; // cmr - charge to mass ratio,
    #ifndef QUILL_NOQED
    double chi;
    #endif
    int trn; // track name
    particle* next;
    particle* previous;
    particle();
    vector3d get_displacement(double);
    void coordinate_advance(vector3d&);
};

class plist
{
    public:
        particle* head;
        particle* start;
        plist();
        void xplus(double);
};

struct cellp
{
    plist pl;
};

enum class maxwell_solver_enum {
    FDTD, NDFX
};

enum class pusher_enum {
    VAY, BORIS
};


enum class moving_window {
    OFF, ON, AUTO
};

#endif /* CONTAINERS_H_ */
