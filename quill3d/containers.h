#ifndef CONTAINERS_H_
#define CONTAINERS_H_

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
    celle();
};
struct cellj
{
    double jx,jy,jz;
    cellj();
};
struct cellb
{
    double bx,by,bz;
    cellb();
};
struct cellbe
{
    double bex,bey,bez;
    cellbe();
};

template <typename T>
class field3d
{
public:
    field3d() : nx(0), ny(0), nz(0), p(nullptr) {}
    field3d(field3d<T>&) = delete;
    field3d(int nx, int ny, int nz, int node_number);
    ~field3d();
    inline T** const & operator[](int i) const { return p[i]; }
    field3d<T> & operator=(field3d<T> &) = delete;
    field3d<T> & operator=(field3d<T> &&);
    int get_nx() { return nx; }
    int get_ny() { return ny; }
    int get_nz() { return nz; }
private:
    int nx, ny, nz;
    T*** p;
    void free();
};

#endif /* CONTAINERS_H_ */
