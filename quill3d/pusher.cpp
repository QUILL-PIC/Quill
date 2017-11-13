#include <cmath>
#include "main.h"
#include "pusher.h"

void pusher_vay(spatial_region::plist::particle* p, vector3d& e_field, vector3d& b_field, double& dt)
{
    // See PoP, J.-L. Vay, 2008
    vector3d e;
    vector3d b;
    double tmp;
    double wx,wy,wz;
    double bw;
    double & cmr = p->cmr;
    double & ux = p->ux;
    double & uy = p->uy;
    double & uz = p->uz;
    double & g = p->g;
    tmp = dt*cmr;
    e.x = e_field.x*tmp;
    e.y = e_field.y*tmp;
    e.z = e_field.z*tmp;
    //
    tmp = 0.5*tmp;
    b.x = b_field.x*tmp;
    b.y = b_field.y*tmp;
    b.z = b_field.z*tmp;
    //
    wx = ux + e.x + (uy*b.z-uz*b.y)/g;
    wy = uy + e.y + (uz*b.x-ux*b.z)/g;
    wz = uz + e.z + (ux*b.y-uy*b.x)/g;
    //
    tmp = b.x*b.x + b.y*b.y + b.z*b.z;
    bw = b.x*wx+b.y*wy+b.z*wz;
    g = 1 + wx*wx + wy*wy + wz*wz - tmp;
    g = sqrt( 0.5*(g+sqrt(g*g+4*(tmp+bw*bw))) );
    //
    tmp = 1/(1+tmp/(g*g));
    bw = bw/(g*g);
    ux = (wx+(wy*b.z-wz*b.y)/g+b.x*bw)*tmp;
    uy = (wy+(wz*b.x-wx*b.z)/g+b.y*bw)*tmp;
    uz = (wz+(wx*b.y-wy*b.x)/g+b.z*bw)*tmp;
}

void pusher_boris(spatial_region::plist::particle* p, vector3d& e_field, vector3d& b_field, double& dt)
{
    // See Birdsall, Langdon, "Plasma Physics via Computer Simulation"
    vector3d e;
    vector3d b;
    double tmp;
    double u1x, u1y, u1z;
    double wx,wy,wz;
    double & cmr = p->cmr;
    double & ux = p->ux;
    double & uy = p->uy;
    double & uz = p->uz;
    double & g = p->g;
    tmp = 0.5 * dt * cmr;

    //e = dt qE/2m
    e.x = e_field.x*tmp;
    e.y = e_field.y*tmp;
    e.z = e_field.z*tmp;
    //

    // u- = u0 + e
    u1x = ux + e.x;
    u1y = uy + e.y;
    u1z = uz + e.z;

    // gamma^2 = 1 + (u-)^2
    g = sqrt(1 + u1x * u1x + u1y * u1y + u1z * u1z);

    tmp /= g;
    // b = dt qB/2mc
    b.x = b_field.x*tmp;
    b.y = b_field.y*tmp;
    b.z = b_field.z*tmp;

    // u' = u- + [u-, b]
    wx = u1x + (u1y * b.z - u1z * b.y);
    wy = u1y + (u1z * b.x - u1x * b.z);
    wz = u1z + (u1x * b.y - u1y * b.x);

    tmp = 2 / (1 + b.x * b.x + b.y * b.y + b.z * b.z);

    // u+ = u- + 2 * [u, b] / (1 + b^2)
    u1x += (wy * b.z - wz * b.y) * tmp;
    u1y += (wz * b.x - wx * b.z) * tmp;
    u1z += (wx * b.y - wy * b.x) * tmp;

    // u1 = u+ + e
    ux = u1x + e.x;
    uy = u1y + e.y;
    uz = u1z + e.z;

    g = sqrt(1 + ux * ux + uy * uy + uz * uz);
}
