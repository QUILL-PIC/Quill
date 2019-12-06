#include <cmath>
#include "pusher.h"

void push_vay(thinparticle &p,
              double ex, double ey, double ez,
              double bx, double by, double bz,
              double dt)
{
    // See PoP, J.-L. Vay, 2008
    double tmp;
    double wx,wy,wz;
    double bw;
    double & cmr = p.cmr;
    double & ux = p.ux;
    double & uy = p.uy;
    double & uz = p.uz;
    double & g = p.g;
    tmp = dt*cmr;
    double ex1 = ex*tmp;
    double ey1 = ey*tmp;
    double ez1 = ez*tmp;
    //
    tmp = 0.5*tmp;
    double bx1 = bx*tmp;
    double by1 = by*tmp;
    double bz1 = bz*tmp;
    //
    wx = ux + ex1 + (uy*bz1-uz*by1)/g;
    wy = uy + ey1 + (uz*bx1-ux*bz1)/g;
    wz = uz + ez1 + (ux*by1-uy*bx1)/g;
    //
    tmp = bx1*bx1 + by1*by1 + bz1*bz1;
    bw = bx1*wx+by1*wy+bz1*wz;
    g = 1 + wx*wx + wy*wy + wz*wz - tmp;
    g = sqrt( 0.5*(g+sqrt(g*g+4*(tmp+bw*bw))) );
    //
    tmp = 1/(1+tmp/(g*g));
    bw = bw/(g*g);
    ux = (wx+(wy*bz1-wz*by1)/g+bx1*bw)*tmp;
    uy = (wy+(wz*bx1-wx*bz1)/g+by1*bw)*tmp;
    uz = (wz+(wx*by1-wy*bx1)/g+bz1*bw)*tmp;
}

void push_boris(thinparticle &p,
                double ex, double ey, double ez,
                double bx, double by, double bz,
                double dt)
{
    // See Birdsall, Langdon, "Plasma Physics via Computer Simulation"
    double tmp;
    double u1x, u1y, u1z;
    double wx,wy,wz;
    double & cmr = p.cmr;
    double & ux = p.ux;
    double & uy = p.uy;
    double & uz = p.uz;
    double & g = p.g;
    tmp = 0.5 * dt * cmr;

    //e = dt qE/2m
    double ex1 = ex*tmp;
    double ey1 = ey*tmp;
    double ez1 = ez*tmp;

    // u- = u0 + e
    u1x = ux + ex1;
    u1y = uy + ey1;
    u1z = uz + ez1;

    // gamma^2 = 1 + (u-)^2
    g = sqrt(1 + u1x * u1x + u1y * u1y + u1z * u1z);

    tmp /= g;
    // b = dt qB/2mc
    double bx1 = bx*tmp;
    double by1 = by*tmp;
    double bz1 = bz*tmp;

    // u' = u- + [u-, b]
    wx = u1x + (u1y * bz1 - u1z * by1);
    wy = u1y + (u1z * bx1 - u1x * bz1);
    wz = u1z + (u1x * by1 - u1y * bx1);

    tmp = 2 / (1 + bx1 * bx1 + by1 * by1 + bz1 * bz1);

    // u+ = u- + 2 * [u, b] / (1 + b^2)
    u1x += (wy * bz1 - wz * by1) * tmp;
    u1y += (wz * bx1 - wx * bz1) * tmp;
    u1z += (wx * by1 - wy * bx1) * tmp;

    // u1 = u+ + e
    ux = u1x + ex1;
    uy = u1y + ey1;
    uz = u1z + ez1;

    g = sqrt(1 + ux * ux + uy * uy + uz * uz);
}
