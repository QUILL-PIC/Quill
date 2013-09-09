#include <cmath>
#include "main.h"

int max(int a, int b)
{
    if(a>b)
        return a;
    else
        return b;
}

void spatial_region::simple_jdep(spatial_region::plist::particle& p, vector3d& d, int_vector3d& w)
{
    // d - displacement, w - walls for deposition, p - initial position of particle
    double xa,xb,ya,yb,za,zb;
    xa = p.x - w.i;
    xb = 1 - xa;
    ya = p.y - w.j;
    yb = 1 - ya;
    za = p.z - w.k;
    zb = 1 - za;
    //
    cj[w.i][w.j+1][w.k+1].jx = cj[w.i][w.j+1][w.k+1].jx + d.x*p.q*(ya*za + 0.5*(d.z*ya+d.y*za) + d.y*d.z/3);
    cj[w.i][w.j][w.k+1].jx = cj[w.i][w.j][w.k+1].jx + d.x*p.q*(yb*za + 0.5*(d.z*yb-d.y*za) - d.y*d.z/3);
    cj[w.i][w.j+1][w.k].jx = cj[w.i][w.j+1][w.k].jx + d.x*p.q*(ya*zb + 0.5*(-d.z*ya+d.y*zb) - d.y*d.z/3);
    cj[w.i][w.j][w.k].jx = cj[w.i][w.j][w.k].jx + d.x*p.q*(yb*zb - 0.5*(d.z*yb+d.y*zb) + d.y*d.z/3);
    //
    cj[w.i+1][w.j][w.k+1].jy = cj[w.i+1][w.j][w.k+1].jy + d.y*p.q*(xa*za + 0.5*(d.x*za+d.z*xa) + d.x*d.z/3);
    cj[w.i][w.j][w.k+1].jy = cj[w.i][w.j][w.k+1].jy + d.y*p.q*(xb*za + 0.5*(-d.x*za+d.z*xb) - d.x*d.z/3);
    cj[w.i+1][w.j][w.k].jy = cj[w.i+1][w.j][w.k].jy + d.y*p.q*(xa*zb + 0.5*(d.x*zb-d.z*xa) - d.x*d.z/3);
    cj[w.i][w.j][w.k].jy = cj[w.i][w.j][w.k].jy + d.y*p.q*(xb*zb - 0.5*(d.x*zb+d.z*xb) + d.x*d.z/3);
    //
    cj[w.i+1][w.j+1][w.k].jz = cj[w.i+1][w.j+1][w.k].jz + d.z*p.q*(xa*ya + 0.5*(d.x*ya+d.y*xa) + d.x*d.y/3);
    cj[w.i][w.j+1][w.k].jz = cj[w.i][w.j+1][w.k].jz + d.z*p.q*(xb*ya + 0.5*(-d.x*ya+d.y*xb) - d.x*d.y/3);
    cj[w.i+1][w.j][w.k].jz = cj[w.i+1][w.j][w.k].jz + d.z*p.q*(xa*yb + 0.5*(d.x*yb-d.y*xa) - d.x*d.y/3);
    cj[w.i][w.j][w.k].jz = cj[w.i][w.j][w.k].jz + d.z*p.q*(xb*yb - 0.5*(d.x*yb+d.y*xb) + d.x*d.y/3);
    //
}

void spatial_region::jdeposition(spatial_region::plist::particle& p, vector3d& d)
{
    // d - displacement
    int i1,i2,j1,j2,k1,k2;
    int_vector3d w;
    vector3d d2,v;
    double path1,path2,path3;
    i1 = floor(p.x);
    j1 = floor(p.y);
    k1 = floor(p.z);
    i2 = floor(p.x+d.x);
    j2 = floor(p.y+d.y);
    k2 = floor(p.z+d.z);
    if(i1==i2)
    {
        if(j1==j2)
        {
            if(k1==k2)
            {
                // deposition without intersection
                w.i = i1;
                w.j = j1;
                w.k = k1;
                simple_jdep(p,d,w);
                p.coordinate_advance(d);
            }
            else
            {
                // z intersection
                d2.z = max(k1,k2) - p.z;
                d2.x = d.x*d2.z/d.z;
                d2.y = d.y*d2.z/d.z;
                w.i = i1;
                w.j = j1;
                w.k = k1;
                simple_jdep(p,d2,w);
                p.coordinate_advance(d2);
                d2.x = d.x - d2.x;
                d2.y = d.y - d2.y;
                d2.z = d.z - d2.z;
                w.k = k2;
                simple_jdep(p,d2,w);
                p.coordinate_advance(d2);

            }
        }
        else if(k1==k2)
        {
            // y intersection
            d2.y = max(j1,j2) - p.y;
            d2.x = d.x*d2.y/d.y;
            d2.z = d.z*d2.y/d.y;
            w.i = i1;
            w.j = j1;
            w.k = k1;
            simple_jdep(p,d2,w);
            p.coordinate_advance(d2);
            d2.x = d.x - d2.x;
            d2.y = d.y - d2.y;
            d2.z = d.z - d2.z;
            w.j = j2;
            simple_jdep(p,d2,w);
            p.coordinate_advance(d2);
        }
        else
        {
            // y & z intersection
            d2.y = max(j1,j2) - p.y;
            d2.z = max(k1,k2) - p.z;
            path1 = d2.y/d.y;
            path2 = d2.z/d.z;
            if(path1>path2)
            {
                // z intersection happens before y intersection
                d2.x = d.x*path2;
                d2.y = d.y*path2;
                d2.z = d.z*path2;
                w.i = i1;
                w.j = j1;
                w.k = k1;
                simple_jdep(p,d2,w);
                p.coordinate_advance(d2);
                //
                d2.x = d.x*(path1-path2);
                d2.y = d.y*(path1-path2);
                d2.z = d.z*(path1-path2);
                w.k = k2;
                simple_jdep(p,d2,w);
                p.coordinate_advance(d2);
                //
                d2.x = d.x*(1-path1);
                d2.y = d.y*(1-path1);
                d2.z = d.z*(1-path1);
                w.j = j2;
                simple_jdep(p,d2,w);
                p.coordinate_advance(d2);

            }
            else
            {
                // y intersection happens before z intersection
                d2.x = d.x*path1;
                d2.y = d.y*path1;
                d2.z = d.z*path1;
                w.i = i1;
                w.j = j1;
                w.k = k1;
                simple_jdep(p,d2,w);
                p.coordinate_advance(d2);
                //
                d2.x = d.x*(path2-path1);
                d2.y = d.y*(path2-path1);
                d2.z = d.z*(path2-path1);
                w.j = j2;
                simple_jdep(p,d2,w);
                p.coordinate_advance(d2);
                //
                d2.x = d.x*(1-path2);
                d2.y = d.y*(1-path2);
                d2.z = d.z*(1-path2);
                w.k = k2;
                simple_jdep(p,d2,w);
                p.coordinate_advance(d2);
            }
        }
    }
    else if(j1==j2)
    {
        if(k1==k2)
        {
            // x intersection
            d2.x = max(i1,i2) - p.x;
            d2.y = d.y*d2.x/d.x;
            d2.z = d.z*d2.x/d.x;
            w.i = i1;
            w.j = j1;
            w.k = k1;
            simple_jdep(p,d2,w);
            p.coordinate_advance(d2);
            d2.x = d.x - d2.x;
            d2.y = d.y - d2.y;
            d2.z = d.z - d2.z;
            w.i = i2;
            simple_jdep(p,d2,w);
            p.coordinate_advance(d2);
        }
        else
        {
            // x & z intersection
            d2.x = max(i1,i2) - p.x;
            d2.z = max(k1,k2) - p.z;
            path1 = d2.x/d.x;
            path2 = d2.z/d.z;
            if(path1>path2)
            {
                // z intersection happens before x intersection
                d2.x = d.x*path2;
                d2.y = d.y*path2;
                d2.z = d.z*path2;
                w.i = i1;
                w.j = j1;
                w.k = k1;
                simple_jdep(p,d2,w);
                p.coordinate_advance(d2);
                //
                d2.x = d.x*(path1-path2);
                d2.y = d.y*(path1-path2);
                d2.z = d.z*(path1-path2);
                w.k = k2;
                simple_jdep(p,d2,w);
                p.coordinate_advance(d2);
                //
                d2.x = d.x*(1-path1);
                d2.y = d.y*(1-path1);
                d2.z = d.z*(1-path1);
                w.i = i2;
                simple_jdep(p,d2,w);
                p.coordinate_advance(d2);

            }
            else
            {
                // x intersection happens before z intersection
                d2.x = d.x*path1;
                d2.y = d.y*path1;
                d2.z = d.z*path1;
                w.i = i1;
                w.j = j1;
                w.k = k1;
                simple_jdep(p,d2,w);
                p.coordinate_advance(d2);
                //
                d2.x = d.x*(path2-path1);
                d2.y = d.y*(path2-path1);
                d2.z = d.z*(path2-path1);
                w.i = i2;
                simple_jdep(p,d2,w);
                p.coordinate_advance(d2);
                //
                d2.x = d.x*(1-path2);
                d2.y = d.y*(1-path2);
                d2.z = d.z*(1-path2);
                w.k = k2;
                simple_jdep(p,d2,w);
                p.coordinate_advance(d2);
            }
        }
    }
    else if(k1==k2)
    {
        // x & y intersection
        d2.y = max(j1,j2) - p.y;
        d2.x = max(i1,i2) - p.x;
        path1 = d2.y/d.y;
        path2 = d2.x/d.x;
        if(path1>path2)
        {
            // x intersection happens before y intersection
            d2.x = d.x*path2;
            d2.y = d.y*path2;
            d2.z = d.z*path2;
            w.i = i1;
            w.j = j1;
            w.k = k1;
            simple_jdep(p,d2,w);
            p.coordinate_advance(d2);
            //
            d2.x = d.x*(path1-path2);
            d2.y = d.y*(path1-path2);
            d2.z = d.z*(path1-path2);
            w.i = i2;
            simple_jdep(p,d2,w);
            p.coordinate_advance(d2);
            //
            d2.x = d.x*(1-path1);
            d2.y = d.y*(1-path1);
            d2.z = d.z*(1-path1);
            w.j = j2;
            simple_jdep(p,d2,w);
            p.coordinate_advance(d2);
        }
        else
        {
            // y intersection happens before x intersection
            d2.x = d.x*path1;
            d2.y = d.y*path1;
            d2.z = d.z*path1;
            w.i = i1;
            w.j = j1;
            w.k = k1;
            simple_jdep(p,d2,w);
            p.coordinate_advance(d2);
            //
            d2.x = d.x*(path2-path1);
            d2.y = d.y*(path2-path1);
            d2.z = d.z*(path2-path1);
            w.j = j2;
            simple_jdep(p,d2,w);
            p.coordinate_advance(d2);
            //
            d2.x = d.x*(1-path2);
            d2.y = d.y*(1-path2);
            d2.z = d.z*(1-path2);
            w.i = i2;
            simple_jdep(p,d2,w);
            p.coordinate_advance(d2);
        }
    }
    else
    {
        // x & y & z intersection
        d2.x = max(i1,i2) - p.x;
        d2.y = max(j1,j2) - p.y;
        d2.z = max(k1,k2) - p.z;
        path1 = d2.x/d.x;
        path2 = d2.y/d.y;
        path3 = d2.z/d.z;
        v = regulate(path1,path2,path3);
        //
        d2.x = d.x*v.z;
        d2.y = d.y*v.z;
        d2.z = d.z*v.z;
        w.i = i1;
        w.j = j1;
        w.k = k1;
        simple_jdep(p,d2,w);
        p.coordinate_advance(d2);
        //
        d2.x = d.x*(v.y - v.z);
        d2.y = d.y*(v.y - v.z);
        d2.z = d.z*(v.y - v.z);
        w.i = floor(p.x+0.5*d2.x);
        w.j = floor(p.y+0.5*d2.y);
        w.k = floor(p.z+0.5*d2.z);
        simple_jdep(p,d2,w);
        p.coordinate_advance(d2);
        //
        d2.x = d.x*(v.x - v.y);
        d2.y = d.y*(v.x - v.y);
        d2.z = d.z*(v.x - v.y);
        w.i = floor(p.x+0.5*d2.x);
        w.j = floor(p.y+0.5*d2.y);
        w.k = floor(p.z+0.5*d2.z);
        simple_jdep(p,d2,w);
        p.coordinate_advance(d2);
        //
        d2.x = d.x*(1 - v.x);
        d2.y = d.y*(1 - v.x);
        d2.z = d.z*(1 - v.x);
        w.i = floor(p.x+0.5*d2.x);
        w.j = floor(p.y+0.5*d2.y);
        w.k = floor(p.z+0.5*d2.z);
        simple_jdep(p,d2,w);
        p.coordinate_advance(d2);
    }
}

void spatial_region::rhodeposition(spatial_region::plist::particle& p)
{
    int i,j,k;
    double xa,xb,ya,yb,za,zb;
    i = floor(p.x);
    j = floor(p.y);
    k = floor(p.z);
    xa = p.x - i;
    xb = 1 - xa;
    ya = p.y - j;
    yb = 1 - ya;
    za = p.z - k;
    zb = 1 - za;
    if (is_inside(i,j,k))
    {
	if (p.cmr==-1)
	{ // electrons
	    cj[i][j][k].jx = cj[i][j][k].jx + p.q*xb*yb*zb;
	    cj[i+1][j][k].jx = cj[i+1][j][k].jx + p.q*xa*yb*zb;
	    cj[i][j+1][k].jx = cj[i][j+1][k].jx + p.q*xb*ya*zb;
	    cj[i][j][k+1].jx = cj[i][j][k+1].jx + p.q*xb*yb*za;
	    cj[i+1][j+1][k].jx = cj[i+1][j+1][k].jx + p.q*xa*ya*zb;
	    cj[i+1][j][k+1].jx = cj[i+1][j][k+1].jx + p.q*xa*yb*za;
	    cj[i][j+1][k+1].jx = cj[i][j+1][k+1].jx + p.q*xb*ya*za;
	    cj[i+1][j+1][k+1].jx = cj[i+1][j+1][k+1].jx + p.q*xa*ya*za;
	}
	else if (p.cmr==1)
	{ // positrons
	    cj[i][j][k].jy = cj[i][j][k].jy + p.q*xb*yb*zb;
	    cj[i+1][j][k].jy = cj[i+1][j][k].jy + p.q*xa*yb*zb;
	    cj[i][j+1][k].jy = cj[i][j+1][k].jy + p.q*xb*ya*zb;
	    cj[i][j][k+1].jy = cj[i][j][k+1].jy + p.q*xb*yb*za;
	    cj[i+1][j+1][k].jy = cj[i+1][j+1][k].jy + p.q*xa*ya*zb;
	    cj[i+1][j][k+1].jy = cj[i+1][j][k+1].jy + p.q*xa*yb*za;
	    cj[i][j+1][k+1].jy = cj[i][j+1][k+1].jy + p.q*xb*ya*za;
	    cj[i+1][j+1][k+1].jy = cj[i+1][j+1][k+1].jy + p.q*xa*ya*za;
	}
	else if (p.cmr==0)
	{ // photons
	    cj[i][j][k].jz = cj[i][j][k].jz + p.q*xb*yb*zb;
	    cj[i+1][j][k].jz = cj[i+1][j][k].jz + p.q*xa*yb*zb;
	    cj[i][j+1][k].jz = cj[i][j+1][k].jz + p.q*xb*ya*zb;
	    cj[i][j][k+1].jz = cj[i][j][k+1].jz + p.q*xb*yb*za;
	    cj[i+1][j+1][k].jz = cj[i+1][j+1][k].jz + p.q*xa*ya*zb;
	    cj[i+1][j][k+1].jz = cj[i+1][j][k+1].jz + p.q*xa*yb*za;
	    cj[i][j+1][k+1].jz = cj[i][j+1][k+1].jz + p.q*xb*ya*za;
	    cj[i+1][j+1][k+1].jz = cj[i+1][j+1][k+1].jz + p.q*xa*ya*za;
	}
	else
	{
	    int n;
	    n = 0;
	    while (n!=n_ion_populations)
	    {
		if (p.cmr==icmr[n])
		{
		    irho[n][i][j][k] += p.q*xb*yb*zb;
		    irho[n][i+1][j][k] += p.q*xa*yb*zb;
		    irho[n][i][j+1][k] += p.q*xb*ya*zb;
		    irho[n][i][j][k+1] += p.q*xb*yb*za;
		    irho[n][i+1][j+1][k] += p.q*xa*ya*zb;
		    irho[n][i+1][j][k+1] += p.q*xa*yb*za;
		    irho[n][i][j+1][k+1] += p.q*xb*ya*za;
		    irho[n][i+1][j+1][k+1] += p.q*xa*ya*za;
		    n = n_ion_populations;
		}
		else
		    n++;
	    }
	}
    }
}

vector3d spatial_region::e_to_particle(double& x, double& y, double& z)
{
    int i,j,k;
    vector3d e;
    double xa,xb,ya,yb,za,zb;
    i = floor(x);
    j = floor(y);
    k = floor(z);
    xa = x - i;
    xb = 1 - xa;
    ya = y - j;
    yb = 1 - ya;
    za = z - k;
    zb = 1 - za;
    e.x = ce[i][j][k].ex*yb*zb + ce[i][j+1][k].ex*ya*zb + ce[i][j][k+1].ex*yb*za + ce[i][j+1][k+1].ex*ya*za;
    e.y = ce[i][j][k].ey*xb*zb + ce[i+1][j][k].ey*xa*zb + ce[i][j][k+1].ey*xb*za + ce[i+1][j][k+1].ey*xa*za;
    e.z = ce[i][j][k].ez*xb*yb + ce[i+1][j][k].ez*xa*yb + ce[i][j+1][k].ez*xb*ya + ce[i+1][j+1][k].ez*xa*ya;
    return e;
}

vector3d spatial_region::b_to_particle(double& x, double& y, double& z)
{
    int i,j,k;
    vector3d b;
    double xa,xb,ya,yb,za,zb;
    i = floor(x);
    j = floor(y);
    k = floor(z);
    xa = x - i;
    xb = 1 - xa;
    ya = y - j;
    yb = 1 - ya;
    za = z - k;
    zb = 1 - za;

    b.x = cbe[i][j][k].bex*yb*zb + cbe[i][j+1][k].bex*ya*zb + cbe[i][j][k+1].bex*yb*za + cbe[i][j+1][k+1].bex*ya*za;
    b.y = cbe[i][j][k].bey*xb*yb + cbe[i+1][j][k].bey*xa*yb + cbe[i][j+1][k].bey*xb*ya + cbe[i+1][j+1][k].bey*xa*ya;
    b.z = cbe[i][j][k].bez*xb*zb + cbe[i+1][j][k].bez*xa*zb + cbe[i][j][k+1].bez*xb*za + cbe[i+1][j][k+1].bez*xa*za;

    return b;
}
