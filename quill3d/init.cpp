#include <cmath>
#include "main.h"

void spatial_region::f_init_cos(double a0y, double a0z, double xsigma, double ysigma, double zsigma, double x0, bool sscos, bool b_sign, double x1, double phase, double y0, double z0, bool append, double phi, double ytarget, double ztarget)
{
    /* Для f_init_gauss xsigma, ysigma и zsigma имели смысл
     * половины масштабов лазерного импульса на уровне 1/e^2 (по
     * интенсивности) */
    /* x0 - расстояние от левой границы области (spatial_region'а) до
     * центра лазерного импульса, x1 - координата по x центра импульса
     * относительно центра всей области моделирования («суммы»
     * spatial_region'ов) */
    /* Интеграл от квадратов полей не вычисляется аналитически в
     * простой форме, поэтому используется \int a^2 dV. Поля задаются
     * таким образом, чтобы \int a^2 dV имел то же значение, что и для
     * f_init_gauss, поэтому вместо xsigma, ysigma и zsigma для
     * описания огибающей используются половины полных размеров
     * импульса xs, ys и zs */
    /* Если sscos==1, то огибающая задаётся в форме "супер-супер
     * косинуса"; для того, чтобы вычисление a_0 по энергии импульса
     * (при задании в конфиг-файле W) было верным, ширина импульса
     * вычисляется так, чтобы интеграл от квадрата его огибающей
     * совпадал с соответствующим интегралом от гауссова импульса */
    double x,y,z,xi;
    double xs, ys, zs;
    if (sscos==0) {
        xs = xsigma*2*sqrt(2*PI)/3;
        ys = ysigma*2*sqrt(2*PI)/3;
        zs = zsigma*2*sqrt(2*PI)/3;
    } else {
        xs = 0.822*xsigma;
        ys = 0.822*ysigma;
        zs = 0.822*zsigma;
    }
    double cosx,cosy,cosz,sinx,siny,sinz,tr_envelope;
    double y12,z12;
    y12 = 0.5*ny*dy;
    z12 = 0.5*nz*dy;
    double r0x,r0y,r0z;
    r0x = sqrt(x1*x1+y0*y0+z0*z0);
    r0y = -y0/r0x;
    r0z = -z0/r0x;
    r0x = -x1/r0x;
    double y0x,y0y,y0z,z0x,z0y,z0z;
    y0z = sqrt( r0x*r0x + r0y*r0y );
    if (y0z!=0)
    {
        y0x = -r0y/y0z;
        y0y = r0x/y0z;
        y0z = 0;
    }
    else
    {
        y0x = 0;
        y0y = 1;
        y0z = 0;
    }
    z0x = r0y*y0z - r0z*y0y;
    z0y = r0z*y0x - r0x*y0z;
    z0z = r0x*y0y - r0y*y0x;
    if (phi!=0) {
        double c,s,yx,yy,yz,zx,zy,zz;
        c = cos(phi);
        s = sin(phi);
        yx = y0x;
        yy = y0y;
        yz = y0z;
        zx = z0x;
        zy = z0y;
        zz = z0z;
        y0x = yx*c + zx*s;
        y0y = yy*c + zy*s;
        y0z = yz*c + zz*s;
        z0x = zx*c - yx*s;
        z0y = zy*c - yy*s;
        z0z = zz*c - yz*s;
    }
    double ex,ey,ez,bx,by,bz;
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int k=0;k<nz;k++)
            {
                x = r0x*(i*dx-x0) + r0y*(j*dy-y12-y0-ytarget) + r0z*(k*dz-z12-z0-ztarget);
                xi = x + phase;
                y = y0x*(i*dx-x0) + y0y*(j*dy-y12-y0-ytarget) + y0z*(k*dz-z12-z0-ztarget);
                z = z0x*(i*dx-x0) + z0y*(j*dy-y12-y0-ytarget) + z0z*(k*dz-z12-z0-ztarget);
                // vector potential envelope = (cosx*cosy*cosz)^2;
                if (x<xs&&x>-xs&&y<ys&&y>-ys&&z<zs&&z>-zs)
                {
                    if (sscos==0) {
                        cosx = cos(PI*x/2/xs);
                        cosy = cos(PI*y/2/ys);
                        cosz = cos(PI*z/2/zs);
                        sinx = sin(PI*x/2/xs);
                        siny = sin(PI*y/2/ys);
                        sinz = sin(PI*z/2/zs);
                        tr_envelope = cosy*cosy*cosz*cosz;
                        ey = a0y*tr_envelope*(cos(xi)*cosx*cosx - PI/xs*sin(xi)*cosx*sinx);
                        ez = a0z*tr_envelope*(sin(xi)*cosx*cosx + PI/xs*cos(xi)*cosx*sinx);
                        bz = ey;
                        by = -ez;
                        ex = a0y*sin(xi)*cosx*cosx*cosz*cosz*PI/ys*cosy*siny - a0z*cos(xi)*cosx*cosx*cosy*cosy*PI/zs*cosz*sinz;
                        bx = a0z*cos(xi)*cosx*cosx*PI/ys*cosz*cosz*cosy*siny + a0y*sin(xi)*cosx*cosx*cosy*cosy*PI/zs*cosz*sinz;
                    } else {
                        cosx = cos(PI*x*x*x*x/2/(xs*xs*xs*xs));
                        cosy = cos(PI*y*y*y*y/2/(ys*ys*ys*ys));
                        cosz = cos(PI*z*z*z*z/2/(zs*zs*zs*zs));
                        sinx = sin(PI*x*x*x*x/2/(xs*xs*xs*xs));
                        siny = sin(PI*y*y*y*y/2/(ys*ys*ys*ys));
                        sinz = sin(PI*z*z*z*z/2/(zs*zs*zs*zs));
                        tr_envelope = cosy*cosy*cosz*cosz;
                        ey = a0y*tr_envelope*(cos(xi)*cosx*cosx - 4*PI*x*x*x/(xs*xs*xs*xs)*sin(xi)*cosx*sinx);
                        ez = a0z*tr_envelope*(sin(xi)*cosx*cosx + 4*PI*x*x*x/(xs*xs*xs*xs)*cos(xi)*cosx*sinx);
                        bz = ey;
                        by = -ez;
                        ex = a0y*sin(xi)*4*PI*y*y*y/(ys*ys*ys*ys)*cosx*cosx*cosz*cosz*cosy*siny - a0z*cos(xi)*4*PI*z*z*z/(zs*zs*zs*zs)*cosx*cosx*cosy*cosy*cosz*sinz;
                        bx = a0z*cos(xi)*4*PI*y*y*y/(ys*ys*ys*ys)*cosx*cosx*cosz*cosz*cosy*siny + a0y*sin(xi)*4*PI*z*z*z/(zs*zs*zs*zs)*cosx*cosx*cosy*cosy*cosz*sinz;
                    }
                    if (b_sign==0)
                    {
                        bx = -bx;
                        by = -by;
                        bz = -bz;
                    }
                    if (append==0)
                    {
                        // прежние значения перезаписываются
                        ce[i][j][k].ex = ex*r0x + ey*y0x + ez*z0x;
                        ce[i][j][k].ey = ex*r0y + ey*y0y + ez*z0y;
                        ce[i][j][k].ez = ex*r0z + ey*y0z + ez*z0z;
                        cb[i][j][k].bx = bx*r0x + by*y0x + bz*z0x;
                        cb[i][j][k].by = bx*r0y + by*y0y + bz*z0y;
                        cb[i][j][k].bz = bx*r0z + by*y0z + bz*z0z;
                    }
                    else
                    {
                        // импульс добавляется к имеющимся полям
                        ce[i][j][k].ex += ex*r0x + ey*y0x + ez*z0x;
                        ce[i][j][k].ey += ex*r0y + ey*y0y + ez*z0y;
                        ce[i][j][k].ez += ex*r0z + ey*y0z + ez*z0z;
                        cb[i][j][k].bx += bx*r0x + by*y0x + bz*z0x;
                        cb[i][j][k].by += bx*r0y + by*y0y + bz*z0y;
                        cb[i][j][k].bz += bx*r0z + by*y0z + bz*z0z;
                    }
                }
                else if (append==0)
                {
                    ce[i][j][k].ex = 0;
                    ce[i][j][k].ey = 0;
                    ce[i][j][k].ez = 0;
                    cb[i][j][k].bx = 0;
                    cb[i][j][k].by = 0;
                    cb[i][j][k].bz = 0;
                }
            }
        }
    }
}

void spatial_region::f_init_focused(double a0y, double a0z, double xsigma, double sigma0, double x0, double x1, bool b_sign, double phase, double y0, double z0, bool append, double phi, bool sscos, double K)
{
    /* sigma0 - поперечный размер в перетяжке (импульс
     * аксиально-симметричный), x0 - положение центра лазерного
     * импульса, x1 - координата по x центра лазерного импульса
     * относительно центра перетяжки */
    /* y0, z0 - смещение центра импульса от оси области,
     * обеспечивающее нужный угол падения. Поворот импульса
     * производится так, что ey остаётся лежать в плоскости xy */
    /* ey задаётся симметричным (если бы импульс был плоским), ez -
     * антисимметричным (по продольной координате) */
    /* Интеграл от квадратов полей не вычисляется аналитически в
     * простой форме, поэтому используется \int a^2 dV. Поля задаются
     * таким образом, чтобы \int a^2 dV имел то же значение, что и для
     * f_init_gauss, поэтому вместо xsigma, ysigma и zsigma для
     * описания огибающей используются половины полных размеров
     * импульса xs и s (половина поперечного размера) */
    // K is a wavevector, K = 1 for first harmonic and 2 for the second.
    double sign = (x1<=0)-(x1>0);
    x1 = sqrt(x1*x1+y0*y0+z0*z0);
    if (x1==0) x1 = dx; // иначе не определён угол поворота импульса
    double xR = K * sigma0*sigma0/2;
    double sigma = sigma0*sqrt(1+x1*x1/xR/xR); // в начальном положении
    // ratio of first and second harmonics beam sizes at initial position
    double a_multiplier = sqrt(1 + x1 * x1 / (sigma0 * sigma0 * sigma0 * sigma0 / 4)) / sqrt(1 + x1 * x1 / (xR * xR));
    double xs = xsigma*2*sqrt(2*PI)/3;
    if ( sscos == 1 )
        xs = 0.822*xsigma;
    double s = sigma/2/sqrt(((double)3)/16-1/PI/PI);
    double s0 = sigma0/2/sqrt(((double)3)/16-1/PI/PI);
    double y12,z12;
    y12 = ny*dy/2;
    z12 = nz*dz/2;
    double r0x,r0y,r0z;
    r0x = sign*sqrt(x1*x1-y0*y0-z0*z0);
    r0x = r0x/x1;
    r0y = -y0/x1;
    r0z = -z0/x1;
    double y0x,y0y,y0z,z0x,z0y,z0z;
    y0z = sqrt( r0x*r0x + r0y*r0y );
    if (y0z!=0)
    {
        y0x = -r0y/y0z;
        y0y = r0x/y0z;
        y0z = 0;
    }
    else
    {
        y0x = 0;
        y0y = 1;
        y0z = 0;
    }
    z0x = r0y*y0z - r0z*y0y;
    z0y = r0z*y0x - r0x*y0z;
    z0z = r0x*y0y - r0y*y0x;
    if (phi!=0) {
        double c,s,yx,yy,yz,zx,zy,zz;
        c = cos(phi);
        s = sin(phi);
        yx = y0x;
        yy = y0y;
        yz = y0z;
        zx = z0x;
        zy = z0y;
        zz = z0z;
        y0x = yx*c + zx*s;
        y0y = yy*c + zy*s;
        y0z = yz*c + zz*s;
        z0x = zx*c - yx*s;
        z0y = zy*c - yy*s;
        z0z = zz*c - yz*s;
    }
    double xi,x,r,cosx,cosr,sinx,sinr;
    double slocal,alocal;
    double ex,ey,ez,bx,by,bz;
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int k=0;k<nz;k++)
            {
                x = r0x*(i*dx-x0) + r0y*(j*dy-y12-y0) + r0z*(k*dz-z12-z0);
                r = sqrt( (i*dx-x0)*(i*dx-x0) + (j*dy-y12-y0)*(j*dy-y12-y0) + (k*dz-z12-z0)*(k*dz-z12-z0) - x*x );
                xi = x + r*r/2*(x-x1)/((x-x1)*(x-x1)+xR*xR) - atan((x-x1)/xR) - atan(x1/xR);
                slocal = s0*sqrt( 1 + (xi-x1+atan(x1/xR))*(xi-x1+atan(x1/xR))/(xR*xR) );
                alocal = s/slocal * a_multiplier;
                if ( xi>-xs && xi<xs && r<slocal )
                {
                    if ( sscos == 1 ) {
                        cosx = cos(PI*xi*xi*xi*xi/2/(xs*xs*xs*xs));
                        cosr = cos(PI*r/2/slocal);
                        sinx = sin(PI*xi*xi*xi*xi/2/(xs*xs*xs*xs));
                        sinr = sin(PI*r/2/slocal);
                        double xxi = K * xi + phase;
                        // envelope = (cosx*cosr)^2;
                        ey = alocal*a0y*cosr*cosr / K * (K * cos(xxi)*cosx*cosx - 4*PI*xi*xi*xi/(xs*xs*xs*xs)*sin(xxi)*cosx*sinx);
                        ez = -alocal*a0z*cosr*cosr / K * (K * sin(xxi)*cosx*cosx + 4*PI*xi*xi*xi/(xs*xs*xs*xs)*cos(xxi)*cosx*sinx);
                        bz = ey;
                        by = -ez;
                        if (r!=0)
                            ex = alocal*a0y / K * ( cosr*sinr*PI/slocal*((i*dx-x0)*y0x+(j*dy-y12-y0)*y0y+(k*dz-z12-z0)*y0z)/r*( K * sin(xxi)*cosx*cosx+4*PI*xi*xi*xi/(xs*xs*xs*xs)*cos(xxi)*cosx*sinx) + ((i*dx-x0)*y0x+(j*dy-y12-y0)*y0y+(k*dz-z12-z0)*y0z)*(x-x1)/((x-x1)*(x-x1)+xR*xR)*cosr*cosr*(cosx*cosx*(-K * cos(xxi))+4*PI*xi*xi*xi/(xs*xs*xs*xs)*cosx*sinx*sin(xxi)) ) - alocal*a0z / K * ( cosr*sinr*PI/slocal*((i*dx-x0)*z0x+(j*dy-y12-y0)*z0y+(k*dz-z12-z0)*z0z)/r*((-K * cos(xxi))*cosx*cosx-4*PI*xi*xi*xi/(xs*xs*xs*xs)*sin(xxi)*cosx*sinx) + ((i*dx-x0)*z0x+(j*dy-y12-y0)*z0y+(k*dz-z12-z0)*z0z)*(x-x1)/((x-x1)*(x-x1)+xR*xR)*cosr*cosr*(cosx*cosx*(-K * sin(xxi))-4*PI*xi*xi*xi/(xs*xs*xs*xs)*cosx*sinx*cos(xxi)) );
                        else
                            ex = 0;
                        if (r!=0)
                            bx = alocal*a0y / K * ( cosr*sinr*PI/slocal*((i*dx-x0)*z0x+(j*dy-y12-y0)*z0y+(k*dz-z12-z0)*z0z)/r*(K * sin(xxi)*cosx*cosx+4*PI*xi*xi*xi/(xs*xs*xs*xs)*cos(xxi)*cosx*sinx) + ((i*dx-x0)*z0x+(j*dy-y12-y0)*z0y+(k*dz-z12-z0)*z0z)*(x-x1)/((x-x1)*(x-x1)+xR*xR)*cosr*cosr*(cosx*cosx*(-K * cos(xxi))+4*PI*xi*xi*xi/(xs*xs*xs*xs)*cosx*sinx*sin(xxi)) ) + alocal*a0z / K * ( cosr*sinr*PI/slocal*((i*dx-x0)*y0x+(j*dy-y12-y0)*y0y+(k*dz-z12-z0)*y0z)/r*((-K * cos(xxi))*cosx*cosx-4*PI*xi*xi*xi/(xs*xs*xs*xs)*sin(xxi)*cosx*sinx) + ((i*dx-x0)*y0x+(j*dy-y12-y0)*y0y+(k*dz-z12-z0)*y0z)*(x-x1)/((x-x1)*(x-x1)+xR*xR)*cosr*cosr*(cosx*cosx*(-K * sin(xxi))-4*PI*xi*xi*xi/(xs*xs*xs*xs)*cosx*sinx*cos(xxi)) );
                        else
                            bx = 0;
                    } else if (-2 * s0 * (x1 - x) / xR <= r && r <= 2 * s0 * (x1 - x) / xR) { // otherwise unphysical wings arise
                        cosx = cos(PI*xi/2/xs);
                        cosr = cos(PI*r/2/slocal);
                        sinx = sin(PI*xi/2/xs);
                        sinr = sin(PI*r/2/slocal);
                        double xxi = K * xi + phase;
                        // envelope = (cosx*cosr)^2;
                        ey = alocal*a0y*cosr*cosr / K * (K * cos(xxi)*cosx*cosx - PI/xs*sin(xxi)*cosx*sinx);
                        ez = -alocal*a0z*cosr*cosr / K * (K * sin(xxi)*cosx*cosx + PI/xs*cos(xxi)*cosx*sinx);
                        bz = ey;
                        by = -ez;
                        if (r!=0)
                            ex = alocal*a0y / K * ( cosr*sinr*PI/slocal*((i*dx-x0)*y0x+(j*dy-y12-y0)*y0y+(k*dz-z12-z0)*y0z)/r*(K * sin(xxi)*cosx*cosx+PI/xs*cos(xxi)*cosx*sinx) + ((i*dx-x0)*y0x+(j*dy-y12-y0)*y0y+(k*dz-z12-z0)*y0z)*(x-x1)/((x-x1)*(x-x1)+xR*xR)*cosr*cosr*(cosx*cosx*(-K * cos(xxi))+PI/xs*cosx*sinx*sin(xxi)) ) - alocal*a0z / K *( cosr*sinr*PI/slocal*((i*dx-x0)*z0x+(j*dy-y12-y0)*z0y+(k*dz-z12-z0)*z0z)/r*((-K * cos(xxi))*cosx*cosx-PI/xs*sin(xxi)*cosx*sinx) + ((i*dx-x0)*z0x+(j*dy-y12-y0)*z0y+(k*dz-z12-z0)*z0z)*(x-x1)/((x-x1)*(x-x1)+xR*xR)*cosr*cosr*(cosx*cosx*(-K * sin(xxi))-PI/xs*cosx*sinx*cos(xxi)) );
                        else
                            ex = 0;
                        if (r!=0)
                            bx = alocal*a0y / K * ( cosr*sinr*PI/slocal*((i*dx-x0)*z0x+(j*dy-y12-y0)*z0y+(k*dz-z12-z0)*z0z)/r*(K * sin(xxi)*cosx*cosx+PI/xs*cos(xxi)*cosx*sinx) + ((i*dx-x0)*z0x+(j*dy-y12-y0)*z0y+(k*dz-z12-z0)*z0z)*(x-x1)/((x-x1)*(x-x1)+xR*xR)*cosr*cosr*(cosx*cosx*(-K * cos(xxi))+PI/xs*cosx*sinx*sin(xxi)) ) + alocal*a0z / K *( cosr*sinr*PI/slocal*((i*dx-x0)*y0x+(j*dy-y12-y0)*y0y+(k*dz-z12-z0)*y0z)/r*((-K * cos(xxi))*cosx*cosx-PI/xs*sin(xxi)*cosx*sinx) + ((i*dx-x0)*y0x+(j*dy-y12-y0)*y0y+(k*dz-z12-z0)*y0z)*(x-x1)/((x-x1)*(x-x1)+xR*xR)*cosr*cosr*(cosx*cosx*(-K * sin(xxi))-PI/xs*cosx*sinx*cos(xxi)) );
                        else
                            bx = 0;
                    } else {
                        ex = 0;
                        ey = 0;
                        ez = 0;
                        bx = 0;
                        by = 0;
                        bz = 0;
                    }
                    if (b_sign==0)
                    {
                        bx = -bx;
                        by = -by;
                        bz = -bz;
                    }
                    if (append==0)
                    {
                        // прежние значения перезаписываются
                        ce[i][j][k].ex = ex*r0x + ey*y0x + ez*z0x;
                        ce[i][j][k].ey = ex*r0y + ey*y0y + ez*z0y;
                        ce[i][j][k].ez = ex*r0z + ey*y0z + ez*z0z;
                        cb[i][j][k].bx = bx*r0x + by*y0x + bz*z0x;
                        cb[i][j][k].by = bx*r0y + by*y0y + bz*z0y;
                        cb[i][j][k].bz = bx*r0z + by*y0z + bz*z0z;
                    }
                    else
                    {
                        // импульс добавляется к имеющимся полям
                        ce[i][j][k].ex += ex*r0x + ey*y0x + ez*z0x;
                        ce[i][j][k].ey += ex*r0y + ey*y0y + ez*z0y;
                        ce[i][j][k].ez += ex*r0z + ey*y0z + ez*z0z;
                        cb[i][j][k].bx += bx*r0x + by*y0x + bz*z0x;
                        cb[i][j][k].by += bx*r0y + by*y0y + bz*z0y;
                        cb[i][j][k].bz += bx*r0z + by*y0z + bz*z0z;
                    }
                }
                else if (append==0)
                {
                    ce[i][j][k].ex = 0;
                    ce[i][j][k].ey = 0;
                    ce[i][j][k].ez = 0;
                    cb[i][j][k].bx = 0;
                    cb[i][j][k].by = 0;
                    cb[i][j][k].bz = 0;
                }
            }
        }
    }
}

void spatial_region::f_init_uniformB(double a0y, double a0z) {
    for(int i = 0; i < nx; ++i) {
        for(int j = 0; j < ny; ++j) {
            for(int k = 0; k < nz; ++k) {
                ce[i][j][k].ex = 0;
                ce[i][j][k].ey = 0;
                ce[i][j][k].ez = 0;
                cb[i][j][k].bx = 0;
                cb[i][j][k].by = a0y;
                cb[i][j][k].bz = a0z;
            }
        }
    }
}

void spatial_region::fill_cell_by_particles(double cmr, int_vector3d& a, int_vector3d& b, double n, double ux0, double dsplmt, double T)
{
    // a = {i,j,k} - cell position, b = {xnpic,ynpic,znpic}, n - density
    double x0;
    double y0;
    double z0;
    double q0;
    spatial_region::plist::particle* tmp_p;
    x0 = 0.5/b.i + dsplmt;
    y0 = 0.5/b.j;
    z0 = 0.5/b.k;
    q0 = 1/double(b.i*b.j*b.k)*n;
    if (cmr<0) q0 = -q0;
    for(int ip=0;ip<b.i;ip++)
    {
        for(int jp=0;jp<b.j;jp++)
        {
            for(int kp=0;kp<b.k;kp++)
            {
                tmp_p = new_particle();
                // initial conditions
                tmp_p->cmr = cmr;
                tmp_p->q = q0;
                tmp_p->x = x0 + double(a.i) + double(ip)/double(b.i);
                tmp_p->y = y0 + double(a.j) + double(jp)/double(b.j);
                tmp_p->z = z0 + double(a.k) + double(kp)/double(b.k);
                if ( T!=0 && cmr!=0 ){
                    double v0 = ux0/sqrt( 1 + ux0*ux0 );
                    double a, b, c, r, gamma, v;
                    // gives gamma with approximately e^(-(x-1)/T) distribution
                    do {
                        gamma = 1 + get_rand() * 5 * T;
                    } while ( get_rand() > exp( -(gamma - 1) / T));
                    // gives random direction (uniform)
                    do {
                        a = 2*get_rand()-1;
                        b = 2*get_rand()-1;
                        c = 2*get_rand()-1;
                        r = sqrt( a*a + b*b + c*c );
                    } while ( r > 1 || r == 0);
                    a = a / r;
                    b = b / r;
                    c = c / r;
                    // (a, b, c) now is the random vector uniformly distributed on a unit
                    // sphere
                    v = sqrt(1 - 1 / (gamma * gamma));
                    a = v * a; // vx'
                    b = v * b; // vy'
                    c = v * c; // vz'
                    // Lorentz transformation to lab reference frame
                    b = b*sqrt( 1 - v0*v0 )/( 1 + a*v0 );
                    c = c*sqrt( 1 - v0*v0 )/( 1 + a*v0 );
                    a = ( a + v0 )/( 1 + a*v0 );
                    double g = 1/sqrt( 1 - a*a - b*b - c*c );
                    tmp_p->ux = a*g;
                    tmp_p->uy = b*g;
                    tmp_p->uz = c*g;
                } else {
                    tmp_p->ux = ux0;
                    tmp_p->uy = 0;
                    tmp_p->uz = 0;
                }
                if (cmr!=0)
                    tmp_p->g = sqrt(1 + (*tmp_p).ux*(*tmp_p).ux + (*tmp_p).uy*(*tmp_p).uy + (*tmp_p).uz*(*tmp_p).uz);
                else
                    tmp_p->g = sqrt((*tmp_p).ux*(*tmp_p).ux + (*tmp_p).uy*(*tmp_p).uy + (*tmp_p).uz*(*tmp_p).uz);
                tmp_p->chi = 0;
                tmp_p->trn = 0;
                place(*tmp_p);
            }
        }
    }
    cp[a.i][a.j][a.k].pl.start = cp[a.i][a.j][a.k].pl.head;
}

void spatial_region::add_beam(double cmr, double n0, double ux0, double xb, double rb, double x0b)
{
    /* Добавляет электронный пучок с полями, вычисленными в
     * приближении бесконечного гамма-фактора пучка. n0 - максимальная
     * концентрация электронов в пучке, нормированная на критическую
     * концентрацию */

    int_vector3d a,b;
    b.i = xnpic;
    b.j = ynpic;
    b.k = znpic;
    double x;
    double r;
    double n;
    double tmp;
    double signum;
    double type;
    signum = (ux0>0) - (ux0<0);
    type = (cmr>0)-(cmr<0);
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int k=0;k<nz;k++)
            {
                a.i = i;
                a.j = j;
                a.k = k;
                x = i*dx-x0b;
                r = sqrt((j-ny/2)*(j-ny/2)*dy*dy+(k-nz/2)*(k-nz/2)*dz*dz);
                if (x>-xb&&x<xb&&r<rb)
                {
                    n = n0*(1-x*x/(xb*xb))*(1-r*r/(rb*rb));
                    fill_cell_by_particles(cmr,a,b,n,ux0);
                    tmp = 0.5*n0*(1-x*x/(xb*xb))*(1-r*r/(2*rb*rb));
                    ce[i][j][k].ey += type*tmp*(j-ny/2)*dy;
                    ce[i][j][k].ez += type*tmp*(k-nz/2)*dz;
                    cb[i][j][k].by -= type*signum*tmp*(k-nz/2)*dz;
                    cb[i][j][k].bz += type*signum*tmp*(j-ny/2)*dy;
                }
                if (x>-xb&&x<xb&&r>rb)
                {
                    tmp = 0.25*n0*(1-x*x/(xb*xb))*rb*rb/(r*r);
                    ce[i][j][k].ey += type*tmp*(j-ny/2)*dy;
                    ce[i][j][k].ez += type*tmp*(k-nz/2)*dz;
                    cb[i][j][k].by -= type*signum*tmp*(k-nz/2)*dz;
                    cb[i][j][k].bz += type*signum*tmp*(j-ny/2)*dy;
                }
            }
        }
    }
}

void spatial_region::film(double x0, double x1, double ne, bool ions, double cmr, double gradwidth, double y0, double y1, double z0, double z1, double T, double vx, bool is_profiled)
{ /* x0 - координата левой границы плёнки, x1 - правой, ne -
     концентрация электронов в плёнке, нормированная на критическую
     концентрацию */
    /* gradwidth - толщина части плёнки с линейным ростом плотности от
     * 0 на левой границе плёнки до ne при x0+gradwidth */
    /* при gradwidth<0 плотность плёнки линейно спадает от ne до 0 на
     * участке с градиентом */
    // if is_profiled == 1 then film has transverse envelope

    int_vector3d a,b;
    b.i = xnpic;
    b.j = ynpic;
    b.k = znpic;
    int i0,i1;
    i0 = x0/dx;
    i1 = x1/dx;
    if (i0<0) i0 = 0;
    if (i1<0) i1 = 0;
    if (i0>nx) i0 = nx;
    if (i1>nx) i1 = nx;
    for(int i=i0;i<i1;i++)
    {
        for(int j=int(y0/dy)+1;j<int(y1/dy)-1;j++)
        {
            for(int k=int(z0/dz)+1;k<int(z1/dz)-1;k++)
            {
                double nes;
                if (is_profiled == 1) {
                    double tmp = (j * dy - ny * dy / 2) / (ny * dy / 2);
                    tmp = cos(0.5 * PI * tmp * tmp);
                    double tr_env = tmp * tmp;
                    tmp =  (k * dz - nz * dz / 2) / (nz * dz / 2);
                    tmp = cos(0.5 * PI * tmp * tmp);
                    tr_env *= tmp * tmp;
                    nes = ne * tr_env;
                } else {
                    nes = ne;
                }
                //
                a.i = i;
                a.j = j;
                a.k = k;
                double ux0 = vx / sqrt(1 - vx * vx);
                if (gradwidth>=0) {
                    if (i>=int((x0+gradwidth)/dx))
                    {
                        fill_cell_by_particles(-1, a, b, nes, ux0, 0, T);
                        if (ions)
                            fill_cell_by_particles(cmr, a, b, nes, ux0, 0, T * cmr);
                    }
                    else
                    {
                        fill_cell_by_particles(-1, a, b, nes*(i*dx-x0)/gradwidth, ux0, 0, T);
                        if (ions)
                            fill_cell_by_particles(cmr, a, b, nes*(i*dx-x0)/gradwidth, ux0, 0, T * cmr);
                    }
                } else {
                    if (i>=int((x0-gradwidth)/dx))
                    {
                        fill_cell_by_particles(-1, a, b, nes, ux0, 0, T);
                        if (ions)
                            fill_cell_by_particles(cmr, a, b, nes, ux0, 0, T * cmr);
                    }
                    else
                    {
                        fill_cell_by_particles(-1, a, b, nes * (1 + (i * dx - x0) / gradwidth), ux0, 0, T);
                        if (ions)
                            fill_cell_by_particles(cmr, a, b, nes * (1 + (i * dx - x0) / gradwidth), ux0, 0, T * cmr);
                    }
                }
            }
        }
    }
}
