/* Распределение частиц с энергией от mineps до maxeps по углам theta
 * и phi, определяющим направление скорости частицы. Угол theta
 * отсчитывается от плоскости yz (theta in [-pi/2,pi/2]), угол phi - от
 * оси y в плоскости yz (phi in [-pi,pi]). */

/* Если maxeps задана меньше или равной mineps, то maxeps полагается
  * равной бесконечности */

/* L определяет расстояние от центра области до плоского детектора,
 * измеряющего N(t) */

real file_number = 12;
real mineps = 0; // MeV
real maxeps = 0; // MeV
real max_fe = 0; // если max_f==0, то используется максимальное значение f
real max_fp = 0; // если max_f==0, то используется максимальное значение f
real max_fg = 0; // если max_f==0, то используется максимальное значение f
int ntheta = 40;
int nphi = 80;
// проекции направления для построения N(t); если хотя бы одно из них задано вне разрешённой области значений, то используется направление, в котором значение dN/d\Omega максимально
real Omega_theta = 0.05*pi;
real Omega_phi = 1*pi;
real rdOmega = 0.1*pi;
real L = 1; // cm
real dts = 0.1; // lambda; шаг по времени для построения N(t)
real level = 0.1; // N(t) строится только в области, где N(t)>level*max(N(t))
real max_epsdndeps = 0e12; // если задана равной нулю, то вычисляется из спектров
real min_epsdndeps = 0; // если задана равной нулю, то полагается равной 1e-5*max_dndeps
string scintillating_screen = "on"; /* если задано значение «on», то
вместо электронного и позитронного спектров рисуется распределение
фотонов на плоском экране, установленом перпендикулярно направлению
Omega (которое проецируется точно в центр экрана); ось z' (нижняя
сторона экрана) параллельна плоскости x,z и при phi=0 сонаправлена с
осью z, а ось y' (левая сторона экрана) при phi=0, theta=pi/2
сонаправлена с осью y */
real l_scsc = 1; /* длина стороны экрана, в расстояниях от центра
области до центра экрана */

real dOmega = rdOmega*rdOmega;
int n_scsc = 50;
real dl_scsc = l_scsc/n_scsc;

//----------------------------------------

real dt,dx,dy,dz,lambda;
int nx,ny,nz;
real deps;
int neps;
int enthp;
real deps_p;
int neps_p;
int enthp_p;
real deps_ph;
int neps_ph;
int enthp_ph;

string results_folder = "../results/";
file fin_param = input( results_folder+"log", comment="" );
string var_name;
while (var_name!="$")
{
    var_name = fin_param;
    if (var_name=="dx") dx = fin_param;
    if (var_name=="dy") dy = fin_param;
    if (var_name=="dz") dz = fin_param;
    if (var_name=="lambda") lambda = fin_param;
    if (var_name=="nx") nx = fin_param;
    if (var_name=="ny") ny = fin_param;
    if (var_name=="nz") nz = fin_param;
    if (var_name=="dt") dt = fin_param;
    if (var_name=="deps") deps = fin_param;
    if (var_name=="neps") neps = fin_param;
    if (var_name=="enthp") enthp = fin_param;
    if (var_name=="deps_p") deps_p = fin_param;
    if (var_name=="neps_p") neps_p = fin_param;
    if (var_name=="enthp_p") enthp_p = fin_param;
    if (var_name=="deps_ph") deps_ph = fin_param;
    if (var_name=="neps_ph") neps_ph = fin_param;
    if (var_name=="enthp_ph") enthp_ph = fin_param;
}

real epsmax = max(neps*deps,neps_p*deps_p,neps_ph*deps_ph);
real xlength,ylength,zlength;
xlength = dx*nx;
ylength = dy*ny;
zlength = dz*nz;

//----------------------------------------

real dtheta = pi/ntheta;
real dphi = 2*pi/nphi;

real[][] fe = new real[ntheta][nphi];
real[][] fp = new real[ntheta][nphi];
real[][] fg = new real[ntheta][nphi];

int i,j;
for( i=0; i<ntheta; i=i+1 )
{
    for( j=0; j<nphi; j=j+1 )
    {
	fe[i][j] = 0;
	fp[i][j] = 0;
	fg[i][j] = 0;
    }
}

//----------------------------------------

real[][] fg_scsc;
real x_Omega,y_Omega,z_Omega,beta,eta,sinbeta,cosbeta,sineta,coseta;
int ys,zs;
real max_fg_scsc=0;
if (scintillating_screen=="on")
{
    fg_scsc = new real[n_scsc][n_scsc];
    for (i=0;i<n_scsc;i=i+1)
    {
	for (j=0;j<n_scsc;j=j+1)
	{
	    fg_scsc[i][j] = 0;
	}
    }
}

//----------------------------------------

file fin=input(results_folder+"phasespace"+format("%g",file_number));
real[] data=fin;

int nparticles = Floor(data.length/8);
real q,ux,uy,uz,utr,g;
int n;
for (n=0;n<nparticles;n=n+1)
{
    q = data[8*n];
    ux = data[8*n+4];
    uy = data[8*n+5];
    uz = data[8*n+6];
    g = data[8*n+7];
    utr = sqrt(uy*uy+uz*uz);
    if ((ux!=0||utr!=0)&&(uy!=0||uz!=0))
    {
	if (g*0.511>mineps&&(maxeps<=mineps||g*0.511<maxeps))
	{
	    i = Floor((atan2(ux,utr)+pi/2)/dtheta);
	    j = Floor((atan2(uz,uy)+pi)/dphi);
	    if (j==nphi) j = 0;
	    fe[i][j] -= q;
	}
    }
}

//----------------------------------------

fin=input(results_folder+"phasespace_p"+format("%g",file_number));
data=fin;

nparticles = Floor(data.length/8);
for (n=0;n<nparticles;n=n+1)
{
    q = data[8*n];
    ux = data[8*n+4];
    uy = data[8*n+5];
    uz = data[8*n+6];
    g = data[8*n+7];
    utr = sqrt(uy*uy+uz*uz);
    if ((ux!=0||utr!=0)&&(uy!=0||uz!=0))
    {
	if (g*0.511>mineps&&(maxeps<=mineps||g*0.511<maxeps))
	{
	    i = Floor((atan2(ux,utr)+pi/2)/dtheta);
	    j = Floor((atan2(uz,uy)+pi)/dphi);
	    if (j==nphi) j = 0;
	    fp[i][j] += q;
	}
    }
}

//----------------------------------------

fin=input(results_folder+"phasespace_ph"+format("%g",file_number));
data=fin;

nparticles = Floor(data.length/8);
for (n=0;n<nparticles;n=n+1)
{
    q = data[8*n];
    ux = data[8*n+4];
    uy = data[8*n+5];
    uz = data[8*n+6];
    g = data[8*n+7];
    utr = sqrt(uy*uy+uz*uz);
    if ((ux!=0||utr!=0)&&(uy!=0||uz!=0))
    {
	if (g*0.511>mineps&&(maxeps<=mineps||g*0.511<maxeps))
	{
	    i = Floor((atan2(ux,utr)+pi/2)/dtheta);
	    j = Floor((atan2(uz,uy)+pi)/dphi);
	    if (j==nphi) j = 0;
	    fg[i][j] += q;
	}
    }
}

//----------------------------------------

real Ne = 0;
real Np = 0;
real Ng = 0;
for( i=0; i<ntheta; i=i+1 )
{
    for( j=0; j<nphi; j=j+1 )
    {
	fe[i][j] = fe[i][j]*dx*dy*dz*1.11485e13*lambda*enthp/(dtheta*dphi);
	fp[i][j] = fp[i][j]*dx*dy*dz*1.11485e13*lambda*enthp_p/(dtheta*dphi);
	fg[i][j] = fg[i][j]*dx*dy*dz*1.11485e13*lambda*enthp_ph/(dtheta*dphi);
	Ne += fe[i][j]*dtheta*dphi;
	Np += fp[i][j]*dtheta*dphi;
	Ng += fg[i][j]*dtheta*dphi;
	fe[i][j] = fe[i][j]/cos((i+0.5)*dtheta-pi/2);
	fp[i][j] = fp[i][j]/cos((i+0.5)*dtheta-pi/2);
	fg[i][j] = fg[i][j]/cos((i+0.5)*dtheta-pi/2);
    }
}
write("N_e = ",format("%g",Ne));
write("N_p = ",format("%g",Np));
write("N_g = ",format("%g",Ng));

//----------------------------------------

fin=input(results_folder+"spectrum"+format("%g",file_number));
data=fin;

real spmax = 0;
real[] epsdndeps = new real[neps];
real[] eps = new real[neps];

for( i=0;i<neps;i=i+1 )
{
    eps[i] = (i+0.5)*deps;
    epsdndeps[i] = eps[i]*data[i];
    if (epsdndeps[i]>spmax) spmax = epsdndeps[i];
}

for( i=0;i<neps;i=i+1 )
{
    if (epsdndeps[i]!=0)
	epsdndeps[i] = log(epsdndeps[i])/log(10);
    else
	epsdndeps[i] = -100;
}

//----------------------------------------

fin=input(results_folder+"spectrum_p"+format("%g",file_number));
data=fin;

real[] epsdndeps_p = new real[neps_p];
real[] eps_p = new real[neps_p];

for( i=0;i<neps_p;i=i+1 )
{
    eps_p[i] = (i+0.5)*deps_p;
    epsdndeps_p[i] = eps_p[i]*data[i];
    if (epsdndeps_p[i]>spmax) spmax = epsdndeps_p[i];
}

for( i=0;i<neps_p;i=i+1 )
{
    if (epsdndeps_p[i]!=0)
	epsdndeps_p[i] = log(epsdndeps_p[i])/log(10);
    else
	epsdndeps_p[i] = -100;
}

//----------------------------------------

fin=input(results_folder+"spectrum_ph"+format("%g",file_number));
data=fin;

real[] epsdndeps_ph = new real[neps_ph];
real[] eps_ph = new real[neps_ph];

for( i=0;i<neps_ph;i=i+1 )
{
    eps_ph[i] = (i+0.5)*deps_ph;
    epsdndeps_ph[i] = eps_ph[i]*data[i];
    if (epsdndeps_ph[i]>spmax) spmax = epsdndeps_ph[i];
}

for( i=0;i<neps_ph;i=i+1 )
{
    if (epsdndeps_ph[i]!=0)
	epsdndeps_ph[i] = log(epsdndeps_ph[i])/log(10);
    else
	epsdndeps_ph[i] = -100;
}

//------------------------------

if (max_fe==0) max_fe = max(fe);
if (max_fp==0) max_fp = max(fp);
if (max_fg==0) max_fg = max(fg);

if (max_epsdndeps==0) max_epsdndeps = spmax;
if (min_epsdndeps==0) min_epsdndeps = 1e-5*max_epsdndeps;

//------------------------------

real tmp;
real dNdOmega_max = 0;

if (Omega_theta>pi/2||Omega_theta<-pi/2||Omega_phi>pi||Omega_phi<-pi)
{
    for( i=1; i<ntheta-1; i=i+1 )
    {
	for( j=1; j<nphi-1; j=j+1 )
	{
	    tmp = (fg[i][j] + fg[i+1][j] + fg[i+1][j+1] + fg[i][j+1] + fg[i-1][j+1] + fg[i-1][j] + fg[i-1][j-1] + fg[i][j-1] + fg[i+1][j-1])/9;
	    if (tmp>dNdOmega_max)
	    {
		dNdOmega_max = tmp;
		Omega_theta = (i+0.5)*dtheta;
		Omega_phi = (j+0.5)*dphi;
	    }
	}
    }

    for( i=1; i<ntheta-1; i=i+1 )
    {
	j=0;
	tmp = (fg[i][j] + fg[i+1][j] + fg[i+1][j+1] + fg[i][j+1] + fg[i-1][j+1] + fg[i-1][j] + fg[i-1][nphi-1] + fg[i][nphi-1] + fg[i+1][nphi-1])/9;
	if (tmp>dNdOmega_max)
	{
	    dNdOmega_max = tmp;
	    Omega_theta = (i+0.5)*dtheta;
	    Omega_phi = (j+0.5)*dphi;
	}
	j=nphi-1;
	tmp = (fg[i][j] + fg[i+1][j] + fg[i+1][0] + fg[i][0] + fg[i-1][0] + fg[i-1][j] + fg[i-1][j-1] + fg[i][j-1] + fg[i+1][j-1])/9;
	if (tmp>dNdOmega_max)
	{
	    dNdOmega_max = tmp;
	    Omega_theta = (i+0.5)*dtheta;
	    Omega_phi = (j+0.5)*dphi;
	}
    }

    i = 0;
    for( j=2; j<nphi-2; j=j+1 )
    {
	tmp = (fg[i][j] + fg[i][j+1] + fg[i][j+2] + fg[i][j-1] + fg[i][j-2])/5;
	if (tmp>dNdOmega_max)
	{
	    dNdOmega_max = tmp;
	    Omega_theta = (i+0.5)*dtheta;
	    Omega_phi = (j+0.5)*dphi;
	}
    }
    i = ntheta-1;
    for( j=2; j<nphi-2; j=j+1 )
    {
	tmp = (fg[i][j] + fg[i][j+1] + fg[i][j+2] + fg[i][j-1] + fg[i][j-2])/5;
	if (tmp>dNdOmega_max)
	{
	    dNdOmega_max = tmp;
	    Omega_theta = (i+0.5)*dtheta;
	    Omega_phi = (j+0.5)*dphi;
	}
    }

    tmp = (fg[0][0] + fg[1][0] + fg[1][1] + fg[1][2] + fg[0][2] + fg[0][1])/6;
    if (tmp>dNdOmega_max)
    {
	dNdOmega_max = tmp;
	Omega_theta = 0.5*dtheta;
	Omega_phi = 0.5*dphi;
    }
    tmp = (fg[ntheta-1][nphi-1] + fg[ntheta-2][nphi-1] + fg[ntheta-2][nphi-2] + fg[ntheta-2][nphi-3] + fg[ntheta-1][nphi-3] + fg[ntheta-1][nphi-2])/6;
    if (tmp>dNdOmega_max)
    {
	dNdOmega_max = tmp;
	Omega_theta = 0.5*dtheta;
	Omega_phi = 0.5*dphi;
    }

    Omega_theta -= pi/2;
    Omega_phi -= pi;
}
write("dN/dOmega_max = ",format("%g",dNdOmega_max));

//----------------------------------------

if (scintillating_screen=="on")
{
    x_Omega = sin(Omega_theta);
    y_Omega = cos(Omega_theta)*cos(Omega_phi);
    z_Omega = cos(Omega_theta)*sin(Omega_phi);
    beta = atan2(z_Omega,x_Omega);
    sinbeta = sin(beta);
    cosbeta = cos(beta);
    eta = atan2(y_Omega*cosbeta,x_Omega);
    sineta = sin(eta);
    coseta = cos(eta);
}

//----------------------------------------

fin=input(results_folder+"phasespace_ph"+format("%g",file_number));
data=fin;
nparticles = Floor(data.length/8);

real theta,phi,N_in_dOmega,r2,psi,sigma4pi,sigma2pi,sigmapi;
N_in_dOmega = 0;
sigma4pi = 0;
sigma2pi = 0;
sigmapi = 0;

real diagonal = sqrt(xlength*xlength+ylength*ylength+zlength*zlength);
real x,y,z,vOmega,rOmega;
int nt = Floor(diagonal/dt);

real[] dNdOmegadt = new real[nt];

for (n=0;n<nt;n=n+1)
{
    dNdOmegadt[n] = 0;
}

for (n=0;n<nparticles;n=n+1)
{
    q = data[8*n];
    x = data[8*n+1]-xlength/2;
    y = data[8*n+2]-ylength/2;
    z = data[8*n+3]-zlength/2;
    ux = data[8*n+4];
    uy = data[8*n+5];
    uz = data[8*n+6];
    g = data[8*n+7];
    utr = sqrt(uy*uy+uz*uz);
    if ((ux!=0||utr!=0)&&(uy!=0||uz!=0))
    {
	if (g*0.511>mineps&&(maxeps<=mineps||g*0.511<maxeps))
	{
	    theta = atan2(ux,utr);
	    phi = atan2(uz,uy);
	    r2 = (sin(theta)-sin(Omega_theta))^2+(cos(theta)*cos(phi)-cos(Omega_theta)*cos(Omega_phi))^2+(cos(theta)*sin(phi)-cos(Omega_theta)*sin(Omega_phi))^2;
	    if (r2!=0)
		psi = 2*asin(sqrt(r2)/2);
	    else
		psi = 0;
	    sigma4pi += q*psi;
	    if (psi<pi/2)
		sigma2pi += q*psi;
	    if (psi<pi/3)
		sigmapi += q*psi;
	    if (r2<dOmega)
	    {
		N_in_dOmega += q;
		vOmega = ux/g*sin(theta) + uy/g*cos(theta)*cos(phi) + uz/g*cos(theta)*sin(phi);
		rOmega = x*sin(theta) + y*cos(theta)*cos(phi) + z*cos(theta)*sin(phi);
		i = Floor( nt/2 + ((L/lambda-rOmega)/vOmega-L/lambda)/dt );
		dNdOmegadt[i] += q*dx*dy*dz*1.11485e13*lambda*enthp/dt/(pi*rdOmega*rdOmega);
	    }
	    if (scintillating_screen=="on")
	    {
		ys = Floor( (l_scsc/2-sineta/coseta+1/((sineta+coseta*(tan(phi)*sinbeta+tan(theta)*cosbeta/cos(phi)))*coseta))/dl_scsc );
		zs = Floor( (l_scsc/2+(cosbeta*tan(phi)-sinbeta*tan(theta)/cos(phi))/(sineta+coseta*(cosbeta*tan(theta)/cos(phi)+tan(phi)*sinbeta)))/dl_scsc );
		if (psi<pi/2&&ys>-1&&ys<n_scsc&&zs>-1&&zs<n_scsc)
		    fg_scsc[zs][ys] += q;
	    }
	}
    }
}

N_in_dOmega = N_in_dOmega*dx*dy*dz*1.11485e13*lambda*enthp_ph;
sigma4pi = sigma4pi*dx*dy*dz*1.11485e13*lambda*enthp_ph/Ng;
sigma2pi = sigma2pi*dx*dy*dz*1.11485e13*lambda*enthp_ph/Ng;
sigmapi = sigmapi*dx*dy*dz*1.11485e13*lambda*enthp_ph/Ng;

write("N in dOmega = / pi rdOmega^2 = ", format("%g",N_in_dOmega/(pi*rdOmega*rdOmega)));
write("Для равномерного распределения sigma4pi/2pi/pi = 90°/57°/39°");
write("sigma4pi = ±"+format("%g",sigma4pi/pi*180)+"°");
write("sigma2pi = ±"+format("%g",sigma2pi/pi*180)+"°");
write("sigmapi = ±"+format("%g",sigmapi/pi*180)+"°");

if (scintillating_screen=="on")
{
    for (i=0;i<n_scsc;i=i+1)
    {
	for (j=0;j<n_scsc;j=j+1)
	{
	    fg_scsc[i][j] = fg_scsc[i][j]*dx*dy*dz*1.11485e13*lambda*enthp/(dl_scsc*dl_scsc);
	}
    }
    max_fg_scsc = max(fg_scsc);
}

//------------------------------

real t1,t2;
t1 = 0;
tmp = max(dNdOmegadt);

for (n=0;n<nt;n=n+1)
{
    if (dNdOmegadt[n]>level*tmp)
    {
	t2 = n*dt;
	if (t1==0) t1 = t2;
    }
}

// для удобства представления результатов область построения расширена на 2*lambda
int nts = Floor((t2-t1+2)/dts);
real[] ts = new real[nts];
real[] dNdOmegadts = new real[nts];
for (n=0;n<nts;n=n+1)
{
    ts[n] = n*dts;
    dNdOmegadts[n] = 0;
}
for (n=0;n<nt;n=n+1)
{
    i = Floor((n*dt-(t1-1))/dts);
    if (i>=0&&i<nts)
	dNdOmegadts[i] += dNdOmegadt[n]*dt/dts;
}

//------------------------------

import graph;
import palette;
defaultpen(linewidth(0.7)+fontsize(10));
pen Tickpen=black;
pen tickpen=gray+0.5*linewidth(currentpen);

if (max_fe==0) max_fe = 1;
if (max_fp==0) max_fp = 1;
if (max_fg==0) max_fg = 1;

range fe_range = Range(0,max_fe);
range fp_range = Range(0,max_fp);
range fg_range = Range(0,max_fg);
range fg_scsc_range = Range(0,max_fg_scsc);

//------------------------------

int npalette = 200;

pen[] palettee = new pen[npalette];

string[] plble = {"0",format("%g",max_fe/2),format("%g",max_fe)};

// pn1 - pne - pn3
int n_intermediate = Floor(npalette/2);
pen pn1 = white;
pen pne = 0.5*red+1.0*blue+0.6*green;
pen pn3 = black;
for(i=0;i<n_intermediate;i=i+1)
{
    palettee[i] = (1-i/(n_intermediate-1))*pn1 + i/(n_intermediate-1)*pne;
}
for(i=n_intermediate;i<npalette;i=i+1)
{
    palettee[i] = (1-(i-n_intermediate)/(npalette-1-n_intermediate))*pne + (i-n_intermediate)/(npalette-1-n_intermediate)*pn3;
}

//------------------------------

pen[] palettep = new pen[npalette];

string[] plblp = {"0",format("%g",max_fp/2),format("%g",max_fp)};

// pn1 - pnp - pn3
pen pnp = 1.0*red+0.6*green+0.0*blue;
for(i=0;i<n_intermediate;i=i+1)
{
    palettep[i] = (1-i/(n_intermediate-1))*pn1 + i/(n_intermediate-1)*pnp;
}
for(i=n_intermediate;i<npalette;i=i+1)
{
    palettep[i] = (1-(i-n_intermediate)/(npalette-1-n_intermediate))*pnp + (i-n_intermediate)/(npalette-1-n_intermediate)*pn3;
}

//------------------------------

pen[] paletteg = new pen[npalette];

string[] plblg = {"0",format("%g",max_fg/2),format("%g",max_fg)};

// pn1 - png - pn3
pen png = 0.4*red+0.8*green+0.0*blue;
for(i=0;i<n_intermediate;i=i+1)
{
    paletteg[i] = (1-i/(n_intermediate-1))*pn1 + i/(n_intermediate-1)*png;
}
for(i=n_intermediate;i<npalette;i=i+1)
{
    paletteg[i] = (1-(i-n_intermediate)/(npalette-1-n_intermediate))*png + (i-n_intermediate)/(npalette-1-n_intermediate)*pn3;
}

//------------------------------

real picturewidth = 4;
real xindent = 1.3; // relative, in picturewidth
real yindent = 2.85; // relative
string[] xtcks = {"$-\pi/2$","0","$\pi/2$"};
string[] ytcks = {"$-\pi$","$-\pi/2$","0","$\pi/2$","$\pi$"};

bounds prange;
if (scintillating_screen=="on")
{
    picture pic1;
    prange=image(pic1,fg_scsc,fg_scsc_range,(0,0),(l_scsc,l_scsc),paletteg);
    palette(pic1,"$d^2N/dz'dy'$",prange,point(pic1,NW)+0.04N*l_scsc,point(pic1,NE)+0.08N*l_scsc,Top,paletteg,PaletteTicks(new string(real x){return format("%g",x);},max_fg_scsc/2,max_fg_scsc/4,Tickpen,tickpen));
    xlimits(pic1,0,l_scsc);
    ylimits(pic1,0,l_scsc);
    xaxis(pic1,"$z'$",BottomTop,LeftTicks(l_scsc/5,l_scsc/10),true);
    yaxis(pic1,"$y'$",LeftRight,RightTicks(l_scsc/5,l_scsc/10),true);
    size(pic1,2*picturewidth*1cm,2*picturewidth*1cm,point(pic1,SW),point(pic1,NE));
    frame frm1=pic1.fit();
    frm1=shift(-0.5*picturewidth*1cm,-picturewidth*1cm)*frm1;
    add(frm1);
}
else
{
    picture pic1;
    prange=image(pic1,fe,fe_range,(-pi/2,-pi),(pi/2,pi),palettee);
    palette(pic1,"$dN_e/d\Omega$",prange,point(pic1,NW)+0.04N*2*pi,point(pic1,NE)+0.08N*2*pi,Top,palettee,PaletteTicks(new string(real x){return plble[Floor(2*x/max_fe)];},max_fe/2,max_fe/4,Tickpen,tickpen));
    xlimits(pic1,-pi/2,pi/2);
    ylimits(pic1,-pi,pi);
    xaxis(pic1,"$\theta$",BottomTop,LeftTicks(new string(real x){return xtcks[Floor(1+2*x/pi)];},pi/2,pi/4),true);
    yaxis(pic1,"$\varphi$",LeftRight,RightTicks(new string(real x){return ytcks[Floor(2+2*x/pi)];},pi/2,pi/4),true);
    size(pic1,picturewidth*1cm,picturewidth*1cm*2,point(pic1,SW),point(pic1,NE));
    frame frm1=pic1.fit();
    add(frm1);

    picture pic2;
    prange=image(pic2,fp,fp_range,(-pi/2,-pi),(pi/2,pi),palettep);
    palette(pic2,"$dN_p/d\Omega$",prange,point(pic2,NW)+0.04N*2*pi,point(pic2,NE)+0.08N*2*pi,Top,palettep,PaletteTicks(new string(real x){return plblp[Floor(2*x/max_fp)];},max_fp/2,max_fp/4,Tickpen,tickpen));
    xlimits(pic2,-pi/2,pi/2);
    ylimits(pic2,-pi,pi);
    xaxis(pic2,"$\theta$",BottomTop,LeftTicks(new string(real x){return xtcks[Floor(1+2*x/pi)];},pi/2,pi/4),true);
    yaxis(pic2,"",LeftRight,RightTicks(new string(real x){return "";},pi/2,pi/4),true);
    size(pic2,picturewidth*1cm,picturewidth*1cm*2,point(pic2,SW),point(pic2,NE));
    frame frm2=pic2.fit();
    frm2=shift(xindent*picturewidth*1cm,0)*frm2;
    add(frm2);
}

picture pic3;
prange=image(pic3,fg,fg_range,(-pi/2,-pi),(pi/2,pi),paletteg);
palette(pic3,"$dN_{ph}/d\Omega$",prange,point(pic3,NW)+0.04N*2*pi,point(pic3,NE)+0.08N*2*pi,Top,paletteg,PaletteTicks(new string(real x){return plblg[Floor(2*x/max_fg)];},max_fg/2,max_fg/4,Tickpen,tickpen));
xlimits(pic3,-pi/2,pi/2);
ylimits(pic3,-pi,pi);
xaxis(pic3,"$\theta$",BottomTop,LeftTicks(new string(real x){return xtcks[Floor(1+2*x/pi)];},pi/2,pi/4),true);
if (scintillating_screen=="on")
    yaxis(pic3,"$\varphi$",LeftRight,RightTicks(new string(real x){return ytcks[Floor(2+2*x/pi)];},pi/2,pi/4),true);
else
    yaxis(pic3,"",LeftRight,RightTicks(new string(real x){return "";},pi/2,pi/4),true);
label(pic3,"$\;\,\mathbf \Omega$",(Omega_theta,Omega_phi),Align,red);
draw(pic3,(Omega_theta,Omega_phi),red+linewidth(4.0)); // dot
size(pic3,picturewidth*1cm,picturewidth*1cm*2,point(pic3,SW),point(pic3,NE));
frame frm3=pic3.fit();
frm3=shift(2*xindent*picturewidth*1cm,0)*frm3;
add(frm3);

picture pic4;
pen p = linejoin(2)+linewidth(1.0);
real pencorrection = 0.7;
scale(pic4,Linear,Log);
draw(pic4,graph(eps,epsdndeps,Straight),p+pencorrection*pne);
draw(pic4,graph(eps_ph,epsdndeps_ph,Straight),p+pencorrection*png);
draw(pic4,graph(eps_p,epsdndeps_p,Straight),p+pencorrection*pnp);
if (mineps!=0)
    draw(pic4,(mineps,log(min_epsdndeps)/log(10))--(mineps,log(max_epsdndeps)/log(10)),p+dashed); // -- - Straight
if (maxeps<epsmax&&maxeps!=0)
    draw(pic4,(maxeps,log(min_epsdndeps)/log(10))--(maxeps,log(max_epsdndeps)/log(10)),p+dashed); // -- - Straight
xlimits(pic4,0,epsmax);
ylimits(pic4,min_epsdndeps,max_epsdndeps,Crop);
xaxis(pic4,"$\varepsilon$, MeV",BottomTop,LeftTicks,true);
yaxis(pic4,"$\varepsilon \times dN/d\varepsilon$",LeftRight,RightTicks,true);
size(pic4,1.5*picturewidth*1cm,1.5*picturewidth*1cm,point(pic4,SW),point(pic4,NE));
frame frm4=pic4.fit();
frm4=shift(2*xindent*picturewidth*1cm-picturewidth*1cm,-yindent*picturewidth*1cm-1.5*picturewidth*log(min_epsdndeps)/(log(max_epsdndeps)-log(min_epsdndeps))*1cm)*frm4;
add(frm4);

picture pic5;
draw(pic5,graph(ts,dNdOmegadts,Straight),p+pencorrection*png);
xlimits(pic5,0,nts*dts);
xaxis(pic5,"$ct/\lambda$",BottomTop,LeftTicks,true);
yaxis(pic5,"$dN/d\Omega dt (\mathbf \Omega),\; c/\lambda$",LeftRight,RightTicks,true);
size(pic5,1.5*picturewidth*1cm,1.5*picturewidth*1cm,point(pic5,SW),point(pic5,NE));
frame frm5=pic5.fit();
frm5=shift(-0.5*picturewidth*1cm,-yindent*picturewidth*1cm)*frm5;
add(frm5);

// Размеры итогового рисунка в сантиметрах
write( (max(currentpicture)-min(currentpicture))*0.03528 ); // US pt -> cm
