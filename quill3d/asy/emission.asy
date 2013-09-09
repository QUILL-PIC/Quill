import graph;
import palette;
defaultpen(linewidth(0.7)+fontsize(10));

// max_inv, max_f, max_dndeps или max_density могут быть заданы равными нулю, в этом
// случае они будут найдены из массивов данных.
//
// Если задать значение epsmax большее, чем максимальная энергия,
// использовавшаяся в PIC-коде, то переменной epsmax будет присвоено
// значение, использовавшееся в PIC-коде. В случае, если epsmax
// меньше, чем значение, использовавшееся в PIC-коде, neps заменяется
// на Floor(epsmax/deps).
//
// Если min_dndeps = 0, то используется значение min_dndeps =
// 1e-5*max_dndeps.

real file_number = 12;
real max_density = 1500; // ncr
real max_f = 5e10; // 1/MeV
real max_inv = 0; // (mc\omega/e)^2
real max_dndeps = 1e11;
real min_dndeps = 0;
real dpsi = 0.05;
real epsmax = 1e10; // MeV
string plane = "xy"; // possible: xy, xz, yz
real psi0 = -3*pi/4; // на графика f psi \in (psi0, psi1)
real psi1 = 10; /* если psi1>2 pi и psi0<0, то psi \in (psi0, 2
pi+psi0), если же psi1>2 pi и psi0>0, то psi \in (psi0, 2 pi) */
if (psi1>2*pi&&psi0<0) psi1=2*pi+psi0;
else if(psi1>2*pi) psi1 = 2*pi;
string particles = "photons"; // possible: electrons, positrons, photons
int npalette = 200;

//----------------------------------------

if (particles=="positrons") particles = "_p";
else if (particles=="photons") particles = "_ph";
else particles = "";
real csign;
if (particles=="") csign = -1;
else csign = 1;

//----------------------------------------

real dx;
real dy;
real dz;
int nx;
int ny;
int nz;
real lambda;
real xlength;
real ylength;
real zlength;
real deps;
int neps;
int enthp;
real deps_p;
int neps_p;
int enthp_p;
real deps_ph;
int neps_ph;
int enthp_ph;

//----------------------------------------

string results_folder = "../results/";
file fin_param = input( results_folder+"log", comment="" );
string var_name;
while (var_name!="$")
{
    var_name = fin_param;
    if (var_name=="dx") dx = fin_param;
    if (var_name=="dy") dy = fin_param;
    if (var_name=="dz") dz = fin_param;
    if (var_name=="nx") nx = fin_param;
    if (var_name=="ny") ny = fin_param;
    if (var_name=="nz") nz = fin_param;
    if (var_name=="lambda") lambda = fin_param;
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
xlength = dx*nx;
ylength = dy*ny;
zlength = dz*nz;
if (particles=="_p")
{
    deps = deps_p;
    neps = neps_p;
    enthp = enthp_p;
}
else if (particles=="_ph")
{
    deps = deps_ph;
    neps = neps_ph;
    enthp = enthp_ph;
}
epsmax = min(neps*deps,epsmax);
neps = Floor(epsmax/deps);

//----------------------------------------

file fin=input(results_folder+"rho"+particles+format("%g",file_number));
real[] data=fin;

int i,j,k;
real length1,length2;
string axis1_name,axis2_name;
real[][] density;
if (plane=="yz")
{
    length1 = ylength;
    length2 = zlength;
    axis1_name = "$y/\lambda$";
    axis2_name = "$z/\lambda$";
    density = new real[ny][nz];
    for( j=0; j<ny; j=j+1 )
    {
	for( k=0; k<nz; k=k+1 )
	{
	    density[j][k] = csign*data[(ny+nz)*nx+j*nz+k];
	}
    }
}
else if (plane=="xz")
{
    length1 = xlength;
    length2 = zlength;
    axis1_name = "$x/\lambda$";
    axis2_name = "$z/\lambda$";
    density = new real[nx][nz];
    for( i=0; i<nx; i=i+1 )
    {
	for( k=0; k<nz; k=k+1 )
	{
	    density[i][k] = csign*data[(ny+nz)*i+ny+k];
	}
    }
}
else
{
    length1 = xlength;
    length2 = ylength;
    axis1_name = "$x/\lambda$";
    axis2_name = "$y/\lambda$";
    density = new real[nx][ny];
    for( i=0; i<nx; i=i+1 )
    {
	for( j=0; j<ny; j=j+1 )
	{
	    density[i][j] = csign*data[(ny+nz)*i+j];
	}
    }
}

if (max_density==0) max_density = max(density);
range density_range = Range(0,max_density);

//----------------------------------------

fin=input(results_folder+"inv"+format("%g",file_number));
data=fin;

real[][] inv;
if (plane=="yz")
{
    inv = new real[ny][nz];
    for( j=0; j<ny; j=j+1 )
    {
	for( k=0; k<nz; k=k+1 )
	{
	    inv[j][k] = data[(ny+nz)*nx+j*nz+k];
	}
    }
}
else if (plane=="xz")
{
    inv = new real[nx][nz];
    for( i=0; i<nx; i=i+1 )
    {
	for( k=0; k<nz; k=k+1 )
	{
	    inv[i][k] = data[(ny+nz)*i+ny+k];
	}
    }
}
else
{
    inv = new real[nx][ny];
    for( i=0; i<nx; i=i+1 )
    {
	for( j=0; j<ny; j=j+1 )
	{
	    inv[i][j] = data[(ny+nz)*i+j];
	}
    }
}

if (max_inv==0) max_inv = max(max(inv),-min(inv));
range inv_range = Range(-max_inv,max_inv);

//----------------------------------------

fin=input(results_folder+"phasespace"+particles+format("%g",file_number));
data=fin;

int n_particles = Floor(data.length/8);
int n_psi = Floor((psi1-psi0)/dpsi);
real q;
real u1;
real u2;
real g;
real psi;
real[][] f = new real[neps][n_psi];
for( i=0; i<neps; i=i+1 )
{
    for( j=0; j<n_psi; j=j+1 )
    {
	f[i][j] = 0;
    }
}
if (plane=="yz")
{
    for( i=0; i<n_particles; i=i+1 )
    {
	q = data[8*i];
	u1 = data[8*i+5];
	u2 = data[8*i+6];
	g = data[8*i+7];
	j = Floor(g*0.511/deps);
	if (u1>0&&u2>=0) psi = atan(u2/u1);
	else if (u1<0) psi = atan(u2/u1) + pi;
	else if (u1!=0) psi = atan(u2/u1) + 2*pi;
	else psi=0;
	psi = psi - psi0;
	if(psi>2*pi) psi = psi - 2*pi;
	k = Floor(psi/dpsi);
	if(j>-1&&j<neps&&k>-1&&k<n_psi) f[j][k] = f[j][k] + csign*q;
    }
}
else if (plane=="xz")
{
    for( i=0; i<n_particles; i=i+1 )
    {
	q = data[8*i];
	u1 = data[8*i+4];
	u2 = data[8*i+6];
	g = data[8*i+7];
	j = Floor(g*0.511/deps);
	if (u1>0&&u2>=0) psi = atan(u2/u1);
	else if (u1<0) psi = atan(u2/u1) + pi;
	else if (u1!=0) psi = atan(u2/u1) + 2*pi;
	else psi=0;
	psi = psi - psi0;
	if(psi>2*pi) psi = psi - 2*pi;
	k = Floor(psi/dpsi);
	if(j>-1&&j<neps&&k>-1&&k<n_psi) f[j][k] = f[j][k] + csign*q;
    }
}
else
{
    for( i=0; i<n_particles; i=i+1 )
    {
	q = data[8*i];
	u1 = data[8*i+4];
	u2 = data[8*i+5];
	g = data[8*i+7];
	j = Floor(g*0.511/deps);
	if (u1>0&&u2>=0) psi = atan(u2/u1);
	else if (u1<0) psi = atan(u2/u1) + pi;
	else if (u1!=0) psi = atan(u2/u1) + 2*pi;
	else psi=0;
	psi = psi - psi0;
	if(psi>2*pi) psi = psi - 2*pi;
	k = Floor(psi/dpsi);
	if(j>-1&&j<neps&&k>-1&&k<n_psi) f[j][k] = f[j][k] + csign*q;
    }
}

real tmp; /* tmp =
n_cr[cm^{-3}]*dx[cm]*dy[cm]*dz[cm]/(deps*dpsi) */
tmp = 1.11485e13*lambda*dx*dy*dz/(deps*dpsi);
for( i=0; i<neps; i=i+1 )
{
    for( j=0; j<n_psi; j=j+1 )
    {
	f[i][j] = f[i][j]*tmp*enthp;
    }
}

if (max_f==0) max_f = max(f);
range f_range = Range(0,max_f);

//----------------------------------------

fin=input(results_folder+"spectrum"+particles+format("%g",file_number));

real[] dndeps = new real[neps];
real[] eps = new real[neps];
data=fin;
for( i=1;i<neps;i=i+1 )
{
    dndeps[i] = data[i];
}
dndeps[0] = dndeps[1];

for( i=0;i<neps;i=i+1 )
{
    eps[i] = i*deps;
}

for( i=0;i<neps;i=i+1 )
{
    if (dndeps[i]>0)
	dndeps[i] = log(dndeps[i])/log(10);
    else
	dndeps[i] = -100;
}

if (max_dndeps==0) max_dndeps = max(dndeps);
if (min_dndeps==0) min_dndeps = 1e-5*max_dndeps;

//----------------------------------------

pen[] palette1 = new pen[npalette];

pen Tickpen=black;
pen tickpen=gray+0.5*linewidth(currentpen);
string[] plbl1 = {"0",format("%g",0.5*max_density),format("%g",max_density)};

// pn1 - pn2 - pn3
int n_intermediate = Floor( npalette/10 );
pen pn1 = white;
pen pn2 = 0.8*red+0.85*green+1.0*blue;
pen pn3 = black;
for(i=0;i<n_intermediate;i=i+1)
{
    palette1[i] = (1-i/(n_intermediate-1))*pn1 + i/(n_intermediate-1)*pn2;
}
for(i=n_intermediate;i<npalette;i=i+1)
{
    palette1[i] = (1-(i-n_intermediate)/(npalette-1-n_intermediate))*pn2 + (i-n_intermediate)/(npalette-1-n_intermediate)*pn3;
}

//----------------------------------------

pen[] palette2 = new pen[npalette];

string[] plbl2 = {format("%g",-max_inv),"0",format("%g",max_inv)};

// pn1 - pn2 - pn3
n_intermediate = Floor( npalette/2 );
pn1 = 1/7*green+blue;
pn2 = white;
pn3 = red;
for(i=0;i<n_intermediate;i=i+1)
{
    palette2[i] = (1-i/(n_intermediate-1))*pn1 + i/(n_intermediate-1)*pn2;
}
for(i=n_intermediate;i<npalette;i=i+1)
{
    palette2[i] = (1-(i-n_intermediate)/(npalette-1-n_intermediate))*pn2 + (i-n_intermediate)/(npalette-1-n_intermediate)*pn3;
}

//------------------------------

pen[] palette3 = new pen[npalette];

string[] plbl3 = {"0",format("%g",max_f/2),format("%g",max_f)};
string[] psi_ticks = {"$-2 \pi$","$-\frac{3\pi}{2}$","$-\pi$","$-\frac{\pi}{2}$","$0$","$\frac{\pi}{2}$","$\pi$","$\frac{3\pi}{2}$","$2\pi$"};

// pn1 - pn2 - pn3
n_intermediate = Floor( npalette/3 );
pn1 = white;
pn2 = 0.4*red+0.8*green;
pn3 = black;
for(i=0;i<n_intermediate;i=i+1)
{
    palette3[i] = (1-i/(n_intermediate-1))*pn1 + i/(n_intermediate-1)*pn2;
}
for(i=n_intermediate;i<npalette;i=i+1)
{
    palette3[i] = (1-(i-n_intermediate)/(npalette-1-n_intermediate))*pn2 + (i-n_intermediate)/(npalette-1-n_intermediate)*pn3;
}

//------------------------------

real picturewidth = 7; // 6.7;
real xindent = 0.06;
real yindent = 0.04;

picture pic1;
bounds prange=image(pic1,density,density_range,(0,0),(length1,length2),palette1);
palette(pic1,"$n/n_{cr}$",prange,point(pic1,NW)+0.04N*length2,point(pic1,NE)+0.08N*length2,Top,palette1,PaletteTicks(new string(real x){return plbl1[Floor(2*x/max_density)];},max_density/2,max_density/4,Tickpen,tickpen));
xlimits(pic1,0,length1);
ylimits(pic1,0,length2);
xaxis(pic1,axis1_name,BottomTop,LeftTicks,true);
yaxis(pic1,axis2_name,LeftRight,RightTicks,true);
size(pic1,picturewidth*1cm,picturewidth*1cm*length2/length1,point(pic1,SW),point(pic1,NE));
pic1 = shift(0,(point(pic1,NE).y-point(pic1,SW).y)*log(min_dndeps)/log(max_dndeps/min_dndeps))*pic1;
frame frm1=pic1.fit();
add(frm1);

picture pic2;
prange=image(pic2,inv,inv_range,(0,0),(length1,length2),palette2);
palette(pic2,"$\left(e/mc\omega\right)^2\left({\mathbf E}^2-{\mathbf B}^2\right)$",prange,point(pic2,NW)+0.04N*length2,point(pic2,NE)+0.08N*length2,Top,palette2,PaletteTicks(new string(real x){return plbl2[Floor(x/max_inv+1)];},max_inv,max_inv/2,Tickpen,tickpen));
xlimits(pic2,0,length1);
ylimits(pic2,0,length2);
xaxis(pic2,axis1_name,BottomTop,LeftTicks,true);
yaxis(pic2,axis2_name,LeftRight,RightTicks,true);
size(pic2,picturewidth*1cm,picturewidth*1cm*length2/length1,point(pic2,SW),point(pic2,NE));
pic2 = shift(0,(point(pic2,NE).y-point(pic2,SW).y)*log(min_dndeps)/log(max_dndeps/min_dndeps))*pic2;
frame frm2=pic2.fit();
frm2=shift(0,min(frm1).y-max(frm1).y-yindent*picturewidth*1cm)*frm2;
add(frm2);

picture pic3;
prange=image(pic3,f,f_range,(0,psi0),(epsmax,psi1),palette3);
palette(pic3,"$d^2N/d\varepsilon d\psi$",prange,point(pic3,NW)+0.04N*(psi1-psi0),point(pic3,NE)+0.08N*(psi1-psi0),Top,palette3,PaletteTicks(new string(real x){return plbl3[Floor(2*x/max_f)];},max_f/2,max_f/4,Tickpen,tickpen));
xlimits(pic3,0,epsmax);
ylimits(pic3,psi0,psi1);
xaxis(pic3,"$\varepsilon$, MeV",BottomTop,LeftTicks,true);
yaxis(pic3,"$\psi$",LeftRight,RightTicks(new string(real x){return psi_ticks[Floor(2*x/pi)+4];},pi/2,pi/4),true);
size(pic3,picturewidth*1cm,picturewidth*1cm*length2/length1,point(pic3,SW),point(pic3,NE));
pic3 = shift(0,-psi0+(point(pic3,NE).y-point(pic3,SW).y)*log(min_dndeps)/log(max_dndeps/min_dndeps))*pic3;
frame frm3=pic3.fit();
frm3=shift(max(frm1).x-min(frm1).x+xindent*picturewidth*1cm,0)*frm3;
add(frm3);

pen p_line = linejoin(2)+linewidth(1.0);
pen p_dndeps = p_line + 0.8*red+0.8*blue;
picture pic4;
scale(pic4,Linear,Log);
draw(pic4,graph(eps,dndeps,Straight),p_dndeps);
xlimits(pic4,0,epsmax);
ylimits(pic4,min_dndeps,max_dndeps,Crop); /* Crop + shift(...)*pic4
работают некорректно, поэтому сдвигаются остальные pictures, а не pic4 */
xaxis(pic4,"$\varepsilon$, MeV",BottomTop,LeftTicks,true);
yaxis(pic4,"$\lg dN/d\varepsilon$",LeftRight,RightTicks,true);
size(pic4,picturewidth*1cm,picturewidth*1cm*length2/length1,point(pic4,SW),point(pic4,NE));
frame frm4=pic4.fit();
frm4=shift(max(frm1).x-min(frm1).x+xindent*picturewidth*1cm, min(frm1).y-max(frm1).y-yindent*picturewidth*1cm)*frm4;
add(frm4);

// Размеры итогового рисунка в сантиметрах
write( (max(currentpicture)-min(currentpicture))*0.03528 ); // US pt -> cm
