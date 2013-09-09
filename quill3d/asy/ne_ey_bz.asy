/* строит зависимость концентрации электронов, а также e_y и b_z от x
 * на оси области моделирования */

real file_number = 5.4;
real x1 = 9; // в длинах волн
real x2 = 11; // если x2=0, то приравнивается равным nx*dx
/* если norm=0, то определяется по максимальному значению
 * автоматически */
real nnorm = 0;
real enorm = 0;
real bnorm = 0;

//----------------------------------------
// чтение параметров из log-файла

real dx;
int nx,ny,nz;

string results_folder = "../results/";
file fin_param = input( results_folder+"log", comment="" );
string var_name;
while (var_name!="$")
{
    var_name = fin_param;
    if (var_name=="dx") dx = fin_param;
    if (var_name=="nx") nx = fin_param;
    if (var_name=="ny") ny = fin_param;
    if (var_name=="nz") nz = fin_param;
}

//----------------------------------------

if (x2==0)
    x2 = nx*dx;
int n1 = Floor(x1/dx);
int n2 = Floor(x2/dx);
int n = n2-n1;

//----------------------------------------

file fin=input(results_folder+"rho"+format("%g",file_number));
real[] data=fin;

real[] ne = new real[n];
real[] x = new real[n];

int i;
for (i=n1;i<n2;i=i+1)
{
    x[i-n1] = i*dx;
    ne[i-n1] = data[Floor(ny/2)+i*(ny+nz)];
}

if (nnorm==0)
{
    nnorm = max(max(ne),-min(ne));
    write("nnorm = ",nnorm);
    if (nnorm==0) nnorm = 1;
}
for (i=n1;i<n2;i=i+1)
{
    ne[i-n1] = ne[i-n1]/nnorm;
}

//----------------------------------------

fin=input(results_folder+"ey"+format("%g",file_number));
data=fin;

real[] ey = new real[n];

for (i=n1;i<n2;i=i+1)
{
    ey[i-n1] = data[Floor(ny/2)+i*(ny+nz)];
}

if (enorm==0)
{
    enorm = max(max(ey),-min(ey));
    write("enorm = ",enorm);
    if (enorm==0) enorm = 1;
}
for (i=n1;i<n2;i=i+1)
{
    ey[i-n1] = ey[i-n1]/enorm;
}

//----------------------------------------

fin=input(results_folder+"bz"+format("%g",file_number));
data=fin;

real[] bz = new real[n];

for (i=n1;i<n2;i=i+1)
{
    bz[i-n1] = data[Floor(ny/2)+i*(ny+nz)];
}

if (bnorm==0)
{
    bnorm = max(max(bz),-min(bz));
    write("bnorm = ",bnorm);
    if (bnorm==0) bnorm = 1;
}
for (i=n1;i<n2;i=i+1)
{
    bz[i-n1] = bz[i-n1]/bnorm;
}

//----------------------------------------

import graph;
defaultpen(linewidth(0.7)+fontsize(10));

//----------------------------------------

real picturewidth = 8;

pen p = linejoin(2)+linewidth(0.7);

picture pic1;
draw(pic1,graph(x,ne,Straight),p+0.5*green);
draw(pic1,graph(x,ey,Straight),p+red);
draw(pic1,graph(x,bz,Straight),p+blue);
xlimits(pic1,x1,x2);
ylimits(pic1,-1,1,Crop);
xaxis(pic1,"$x/\lambda$",BottomTop,LeftTicks,true);
yaxis(pic1,"$n_e,\;E_y,\;B_z$",LeftRight,RightTicks,true);
size(pic1,picturewidth*1cm,picturewidth*1cm,point(pic1,SW),point(pic1,NE));
frame frm1=pic1.fit();
add(frm1);

// Размеры итогового рисунка в сантиметрах
write( (max(currentpicture)-min(currentpicture))*0.03528 ); // US pt -> cm
