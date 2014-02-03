// Строит зависимости rho и w от x на оси области моделирования

real file_number = 16.12;
real x1 = 10; // lambda
real x2 = 15; // если x2=0, то приравнивается равным nx*dx
/* если norm=0, то определяется по максимальному значению
 * автоматически */
real enorm = 0;
real pnorm = 0;
real gnorm = 0;
real wnorm = 0;

//----------------------------------------

string f1name,f2name,f3name,f4name;
f1name = "rho_p"; // red
f2name = "rho_ph"; // blue
f3name = "rho"; // green
f4name = "w"; // orange

//----------------------------------------

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

int i;

//----------------------------------------

file fin=input(results_folder+f1name+format("%g",file_number));
real[] data=fin;

real[] f1 = new real[n];
real[] x = new real[n];

for (i=n1;i<n2;i=i+1)
{
    x[i-n1] = i*dx;
    f1[i-n1] = data[Floor(ny/2)+i*(ny+nz)];
}
if (pnorm==0)
{
    pnorm = max(f1);
    write("pnorm = ",pnorm);
    if (pnorm==0) pnorm = 1;
}
for (i=n1;i<n2;i=i+1)
{
    f1[i-n1] = f1[i-n1]/pnorm;
}

//----------------------------------------

fin=input(results_folder+f2name+format("%g",file_number));
data=fin;

real[] f2 = new real[n];

for (i=n1;i<n2;i=i+1)
{
    f2[i-n1] = data[Floor(ny+nz/2)+i*(ny+nz)];
}
if (gnorm==0)
{
    gnorm = max(f2);
    write("gnorm = ",gnorm);
    if (gnorm==0) gnorm = 1;
}
for (i=n1;i<n2;i=i+1)
for (i=n1;i<n2;i=i+1)
{
    f2[i-n1] = f2[i-n1]/gnorm;
}

//----------------------------------------

fin=input(results_folder+f3name+format("%g",file_number));
data=fin;

real[] f3 = new real[n];

for (i=n1;i<n2;i=i+1)
{
    f3[i-n1] = -data[Floor(ny+nz/2)+i*(ny+nz)];
}
if (enorm==0)
{
    enorm = max(f3);
    write("enorm = ",enorm);
    if (enorm==0) enorm = 1;
}
for (i=n1;i<n2;i=i+1)
{
    f3[i-n1] = f3[i-n1]/enorm;
}

//----------------------------------------

fin=input(results_folder+f4name+format("%g",file_number));
data=fin;

real[] f4 = new real[n];

for (i=n1;i<n2;i=i+1)
{
    f4[i-n1] = data[Floor(ny+nz/2)+i*(ny+nz)];
}
if (wnorm==0)
{
    wnorm = max(f4);
    write("wnorm = ",wnorm);
    if (wnorm==0) wnorm = 1;
}
for (i=n1;i<n2;i=i+1)
{
    f4[i-n1] = f4[i-n1]/wnorm;
}

//----------------------------------------

import graph;
import palette;
defaultpen(linewidth(0.7)+fontsize(10));

//----------------------------------------

real picturewidth = 10;

pen p = linejoin(2)+linewidth(0.5);

picture pic1;
if (f1name!="") draw(pic1,graph(x,f1,Straight),p+red);
//if (f2name!="") draw(pic1,graph(x,f2,Straight),p+blue);
if (f3name!="") draw(pic1,graph(x,f3,Straight),p+0.7*green);
//if (f4name!="") draw(pic1,graph(x,f4,Straight),p+red+0.3*green);
xlimits(pic1,x1,x2);
xaxis(pic1,"$x/\lambda$",BottomTop,LeftTicks,true);
yaxis(pic1,"$n$",LeftRight,RightTicks,true);
size(pic1,picturewidth*1cm,picturewidth*1cm,point(pic1,SW),point(pic1,NE));
frame frm1=pic1.fit();
add(frm1);

// Размеры итогового рисунка в сантиметрах
write( (max(currentpicture)-min(currentpicture))*0.03528 ); // US pt -> cm
