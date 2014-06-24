// Строит зависимости f1, f2, f3 от x на оси области моделирования

// qwe! в случае вывода не всех компонент поля (yz, e.g.) даёт ошибку

string f1name,f2name,f3name; // possible: ex, ey, ez, bx, by, bz
f1name = "bx"; // red
f2name = "by"; // green
f3name = "bz"; // blue
real flim = 10;

//----------------------------------------

real file_number = 16;

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

int i;

//----------------------------------------

file fin=input(results_folder+f1name+format("%g",file_number));
real[] data=fin;

real[] f1 = new real[nx];
real[] x = new real[nx];

for (i=0;i<nx;i=i+1)
{
    x[i] = i*dx;
    f1[i] = data[Floor(ny/2)+i*(ny+nz)];
}

//----------------------------------------

fin=input(results_folder+f2name+format("%g",file_number));
data=fin;

real[] f2 = new real[nx];

for (i=0;i<nx;i=i+1)
{
    f2[i] = data[Floor(ny+nz/2)+i*(ny+nz)];
}

//----------------------------------------

fin=input(results_folder+f3name+format("%g",file_number));
data=fin;

real[] f3 = new real[nx];

for (i=0;i<nx;i=i+1)
{
    f3[i] = data[Floor(ny+nz/2)+i*(ny+nz)];
}

//----------------------------------------

import graph;
import palette;
defaultpen(linewidth(0.7)+fontsize(10));

//----------------------------------------

real picturewidth = 10;

pen p = linejoin(2)+linewidth(0.5);

picture pic1;
draw(pic1,graph(x,f1,Straight),p+red);
draw(pic1,graph(x,f2,Straight),p+0.7*green);
draw(pic1,graph(x,f3,Straight),p+blue);
xlimits(pic1,0,dx*nx);
ylimits(pic1,-flim,flim,Crop);
//qwe; ylimits(pic1,min(f3),max(f3),Crop);
/*qwe real tmpmax = 0;
real tmp;
for (i=0;i<nx;i=i+1) {
    tmp = f1[i]*f1[i]+f2[i]*f2[i]+f3[i]*f3[i];
    if (tmpmax<tmp)
	tmpmax = tmp;
}
write(tmpmax);
*/
xaxis(pic1,"$x/\lambda$",BottomTop,LeftTicks,true);
yaxis(pic1,"$eE/mc\omega$",LeftRight,RightTicks,true);
size(pic1,picturewidth*1cm,picturewidth*1cm,point(pic1,SW),point(pic1,NE));
frame frm1=pic1.fit();
add(frm1);

// Размеры итогового рисунка в сантиметрах
write( (max(currentpicture)-min(currentpicture))*0.03528 ); // US pt -> cm
