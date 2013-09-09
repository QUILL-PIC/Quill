import graph;
defaultpen(linewidth(0.7)+fontsize(10));

// Если Nmax задан равными нулю, то его значение вычисляется из
// данных, если Nmin задан меньшим нуля (или Nmin=0 и scale="log"), то
// Nmin вычисляется как минимальное положительное значение.
// Если t1<=t0, то строится график от t0 до максимального значения t

real Nmax = 0;
real Nmin = 0;

real t0 = 0;
real t1 = 0;

string scale = "log"; // possible: linear, log
string particles = "epph"; // possible: e, p, ph, ep, eph, pph, epph

//----------------------------------------

real dt;

string results_folder = "../results/";
file fin_param = input(results_folder+"log", comment="");
string var_name="";
while (var_name!="$")
{
    var_name = fin_param;
    if (var_name=="dt") dt = fin_param;
}

//----------------------------------------

file fin = input(results_folder+"N");
real[] data = fin;

int n = Floor(data.length/3);
int n0 = Floor(t0/dt);
int n1 = Floor(t1/dt);

if (n1>n||n1<=n0)
{
    n1 = n;
    t1 = (n-1)*dt;
}
n = n1 - n0;

real[] t = new real[n];
real[] N_e = new real[n];
real[] N_p = new real[n];
real[] N_ph = new real[n];

for (int i=0;i<n;i=i+1)
{
    t[i] = (i+n0)*dt;
    N_e[i] = data[3*i];
    N_p[i] = data[3*i+1];
    N_ph[i] = data[3*i+2];
}

//----------------------------------------

write("N_e = ",format("%g",N_e[n-1]));
write("N_p = ",format("%g",N_p[n-1]));
write("N_ph = ",format("%g",N_ph[n-1]));

//----------------------------------------

if (Nmax==0)
{
    if (particles=="e"||particles=="ep"||particles=="eph"||particles=="epph")
	Nmax = max(N_e);
    if (particles=="p"||particles=="ep"||particles=="pph"||particles=="epph")
	Nmax = max(Nmax,max(N_p));
    if (particles=="ph"||particles=="eph"||particles=="pph"||particles=="epph")
	Nmax = max(Nmax,max(N_ph));
}
if (Nmin<0||(Nmin==0&&scale=="log"))
{
    Nmin = Nmax;
    if (particles=="e"||particles=="ep"||particles=="eph"||particles=="epph")
    {
	for (int i=0;i<n;i=i+1)
	    if (N_e[i]<Nmin&&N_e[i]>0) Nmin = N_e[i];
    }
    if (particles=="p"||particles=="ep"||particles=="pph"||particles=="epph")
    {
	for (int i=0;i<n;i=i+1)
	    if (N_p[i]<Nmin&&N_p[i]>0) Nmin = N_p[i];
    }
    if (particles=="ph"||particles=="eph"||particles=="pph"||particles=="epph")
    {
	for (int i=0;i<n;i=i+1)
	    if (N_ph[i]<Nmin&&N_ph[i]>0) Nmin = N_ph[i];
    }
}

//----------------------------------------

pen p = linejoin(2)+linewidth(1.0);

picture pic1;
if (scale=="log")
{
    scale(pic1,Linear,Log);
    real m;
    m = log(Nmin)/log(10) - 1;
    for (int i=0;i<n;i=i+1)
    {
	if (N_e[i]!=0) N_e[i] = log(N_e[i])/log(10); else N_e[i] = m;
	if (N_p[i]!=0) N_p[i] = log(N_p[i])/log(10); else N_p[i] = m;
	if (N_ph[i]!=0) N_ph[i] = log(N_ph[i])/log(10); else N_ph[i] = m;
    }
}
if (particles=="e"||particles=="ep"||particles=="eph"||particles=="epph")
{
    draw(pic1,graph(t,N_e,Straight),p+blue);
}
if (particles=="p"||particles=="ep"||particles=="pph"||particles=="epph")
    draw(pic1,graph(t,N_p,Straight),p+red);
if (particles=="ph"||particles=="eph"||particles=="pph"||particles=="epph")
    draw(pic1,graph(t,N_ph,Straight),p+0.7*green);
xlimits(pic1,t0,t1);
ylimits(pic1,Nmin,Nmax,Crop);
xaxis(pic1,"$ct/\lambda$",BottomTop,LeftTicks,false);
yaxis(pic1,"$N$",LeftRight,RightTicks,true);
size(pic1,8cm,6cm,point(pic1,SW),point(pic1,NE));
frame frm1=pic1.fit();
add(frm1);

// Размеры итогового рисунка в сантиметрах
write( (max(currentpicture)-min(currentpicture))*0.03528 ); // US pt -> cm
