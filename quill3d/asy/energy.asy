import graph;
defaultpen(linewidth(0.7)+fontsize(10));

// Если energymax задан равными нулю, то его значение вычисляется из
// данных, если energymin задан меньшим нуля (или energymin=0 и scale="log"), то
// energymin вычисляется как минимальное положительное значение.
// Если t1<=t0, то строится график от t0 до максимального значения t

real energymin = 0e2;
real energymax = 0e4;

real t0 = 0; // lambda/c
real t1 = 0;

string scale = "linear"; // possible: linear, log
string species = "ofepgi"; // possible: o (overall), f (fields), e, p, g (gamma-quants), i (ions), ep, eg, etc.

//----------------------------------------

bool overall,fields,electrons,positrons,gamma,ions;
overall = find(species,"o") != -1;
fields = find(species,"f") != -1;
electrons = find(species,"e") != -1;
positrons = find(species,"p") != -1;
gamma = find(species,"g") != -1;
ions = find(species,"i") != -1;

//----------------------------------------

real dt;
int n_ion_populations;

string results_folder = "../results/";
file fin_param = input(results_folder+"log", comment="");
string var_name="";
while (var_name!="$")
{
    var_name = fin_param;
    if (var_name=="dt") dt = fin_param;
    if (var_name=="n_ion_populations") n_ion_populations = fin_param;
}

int ncolumns = 4 + n_ion_populations;

//----------------------------------------

file fin = input(results_folder+"energy");
real[] data = fin;

int n = Floor(data.length/(ncolumns));
int n0 = Floor(t0/dt);
int n1 = Floor(t1/dt);

if (n1>n||n1<=n0)
{
    n1 = n;
    t1 = (n-1)*dt;
}
n = n1 - n0;

real[] t = new real[n];
real[] energy_o = new real[n];
real[] energy_f = new real[n];
real[] energy_e = new real[n];
real[] energy_p = new real[n];
real[] energy_g = new real[n];
real[][] energy_i = new real[n_ion_populations][n];

for (int i=0;i<n;i=i+1)
{
    t[i] = (i+n0)*dt;
    energy_f[i] = data[ncolumns*i];
    energy_e[i] = data[ncolumns*i+1];
    energy_p[i] = data[ncolumns*i+2];
    energy_g[i] = data[ncolumns*i+3];
    energy_o[i] = 0;
    for (int j=0;j<n_ion_populations;j=j+1)
    {
	energy_i[j][i] = data[ncolumns*i+4+j];
	energy_o[i] += energy_i[j][i];
    }
    energy_o[i] += energy_f[i]+energy_e[i]+energy_p[i]+energy_g[i];
}

write("overall_energy [J; MeV] = "+format("%g",energy_o[n-1])+"; "+format("%g",energy_o[n-1]*6.24e12));
write("f_energy [J; MeV] = "+format("%g",energy_f[n-1])+"; "+format("%g",energy_f[n-1]*6.24e12));
write("e_energy [J; MeV] = "+format("%g",energy_e[n-1])+"; "+format("%g",energy_e[n-1]*6.24e12));
write("p_energy [J; MeV] = "+format("%g",energy_p[n-1])+"; "+format("%g",energy_p[n-1]*6.24e12));
write("ph_energy [J; MeV] = "+format("%g",energy_g[n-1])+"; "+format("%g",energy_g[n-1]*6.24e12));
for (int j=0;j<n_ion_populations;j=j+1)
    write("i_energy [J; MeV] = "+format("%g",energy_i[j][n-1])+"; "+format("%g",energy_i[j][n-1]*6.24e12));

//----------------------------------------

if (energymax==0)
{
    if (overall==true)
    {
	energymax = max(energy_o);
	energymax = 1.1*energymax;
    }
    else
    {
	if (electrons==true)
	    energymax = max(energy_e);
	if (positrons==true)
	    energymax = max(energymax,max(energy_p));
	if (gamma==true)
	    energymax = max(energymax,max(energy_g));
	if (fields==true)
	    energymax = max(energymax,max(energy_f));
	if (ions==true)
	{
	    for (int j=0;j<n_ion_populations;j=j+1)
		energymax = max(energymax,max(energy_i[j]));
	}
    }
}
if (energymin<0||(energymin==0&&scale=="log"))
{
    energymin = energymax;
    if (electrons==true)
    {
	for (int i=0;i<n;i=i+1)
	    if (energy_e[i]<energymin&&energy_e[i]>0) energymin = energy_e[i];
    }
    if (positrons==true)
    {
	for (int i=0;i<n;i=i+1)
	    if (energy_p[i]<energymin&&energy_p[i]>0) energymin = energy_p[i];
    }
    if (gamma==true)
    {
	for (int i=0;i<n;i=i+1)
	    if (energy_g[i]<energymin&&energy_g[i]>0) energymin = energy_g[i];
    }
    if (fields==true)
    {
	for (int i=0;i<n;i=i+1)
	    if (energy_f[i]<energymin&&energy_f[i]>0) energymin = energy_f[i];
    }
    if (ions==true)
    {
	for (int j=0;j<n_ion_populations;j=j+1)
	{
	    for (int i=0;i<n;i=i+1)
		if (energy_i[j][i]<energymin&&energy_i[j][i]>0) energymin = energy_i[j][i];
	}
    }
}

//----------------------------------------

pen p = linejoin(2)+linewidth(1.0);

picture pic1;
if (scale=="log")
{
    scale(pic1,Linear,Log);
    real m;
    m = log(energymin)/log(10) - 1;
    for (int i=0;i<n;i=i+1)
    {
	if (energy_e[i]>0) energy_e[i] = log(energy_e[i])/log(10); else energy_e[i] = m;
	if (energy_p[i]>0) energy_p[i] = log(energy_p[i])/log(10); else energy_p[i] = m;
	if (energy_g[i]>0) energy_g[i] = log(energy_g[i])/log(10); else energy_g[i] = m;
	if (energy_f[i]>0) energy_f[i] = log(energy_f[i])/log(10); else energy_f[i] = m;
	if (energy_o[i]>0) energy_o[i] = log(energy_o[i])/log(10); else energy_o[i] = m;
	for (int j=0;j<n_ion_populations;j=j+1)
	    if (energy_i[j][i]>0) energy_i[j][i] = log(energy_i[j][i])/log(10); else energy_i[j][i] = m;
    }
}
if (electrons==true)
    draw(pic1,graph(t,energy_e,Straight),p+blue+0.1*(red+green));
if (positrons==true)
    draw(pic1,graph(t,energy_p,Straight),p+red);
if (gamma==true)
    draw(pic1,graph(t,energy_g,Straight),p+0.3*green);
if (fields==true)
    draw(pic1,graph(t,energy_f,Straight),p);
if (overall==true)
    draw(pic1,graph(t,energy_o,Straight),p+dashed);
if (ions==true)
{
    srand(seconds());
    real r = rand()/randMax;
    real g = rand()/randMax;
    real b = 0.5*rand()/randMax;
    for (int j=0;j<n_ion_populations;j=j+1)
    {
	draw(pic1,graph(t,energy_i[j],Straight),p+r*red+g*green+b*blue);
    }
}
xlimits(pic1,t0,t1);
ylimits(pic1,energymin,energymax,Crop);
xaxis(pic1,"$ct/\lambda$",BottomTop,LeftTicks,false);
yaxis(pic1,"$W$, J",LeftRight,RightTicks,true);
size(pic1,8cm,6cm,point(pic1,SW),point(pic1,NE));
frame frm1=pic1.fit();
add(frm1);

// Размеры итогового рисунка в сантиметрах
write( (max(currentpicture)-min(currentpicture))*0.03528 ); // US pt -> cm
