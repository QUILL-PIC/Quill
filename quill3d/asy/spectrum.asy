import graph;
defaultpen(linewidth(0.7)+fontsize(10));

// Если задать значение epsmax большее, чем максимальная энергия,
// использовавшаяся в PIC-коде, или равное нулю, то переменной epsmax
// будет присвоено значение, использовавшееся в PIC-коде (максимальное
// по сортам частиц). В случае, если epsmax меньше, чем значение,
// использовавшееся в PIC-коде, neps заменяется на Floor(epsmax/deps).
//
// Если sp_type="energy", то строятся не распределение частиц по энергиям, а
// распределение полной энергии частиц по энергиям eps*dN/deps
//
// Если min_dndeps=0 и scale="log", то min_dndeps вычисляется как 1e-5*max_dndeps

real file_number = 12;
real max_dndeps = 1e14;
real min_dndeps = 1e0;
real max_dndeps_energy = 0e10;
real min_dndeps_energy = 0;
real epsmax = 1e4; // MeV

string particles = "i"; // possible: e, p, g, i, ep, igp, etc.
string scale = "log"; // possible: linear, log
string sp_type = "simple"; // possible: simple, energy

//----------------------------------------

bool electrons,positrons,gamma,ions;
electrons = find(particles,"e") != -1;
positrons = find(particles,"p") != -1;
gamma = find(particles,"g") != -1;
ions = find(particles,"i") != -1;

//----------------------------------------

real type;
if (sp_type=="energy")
{
    max_dndeps = max_dndeps_energy;
    min_dndeps = min_dndeps_energy;
    type = 1;
}
else
{
    type = 0;
}

real overall_en_e=0;
real overall_en_p=0;
real overall_en_g=0;
real Ne = 0;
real Np = 0;
real Ng = 0;

real deps_e;
int neps_e;
real deps_p;
int neps_p;
real deps_g;
int neps_g;

int n_ion_populations;
real deps_i;
int neps_i;

//----------------------------------------

string results_folder = "../results/";
file fin_param = input( results_folder+"log", comment="" );
string var_name;
while (var_name!="$")
{
    var_name = fin_param;
    if (var_name=="deps") deps_e = fin_param;
    if (var_name=="neps") neps_e = fin_param;
    if (var_name=="deps_p") deps_p = fin_param;
    if (var_name=="neps_p") neps_p = fin_param;
    if (var_name=="deps_ph") deps_g = fin_param;
    if (var_name=="neps_ph") neps_g = fin_param;
    if (var_name=="n_ion_populations") n_ion_populations = fin_param;
    if (var_name=="deps_i") deps_i = fin_param;
    if (var_name=="neps_i") neps_i = fin_param;
}

if (!electrons)
    neps_e = 0;
if (!positrons)
    neps_p = 0;
if (!gamma)
    neps_g = 0;
if (!ions)
    neps_i = 0;

if (epsmax>max(deps_e*neps_e,deps_p*neps_p,deps_g*neps_g,deps_i*neps_i)||epsmax==0)
    epsmax = max(deps_e*neps_e,deps_p*neps_p,deps_g*neps_g,deps_i*neps_i);
else
{
    if(epsmax<deps_e*neps_e)
	neps_e = Floor(epsmax/deps_e);
    if(epsmax<deps_p*neps_p)
	neps_p = Floor(epsmax/deps_p);
    if(epsmax<deps_g*neps_g)
	neps_g = Floor(epsmax/deps_g);
    if(epsmax<deps_i*neps_i)
	neps_i = Floor(epsmax/deps_i);
}

//----------------------------------------

file fin;
real[] data;
int i;

real[] dndeps_e = new real[neps_e];
real[] eps_e = new real[neps_e];
real[] dndeps_p = new real[neps_p];
real[] eps_p = new real[neps_p];
real[] dndeps_g = new real[neps_g];
real[] eps_g = new real[neps_g];

real[][] dndeps_i = new real[n_ion_populations][neps_i];
real[] eps_i = new real[neps_i];

real max_dndeps_e = 0;
real max_dndeps_p = 0;
real max_dndeps_g = 0;
real max_dndeps_i = 0;

//----------------------------------------

if (electrons)
{
    fin=input(results_folder+"spectrum"+format("%g",file_number));

    data=fin;

    for( i=0;i<neps_e;i=i+1 )
    {
	eps_e[i] = (i+0.5)*deps_e;
    }

    for( i=0;i<neps_e;i=i+1 )
    {
	dndeps_e[i] = (1-type+type*eps_e[i])*data[i];
	overall_en_e += eps_e[i]*data[i]*deps_e;
	Ne += data[i]*deps_e;
    }
    //dndeps_e[0] = dndeps_e[1];

    max_dndeps_e = max(dndeps_e);

    write("overall electron energy [MeV] = "+format("%g",overall_en_e));
    write("Ne = "+format("%g",Ne));
    if (Ne!=0)
	write("<eps_e>[MeV] = "+format("%g",overall_en_e/Ne));
}

//----------------------------------------

if (positrons)
{
    fin=input(results_folder+"spectrum_p"+format("%g",file_number));

    data=fin;

    for( i=0;i<neps_p;i=i+1 )
    {
	eps_p[i] = (i+0.5)*deps_p;
    }

    for( i=0;i<neps_p;i=i+1 )
    {
	dndeps_p[i] = (1-type+type*eps_p[i])*data[i];
	overall_en_p += eps_p[i]*data[i]*deps_p;
	Np += data[i]*deps_p;
    }
    //dndeps_p[0] = dndeps_p[1];

    max_dndeps_p = max(dndeps_p);

    write("overall positron energy [MeV] = "+format("%g",overall_en_p));
    write("Np = "+format("%g",Np));
    if (Np!=0)
	write("<eps_p>[MeV] = "+format("%g",overall_en_p/Np));
}

//----------------------------------------

if (gamma)
{
    fin=input(results_folder+"spectrum_ph"+format("%g",file_number));

    data=fin;

    for( i=0;i<neps_g;i=i+1 )
    {
	eps_g[i] = (i+0.5)*deps_g;
    }

    for( i=0;i<neps_g;i=i+1 )
    {
	dndeps_g[i] = (1-type+type*eps_g[i])*data[i];
	overall_en_g += eps_g[i]*data[i]*deps_g;
	Ng += data[i]*deps_g;
    }
    //dndeps_g[0] = dndeps_g[1];

    max_dndeps_g = max(dndeps_g);

    write("overall photon energy [MeV] = "+format("%g",overall_en_g));
    write("Ng = "+format("%g",Ng));
    if (Ng!=0)
	write("<eps_g>[MeV] = "+format("%g",overall_en_g/Ng));
}

if (max_dndeps==0) max_dndeps = max(max_dndeps_e,max_dndeps_p,max_dndeps_g);
if (min_dndeps==0&&scale=="log")
    min_dndeps = 1e-5*max_dndeps;

//----------------------------------------

real[] icmr = new real[n_ion_populations];
int ii=0;
fin_param = input( results_folder+"log", comment="" );
var_name="";
while (var_name!="$")
{
    var_name = fin_param;
    if (var_name=="icmr")
    {
	icmr[ii] = fin_param;
	ii += 1;
    }
}
if (ions)
{
    for( i=0;i<neps_i;i=i+1 )
    {
	eps_i[i] = (i+0.5)*deps_i;
    }
    for (ii=0;ii<n_ion_populations;ii+=1)
    {
	fin=input(results_folder+"spectrum_"+format("%g",icmr[ii])+"_"+format("%g",file_number));
	data=fin;
	for( i=0;i<neps_i;i=i+1 )
	{
	    dndeps_i[ii][i] = (1-type+type*eps_i[i])*data[i];
	}
	//dndeps_i[ii][0] = dndeps_i[ii][1];
	max_dndeps_i = max(max_dndeps_i,max(dndeps_i[ii]));
	write("icmr = "+format("%g",icmr[ii]));
    }
}

if (max_dndeps==0)
    max_dndeps = max(max_dndeps_e,max_dndeps_p,max_dndeps_g,max_dndeps_i);
if (min_dndeps==0&&scale=="log")
    min_dndeps = 1e-5*max_dndeps;

//----------------------------------------

pen p = linejoin(2)+linewidth(1.0);

picture pic1;
if (scale=="log")
{
    real m;
    m = log(min_dndeps)/log(10)-1;
    scale(pic1,Linear,Log);
    if (electrons)
    {
	for (i=0;i<neps_e;i=i+1) {
	    if (dndeps_e[i]!=0) dndeps_e[i] = log(dndeps_e[i])/log(10);
	    else dndeps_e[i]=m;
	}
    }
    if (positrons)
    {
	for (i=0;i<neps_p;i=i+1) {
	    if (dndeps_p[i]!=0) dndeps_p[i] = log(dndeps_p[i])/log(10);
	    else dndeps_p[i]=m;
	}
    }
    if (gamma)
    {
	for (i=0;i<neps_g;i=i+1) {
	    if (dndeps_g[i]!=0) dndeps_g[i] = log(dndeps_g[i])/log(10);
	    else dndeps_g[i]=m;
	}
    }
    if (ions)
    {
	for (ii=0;ii<n_ion_populations;ii+=1)
	{
	    for (i=0;i<neps_i;i=i+1)
	    {
		if (dndeps_i[ii][i]!=0) dndeps_i[ii][i] = log(dndeps_i[ii][i])/log(10);
		else dndeps_i[ii][i]=m;
	    }
	}
    }
}
if (electrons)
    draw(pic1,graph(eps_e,dndeps_e,Straight),p+blue);
if (positrons)
    draw(pic1,graph(eps_p,dndeps_p,Straight),p+red);
if (gamma)
    draw(pic1,graph(eps_g,dndeps_g,Straight),p+0.7*green);
if (ions)
{
    real r,g,b;
    srand(seconds());
    for (ii=0;ii<n_ion_populations;ii+=1)
    {
	r = rand()/randMax;
	b = rand()/randMax;
	g = 1 - 0.5*(r+b);
	draw(pic1,graph(eps_i,dndeps_i[ii],Straight),p+r*red+g*green+b*blue);
    }
}
xlimits(pic1,0,epsmax);
ylimits(pic1,min_dndeps,max_dndeps,Crop);
xaxis(pic1,"$\varepsilon$, MeV",BottomTop,LeftTicks,true);
if (sp_type=="energy")
    yaxis(pic1,"$\varepsilon \times dN/d\varepsilon$",LeftRight,RightTicks,true);
else
    yaxis(pic1,"$dN/d\varepsilon$, $1/$MeV",LeftRight,RightTicks,true);
size(pic1,8cm,6cm,point(pic1,SW),point(pic1,NE));
frame frm1=pic1.fit();
add(frm1);

// Размеры итогового рисунка в сантиметрах
write( (max(currentpicture)-min(currentpicture))*0.03528 ); // US pt -> cm
