real start = 1;
real stop = 5.1;
real step = 0.1;
real file_name_accuracy = 10;

string script_name = getstring("","phasespace","script to draw: ",store=false);
write(script_name);
system("./make_script.sh "+script_name);

int n;
int n0;
int n1;
n0 = floor(start/step);
n1 = floor(stop/step)+1;
real current_number;
for(n=n0;n<n1;n=n+1)
{
    current_number = floor(n*step*file_name_accuracy)/file_name_accuracy;
    write(current_number);
    file fout=output("number");
    write(fout,current_number,endl);
    system("asy "+script_name+".asy");
    system("convert -density 100 "+script_name+".eps tmp.png"); 
    if (n==n0) system("./make_background.sh");
    system("convert background.png tmp.png -composite frame"+format("%05d",n)+".png");
}
system("rm background.png tmp.png number "+script_name+".asy "+script_name+".eps");
system("./moviemaker.sh "+script_name);
