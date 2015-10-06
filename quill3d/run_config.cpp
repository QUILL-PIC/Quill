/* run_config.cpp; this program runs quill through MPI processes for various
initial parameters (i.e., a0, n0).
Change the section marked *changeme*, set parameters there and set df (data
folder, used for naming of result folders). Units of ne, a0, etc are the same as
in initial conf file.
Compiling:
$ mpic++ run_config.cpp
Usage:
First, initial config file ("conf") should be prepared. Use
either
$ ./parse.sh .suffix > conf
or 
$ make conf FILE='.suffix'
for human-readable-config "../quill3d-conf/quill.conf.suffix".
Second, provide hostfile containing hostnames for quill execution, for
example, "hfile" file that contains:
n10
n11
n12
Notice that you should have ssh access to the nodes n10-n12 with rsa keys, and
nodes should be free from other heavy computations (check this by $ gstat -a).
You shoud run the program with the option -np equal to the number of parameter sets,
i.e. a0_size * ne_size (do not forget to use screen also):
$ mpirun --quiet --hostfile hfile -np 12 ./a.out
However, it is preferable to use *torque*. For this purpose, check file job.sh,
write paths in it and run (submit):
$ qsub job.sh
Look for quill's output in job.sh.o* and job.sh.e*. Use qstat to check status.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>  
#include <cstdlib>

#include "mpi.h"

//struct of parameters
struct parameters{
	double a0,ne;
	std::string data_folder;
};

bool set_parameter(std::string& content, const std::string& param_name,double& param_value);
bool set_parameter(std::string& content, const std::string& param_name,std::string& param_value);

int main()
{
	int np,id;

  /*changeme*/
	const int arr_size=3,a0_size=5,ne_size=3;
	std::string str[arr_size]={"ne","theta","data_folder"};
	double a0[a0_size]={0, 5, 10, 20, 35};
	double ne[ne_size]={9e22, 18e22, 36e22};
  std::string df = "results1";
  /*changeme*/

	
	//make variants of parameters
	parameters arr_struct_parameters[a0_size*ne_size];
	int k=0;
	for(int i=0;i<a0_size;i++){
		for(int j=0;j<ne_size;j++){
			arr_struct_parameters[k].a0=a0[i];
			arr_struct_parameters[k].ne=ne[j];
			k++;
		}
	}

	//code with mpi
	MPI::Init();
	np=MPI::COMM_WORLD.Get_size();
	id=MPI::COMM_WORLD.Get_rank();


	std::ifstream input_stream("conf");
	if( !input_stream.is_open() )
	{
		std::cout << "No configuration file was found" << std::endl;
		return 0;
	}

	std::stringstream parser;
	parser << input_stream.rdbuf();
	input_stream.close();

	std::string file_content = parser.str();

	//make param data_folder
	std::ostringstream param_d_f;
	param_d_f<<df+"_"+str[0]<<arr_struct_parameters[id].ne<<"-"+str[1]<<arr_struct_parameters[id].a0;
	arr_struct_parameters[id].data_folder=param_d_f.str();

	set_parameter(file_content,str[0],arr_struct_parameters[id].ne);
	set_parameter(file_content,str[1],arr_struct_parameters[id].a0);
	set_parameter(file_content,"data_folder",arr_struct_parameters[id].data_folder);

	//make filename
	std::ostringstream ss;
	std::string result;
	ss<<"file"<<id<<".conf";
	result=ss.str();

	std::ofstream output_stream(result.c_str());
	output_stream << file_content;
	output_stream.close();

	//make dir
	std::ostringstream new_dir;
	new_dir<<"mkdir "+df+"_"+str[0]<<arr_struct_parameters[id].ne<<"-"+str[1]<<arr_struct_parameters[id].a0;
	std::string command_mkdir=new_dir.str();
	system(command_mkdir.c_str());

	std::ostringstream running;
	running<<"cat "<<result<<"| ./quill";
	std::string command=running.str();
	system(command.c_str());


	MPI::Finalize();
	return 0;
}

bool set_parameter(std::string& content, const std::string& param_name,double& param_value)
{
	using std::string;
	
	string::size_type start_pos = string::npos;
	start_pos = content.find(param_name)+1; // start_pos pointes to first '\n'
	if(start_pos == string::npos)
		return false;

	start_pos +=param_name.size() + 0; // now start_pos pointes to param_name
									  // param_name.size + "\n".size()
	std::ostringstream o;
	std::string for_param_value;
	o<<param_value;
	for_param_value=o.str();


	string::size_type stop_pos = string::npos;
	stop_pos = content.find('\n',start_pos); // +1 - to avoid '\n' after param_name

	content.replace(content.begin() + start_pos, content.begin() + stop_pos, for_param_value.c_str());

	return true;
}

bool set_parameter(std::string& content, const std::string& param_name,std::string& param_value)
{
	using std::string;
	
	string::size_type start_pos = string::npos;
	start_pos = content.find(param_name)+1; // start_pos pointes to first '\n'
	if(start_pos == string::npos)
		return false;

	start_pos+=param_name.size()+2; //param_name.size+"\n".size()+"#".size()+"\n".size()

	std::ostringstream o;
	std::string for_param_value;
	o<<param_value;
	for_param_value=o.str();


	string::size_type stop_pos = string::npos;
	stop_pos = content.find('\n',start_pos); // +1 - to avoid '\n' after param_name

	content.replace(content.begin() + start_pos, content.begin() + stop_pos, for_param_value.c_str());

	return true;
}
