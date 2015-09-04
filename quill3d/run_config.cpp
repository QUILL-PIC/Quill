/* run_config.cpp; this program runs quill through MPI processes for various
initial parameters (i.e., a0, n0).
Compiling:
$ mpic++ run_config.cpp
Usage:
First, initial config file ("conf") should be prepared. Use
either
$ ./parse.sh .suffix > conf
or 
$ make conf FILE='.suffix'
for human-readable-config "../quill3d-conf/quill.conf.suffix".
Second, provide "myhostfile" containing hostnames for quill execution, for
example, "myhostfile" file can contain:
n10
n11
n12
Notice that you should have ssh access to the nodes n10-n12 with rsa keys, and
nodes should be free from other heavy computations.
Then, run the program with the option -np equal to the number of parameter sets,
i.e. a0_size*ne_size (do not forget to use screen):
$ mpirun --hostfile myhostfile -np 12 ./a.out
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

	const int arr_size=3,a0_size=1,ne_size=31;
	std::string str[arr_size]={"ne","a0","data_folder"};
	double a0[a0_size]={15};
	double ne[ne_size]={0.01689643,  0.01841711,  0.02007465,  0.02188137,  0.02385069,
          0.02599726,  0.02833701,  0.03088734,  0.0336672 ,  0.03669725,
          0.04      ,  0.0436    ,  0.047524  ,  0.05180116,  0.05646326,
          0.06154496,  0.067084  ,  0.07312156,  0.07970251,  0.08687573,
          0.09469455,  0.10321706,  0.11250659,  0.12263218,  0.13366908,
          0.1456993 ,  0.15881224,  0.17310534,  0.18868482,  0.20566645,
          0.22417643};

	
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
	param_d_f<<"results_ne"<<arr_struct_parameters[id].ne<<"-a0"<<arr_struct_parameters[id].a0;
	arr_struct_parameters[id].data_folder=param_d_f.str();

	set_parameter(file_content,"ne",arr_struct_parameters[id].ne);
	set_parameter(file_content,"a0",arr_struct_parameters[id].a0);
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
	new_dir<<"mkdir results_ne"<<arr_struct_parameters[id].ne<<"-a0"<<arr_struct_parameters[id].a0;
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
