// --------------------------- Parallelizator v 1.1 ---------------------------
//
// Parallelizator is a tool which runs Quill codes on several threads in parallel.
// Parameters of Quill codes can be set in 2D matrix in settings file parallel.ini
//
// To run several quill codes in parallel, one needs perform the following steps:
// 1. Build Parallelizator executable
// 2. Open Quill/quill3d-conf/quill.conf.parallel/ and place a config file here which will be used as template
// 3. In parallel.ini, set path to your template in "config_template_path" variable
// 4. Make sure that location which is in "temp_config_folder" (in parallel.ini) exists
// 5. Choose 2 parameters which should be changed during Quill runs, make sure that they exist in your config file
//    If only one parameter needs to be changed, you may take anything as the second one and assign a constant value to it, see below
// 6. Write both parameters in the format listed below (also see example in parallel.ini):
//
//    <name of parameter> <list of values separated with whitespace> dimension:none|<parameter dimension>
//
//    List should contain at least one value; if exactly one, this should be used in all runs
//    Otherwise it will be subsequently replaced with other values from the list
// 7. Run Parallelizator using the following command:
//
//    ./parallel -threads <number of threads>
//
//    where <number of threads> is the desired number of simultaneously running Quills.
//    If number of tasks (which is (number of values in X-list) * (number of values in Y-list)) is greater than number of threads,
//    a task queue is formed, and thread who first gets free will take next task.
// 8. Enjoy!
//
// ----------------------------------------------------------------------------


// Changelog

// 2014-08-08 
// Dmitry Serebryakov
// v.1.0 - initial version

// 2014-08-09 
// Eugene Nerush
// v.1.0.1 - added support for filmwidth parameter, added command line -help option

// 2014-08-11 
// Dmitry Serebryakov
// v.1.1 - any parameter from config template can be used in parallel.ini (not only ne and a0)


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
#include <pthread.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#define N_THREADS_DEFAULT 6
#define N_THREADS_MAX     64
#define N_TASKS_DEFAULT   80

pthread_t threads[N_THREADS_MAX];
pthread_mutex_t print_mutex;
pthread_mutex_t task_mutex;

//forward declarations
class Parallelizator;
Parallelizator* TheParallelizator;
void* threadFun(void*);
class TaskBase;

bool option_1d = 0;

using namespace std;
typedef queue<TaskBase*> queueBase;

class Param
{
public:
	string name;
	vector<string> value;
	string dimension;
	string comment;
	
	bool isInitialized() {if (name.empty() || value.empty()) return false; return true;};
};

class Parallelizator
{
	queueBase* task_queue;
	int num_threads;
	int num_tasks;
	bool test_mode;
	Param param_x, param_y;
	
public:
	Parallelizator() : temp_config_folder("../quill3d-conf/quill.conf.parallel/"), config_template_path("../quill3d-conf/quill.conf.example") {task_queue = new queueBase;};
	~Parallelizator() {delete task_queue;};
	queueBase* getTaskQueue() {return task_queue;};
	void init(int argc, char** argv);
	void populateTaskQueue();
	void performTasks();
	static void logError(string message);
	static void logErrorAndExit(string message);

	string temp_config_folder;
	string config_template_path;
};

class TaskBase
{
protected:
	int task_id;
public:
	int getTaskId() {return task_id;}
	TaskBase(int ti) : task_id(ti) {};
	virtual int perform(int) = 0;
};

class QuillTask: public TaskBase
{
	string config_name;

public:
	QuillTask(int ti, string cn) : TaskBase(ti), config_name(cn) {};
	int perform(int tid);
};

class TestTask: public TaskBase
{
public:
	TestTask(int ti) : TaskBase(ti) {};
	int perform(int tid)
	{
		pthread_mutex_lock(&print_mutex);
		cout << "Begin TestTask # " << task_id << ", thread id " << tid << endl;
		pthread_mutex_unlock(&print_mutex);	

		// Just simulation of task execution
		double time = ((double) rand()*10.0 / (RAND_MAX)); 
		sleep((long)time);
		
		pthread_mutex_lock(&print_mutex);
		cout << "Completed TestTask # " << task_id << ", thread id " << tid << endl;
		pthread_mutex_unlock(&print_mutex);	
		
		return 0;
	}
};

// Performs a QuillTask - runs Quill using a given config file

int QuillTask::perform(int tid)
{
	if (TheParallelizator == NULL)
		Parallelizator::logErrorAndExit("Parallelizator is NULL");

	string command;
	if (TheParallelizator->temp_config_folder.compare("../quill3d-conf/quill.conf.parallel/") != 0)
		command = "./parse.sh " + TheParallelizator->temp_config_folder + config_name + " | ./quill";
	else
		command = "./parse.sh .parallel/" + config_name + " | ./quill";

	pthread_mutex_lock(&print_mutex);
	cout << "Begin QuillTask # " << task_id << ", thread id " << tid << ", config_name " << config_name << endl;
	cout << command << endl;
	pthread_mutex_unlock(&print_mutex);

	if (system(("mkdir -p results_" + config_name).c_str()))
	{
		pthread_mutex_lock(&print_mutex);
		cout << "Error creating directory" << endl;
		pthread_mutex_unlock(&print_mutex);	
	}
//	int rc = 0;
	int rc = system(command.c_str());
	
	pthread_mutex_lock(&print_mutex);
	cout << "QuillTask # " << task_id << " ended with status " << rc << ", thread id " << tid << endl;
	pthread_mutex_unlock(&print_mutex);	
	return rc;
}

// Performs initialization of Parallelizatior class.

void Parallelizator::init(int argc, char** argv)
{
	cout << argc << endl;
	test_mode = false; 
	num_threads = -1;  //uninitialized
	if (argc == 1)
	{
	    cout << "usage: parallel [-1d] [-threads value] [-test] [-help]" << endl;
	    // the program should stop here, exit() may be not the
	    // best for it
	    exit(0);
	}
	for (int i=1; i<argc; i++)
	{
		if (strcmp(argv[i], "-threads") == 0)
		{
			if (i+1 < argc)
			{
				num_threads = atoi(argv[i+1]);
				if (num_threads > N_THREADS_MAX)
				{
					logError("Max number of threads exceeded");
					num_threads = N_THREADS_MAX;
				}
			}
			else
				logError("option -threads missing an argument");
		} 
		else if (strcmp(argv[i], "-test") == 0)
		{
			test_mode = true; //don't run Quill at all
			num_tasks = N_TASKS_DEFAULT;
		} 
		else if (strcmp(argv[i], "-1d") == 0) {
		    option_1d = 1;
		}
		else if (strcmp(argv[i], "-help") == 0)
		{
		    cout << "usage: parallel [-1d] [-threads value] [-test] [-help]\n";
		    cout << "value - the number of threads\n";
		    cout << "parameter matrix is in parallel.ini\n";
		    cout << "template config is ../quill3d-conf/quill.conf.parallel/config_template" << endl;
		    cout << "if option -1d is set then parallel.ini is treated as vector of pairs of parameters, not as matrix" << endl;
		    exit(0);
		}
	}
	if (num_threads == -1)
		num_threads = N_THREADS_DEFAULT;
		
	cout << "Welcome to Parallelizator running on " << num_threads << " threads" << endl;
	
	ifstream settings_file("parallel.ini");
	if (!settings_file.is_open() && !test_mode)
		logErrorAndExit("unable to open settings file \"parallel.ini\"");

	string buffer;
	bool settings_initialized = false;
	while (getline(settings_file, buffer))
	{
		cout << "Read line: " << buffer << endl;
		if (buffer.find_first_of('#') == 0 || buffer.empty())
			//ignoring comments and empty lines
			continue;
			
		if (buffer.find("---END OF SETTINGS---") == 0)
		{
			settings_initialized = true;
			cout << "Settings initialized, starting reading parameters" << endl;
			continue;
		}

		if (buffer.find("config_template_path=") == 0)
		{
			config_template_path = buffer.substr(strlen("config_template_path="));
			cout << "Initialized: config template path = " << config_template_path << endl;
		}
		if (buffer.find("temp_config_folder=") == 0)
		{
			temp_config_folder = buffer.substr(strlen("temp_config_folder="));
			cout << "Initialized: temp config folder = " << temp_config_folder << endl;
		}
		
		if (settings_initialized)
		{
			if (!param_x.isInitialized())
			{
				istringstream iss(buffer);
				iss >> param_x.name;
				string value;
				while (iss >> value)
				{
					stringstream ss;
					ss << value;
					if (ss.str().find("dimension:") == 0)
					{
						string dim = ss.str().substr(10);
						param_x.dimension = (dim.compare("none") == 0) ? "" : dim;
					}
					else
					{
						param_x.value.push_back(ss.str());
					}
				}
			}
			else if (!param_y.isInitialized())
			{
				istringstream iss(buffer);
				iss >> param_y.name;
				string value;
				while (iss >> value)
				{
					stringstream ss;
					ss << value;
					if (ss.str().find("dimension:") == 0)
					{
						string dim = ss.str().substr(10);
						param_y.dimension = (dim.compare("none") == 0) ? "" : dim;
					}
					else
					{
						param_y.value.push_back(ss.str());
					}
				}
			}
			else logError("Too many parameters in settings file");
		}
	}
	
	if (!param_x.isInitialized() || !param_y.isInitialized())
	{
		stringstream ss;
		ss << "bad info in parallel.ini, param_x.name = " << param_x.name << ", param_y.name = " << param_y.name;
		settings_file.close();
		logErrorAndExit(ss.str());
	}
	settings_file.close();
}

// Populates task queue with "tasks" which refer to config files.

void Parallelizator::populateTaskQueue()
{
	if (test_mode)
	{
		cout << "Test mode detected. Creating test tasks" << endl;
		for (int j=0; j<num_tasks; j++)
		{
			TestTask* pTask = new TestTask(j);
			task_queue->push(pTask);
		}
	}
	else
	{
		ifstream config_template(config_template_path.c_str());
		if (!config_template.is_open())
			logErrorAndExit("unable to open config template " + config_template_path);
			
		string buffer;
		vector<string> config_data;
		bool param_x_found = false;
		bool param_y_found = false;
		while (getline(config_template, buffer))
		{
			if (buffer.find(param_x.name) == 0)
			{
				param_x_found = true;
				int commentFrom = buffer.find_first_of('#');
				if (commentFrom != string::npos)
					param_x.comment = buffer.substr(commentFrom);
				buffer = "$PARAM_X$";
			}
			else if (buffer.find(param_y.name) == 0)
			{
				param_y_found = true;
				int commentFrom = buffer.find_first_of('#');
				if (commentFrom != string::npos)
					param_y.comment = buffer.substr(commentFrom);
				buffer = "$PARAM_Y$";
			}
			config_data.push_back(buffer);
		}
		config_template.close();
		
		if (!(param_x_found && param_y_found))
			logErrorAndExit("Config template doesn't match settings file, unable to continue");
		
		int temp_task_id = 0;
			
		vector<string>::iterator yit = param_y.value.begin( );
		for (vector<string>::iterator it_x = param_x.value.begin(); it_x != param_x.value.end(); ++it_x)
		{
		    vector<string> yvalues;
		    if ( option_1d == 1 ) {
			vector<string> tmp(1);
			tmp[0] = *yit;
			yvalues = tmp;
			yit++;
		    }
		    else
			yvalues = param_y.value;
		    for ( vector<string>::iterator it_y = yvalues.begin( ); it_y != yvalues.end ( ); ++it_y )
		    {
			    stringstream file_name;
			    file_name << temp_config_folder << param_x.name << *it_x << "-" << param_y.name << *it_y;
			    ofstream my_config(file_name.str().c_str());
			    if (!my_config.is_open())
				    logErrorAndExit("Failed to create config file");
			    
			    for (vector<string>::iterator it_data = config_data.begin(); it_data != config_data.end(); it_data++)
			    {
				    if ((*it_data).compare("$PARAM_X$") == 0)
				    {
					    //in this code, ne is always param_y; check anyway param_x because this can be changed in the future
					    my_config << param_x.name << " = " << *it_x << " " << param_x.dimension << " " << param_x.comment << endl;
				    }
				    else if ((*it_data).compare("$PARAM_Y$") == 0)
				    {
					    my_config << param_y.name << " = " << *it_y << " " << param_y.dimension << " " << param_y.comment << endl;
				    }
				    else if ((*it_data).find("data_folder = ") != string::npos)
					    //data folder will be added anyway, ignore this
					    continue;
				    else
					    my_config << *it_data << endl;
			    }

			    stringstream config_name;
			    config_name << param_x.name << *it_x << "-" << param_y.name << *it_y;
			    my_config << "data_folder = results_" << config_name.str();
			    my_config.close();

			    QuillTask* pTask = new QuillTask(++temp_task_id, config_name.str());
			    task_queue->push(pTask);
		    }
		}
	}
	
//commented code is for logging content of Task queue
//	queueBase queue2(*task_queue);
//	while (!queue2.empty())
//	{
//		cout << queue2.front()->getTaskId() << " ";
//		queue2.pop();
//	}
//	for (queue<TaskBase*>::iterator it = task_queue.begin(); it != task_queue.end(); ++it)
//	{
//		cout << "Task #" << (*it)->getTaskId() << endl;
//	}	

}

// Creates threads that perform tasks from the queue.

void Parallelizator::performTasks()
{
	int retVal;
	cout << "Performing parallelization of " << task_queue->size() << " tasks on " << num_threads << " threads" << endl;
	
	for (int i=0; i<num_threads; i++)
	{
		pthread_mutex_lock(&print_mutex);
		cout << "Starting thread #" << i << endl;
		pthread_mutex_unlock(&print_mutex);
		
		int* thread_id = new int(i);
		retVal = pthread_create(&threads[i], NULL, threadFun, (void*)thread_id);
		if (retVal != 0)
		{
			pthread_mutex_lock(&print_mutex);
			cout << "Thread #" << i << " creation failed with error code " << retVal << endl;
			pthread_mutex_unlock(&print_mutex);
			delete thread_id;
		}
	}
	//Waiting for all threads to complete
	for (int i=0; i<num_threads; i++)
	{
		pthread_join(threads[i], NULL);
	}
}

// Logs error to console and exits

void Parallelizator::logErrorAndExit(string message)
{
	logError(message);
	pthread_mutex_destroy(&print_mutex);
	pthread_mutex_destroy(&task_mutex);
	exit(1);
}

// Logs error to console (threadsafe)

void Parallelizator::logError(string message)
{
	pthread_mutex_lock(&print_mutex);
	cout << "ERROR: " << message << endl;
	pthread_mutex_unlock(&print_mutex);
}

// Staring point of each task executor thread.

void* threadFun(void* tid)
{
	pthread_mutex_lock(&print_mutex);
	cout << "This is thread # " << *((int*)tid) << endl;
	pthread_mutex_unlock(&print_mutex);
	
	queueBase* local_task_queue;
	while (1)
	{
		pthread_mutex_lock(&task_mutex);

		local_task_queue = TheParallelizator->getTaskQueue();
		if (local_task_queue->empty())
		{
			pthread_mutex_unlock(&task_mutex);
			break;
		}
		TaskBase* currentTask = local_task_queue->front();
		local_task_queue->pop();
		pthread_mutex_unlock(&task_mutex);
		
		//Performing task 
		int result = currentTask->perform(*((int*)tid));
		if (result)
			Parallelizator::logError("Task finished with unsuccessful status");
		delete currentTask;
	}
		
	pthread_mutex_lock(&print_mutex);
	cout << "Thread # " << *((int*)tid) << " exiting" << endl;
	pthread_mutex_unlock(&print_mutex);
	delete (int*)tid;
	pthread_exit(NULL);
}

int main(int argc, char** argv)
{
	pthread_mutex_init(&print_mutex, NULL);
	pthread_mutex_init(&task_mutex, NULL);
	srand(time(NULL));
	
	TheParallelizator = new Parallelizator();
	if (TheParallelizator != NULL)
	{
		TheParallelizator->init(argc, argv);
		TheParallelizator->populateTaskQueue();
		TheParallelizator->performTasks();
		delete TheParallelizator;
	}
	
	cout << "Exiting..." << endl;
	pthread_mutex_destroy(&print_mutex);
	pthread_mutex_destroy(&task_mutex);
	return 0;
}
