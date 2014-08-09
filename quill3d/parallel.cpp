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

using namespace std;
typedef queue<TaskBase*> queueBase;

class Parallelizator
{
	queueBase* task_queue;
	int num_threads;
	int num_tasks;
	bool test_mode;
	string param_x_name;
	string param_y_name;
	vector<string> param_x_value;
	vector<string> param_y_value;
	
public:
	Parallelizator() {task_queue = new queueBase;};
	~Parallelizator() {delete task_queue;};
	queueBase* getTaskQueue() {return task_queue;};
	void init(int argc, char** argv);
	void populateTaskQueue();
	void performTasks();
	static void logError(string message);
	static void logErrorAndExit(string message);
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
		long long time = rand()*80000000/RAND_MAX; //awful and not working somehow
		usleep(time);
		
		pthread_mutex_lock(&print_mutex);
		cout << "Completed TestTask # " << task_id << ", thread id " << tid << endl;
		pthread_mutex_unlock(&print_mutex);	
		
		return 0;
	}
};

// Performs a QuillTask - runs Quill using a given config file

int QuillTask::perform(int tid)
{
	string command = "./parse.sh .parallel/" + config_name + " | ./quill";

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
	    cout << "usage: parallel [-threads value] [-test] [-help]" << endl;
	    // the program should stop here, exit() may be not the
	    // best for it
	    exit(0);
	}
	for (int i=1; i<argc; i++)
	{
		cout << argv[i] << endl;
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
		else if (strcmp(argv[i], "-help") == 0)
		{
		    cout << "usage: parallel [-threads value] [-test] [-help]\n";
		    cout << "value - the number of threads\n";
		    cout << "parameter matrix is parallel.ini\n";
		    cout << "template config is ../quill3d-conf/quill.conf.parallel/config" << endl;
		    exit(0);
		}
	}
	if (num_threads == -1)
		num_threads = N_THREADS_DEFAULT;
		
	cout << "Welcome to Parallelizator v.1.0. Running on " << num_threads << " threads" << endl;
	
	ifstream settings_file("parallel.ini");
	if (!settings_file.is_open() && !test_mode)
		logErrorAndExit("unable to open settings file \"parallel.ini\"");

	string buffer;
	while (getline(settings_file, buffer))
	{
		cout << "Read line: " << buffer << endl;
		if (buffer.find_first_of('#') == 0)
			//ignoring comments
			continue;
			
		if (buffer.find("a0") == 0)
		{
			istringstream iss(buffer);
			iss >> param_x_name;
			int value;
			while (iss >> value)
			{
				stringstream ss;
				ss << value;
				param_x_value.push_back(ss.str());
			}
		}
		if (buffer.find("ne") == 0)
		{
			istringstream iss(buffer);
			iss >> param_y_name;
			int value;
			while (iss >> value)
			{
				stringstream ss;
				ss << value;
				param_y_value.push_back(ss.str());
			}
		}
		else if (buffer.find("filmwidth") == 0)
		{
			istringstream iss(buffer);
			iss >> param_y_name;
			int value;
			while (iss >> value)
			{
				stringstream ss;
				ss << value;
				param_y_value.push_back(ss.str());
			}
		}
	}
	
	if (param_x_name.empty() || param_y_name.empty() || param_x_value.empty() || param_y_value.empty())
	{
		stringstream ss;
		ss << "bad info in parallel.ini, param_x_name = " << param_x_name << ", param_y_name = " << param_y_name;
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
		ifstream config_template("../quill3d-conf/quill.conf.parallel/config"); //will be customizeable in the future release
		if (!config_template.is_open())
			logErrorAndExit("unable to open config template \"../quill3d-conf/quill.conf.parallel\"");
			
		string buffer;
		vector<string> config_data;
		bool param_x_found = false;
		bool param_y_found = false;
		while (getline(config_template, buffer))
		{
			if (buffer.find(param_x_name) == 0)
			{
				param_x_found = true;
				buffer = "$PARAM_X$";
			}
			else if (buffer.find(param_y_name) == 0)
			{
				param_y_found = true;
				buffer = "$PARAM_Y$";
			}
			config_data.push_back(buffer);
//			cout << buffer << endl;
		}
		config_template.close();
		
		if (!(param_x_found && param_y_found))
			logErrorAndExit("Config template doesn't match settings file, unable to continue");
		
		int temp_task_id = 0;
			
		for (vector<string>::iterator it_x = param_x_value.begin(); it_x != param_x_value.end(); ++it_x)
		{
			for (vector<string>::iterator it_y = param_y_value.begin(); it_y != param_y_value.end(); ++it_y)
			{
				stringstream file_name;
				file_name << "../quill3d-conf/quill.conf.parallel/" << param_x_name << *it_x << "-" << param_y_name << *it_y;
//				cout << "  DEBUG: file_name " << file_name.str();
				ofstream my_config(file_name.str().c_str());
				if (!my_config.is_open())
					logErrorAndExit("Failed to create config file");
				
				for (vector<string>::iterator it_data = config_data.begin(); it_data != config_data.end(); it_data++)
				{
					if ((*it_data).compare("$PARAM_X$") == 0)
					{
						//in this code, ne is always param_y; check anyway param_x because this can be changed in the future
						my_config << param_x_name << " = " << *it_x << (param_x_name.find("ne") != string::npos ? " ncr" : "") << endl;
					}
					else if ((*it_data).compare("$PARAM_Y$") == 0)
					{
						my_config << param_y_name << " = " << *it_y << (param_y_name.find("ne") != string::npos ? " ncr" : "") << endl;
					}
					else
						my_config << *it_data << endl;
				}

				stringstream config_name;
				config_name << param_x_name << *it_x << "-" << param_y_name << *it_y;
				my_config << "data_folder = results_" << config_name.str();
				my_config.close();

				QuillTask* pTask = new QuillTask(++temp_task_id, config_name.str());
				task_queue->push(pTask);
			}
		}
	}
	
// commented code is for logging content of Task queue
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

// Logs error to console

void Parallelizator::logError(string message)
{
	cout << "ERROR: " << message << endl;
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
