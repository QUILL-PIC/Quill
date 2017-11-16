#!/usr/bin/env python3

import collections
import hashlib
import math
import os
import re


# -- How to run multiple QUILLs with SLURM --
#
# 1. Create config template based on the existing Quill config file:
#      - put a variable ($<param>) instead of parameter value for every Quill parameter that needs to be changed;
#          variable name can be any but must start with $
#          (e.g. t_end = $my_t_end)
#      - put $data_folder after "data_folder = " in the config file
# 2. Modify config strings below (path to template, data folder, etc.) according to your needs
# 3. For variables that don't need to be changed, add corresponding entries to dict d
#      (e.g. d['t_end'] = 120)
# 4. For each variable that needs to be changed, add an entry to dict d_change
#      The key is variable name and the value is a list of possible values
#      (e.g. 'mcr' : [1, 2, 4])
#      In this case, each possible combination of variables will go to a different temporary config file.
#
#      One can also lock some of the variables so they will be changed only together.
#      For that, a tuple containing locked variables can be set as a key in d_change,
#      with value containing the list of possible combinations.
#      (e.g. ('x0film', 'filmwidth') : [[13.0, 1.0], [13.2, 0.8], [13.4, 0.6]] )
# 5. Update the slurm script as needed (set partition, requested number of cores, etc.)
# 6. Execute the script


# -- Configuration settings --

# Prefix to job name displayed in Slurm (squeue, sview, etc.)
job_name_prefix = 'quill'

# Prefix of config template
conf_prefix = '../quill3d-conf/quill.conf.'

# Path to config template
template_file = conf_prefix + 'stepanov-template'

# Folder for modified (temporary) config files
tmp_config_folder = conf_prefix + 'tmp_config_stepanov/'

# Slurm script
slurm_script = 'run_slurm.sh'

# Specifying QOS for tasks (e.g. node limits). Empty string if QOS is not used
qos_string = ''
# qos_string = '--qos=limit_10_nodes '

# Prefix for output data folders
df_prefix = 'res_stepanov'


d = {}

d['t_end'] = 120
d['mcr'] = 0.054
d['vx'] = 0.1
d['tflow'] = 2e-3
d['a0'] = 0

d_change = {'mcr' : [0.054, 0.0054, 1, 10], 'tflow' : [2e-3, 4e-3, 1e-3, 5e-4], 'vx' : [0.1, 0.2], 'a0' : [0]}
# d_change = {'a0' : [32, 55, 96], 'theta' : [12], ('dy', 'dz') : [(0.0125, 0.1)], 'nfilm' : [100, 300, 500]}
# d_change = {'a0' : [72], 'theta' : [theta], ('dy', 'dz') : [(0.025, 0.1)]}


def create_conf(template, d):
    new = []
    name = os.path.join(tmp_config_folder, template[len(conf_prefix):])
    with open(template, 'r') as f:
        for l in f:
            for s in d:
                l = re.sub(r'\$%s\b' % s, str(d[s]), l)
            new.append(l)
    name_suffix = ''
    for s in d:
        name_suffix += '_' + s + str(d[s])
    hasher = hashlib.md5()
    hasher.update(name_suffix.encode('utf-8'))
    name += '_' + hasher.hexdigest()
    print('Printing to %s' % name)
    if not os.path.isdir(tmp_config_folder):
        print('Creating %s' % tmp_config_folder)
        os.makedirs(tmp_config_folder)
    with open(name, 'w') as f:
        for l in new:
            f.write(l)
    return name


def run_confs(template, d, d_change, folder='res'):
    if len(d_change) == 0:
        run_conf(template, d, folder)
    else:
        od_change = collections.OrderedDict(d_change)
        k = list(od_change.keys())[0]
        for val in d_change[k]:
            if isinstance(k, tuple):
                if len(val) != len(k):
                    raise ValueError('Different length of {0} and {1}'.format(k, val))
                new_folder = folder
                for i in range(len(k)):
                    new_folder += '_' + k[i] + str(val[i])
                    d[k[i]] = val[i]
                new_d_change = d_change.copy()
                new_d_change.pop(k)
                run_confs(template, d, new_d_change, new_folder)
            else:
                new_folder = folder + '_' + k + str(val)
                d[k] = val
                new_d_change = d_change.copy()
                new_d_change.pop(k)
                run_confs(template, d, new_d_change, new_folder)


def run_conf(template, d, folder=None):
    name = job_name_prefix
    if folder is not None:
        d['data_folder'] = folder
        arr = folder.split('_')
        if len(arr) > 1:
            arr = arr[1:]
        name += '_' + '_'.join(arr)
    conf = create_conf(template, d)

    command = 'sbatch {0}--job-name="{1}" {2} {3}'.format(qos_string, name, slurm_script, conf[len(conf_prefix)-1:])
    print('Running: ' + command)
    os.system(command)


run_confs(template_file, d, d_change, df_prefix)

