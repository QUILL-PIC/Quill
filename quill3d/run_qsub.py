import os
import re

conf_prefix = '../quill3d-conf/quill.conf.'
qsub_folder = conf_prefix + 'qsub/'
qsub_script = 'run_qsub.sh'

def create_conf(template, d):
    new = []
    name = os.path.join(qsub_folder, template[len(conf_prefix):])
    with open(template, "r") as f:
        for l in f:
            for s in d:
                l = re.sub(r'\$%s\b' % s, str(d[s]), l)
            new.append(l)
    for s in d:
        name += '_' + s + str(d[s])
    print("Printing to %s" % name)
    if not os.path.isdir(qsub_folder):
        print("Creating %s" % qsub_folder)
        os.makedirs(qsub_folder)
    with open(name, "w") as f:
        for l in new:
            f.write(l)
    return name


def run_confs(template, d, d_change, folder='res'):
    if len(d_change) == 0:
        run_conf(template, d, folder)
    else:
        k = sorted(list(d_change.keys()))[0]
        for val in d_change[k]:
            new_folder = folder + '_' + k + str(val)
            d[k] = val
            new_d_change = d_change.copy()
            new_d_change.pop(k)
            run_confs(template, d, new_d_change, new_folder)


def run_conf(template, d, folder=None):
    if folder is not None:
        d['data_folder'] = folder
        arr = folder.split('_')
        if len(arr) > 1:
            arr = arr[1:]
        name = 'quill_' + '_'.join(arr) 
    else:
        name = 'quill'
    conf = create_conf(template, d)
    
    command = 'qsub -N %s %s -F %s' % (name, qsub_script, conf[len(conf_prefix)-1:])
    print("Running:" + command)
    os.system(command)


tmp = conf_prefix + 'bubble-beam-template'
d = {}

d['dt']= 0.05
d['dx'] = 0.11160714285714286
d['dtr'] = 0.15625
d['xnpic'] = 1
d['trnpic'] = 1
d['t_end'] = 150
d['mw_trans'] = 10
d['xlength'] = 35
d['trlength'] = 35
d['t_end'] = 150
d['epsb'] = 20000
d['rc'] = 3
d['Nb'] = '2e10'

d_change = {'rc': [1, 2, 3]}

run_confs(tmp, d, d_change, 'res_test')
