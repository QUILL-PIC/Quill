import os
import subprocess
FAIL = '\033[91m'
ENDC = '\033[0m'

conf_folder = os.path.join('..', 'quill3d-conf')
prefix = 'quill.conf'
for root, subs, files in os.walk(conf_folder):
    for f in files:
        if f.startswith(prefix) and not os.path.isdir(f):
            file_name = os.path.join(root, f)
            print('Testing: %s' % file_name)
            suffix = f[len(prefix):]
            try:
                s1 = subprocess.check_output(['./parse.sh', suffix])
            except Exception as err:
                print('parse.sh failed: %s' % err)
                continue
            try:
                s2 = subprocess.check_output(['./parse.py', suffix])
            except Exception as err:
                print('parse.py failed: %s' % err)
                continue

            if s1 == s2:
                print('Equal')
            else:
                s1_arr = s1.decode('utf-8').splitlines()
                s2_arr = s2.decode('utf-8').splitlines()
                swap = False
                if len(s1_arr) > len(s2_arr):
                    s1_arr, s2_arr = s2_arr, s1_arr
                    swap = True
                min_len = len(s1_arr)
                max_len = len(s2_arr)
                for i in range(min_len):
                    s1 = s1_arr[i] if not swap else s2_arr[i]
                    s2 = s2_arr[i] if not swap else s1_arr[i]
                    if s1 != s2:
                        print(FAIL + '%s <-> %s' % (s1, s2) + ENDC)
                    else:
                        print(s1)
                for i in range(min_len, max_len):
                    s1 = ''
                    s2 = s2_arr[i]
                    if swap:
                        s1, s2 = s2, s1
                    print(FAIL + '%s <-> %s' % (s1, s2) + ENDC)
