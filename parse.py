#!/usr/bin/env python
import os
import sys
import re


# regex for numbers
number_regex = r'([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
# regex for array of numbers in format '[1;2;3]'
array_regex = r'(\[' + number_regex + r'(;' + number_regex + r')*\])'
# value can be either a number or an array
value_regex = '(' + number_regex + '|' + array_regex + ')'
# basic dimension is word, possibly with power, e.g. cm, cm^{-3}, cm^2.
simple_dimension_regex = r'\w+(\^({[+-]?\d+}|\d+))?'
# dimension can contain multiple simple dimensions divided by / and *, or percent sign
dimension_regex = '({0}([/*]{0})*|%)'.format(simple_dimension_regex)
# value with the dimension are separated by space
value_dimension_regex = value_regex + '\s+' + dimension_regex
# text only values can contain word symbols, digits, dots, slashes, ampersands and hyphens
text_only_regex = r'([\w\d&-./]+)'
# Python expressions that are to be evaluated during preprocessing
eval_regex = r'eval\{(.+)\}'
# Python statements that are to be executed during preprocessing
exec_regex = r'exec\{(.+)\}'

conf = sys.argv[1]

if not os.path.isfile(conf):
    raise Exception('Config file [%s] not found' % conf)

with open(conf, "r") as f:
    buf = ''
    for l in f:
        if '#' in l:
            l = l[:l.find('#')]
        l = l.strip()
        if l == '':
            continue
        l = re.sub(eval_regex, lambda m: str(eval(m.group(1))), l)
        if re.match('^' + exec_regex + '$', l):
            command = re.findall(exec_regex, l)[0]
            exec(command)
        else:
            buf += l + '\n'
    for l in buf.splitlines():
        key_value = l.split('=')
        if len(key_value) != 2:
            raise Exception('Line [%s] is invalid' % l)

        key = key_value[0].strip()
        value = key_value[1].strip()
        print(key)
        if re.match('^' + value_regex + '$', value):
            print(value)
            print('#')
        elif re.match('^' + text_only_regex + '$', value):
            print('#')
            print(value)
        elif re.match('^' + value_dimension_regex + '$', value):
            v = re.findall(value_regex, value)[0][0]
            print(v)
            d = value[len(v):].strip()
            print(d)
        else:
            raise Exception('Syntax of value [%s] of key [%s] is incorrect' % (value, key))
        print('$')
    print('$')
