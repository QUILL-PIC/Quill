#!/usr/bin/env python
import sys
import re

prefix = '../quill3d-conf/quill.conf'

# regex for numbers
number_regex = r'([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
# regex for array of numbers in format '[1;2;3]'
array_regex = r'(\[' + number_regex + r'(;' + number_regex + r')*\])'
# value can be either a number or an array
value_regex = '(' + number_regex + '|' + array_regex + ')'
# dimension can contain word symbols or percent sign, as well as power e.g. "^{-3}"
dimension_regex = r'(\w+|%)(\^{?[+-]?\d+}?)?'
# value with the dimension are separated by space
value_dimension_regex = value_regex + '\s+' + dimension_regex
# text only values can contain word symbols, digits, dots, slashes, ampersands and hyphens
text_only_regex = r'([\w\d&-./]+)'

conf = prefix + sys.argv[1]
with open(conf, "r") as f:
    for l in f:
        if '#' in l:
            l = l[:l.find('#')]
        l = l.strip()
        if l == '':
            continue
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
