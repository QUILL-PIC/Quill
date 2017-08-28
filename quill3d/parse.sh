#!/bin/bash

cat ../quill3d-conf/quill.conf"$1" | sed -e 's/#\(.*\)$//;s/\s//g;/^$/d;s/$/\n$/' | sed -n 's/\$/&/p;s/\[*\]/&\n/p;t;s/=[+-]\?[0-9]*\(.[0-9]\+\)\?\(e[+-]\?[0-9]\+\)\?/&\n/p' | sed -e 's/^$/#/'| sed -e 's/=[0-9]/\n&/;s/=//;/^$/d' | sed -e 's/^$/#/;$a\$' |sed -e 's/\[/\n&/' | sed -e 's/[+-]/\n&/' >> test
python ./edit_sed_cout.py
cat temp
rm temp

# Содержание строк после '#' удаляется
# Пробелы удаляются
# Пустые строки удаляются
# После каждой строки добавляется строка с символом '$'
# После '=[число]' (в качестве числа может использоваться в том числе
# запись вида +1.2e-5) вставляется переход на новую
# строку
# Вместо '=' вставляется переход на новую строку
# В пустые строки добавляется символ "#"
# В конец потока добавляется символ '$'

# Таким образом, строка вида 'имя = [значение] [размерность]'
# превращается в четыре строки:
# имя
# [значение]
# [размерность]
# $
# Если [значение] и/или [размерность] не указаны, то вместо них
# пишется символ "#"

#This script parses input data stream for Quill's initialization.
#First part deletes all string after #, spacing, empty strings.
#Sed with key '-n' searches strings with $ and prints its, then finds expression in [] or +1.2e-5 (prints one of them)
#Next in place empty stirng paste #
#Thereafter = inserts new line, after = and empty line deletes
#Append $ to the end stream
#Insetrs new line after [
