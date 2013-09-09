#! /bin/bash

touch number && touch "$1".asy
cat ../"$1".asy | sed -e 's/real\s\+file_number.*$/real file_number = input(\"number\");/;s/string\s\+results_folder.*$/string results_folder = \"..\/..\/results\/\";/' > "$1".asy
