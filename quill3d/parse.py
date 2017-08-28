import os
import re

file_path = os.path.dirname(os.path.realpath(__file__))
file_path +="/parse_temp"

def parse_fun(file_path):

	export_file=open(file_path,"r")
	temp_file=open("temp","w")

	for line in export_file:
		n=line.find('#') #delete comments
		line_new=line[:n]
		result1=re.split(r'=',line_new) #first parse line
		if result1[0]!='':
			st=result1[1]
			result2=re.split(r' ',st) #parse expression after "="
			print(result1[0])
			if result2[1]!='' and result2[1].isalpha(): #chech element is string or number
				print('#')
				print(result2[1])
				print('$')
			else:
				print(result2[1])
				if len(result2)>2:
					if result2[2]!='':
						print(result2[2])
						print('$')
					else:
						print('#')
						print('$')
				else:
					print('#')
					print('$')

	print('$')

	export_file.close()
	temp_file.close()

def main():
	parse_fun(file_path)

if __name__ == "__main__":
    main()