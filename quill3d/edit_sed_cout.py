import os

file_path = os.path.dirname(os.path.realpath(__file__))
file_path +="/test"

def add_grid_symbol(file_path):
	""" Replaces standart design name with a variable """

	export_file=open(file_path,"r")
	temp_file=open("temp","w")

	file=export_file.readlines()
	length=len(file)

	for i in range(length):
		if file[i].find("$")!=-1:
			if file[i+3].find("$")!=-1 and file[i+1].find("data_folder")==-1 and file[i+1].find("output_mode")==1 and file[i+1].find("type_of_interpolation")==-1 and file[i+1].find("polarization")==-1 and file[i+1].find("e_components_for_output")==-1 and file[i+1].find("b_components_for_output")==-1 and file[i+1].find("particles_for_output")==-1 and file[i+1].find("mwindow")==-1 and file[i+1].find("tr_init")==-1:
				file.insert(i+3,"#\n")
				length+=1

	for i in range(length):
		if file[i].find("$")!=-1:
			if file[i+1].find("data_folder")!=-1:
				file.insert(i+2,"#\n")
			if file[i+1].find("output_mode")!=-1:
				file.insert(i+2,"#\n")
			 #if file[i+1].find("type_of_interpolation")!=-1:
			 #	file.insert(i+2,"#\n")
			# if file[i+1].find("polarization")!=-1:
			# 	file.insert(i+2,"#\n")
			# if file[i+1].find("e_components_for_output")!=-1:
			# 	file.insert(i+2,"#\n")
			# if file[i+1].find("b_components_for_output")!=-1:
			# 	file.insert(i+2,"#\n")
			# if file[i+1].find("particles_for_output")!=-1:
			# 	file.insert(i+2,"#\n")
			# if file[i+1].find("mwindow")!=-1:
			# 	file.insert(i+2,"#\n")
			# if file[i+1].find("tr_init")!=-1:
			# 	file.insert(i+2,"#\n")

	for line in file:
		temp_file.write(line)
		#print(line)

	os.remove(file_path)

	export_file.close()
	temp_file.close()

def main():

	add_grid_symbol(file_path)

if __name__ == "__main__":
    main()