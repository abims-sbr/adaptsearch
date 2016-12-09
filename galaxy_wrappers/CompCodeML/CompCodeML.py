#!/usr/bin/python
#Build directories (depending on the model 0, 1, 2) pour CODEML tool
#Written by mmonsoor ABIMS

import string, os, sys,shutil
script_path = os.path.dirname(sys.argv[0])

tree_file=sys.argv[1]
model=sys.argv[2]
concat_nuc_file=sys.argv[3]
model_list=["0","1","2","all"]
#Build the directories for codeml
if model in model_list:
	if model=="0":
		print("Creating directories for model 0")
		root_model0="model0"
		os.mkdir("model0")
		folders_model0=['ns0','ns2','ns2/wfixed','ns2/wfree']
		for folder in folders_model0:
			os.mkdir(os.path.join(root_model0,folder))
			if folder!="ns2":
				path_tree="%s/tree.txt"%(os.path.join(root_model0,folder))
				print(path_tree)
				#shutil.copy2(sys.argv[1],os.path.join(root_model0,folder))
				shutil.copy2(sys.argv[1],path_tree)
				#Copy the concat_nuc fasta file in each directory
				path_concat="%s/concat.fasta"%(os.path.join(root_model0,folder))
				print(path_concat)
				#shutil.copy2(sys.argv[3],os.path.join(root_model0,folder))
				shutil.copy2(sys.argv[3],path_concat)

		#Creates the ctl files for the codeml script
		os.system("bash %s/correctcodeml.shell.sh %s/" %(script_path,os.path.abspath(root_model0)))
		print (script_path)
		#Executes codeml
		for folder in folders_model0:		
			if folder!="ns2":
				print("codeml %s/codeml.ctl"%( os.path.abspath(os.path.join(root_model0,folder)))    )
				os.system("cd %s && codeml %s/codeml.ctl"%(os.path.abspath(os.path.join(root_model0,folder)),os.path.abspath(os.path.join(root_model0,folder))))
	elif model=="1":
		print("Creating directories for model 1")
		root_model1="model1"
		os.mkdir("model1")
		folders_model1=['ns0']
		for folder in folders_model1:
			os.mkdir(os.path.join(root_model1,folder))
			shutil.copy2(sys.argv[1],os.path.join(root_model1,folder))
			#Copy the concat_nuc fasta file in each directory
			shutil.copy2(sys.argv[3],os.path.join(root_model1,folder))
		#Creates the ctl files for the codeml script
		os.system("bash %s/correctcodeml.shell.sh %s/" %(script_path,root_model1))
		
	elif model=="2":
		print("Creating directories for model 2")
		root_model2="model2"
		os.mkdir("model2")
		folders_model2=['ns2','ns2/wfixed','ns2/wfree']
		for folder in folders_model2:
			os.mkdir(os.path.join(root_model2,folder))
			if folder!="ns2":
				shutil.copy2(sys.argv[1],os.path.join(root_model2,folder))
				#Copy the concat_nuc fasta file in each directory
				shutil.copy2(sys.argv[3],os.path.join(root_model2,folder))
		#Creates the ctl files for the codeml script
		os.system("bash correctcodeml.shell.sh %s/" %(root_model2))
	elif model=="all":
		print("Creating directories for all models")
		print("Creating directories for model 0")
		root_model0="model0"
		os.mkdir("model0")
		folders_model0=['ns0','ns2','ns2/wfixed','ns2/wfree']
		for folder in folders_model0:
			os.mkdir(os.path.join(root_model0,folder))
			if folder!="ns2":
				path_tree="%s/tree.txt"%(os.path.join(root_model0,folder))
				print(path_tree)
				#shutil.copy2(sys.argv[1],os.path.join(root_model0,folder))
				shutil.copy2(sys.argv[1],path_tree)
				#Copy the concat_nuc fasta file in each directory
				path_concat="%s/concat.fasta"%(os.path.join(root_model0,folder))
				print(path_concat)
				#shutil.copy2(sys.argv[3],os.path.join(root_model0,folder))
				shutil.copy2(sys.argv[3],path_concat)
		#Creates the ctl files for the codeml script
		os.system("bash %s/correctcodeml.shell.sh %s/" %(script_path,os.path.abspath(root_model0)   ))
		#Executes codeml
		for folder in folders_model0:		
			if folder!="ns2":
				print("codeml %s/codeml.ctl"%( os.path.abspath(os.path.join(root_model0,folder)))    )
				os.system("cd %s && codeml %s/codeml.ctl"%(os.path.abspath(os.path.join(root_model0,folder)),os.path.abspath(os.path.join(root_model0,folder))))

		print("Creating directories for model 1")
		root_model1="model1"
		os.mkdir("model1")
		folders_model1=['ns0']
		for folder in folders_model1:
			os.mkdir(os.path.join(root_model1,folder))
			path_tree="%s/tree.txt"%(os.path.join(root_model1,folder))
			shutil.copy2(sys.argv[1],path_tree)
			#Copy the concat_nuc fasta file in each directory
			path_concat="%s/concat.fasta"%(os.path.join(root_model1,folder))
			shutil.copy2(sys.argv[3],path_concat)
		#Creates the ctl files for the codeml script
		os.system("bash %s/correctcodeml.shell.sh %s/" %( script_path,os.path.abspath(root_model1)   ))
		#Executes codeml
		for folder in folders_model1:		
			print("codeml %s/codeml.ctl"%( os.path.abspath(os.path.join(root_model1,folder)))    )
			os.system("cd %s && codeml %s/codeml.ctl"%(os.path.abspath(os.path.join(root_model1,folder)),os.path.abspath(os.path.join(root_model1,folder))))


else:
	stop("Error, you must choose a model between 0, 1 or 2")

#Zip outputs files
os.system("zip -r output.zip model0/ model1/")




