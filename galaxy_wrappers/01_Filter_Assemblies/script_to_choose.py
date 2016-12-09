import sys, string, os, itertools, re, zipfile

# coding: utf8

i=1
j=i+1
L1 = []
fas = "^.*fas$"
fasta = "^.*fasta$"
script_path = os.path.dirname(sys.argv[0])

##Run scripts
if sys.argv[1]=="velvet" :
    zfile = zipfile.ZipFile(sys.argv[2])
    for name in zfile.namelist() :
	zfile.extract(name)
        suffix=name[:2]
        name2 = "01_%s" %name
        name3="02_%s" %name
        name4="%s%s" %(suffix,name)
        name5="03_%s" % name
        name6="%s" %(suffix)
        name7="%s%s"%(suffix,name)
	#Fasta sequences on one line
        os.system("zcat -f < '%s' | fasta_formatter -w 0 -o '%s'" % (name,name2))
        
        #Format the name of the sequences with good name
        suffix=name[:2]
        #if re.match(fas, name) :
	#os.system("python /w/galaxy/galaxy4misharl/galaxy-dist/tools/abims/julie/oasearch/filter_assembly/remove_redondancy_from_oases_output_v3.1.py %s %s" %(name2, name3))
	os.system("python %s/remove_redondancy_from_oases_output_v3.1.py %s %s" %(script_path,name2, name3))
        #elif re.match(fasta, name) :
        os.system("sed -e 's/Locus_/%s/g' -e 's/_Confidence_/_/g' -e 's/_Transcript_/_/g' -e 's/_Length_/_/g' %s > %s" % (suffix,name3,name4))
        #os.system("python /w/galaxy/galaxy4misharl/galaxy-dist/tools/abims/julie/oasearch/filter_assembly/remove_redondancy_from_oases_output_v3.0.py %s %s" %(name, name2))
        #Pierre guillaume find_orf script for keeping the longuest ORF
        os.system("python %s/find_orf.py %s %s" %(script_path,name4,name5))

	#Apply cap3
        os.system("cap3  %s -p %s -o %s"%(name5,sys.argv[4],sys.argv[5]))
        #Fasta sequences on one line
        #Il faudrait faire un merge des singlets et contigs! TODO
        os.system("zcat -f < '%s.cap.singlets' | fasta_formatter -w 0 -o '%s'" % (name5,name6))
        #Apply pgbrun script filter script TODO length parameter
        os.system("python %s/filter.py %s %s %s" %(script_path,name6,sys.argv[3],name7))  
        L1.append(name7)
    f = zipfile.ZipFile("sequences_filtered.zip", "w")
    for files in L1 :
	f.write("./%s" %files)

##For trinity files
else:
    zfile = zipfile.ZipFile(sys.argv[2])
    for name in zfile.namelist() :
        zfile.extract(name)
        suffix=name[:2]
        name2 = "01%s" %name  
        name3="02%s" %name
        name4="03%s" %name
        name5="04%s" %name
        name6="05%s"% name
        name7="%s"%(suffix)
        name8="%s%s"%(suffix,name)
        #Fasta sequences on one line
        os.system("zcat -f < '%s' | fasta_formatter -w 0 -o '%s'" % (name,name2))
        #replace white space by "_" 
        os.system("sed -e 's/ /_/g' %s > %s " % (name2,name3))
        #Format the name of the sequences with good name
        os.system("python %s/02_format_fasta_name_TRINITY.py %s %s %s" %(script_path,name3, name4,suffix))
        #Apply first script to avoid reductant sequences
        os.system("python %s/01_Choose_One_variants_per_locus_TRINITY_v1.0.py %s %s" %(script_path,name4, name5))
	#Pierre guillaume find_orf script for keeping the longuest ORF
        os.system("python %s/find_orf.py %s %s" %(script_path,name5,name6))
        #Apply cap3
        os.system("cap3  %s -p %s -o %s"%(name6,sys.argv[4],sys.argv[5]))
        #Fasta sequences on one line
        #Il faudrait faire un merge des singlets et contigs! TODO
        os.system("zcat -f < '%s.cap.singlets' | fasta_formatter -w 0 -o '%s'" % (name6,name7))
        #Apply pgbrun script filter script TODO length parameter
        os.system("python %s/filter.py %s %s %s" %(script_path,name7,sys.argv[3],name8)) 
        L1.append(name8)
    f = zipfile.ZipFile("sequences_filtered.zip", "w")
    for files in L1 :
        f.write("./%s" %files)
