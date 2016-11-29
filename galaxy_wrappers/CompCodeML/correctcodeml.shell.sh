#! /bin/bash
#creates a codeml.ctl file in each subdirectory containing a fasta or a phylip file
#bash correctcodeml.sh folder/

current=`pwd`
cd $1
wholepath=`pwd`

for i in `find $wholepath -name *.fasta`
do
  path=`dirname $i`
  python /w/galaxy/galaxy4misharl/galaxy-dist/tools/abims/julie/oasearch/CompCodeML/correctcodeml.py $path'/' 0 $1 && echo 'codeml.ctl generated from '$path
done

for i in `find $wholepath -name *.phy`
do
  path=`dirname $i`
  python /w/galaxy/galaxy4misharl/galaxy-dist/tools/abims/julie/oasearch/CompCodeML/correctcodeml.py $path'/' 0 $1 && echo 'codeml.ctl generated from '$path
done

for i in `find $wholepath -name *.dat`
do
  path=`dirname $i`
  python /w/galaxy/galaxy4misharl/galaxy-dist/tools/abims/julie/oasearch/CompCodeML/correctcodeml.py $path'/' 0 $1 && echo 'codeml.ctl generated from '$path
done

#I hate shell
