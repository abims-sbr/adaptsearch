#!/bin/bash

TARGET=$1
#TARGET="toolshed"
#TARGET="testtoolshed"
#TARGET="local"


realpath() {
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

ROOTPATH=$(dirname $(realpath "$0"))

cd $ROOTPATH/01_Filter_Assemblies
planemo shed_update -t $TARGET

cd $ROOTPATH/02_Pairwise
planemo shed_update -t $TARGET

cd $ROOTPATH/03_POGs
planemo shed_update -t $TARGET

cd $ROOTPATH/03b_Orthogroups_Tool
planemo shed_update -t $TARGET

cd $ROOTPATH/04_BlastAlign
planemo shed_update -t $TARGET

cd $ROOTPATH/05_CDS_search
planemo shed_update -t $TARGET

cd $ROOTPATH/06_ConcatPhyl
planemo shed_update -t $TARGET

cd $ROOTPATH/suite_adaptsearch
planemo shed_update -t $TARGET



