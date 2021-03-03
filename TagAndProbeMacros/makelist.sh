#!/bin/bash
basedir="/eos/cms/store/group/phys_egamma/arun/TriggerEff_Run3/"
if [[ "$#" -lt 1 ]]; 
then
    echo " "
    echo "$basedir"
    echo "====>>  no argument is passed. Pass the remaining path for ntuples"
    echo " for example : ./makelist.sh TrigEff_HWW_EGamma_2018_forpublic/EGamma/Run2018B/190503_221754/0000/"
    exit
fi
vardir=$1
#vardir="TrigEff_HWW_EGamma_2018_forpublic/EGamma/Run2018B/190503_221754/0000/"
#ls -ltrh $basedir$vardir | awk -v a="$basedir$vardir" '{print a,$9}' | sed 's/ //g' > list.txt
ls -ltrh $basedir$vardir | awk -v a="$basedir$vardir" '{print a,$9}' | sed 's/ //g'
