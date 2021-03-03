#!/bin/bash

#now="$(date)"
#now="$(date+'%d_%m_%Y')"
#echo $now
#AllFiles="AllRootFiles.txt"
HistFiles_Ele="histNames_Ele.txt"
root -l -b -q "plotEfficiency_multi_Ele.C(\"$HistFiles_Ele\")"  
#rm -rf TriggerResults_2018
#dirname="TriggerResults_2018"+"$now"
#mkdir "$dirname"
#mkdir TriggerResults_2018
#mv results_Ele TriggerResults_2018
#mv results_Mu12 TriggerResults_2018

