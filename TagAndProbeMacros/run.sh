HistFiles_Ele_PT="histNames_Ele_PT.txt"
HistFiles_Ele_Phi="histNames_Ele_Phi.txt"
HistFiles_Ele_ETA="histNames_Ele_ETA.txt"
HistFiles_Ele_PU="histNames_Ele_PU.txt"

#g++ Run_TnP.cxx -o tnp_PT -std=c++0x `root-config --libs --cflags` -include TagAndProbe_Ele_PT.C
#./tnp_PT listDY2021.txt DY2021
#./tnp_PT listDY2018.txt DY2018

#root -l -b -q "plotEfficiency_Ele_PT.C(\"$HistFiles_Ele_PT\")"
#root -l -b -q "plotEfficiency_Ele_PT_ID.C(\"$HistFiles_Ele_PT\")"
#root -l -b -q "plotEfficiency_multi_Ele.C(\"$HistFiles_Ele_PT\")"


#g++ Run_TnP.cxx -o tnp_ETA -std=c++0x `root-config --libs --cflags` -include TagAndProbe_Ele_ETA.C
#./tnp_ETA listDY2021.txt DY2021
#./tnp_ETA listZprime.txt DY2023
#./tnp_ETA listDY2018.txt DY2018

#root -l -b -q "plotEfficiency_multi_Ele.C(\"$HistFiles_Ele_ETA\")"
#root -l -b -q "plotEfficiency_Ele_ETA.C(\"$HistFiles_Ele_ETA\")"

#g++ Run_TnP.cxx -o tnp_PU -std=c++0x `root-config --libs --cflags` -include TagAndProbe_Ele_PU.C
#./tnp_PU listDY2021.txt DY2021
#./tnp_PU listZprime.txt DY2023
#./tnp_PU listDY2018.txt DY2018


#root -l -b -q "plotEfficiency_multi_Ele.C(\"$HistFiles_Ele_PU\")"
#root -l -b -q "plotEfficiency_Ele_PU.C(\"$HistFiles_Ele_PU\")"

g++ Run_TnP.cxx -o tnp_PHI -std=c++0x `root-config --libs --cflags` -include TagAndProbe_Ele_Phi.C
./tnp_PHI listDY2021.txt DY2021
./tnp_PHI listDY2018.txt DY2018

#root -l -b -q "plotEfficiency_Ele_PT.C(\"$HistFiles_Ele_PT\")"
#root -l -b -q "plotEfficiency_Ele_PT_ID.C(\"$HistFiles_Ele_PT\")"
#root -l -b -q "plotEfficiency_multi_Ele.C(\"$HistFiles_Ele_PT\")"
root -l -b -q "plotEfficiency_multi_Ele.C(\"$HistFiles_Ele_Phi\")"



#mv *.png ~/www/EGMHLT/EleTriggers/Winter20
