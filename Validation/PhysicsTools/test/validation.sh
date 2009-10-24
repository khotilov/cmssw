#!/bin/bash
rm -f plot.C
rm -f -r test
rm -f -r  JETS 
rm -f -r  MUONS
rm -f -r  ELECTRONS
rm -f -r  PHOTONS
rm -f -r  TAUS
rm -f -r  CALOMET
rm -f -r  TCMET
rm -f *.png *.eps
################################
mkdir JETS
mkdir MUONS
mkdir ELECTRONS
mkdir PHOTONS
mkdir TAUS
mkdir CALOMET
mkdir TCMET
##################################
cp captions_jets.jets captions.txt
root -l -b -q plot_jets.C
cp plot_jets.C plot.C
../Tools/indexGen.py
python merge_jetAK5.py
cp test/merged.html test/index.html
rm -f test/merged.html
cp test/* JETS
rm -f -r test
rm -f *.png *.eps
#############
cp captions_muons.muons captions.txt
root -l -b -q plot_muons.C
cp plot_muons.C plot.C
../Tools/indexGen.py
python merge_muon.py
cp test/merged.html test/index.html
rm -f test/merged.html
cp test/* MUONS
rm -f -r test
rm -f *.png *.eps
#############
cp captions_electrons.electrons captions.txt
root -l -b -q plot_electrons.C
cp plot_electrons.C plot.C
../Tools/indexGen.py
python merge_ele.py
cp test/merged.html test/index.html
rm -f test/merged.html
cp test/* ELECTRONS
rm -f -r test
rm -f *.png *.eps
#############
cp captions_photons.photons captions.txt
root -l -b -q plot_photons.C
cp plot_photons.C plot.C
../Tools/indexGen.py
python merge_pho.py
cp test/merged.html test/index.html
rm -f test/merged.html
cp plot_photons.C plot.C
cp test/* PHOTONS
rm -f -r test
rm -f *.png *.eps
############
cp captions_taus.taus captions.txt
root -l -b -q plot_taus.C
cp plot_taus.C plot.C
../Tools/indexGen.py
python merge_tau.py
cp test/merged.html test/index.html
rm -f test/merged.html
cp test/* TAUS
rm -f -r test
rm -f *.png *.eps
############
cp captions_calomets.calomets captions.txt
root -l -b -q plot_calomets.C
cp plot_calomets.C plot.C
../Tools/indexGen.py
python merge_calomet.py
cp test/merged.html test/index.html
rm -f test/merged.html
cp test/* CALOMET
rm -f -r test
rm -f *.png *.eps
##############
cp captions_tcmets.tcmets  captions.txt
root -l -b -q plot_tcmets.C
cp plot_tcmets.C plot.C
../Tools/indexGen.py
python merge_tcmet.py
cp test/merged.html test/index.html
rm -f test/merged.html
cp test/* TCMET
rm -f -r test
rm -f *.png *.eps
##############


#cp JETS/* /afs/fnal.gov/files/home/room3/malik/public_html/PATVALIDATION/PATVALIDATION_Lx345_Lx432_Lx64_432/JETS/IC5/
#cp MUONS/* /afs/fnal.gov/files/home/room3/malik/public_html/PATVALIDATION/PATVALIDATION_Lx345_Lx432_Lx64_432/MUONS/
#cp ELECTRONS/* /afs/fnal.gov/files/home/room3/malik/public_html/PATVALIDATION/PATVALIDATION_Lx345_Lx432_Lx64_432/ELECTRONS/
#cp PHOTONS/* /afs/fnal.gov/files/home/room3/malik/public_html/PATVALIDATION/PATVALIDATION_Lx345_Lx432_Lx64_432/PHOTONS/
#cp TAUS/* /afs/fnal.gov/files/home/room3/malik/public_html/PATVALIDATION/PATVALIDATION_Lx345_Lx432_Lx64_432/TAUS/
#cp TAUS/* /afs/fnal.gov/files/home/room3/malik/public_html/PATVALIDATION/PATVALIDATION_Lx345_Lx432_Lx64_432/MET/

#cp /uscms_data/d2/malik/DONOTDELETE_VALIDATION/PATVALIDATION310pre6/JETS/IC5/index.html /afs/fnal.gov/files/home/room3/malik/public_html/PATVALIDATION/PATVALIDATION_Lx345_Lx432_Lx64_432/JETS/IC5/
#cp /uscms_data/d2/malik/DONOTDELETE_VALIDATION/PATVALIDATION310pre6/MUONS/index.html /afs/fnal.gov/files/home/room3/malik/public_html/PATVALIDATION/PATVALIDATION_Lx345_Lx432_Lx64_432/MUONS/
#cp /uscms_data/d2/malik/DONOTDELETE_VALIDATION/PATVALIDATION310pre6/ELECTRONS/index.html /afs/fnal.gov/files/home/room3/malik/public_html/PATVALIDATION/PATVALIDATION_Lx345_Lx432_Lx64_432/ELECTRONS/
#cp /uscms_data/d2/malik/DONOTDELETE_VALIDATION/PATVALIDATION310pre6/PHOTONS/index.html /afs/fnal.gov/files/home/room3/malik/public_html/PATVALIDATION/PATVALIDATION_Lx345_Lx432_Lx64_432/PHOTONS
#cp /uscms_data/d2/malik/DONOTDELETE_VALIDATION/PATVALIDATION310pre6/TAUS/index.html /afs/fnal.gov/files/home/room3/malik/public_html/PATVALIDATION/PATVALIDATION_Lx345_Lx432_Lx64_432/TAUS
#cp /uscms_data/d2/malik/DONOTDELETE_VALIDATION/PATVALIDATION310pre6/TAUS/index.html /afs/fnal.gov/files/home/room3/malik/public_html/PATVALIDATION/PATVALIDATION_Lx345_Lx432_Lx64_432/MET


