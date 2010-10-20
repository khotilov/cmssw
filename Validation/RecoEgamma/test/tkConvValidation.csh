#!/bin/csh

#This script can be used to generate a web page to compare histograms from 
#two input root files produced using the EDAnalyzers in RecoEgamma/Examples,
#by running one of:
#
#  
#  
#  "Validation/RecoEgamma/test/PhotonValidator_cfg.py
#
# The default list of histograms (configurable) is based on version VXX-XX-XX
# of Validation/RecoEgamma
#
#Two files are created by this script: validation.C and validation.html.
#validation.C should be run inside root to greate a set of gif images
#which can then be viewed in a web browser using validation.html.

#=============BEGIN CONFIGURATION=================
setenv TYPE TrackBasedConversions
setenv CMSSWver1 3_8_5
setenv CMSSWver2 3_9_0
setenv OLDRELEASE 385
setenv NEWRELEASE 390
setenv OLDPRERELEASE 
setenv NEWPRERELEASE 


setenv OLDRELEASE ${OLDRELEASE}${OLDPRERELEASE}
setenv NEWRELEASE ${NEWRELEASE}${NEWPRERELEASE}


#setenv WorkDir1   /afs/cern.ch/user/n/nancy/scratch0/CMSSW/test/CMSSW_${CMSSWver1}/src/Validation/RecoEgamma/test
#setenv WorkDir2   /afs/cern.ch/user/n/nancy/scratch0/CMSSW/test/CMSSW_${CMSSWver2}_${NEWPRERELEASE}/src/Validation/RecoEgamma/test



#setenv WorkDir1   /afs/cern.ch/user/n/nancy/scratch0/CMSSW/test/CMSSW_${CMSSWver1}_${OLDPRERELEASE}/src/Validation/RecoEgamma/test
#setenv WorkDir2   /afs/cern.ch/user/n/nancy/scratch0/CMSSW/test/CMSSW_${CMSSWver2}_${NEWPRERELEASE}/src/Validation/RecoEgamma/test


#setenv WorkDir1   /afs/cern.ch/user/n/nancy/scratch0/CMSSW/test/CMSSW_${CMSSWver1}_${OLDPRERELEASE}/src/Validation/RecoEgamma/test
#setenv WorkDir2   /afs/cern.ch/user/n/nancy/scratch0/CMSSW/test/CMSSW_${CMSSWver2}/src/Validation/RecoEgamma/test


#setenv WorkDir1   /afs/cern.ch/user/n/nancy/scratch0/CMSSW/test/CMSSW_${CMSSWver1}/src/Validation/RecoEgamma/test
#setenv WorkDir2   /afs/cern.ch/user/n/nancy/scratch0/CMSSW/test/CMSSW_${CMSSWver2}/src/Validation/RecoEgamma/test


setenv WorkDir1    /data/pccmsnd1/b/test/CMSSW_${CMSSWver1}/src/Validation/RecoEgamma/test
#setenv WorkDir2    /afs/cern.ch/user/n/nancy/scratch0/CMSSW/test/CMSSW_${CMSSWver2}_${NEWPRERELEASE}/src/Validation/RecoEgamma/test
setenv WorkDir2    /afs/cern.ch/user/n/nancy/scratch0/CMSSW/test/CMSSW_${CMSSWver2}/src/Validation/RecoEgamma/test


#Name of sample (affects output directory name and htmldescription only) 

setenv SAMPLE QCD_Pt_80_120STARTUP
#setenv SAMPLE QCD_Pt_20_30STARTUP
#TYPE must be one ofPixelMatchGsfElectron, Photon 

#==============END BASIC CONFIGURATION==================


#Input root trees for the two cases to be compared 


if ($SAMPLE == QCD_Pt_20_30STARTUP) then 

setenv OLDFILE ${WorkDir1}/ConversionValidationRelVal${OLDRELEASE}_QCD_Pt_20_30.root
setenv NEWFILE ${WorkDir2}/ConversionValidationRelVal${NEWRELEASE}_QCD_Pt_20_30.root

else if ($SAMPLE == QCD_Pt_80_120STARTUP) then 

setenv OLDFILE ${WorkDir1}/ConversionValidationRelVal${OLDRELEASE}_QCD_Pt_80_120.root
setenv NEWFILE ${WorkDir2}/ConversionValidationRelVal${NEWRELEASE}_QCD_Pt_80_120.root



endif



#Location of output.  The default will put your output in:
#http://cmsdoc.cern.ch/Physics/egamma/www/validation/

setenv CURRENTDIR $PWD
setenv OUTPATH /afs/cern.ch/cms/Physics/egamma/www/validation
cd $OUTPATH
if (! -d $NEWRELEASE) then
  mkdir $NEWRELEASE
endif
setenv OUTPATH $OUTPATH/$NEWRELEASE
cd $OUTPATH

if (! -d ${TYPE}) then
  mkdir ${TYPE}
endif
setenv OUTPATH $OUTPATH/${TYPE}
cd  $OUTPATH

if (! -d vs${OLDRELEASE}) then
  mkdir vs${OLDRELEASE}
endif
setenv OUTPATH $OUTPATH/vs${OLDRELEASE}


setenv OUTDIR $OUTPATH/${SAMPLE}
if (! -d $OUTDIR) then
  cd $OUTPATH
  mkdir $OUTDIR
  cd $OUTDIR
  mkdir gifs
endif
cd $OUTDIR
echo "cd $OUTDIR"

#The list of histograms to be compared for each TYPE can be configured below:


if ( $TYPE == TrackBasedConversions ) then

cat > efficiency <<EOF
  convEffVsEtaTwoTracks
  convEffVsPhiTwoTracks
  convEffVsRTwoTracks
  convEffVsZTwoTracks
  convEffVsEtTwoTracks
  convFakeRateVsEtaTwoTracks
  convFakeRateVsPhiTwoTracks
  convFakeRateVsRTwoTracks
  convFakeRateVsZTwoTracks
  convFakeRateVsEtTwoTracks
EOF



cat > scaledhistos <<EOF
convEta
convPhi
convZ
convVtxdEta
convVtxdPhi
convVtxdR
convVtxdZ
convVtxdX
convVtxdY
convVtxdR_barrel
convVtxdZ_barrel
convVtxdX_barrel
convVtxdY_barrel
convVtxdR_endcap
convVtxdZ_endcap
convVtxdX_endcap
convVtxdY_endcap
hDCotTracksAll
hDCotTracksBarrel
hDCotTracksEndcap
convPtResAll
convPtResBarrel
convPtResEndcap
hDPhiTracksAtVtxAll
hDPhiTracksAtVtxBarrel
hDPhiTracksAtVtxEndcap
hInvMassAll_AllTracks
hInvMassBarrel_AllTracks
hInvMassEndcap_AllTracks
hTkD0All
hTkPtPullAll
hTkPtPullBarrel
hTkPtPullEndcap

EOF

cat > logYScaledHistos <<EOF
convRplot
convPt
hDistMinAppTracksAll
hDistMinAppTracksBarrel
hDistMinAppTracksEndcap
vtxChi2ProbAll
vtxChi2ProbBarrel
vtxChi2ProbEndcap

EOF

cat > projections <<EOF
nHitsVsEtaAllTracks


EOF

cat > unscaledhistos <<EOF

EOF


endif

#=================END CONFIGURATION=====================

if (-e validation.C) rm validation.C
touch validation.C
cat > begin.C <<EOF
{
TFile *file_old = TFile::Open("$OLDFILE");
TFile *file_new = TFile::Open("$NEWFILE");

EOF
cat begin.C >>& validation.C
rm begin.C

setenv N 1






foreach i (`cat scaledhistos`)
  cat > temp$N.C <<EOF
TCanvas *c$i = new TCanvas("c$i","c$i",430, 10, 700, 500);
c$i->SetFillColor(10);
file_new->cd("DQMData/EgammaV/ConversionValidator/ConversionInfo");
Double_t mnew=$i->GetMaximum();
Double_t nnew=$i->GetEntries();
file_old->cd("DQMData/EgammaV/ConversionValidator/ConversionInfo");
Double_t mold=$i->GetMaximum();
Double_t nold=$i->GetEntries();
$i->SetStats(0);
$i->SetMinimum(0.);
if ( $i==convVtxdR || $i==convVtxdX  ||$i==convVtxdY  ) {
$i->SetMaximum (1300);
} else if ( $i==convVtxdR_barrel ||  $i==convVtxdX_barrel ||  $i==convVtxdY_barrel ) {
$i->SetMaximum (1000);
} else if ( $i==convVtxdR_endcap  ||  $i==convVtxdX_endcap ||  $i==convVtxdY_endcap ) {
$i->SetMaximum (800);
} else if ( $i==convVtxdZ_endcap ) {
$i->SetMaximum (400);
} else if ( $i==hDPhiTracksAtVtxAll ) {
$i->SetMaximum (2000);
} else if ( $i==hDPhiTracksAtVtxBarrel ||  $i==hDPhiTracksAtVtxEndcap ) {
$i->SetMaximum (1000);
}

$i->SetLineColor(kPink+8);
$i->SetFillColor(kPink+8);
//$i->SetLineWidth(3);
$i->Draw();
file_new->cd("DQMData/EgammaV/ConversionValidator/ConversionInfo");
Double_t nnew=$i->GetEntries();
$i->SetStats(0);
$i->SetLineColor(kBlack);
$i->SetMarkerColor(kBlack);
$i->SetMarkerStyle(20);
$i->SetMarkerSize(1);
//$i->SetLineWidth(1);
$i->Scale(nold/nnew);
$i->Draw("e1same");
c$i->SaveAs("gifs/$i.gif");

EOF
  setenv N `expr $N + 1`
end




foreach i (`cat logYScaledHistos`)
  cat > temp$N.C <<EOF
TCanvas *c$i = new TCanvas("c$i","c$i",430, 10, 700, 500);
c$i->SetFillColor(10);
c$i->SetLogy(1);
file_new->cd("DQMData/EgammaV/ConversionValidator/ConversionInfo");
Double_t mnew=$i->GetMaximum();
Double_t nnew=$i->GetEntries();
file_old->cd("DQMData/EgammaV/ConversionValidator/ConversionInfo");
Double_t mold=$i->GetMaximum();
Double_t nold=$i->GetEntries();
$i->SetStats(0);
$i->SetMinimum(3.);
//$i->GetXaxis()->SetRangeUser(0.,60.);
$i->SetLineColor(kPink+8);
$i->SetFillColor(kPink+8);
$i->Draw();
file_new->cd("DQMData/EgammaV/ConversionValidator/ConversionInfo");
Double_t nnew=$i->GetEntries();
$i->SetStats(0);
$i->SetLineColor(kBlack);
$i->SetMarkerColor(kBlack);
$i->SetMarkerStyle(20);
$i->SetMarkerSize(1);
if ( $i==convRplot) {
$i->SetMarkerSize(0.08);
}
$i->Scale(nold/nnew);
$i->Draw("e1same");
c$i->SaveAs("gifs/$i.gif");

EOF
  setenv N `expr $N + 1`
end


foreach i (`cat efficiency`)
  cat > temp$N.C <<EOF

TCanvas *c$i = new TCanvas("c$i","c$i",430, 10, 700, 500);
c$i->SetFillColor(10);
file_old->cd("DQMData/EgammaV/ConversionValidator/EfficienciesAndFakeRate");
$i->SetStats(0);
if ( $i==convFakeRateVsRTwoTracks ) {
$i->GetYaxis()->SetRangeUser(0.,1.);
} else {
$i->GetYaxis()->SetRangeUser(0.,0.5);
}
$i->SetLineColor(kPink+8);
$i->SetMarkerColor(kPink+8);
$i->SetMarkerStyle(20);
$i->SetMarkerSize(1);
$i->SetLineWidth(1);
$i->Draw("e1");
file_new->cd("DQMData/EgammaV/ConversionValidator/EfficienciesAndFakeRate");
$i->SetStats(0);
$i->SetMinimum(0.);
$i->SetMaximum(1.1);
$i->SetLineColor(kBlack);
$i->SetMarkerColor(kBlack);
$i->SetMarkerStyle(20);
$i->SetMarkerSize(1);
$i->SetLineWidth(1);
$i->Draw("e1same");
c$i->SaveAs("gifs/$i.gif");

EOF
  setenv N `expr $N + 1`
end


foreach i (`cat projections`)
  cat > temp$N.C <<EOF

TCanvas *c$i = new TCanvas("c$i","c$i",430, 10, 700, 500);
c$i->SetFillColor(10);
file_old->cd("DQMData/EgammaV/ConversionValidator/ConversionInfo");
TH1D *tmp1$i= $i->ProjectionY();
Double_t nold=tmp1$i->GetEntries();
Double_t mold=tmp1$i->GetMaximum();
file_new->cd("DQMData/EgammaV/ConversionValidator/ConversionInfo");
TH1D *tmp2$i= $i->ProjectionY();
Double_t nnew=tmp2$i->GetEntries();
Double_t mnew=tmp2$i->GetMaximum();
tmp1$i->SetStats(0);
tmp1$i->SetMinimum(0.);
tmp1$i->SetLineColor(kPink+8);
tmp1$i->SetFillColor(kPink+8);
tmp1$i->SetLineWidth(3);
tmp1$i->Draw();
tmp2$i->SetStats(0);
tmp2$i->SetLineColor(kBlack);
tmp2$i->SetMarkerColor(kBlack);
tmp2$i->SetMarkerStyle(20);
tmp2$i->SetMarkerSize(1);
tmp2$i->SetLineWidth(1);
tmp2$i->Scale(nold/nnew);
tmp2$i->Draw("e1same");
c$i->SaveAs("gifs/$i.gif");


EOF
  setenv N `expr $N + 1`
end




setenv NTOT `expr $N - 1`
setenv N 1
while ( $N <= $NTOT )
  cat temp$N.C >>& validation.C
  rm temp$N.C
  setenv N `expr $N + 1`
end

cat > end.C <<EOF
}
EOF
cat end.C >>& validation.C
rm end.C


if ( $TYPE == PixelMatchGsfElectron ) then
  setenv ANALYZER PixelMatchGsfElectronAnalyzer
  setenv CFG read_gsfElectrons
else if ( $TYPE == TrackBasedConversions ) then
  setenv ANALYZER ConversionValidator
  setenv CFG ConversionValidator_cfg
endif

if (-e validation.html) rm validation.html
if (-e tkConvValidationPlotsTemplate.html) rm tkConvValidationPlotsTemplate.html
cp ${CURRENTDIR}/tkConvValidationPlotsTemplate.html tkConvValidationPlotsTemplate.html
touch validation.html
cat > begin.html <<EOF
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">
<html>
<head>
<meta http-equiv="content-type" content="text/html; charset=UTF-8" />
<title>$NEWRELEASE vs $OLDRELEASE $TYPE validation</title>
</head>

<h1>$NEWRELEASE vs $OLDRELEASE $TYPE validation
<br>
 $SAMPLE 
</h1>

<p>The following plots were made using <a href="http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/Validation/RecoEgamma/src/$ANALYZER.cc">Validation/RecoEgamma/src/$ANALYZER</a>, 
using <a href="http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/Validation/RecoEgamma/test/$CFG.py">Validation/RecoEgamma/test/$CFG.py</a>
<p>The script used to make the plots is <a href="validation.C">here</a>.
<br>
In all plots below, $OLDRELEASE is in purple , $NEWRELEASE in black. 
<br>
Click on the plots to see them enlarged.
<br>
Responsible: N. Marinelli
<br>
<br>



EOF
cat begin.html >>& validation.html
rm begin.html
cat  tkConvValidationPlotsTemplate.html >>& validation.html
rm   tkConvValidationPlotsTemplate.html 

rm efficiency
rm scaledhistos
rm unscaledhistos
rm logYScaledHistos
rm projections

root -b -l -q validation.C
cd $CURRENTDIR

python makeWebpage.py ${OUTDIR}/gifs  ${OLDRELEASE}  ${NEWRELEASE} ${SAMPLE}
echo "Then you can view your validation plots here:"
echo "http://cmsdoc.cern.ch/Physics/egamma/www/$OUTPATH/validation.html"
