#include "RecoBTag/Analysis/interface/TrackIPTagPlotter.h"

TrackIPTagPlotter::TrackIPTagPlotter(const TString & tagName,
	const EtaPtBin & etaPtBin, const edm::ParameterSet& pSet, bool update) :
	BaseTagInfoPlotter(tagName, etaPtBin)
{
  nBinEffPur_  = pSet.getParameter<int>("nBinEffPur");
  startEffPur_ = pSet.getParameter<double>("startEffPur");
  endEffPur_   = pSet.getParameter<double>("endEffPur");

  finalized = false;
  if (update){
  TString dir= "TrackIP"+theExtensionString;
  gFile->cd();
  gFile->cd(dir);
  }

  trkNbr3D = new FlavourHistorgrams<int>
	("selTrksNbr_3D" + theExtensionString, "Number of selected tracks for 3D IPS" + theExtensionString, 31, -0.5, 30.5,
	false, true, true, "b", update);

  trkNbr2D = new FlavourHistorgrams<int>
	("selTrksNbr_2D" + theExtensionString, "Number of selected tracks for 2D IPS" + theExtensionString, 31, -0.5, 30.5,
	false, true, true, "b", update);

  lowerIPSBound = -35.0;
  upperIPSBound = +35.0;

  tkcntHistosSig3D[4] = new FlavourHistorgrams<double>
       ("ips_3D" + theExtensionString, "3D Significance of impact parameter",
	50, lowerIPSBound, upperIPSBound, false, true, true, "b", update) ;

  tkcntHistosSig3D[0] = new FlavourHistorgrams<double>
       ("ips1_3D" + theExtensionString, "3D Significance of impact parameter 1st trk",
	50, lowerIPSBound, upperIPSBound, false, true, true, "b", update) ;

  tkcntHistosSig3D[1] = new FlavourHistorgrams<double>
       ("ips2_3D" + theExtensionString, "3D Significance of impact parameter 2nd trk",
	50, lowerIPSBound, upperIPSBound, false, true, true, "b", update) ;

  tkcntHistosSig3D[2] = new FlavourHistorgrams<double>
       ("ips3_3D" + theExtensionString, "3D Significance of impact parameter 3rd trk",
	50, lowerIPSBound, upperIPSBound, false, true, true, "b", update) ;

  tkcntHistosSig3D[3] = new FlavourHistorgrams<double>
       ("ips4_3D" + theExtensionString, "3D Significance of impact parameter 4th trk",
	50, lowerIPSBound, upperIPSBound, false, true, true, "b", update) ;

  tkcntHistosSig2D[4] = new FlavourHistorgrams<double>
       ("ips_2D" + theExtensionString, "2D Significance of impact parameter",
	50, lowerIPSBound, upperIPSBound, false, true, true, "b", update) ;

  tkcntHistosSig2D[0] = new FlavourHistorgrams<double>
       ("ips1_2D" + theExtensionString, "2D Significance of impact parameter 1st trk",
	50, lowerIPSBound, upperIPSBound, false, true, true, "b", update) ;

  tkcntHistosSig2D[1] = new FlavourHistorgrams<double>
       ("ips2_2D" + theExtensionString, "2D Significance of impact parameter 2nd trk",
	50, lowerIPSBound, upperIPSBound, false, true, true, "b", update) ;

  tkcntHistosSig2D[2] = new FlavourHistorgrams<double>
       ("ips3_2D" + theExtensionString, "2D Significance of impact parameter 3rd trk",
	50, lowerIPSBound, upperIPSBound, false, true, true, "b", update) ;

  tkcntHistosSig2D[3] = new FlavourHistorgrams<double>
       ("ips4" + theExtensionString, "2D Significance of impact parameter 4th trk",
	50, lowerIPSBound, upperIPSBound, false, true, true, "b", update) ;

  tkcntHistosProb3D[4] = new FlavourHistorgrams<float>
       ("prob_3D" + theExtensionString, "3D Probability of impact parameter",
	50, -1.1, 1.1, false, true, true, "b", update) ;

  tkcntHistosProb3D[0] = new FlavourHistorgrams<float>
       ("prob1_3D" + theExtensionString, "3D Probability of impact parameter 1st trk",
	50, -1.1, 1.1, false, true, true, "b", update) ;

  tkcntHistosProb3D[1] = new FlavourHistorgrams<float>
       ("prob2_3D" + theExtensionString, "3D Probability of impact parameter 2nd trk",
	50, -1.1, 1.1, false, true, true, "b", update) ;

  tkcntHistosProb3D[2] = new FlavourHistorgrams<float>
       ("prob3_3D" + theExtensionString, "3D Probability of impact parameter 3rd trk",
	50, -1.1, 1.1, false, true, true, "b", update) ;

  tkcntHistosProb3D[3] = new FlavourHistorgrams<float>
       ("prob4_3D" + theExtensionString, "3D Probability of impact parameter 4th trk",
	50, -1.1, 1.1, false, true, true, "b", update) ;

  tkcntHistosProb2D[4] = new FlavourHistorgrams<float>
       ("prob_2D" + theExtensionString, "2D Probability of impact parameter",
	50, -1.1, 1.1, false, true, true, "b", update) ;

  tkcntHistosProb2D[0] = new FlavourHistorgrams<float>
       ("prob1_2D" + theExtensionString, "2D Probability of impact parameter 1st trk",
	50, -1.1, 1.1, false, true, true, "b", update) ;

  tkcntHistosProb2D[1] = new FlavourHistorgrams<float>
       ("prob2_2D" + theExtensionString, "2D Probability of impact parameter 2nd trk",
	50, -1.1, 1.1, false, true, true, "b", update) ;

  tkcntHistosProb2D[2] = new FlavourHistorgrams<float>
       ("prob3_2D" + theExtensionString, "2D Probability of impact parameter 3rd trk",
	50, -1.1, 1.1, false, true, true, "b", update) ;

  tkcntHistosProb2D[3] = new FlavourHistorgrams<float>
       ("prob4" + theExtensionString, "2D Probability of impact parameter 4th trk",
	50, -1.1, 1.1, false, true, true, "b", update) ;

#if 0
  decayLengthValuHisto = new FlavourHistorgrams<double>
       ("decLen" + theExtensionString, "Decay Length",
	50, -5.0, 5.0, false, true, true, "b", update);
  decayLengthSignHisto = new FlavourHistorgrams<double>
       ("decLenSign" + theExtensionString, "Decay Length significance",
	50, -2000.0, 2000.0, false, true, true, "b", update);
#endif
  jetDistanceValuHisto = new FlavourHistorgrams<double>
       ("jetDist" + theExtensionString, "JetDistance",
	50, -0.1, 0.1, false, true, true, "b", update);
  jetDistanceSignHisto = new FlavourHistorgrams<double>
       ("jetDistSign" + theExtensionString, "JetDistance significance",
	50, -100.0, 100.0, false, true, true, "b", update);

}


TrackIPTagPlotter::~TrackIPTagPlotter ()
{

  delete trkNbr3D;
  delete trkNbr2D;
#if 0
  delete decayLengthValuHisto;
  delete decayLengthSignHisto;
#endif
  delete jetDistanceValuHisto;
  delete jetDistanceSignHisto;

  for(int n=0; n <= 4; n++) {
    delete tkcntHistosSig2D[n];
    delete tkcntHistosSig3D[n];
    delete tkcntHistosProb2D[n];
    delete tkcntHistosProb3D[n];
  }
  if (finalized) {
    for(int n=0; n < 4; n++) delete effPurFromHistos[n];
  }
}


void TrackIPTagPlotter::analyzeTag (const reco::BaseTagInfo * baseTagInfo,
	const JetFlavour & jetFlavour)
{

  const reco::TrackIPTagInfo * tagInfo = 
	dynamic_cast<const reco::TrackIPTagInfo *>(baseTagInfo);

  if (!tagInfo) {
    throw cms::Exception("Configuration")
      << "BTagPerformanceAnalyzer: Extended TagInfo attached to selected JetTag not of type TrackIPTagInfo. " << endl
      << "Select a different JetTag collection or change the algorithm parameter to TrackIPTags in the configuration\n";
  }

  int jetFlav = jetFlavour.flavour();

  vector<reco::TrackIPTagInfo::TrackIPData> ip = tagInfo->impactParameterData();

  vector<float> prob2d, prob3d;
  if (tagInfo->hasProbabilities()) {
    prob2d = tagInfo->probabilities(0);	
    prob3d = tagInfo->probabilities(1);	
  }

  trkNbr3D->fill(jetFlav, ip.size());
  trkNbr2D->fill(jetFlav, ip.size());

  for(unsigned int n=0; n < ip.size() && n < 4; n++) {
    tkcntHistosSig2D[n]->fill(jetFlav, ip[n].ip2d.significance());
    tkcntHistosProb2D[n]->fill(jetFlav, prob2d[n]);
  }
  for(unsigned int n=ip.size(); n < 4; n++)
    tkcntHistosSig2D[n]->fill(jetFlav, lowerIPSBound-1.0);

  for(unsigned int n=0; n < ip.size() && n < 4; n++) {
    tkcntHistosSig3D[n]->fill(jetFlav, ip[n].ip3d.significance());
    tkcntHistosProb3D[n]->fill(jetFlav, prob3d[n]);
  }
  for(unsigned int n=ip.size(); n < 4; n++)
    tkcntHistosSig3D[n]->fill(jetFlav, lowerIPSBound-1.0);

  for(unsigned int n=0; n < ip.size(); n++) {
    tkcntHistosSig2D[4]->fill(jetFlav, ip[n].ip2d.significance());
    tkcntHistosProb2D[4]->fill(jetFlav, prob2d[n]);
  }
  for(unsigned int n=0; n < ip.size(); n++) {
    tkcntHistosSig3D[4]->fill(jetFlav, ip[n].ip3d.significance());
    tkcntHistosProb3D[4]->fill(jetFlav, prob3d[n]);
  }
#if 0
  for(unsigned int n=0; n < ip.size(); n++) {
    decayLengthValuHisto->fill(jetFlav, ip[n].decayLen.value());
    decayLengthSignHisto->fill(jetFlav, ip[n].decayLen.significance());
  }
#endif
  for(unsigned int n=0; n < ip.size(); n++) {
    jetDistanceValuHisto->fill(jetFlav, ip[n].distanceToJetAxis);
//    jetDistanceSignHisto->fill(jetFlav, ip[n].distanceToJetAxis.significance());
  }

}

void TrackIPTagPlotter::finalize ()
{
  //
  // final processing:
  // produce the misid. vs. eff histograms
  //
  effPurFromHistos[0] = new EffPurFromHistos (tkcntHistosSig3D[1],
		nBinEffPur_, startEffPur_,
		endEffPur_);
  effPurFromHistos[1] = new EffPurFromHistos (tkcntHistosSig3D[2],
		nBinEffPur_, startEffPur_,
		endEffPur_);
  effPurFromHistos[2] = new EffPurFromHistos (tkcntHistosSig2D[1],
		nBinEffPur_, startEffPur_,
		endEffPur_);
  effPurFromHistos[3] = new EffPurFromHistos (tkcntHistosSig2D[2],
		nBinEffPur_, startEffPur_,
		endEffPur_);
  for(int n=0; n < 4; n++) effPurFromHistos[n]->compute();
  finalized = true;
}

void TrackIPTagPlotter::psPlot(const TString & name)
{
  TString cName = "TrackIPPlots"+ theExtensionString;
  setTDRStyle()->cd();
  TCanvas canvas(cName, "TrackIPPlots"+ theExtensionString, 600, 900);
  canvas.UseCurrentStyle();
  canvas.Divide(2,3);
  canvas.Print(name + cName + ".ps[");

  canvas.cd(1);
  trkNbr3D->plot();
  canvas.cd(2);
  tkcntHistosSig3D[4]->plot();
  for(int n=0; n < 4; n++) {
    canvas.cd(3+n);
    tkcntHistosSig3D[n]->plot();
  }

  canvas.Print(name + cName + ".ps");
  canvas.Clear();
  canvas.Divide(2,3);

  canvas.cd(1);
  trkNbr3D->plot();
  canvas.cd(2);
  tkcntHistosProb3D[4]->plot();
  for(int n=0; n < 4; n++) {
    canvas.cd(3+n);
    tkcntHistosProb3D[n]->plot();
  }

  canvas.Print(name + cName + ".ps");
  canvas.Clear();
  canvas.Divide(2,3);
  canvas.cd(1);
  trkNbr2D->plot();
  canvas.cd(2);
  tkcntHistosSig2D[4]->plot();
  for(int n=0; n < 4; n++) {
    canvas.cd(3+n);
    tkcntHistosSig2D[n]->plot();
  }

  canvas.Print(name + cName + ".ps");
  canvas.Clear();
  canvas.Divide(2,3);
  canvas.cd(1);
  trkNbr2D->plot();
  canvas.cd(2);
  tkcntHistosProb2D[4]->plot();
  for(int n=0; n < 4; n++) {
    canvas.cd(3+n);
    tkcntHistosProb2D[n]->plot();
  }

  if (finalized) {
    for(int n=0; n < 2; n++) {
      canvas.Print(name + cName + ".ps");
      canvas.Clear();
      canvas.Divide(2,3);
      canvas.cd(1);
      effPurFromHistos[0+n]->discriminatorNoCutEffic()->plot();
      canvas.cd(2);
      effPurFromHistos[0+n]->discriminatorCutEfficScan()->plot();
      canvas.cd(3);
      effPurFromHistos[0+n]->plot();
      canvas.cd(4);
      effPurFromHistos[1+n]->discriminatorNoCutEffic()->plot();
      canvas.cd(5);
      effPurFromHistos[1+n]->discriminatorCutEfficScan()->plot();
      canvas.cd(6);
      effPurFromHistos[1+n]->plot();
    }
  }

  canvas.Print(name + cName + ".ps");
  canvas.Clear();
  canvas.Divide(2,3);
  canvas.cd(1);
  jetDistanceValuHisto->plot();
  canvas.cd(2);
  jetDistanceSignHisto->plot();
#if 0
  canvas.cd(3);
  decayLengthValuHisto->plot();
  canvas.cd(4);
  decayLengthSignHisto->plot();
#endif


  canvas.Print(name + cName + ".ps");
  canvas.Print(name + cName + ".ps]");
}

void TrackIPTagPlotter::write()
{
  TString dir= "TrackIP"+theExtensionString;
  gFile->cd();
  gFile->mkdir(dir);
  gFile->cd(dir);
  trkNbr2D->write();
  trkNbr3D->write();
#if 0
  decayLengthValuHisto->write();
  decayLengthSignHisto->write();
#endif
  jetDistanceValuHisto->write();
  jetDistanceSignHisto->write();
  for(int n=0; n <= 4; n++) {
    tkcntHistosSig2D[n]->write();
    tkcntHistosSig3D[n]->write();
    tkcntHistosProb2D[n]->write();
    tkcntHistosProb3D[n]->write();
  }
  if (finalized) {
    for(int n=0; n < 4; n++) effPurFromHistos[n]->write();
  }
  gFile->cd();
}

void TrackIPTagPlotter::epsPlot(const TString & name)
{
  trkNbr2D->epsPlot(name);
  trkNbr3D->epsPlot(name);
  decayLengthValuHisto->epsPlot(name);
  decayLengthSignHisto->epsPlot(name);
  jetDistanceValuHisto->epsPlot(name);
  jetDistanceSignHisto->epsPlot(name);
  for(int n=0; n <= 4; n++) {
    tkcntHistosSig2D[n]->epsPlot(name);
    tkcntHistosSig3D[n]->epsPlot(name);
    tkcntHistosProb2D[n]->epsPlot(name);
    tkcntHistosProb3D[n]->epsPlot(name);
  }
  if (finalized) {
    for(int n=0; n < 4; n++) effPurFromHistos[n]->epsPlot(name);
  }
}
