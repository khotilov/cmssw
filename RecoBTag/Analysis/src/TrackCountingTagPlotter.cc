#include "RecoBTag/Analysis/interface/TrackCountingTagPlotter.h"

TrackCountingTagPlotter::TrackCountingTagPlotter(JetTagPlotter *jetTagPlotter,
	bool update) :
    BaseBTagPlotter(jetTagPlotter->etaPtBin()), jetTagPlotter_(jetTagPlotter)
{
  finalized = false;
  if (update){
  TString dir= "TrackCounting"+theExtensionString;
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

  tkcntHistosSig3D[4] = new FlavourHistorgrams<double>
       ("ips_3D" + theExtensionString, "3D Significance of impact parameter",
	50, lowerIPSBound, 35.0, false, true, true, "b", update) ;

  tkcntHistosSig3D[0] = new FlavourHistorgrams<double>
       ("ips1_3D" + theExtensionString, "3D Significance of impact parameter 1st trk",
	50, lowerIPSBound, 35.0, false, true, true, "b", update) ;

  tkcntHistosSig3D[1] = new FlavourHistorgrams<double>
       ("ips2_3D" + theExtensionString, "3D Significance of impact parameter 2nd trk",
	50, lowerIPSBound, 35.0, false, true, true, "b", update) ;

  tkcntHistosSig3D[2] = new FlavourHistorgrams<double>
       ("ips3_3D" + theExtensionString, "3D Significance of impact parameter 3rd trk",
	50, lowerIPSBound, 35.0, false, true, true, "b", update) ;

  tkcntHistosSig3D[3] = new FlavourHistorgrams<double>
       ("ips4_3D" + theExtensionString, "3D Significance of impact parameter 4th trk",
	50, lowerIPSBound, 35.0, false, true, true, "b", update) ;

  tkcntHistosSig2D[4] = new FlavourHistorgrams<double>
       ("ips_2D" + theExtensionString, "2D Significance of impact parameter",
	50, lowerIPSBound, 35.0, false, true, true, "b", update) ;

  tkcntHistosSig2D[0] = new FlavourHistorgrams<double>
       ("ips1_2D" + theExtensionString, "2D Significance of impact parameter 1st trk",
	50, lowerIPSBound, 35.0, false, true, true, "b", update) ;

  tkcntHistosSig2D[1] = new FlavourHistorgrams<double>
       ("ips2_2D" + theExtensionString, "2D Significance of impact parameter 2nd trk",
	50, lowerIPSBound, 35.0, false, true, true, "b", update) ;

  tkcntHistosSig2D[2] = new FlavourHistorgrams<double>
       ("ips3_2D" + theExtensionString, "2D Significance of impact parameter 3rd trk",
	50, lowerIPSBound, 35.0, false, true, true, "b", update) ;

  tkcntHistosSig2D[3] = new FlavourHistorgrams<double>
       ("ips4" + theExtensionString, "2D Significance of impact parameter 4th trk",
	50, lowerIPSBound, 35.0, false, true, true, "b", update) ;

}


TrackCountingTagPlotter::~TrackCountingTagPlotter ()
{

  delete jetTagPlotter_;
  delete trkNbr3D;
  delete trkNbr2D;

  for(int n=0; n <= 4; n++) {
    delete tkcntHistosSig2D[n];
    delete tkcntHistosSig3D[n];
  }
  if (finalized) {
    for(int n=0; n < 4; n++) delete effPurFromHistos[n];
  }
}


void TrackCountingTagPlotter::analyzeTag (const reco::TrackCountingTagInfo & tagInfo,
	const reco::JetTag & jetTag, const JetFlavour & jetFlavour)
{

  int jetFlav = jetFlavour.flavour();

  jetTagPlotter_->analyzeJetTag(jetTag, jetFlavour);

  trkNbr3D->fill(jetFlav, tagInfo.selectedTracks(0));
  trkNbr2D->fill(jetFlav, tagInfo.selectedTracks(1));

  for(int n=0; n < tagInfo.selectedTracks(1) && n < 4; n++)
    tkcntHistosSig2D[n]->fill(jetFlav, tagInfo.significance(n,1));
  for(int n=tagInfo.selectedTracks(1); n < 4; n++)
    tkcntHistosSig2D[n]->fill(jetFlav, lowerIPSBound-1.0);

  for(int n=0; n < tagInfo.selectedTracks(0) && n < 4; n++)
    tkcntHistosSig3D[n]->fill(jetFlav, tagInfo.significance(n,0));
  for(int n=tagInfo.selectedTracks(0); n < 4; n++)
    tkcntHistosSig3D[n]->fill(jetFlav, lowerIPSBound-1.0);

  for(int n=0; n < tagInfo.selectedTracks(1); n++)
    tkcntHistosSig2D[4]->fill(jetFlav, tagInfo.significance(n,1));
  for(int n=0; n < tagInfo.selectedTracks(0); n++)
    tkcntHistosSig3D[4]->fill(jetFlav, tagInfo.significance(n,0));
}

void TrackCountingTagPlotter::finalize ()
{
  jetTagPlotter_->finalize();
  //
  // final processing:
  // produce the misid. vs. eff histograms
  //
  effPurFromHistos[0] = new EffPurFromHistos (tkcntHistosSig3D[1],
		jetTagPlotter_->nBinEffPur(), jetTagPlotter_->startEffPur(),
		jetTagPlotter_->endEffPur());
  effPurFromHistos[1] = new EffPurFromHistos (tkcntHistosSig3D[2],
		jetTagPlotter_->nBinEffPur(), jetTagPlotter_->startEffPur(),
		jetTagPlotter_->endEffPur());
  effPurFromHistos[2] = new EffPurFromHistos (tkcntHistosSig2D[1],
		jetTagPlotter_->nBinEffPur(), jetTagPlotter_->startEffPur(),
		jetTagPlotter_->endEffPur());
  effPurFromHistos[3] = new EffPurFromHistos (tkcntHistosSig2D[2],
		jetTagPlotter_->nBinEffPur(), jetTagPlotter_->startEffPur(),
		jetTagPlotter_->endEffPur());
  for(int n=0; n < 4; n++) effPurFromHistos[n]->compute();
  finalized = true;
}

void TrackCountingTagPlotter::psPlot(const TString & name)
{
  jetTagPlotter_->psPlot(name);

  TString cName = "TrackCountingPlots"+ theExtensionString;
  setTDRStyle()->cd();
  TCanvas canvas(cName, "TrackCountingPlots"+ theExtensionString, 600, 900);
  canvas.UseCurrentStyle();
  canvas.Divide(2,3);
  canvas.Print(name + cName + ".ps[");

  canvas.cd(1);
  trkNbr3D->plot((TPad*) canvas.GetPrimitive(cName+"_1"));
  canvas.cd(2);
  tkcntHistosSig3D[4]->plot((TPad*) canvas.GetPrimitive(cName+"_2"));
  for(int n=0; n < 4; n++) {
    canvas.cd(3+n);
    tkcntHistosSig3D[n]->plot((TPad*) canvas.GetPrimitive(cName+"_"+itos(n+3)));
  }

  canvas.Print(name + cName + ".ps");
  canvas.Clear();
  canvas.Divide(2,3);

  canvas.cd(1);
  trkNbr2D->plot((TPad*) canvas.GetPrimitive(cName+"_1"));
  canvas.cd(2);
  tkcntHistosSig2D[4]->plot((TPad*) canvas.GetPrimitive(cName+"_2"));
  for(int n=0; n < 4; n++) {
    canvas.cd(2+n);
    tkcntHistosSig2D[n]->plot((TPad*) canvas.GetPrimitive(cName+"_"+itos(n+3)));
  }

  if (finalized) {
    for(int n=0; n < 2; n++) {
      canvas.Print(name + cName + ".ps");
      canvas.Clear();
      canvas.Divide(2,3);
      canvas.cd(1);
      effPurFromHistos[0+n]->discriminatorNoCutEffic()->plot((TPad*) canvas.GetPrimitive(cName+"_1"));
      canvas.cd(2);
      effPurFromHistos[0+n]->discriminatorCutEfficScan()->plot((TPad*) canvas.GetPrimitive(cName+"_2"));
      canvas.cd(3);
      effPurFromHistos[0+n]->plot((TPad*) canvas.GetPrimitive(cName+"_3"));
      canvas.cd(4);
      effPurFromHistos[1+n]->discriminatorNoCutEffic()->plot((TPad*) canvas.GetPrimitive(cName+"_4"));
      canvas.cd(5);
      effPurFromHistos[1+n]->discriminatorCutEfficScan()->plot((TPad*) canvas.GetPrimitive(cName+"_5"));
      canvas.cd(6);
      effPurFromHistos[1+n]->plot((TPad*) canvas.GetPrimitive(cName+"_6"));
    }
  }

  canvas.Print(name + cName + ".ps");
  canvas.Print(name + cName + ".ps]");
}

void TrackCountingTagPlotter::write()
{
  jetTagPlotter_->write();

  TString dir= "TrackCounting"+theExtensionString;
  gFile->cd();
  gFile->mkdir(dir);
  gFile->cd(dir);
  trkNbr2D->write();
  trkNbr3D->write();
  for(int n=0; n <= 4; n++) {
    tkcntHistosSig2D[n]->write();
    tkcntHistosSig3D[n]->write();
  }
  if (finalized) {
    for(int n=0; n < 4; n++) effPurFromHistos[n]->write();
  }
  gFile->cd();
}

void TrackCountingTagPlotter::epsPlot(const TString & name)
{
  jetTagPlotter_->epsPlot(name);
  trkNbr2D->epsPlot(name);
  trkNbr3D->epsPlot(name);
  for(int n=0; n <= 4; n++) {
    tkcntHistosSig2D[n]->epsPlot(name);
    tkcntHistosSig3D[n]->epsPlot(name);
  }
  if (finalized) {
    for(int n=0; n < 4; n++) effPurFromHistos[n]->epsPlot(name);
  }
}
