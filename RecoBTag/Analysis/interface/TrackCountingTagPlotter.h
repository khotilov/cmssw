#ifndef TrackCountingTagPlotter_H
#define TrackCountingTagPlotter_H

#include "RecoBTag/Analysis/interface/BaseBTagPlotter.h"
#include "RecoBTag/Analysis/interface/JetTagPlotter.h"
#include "DataFormats/BTauReco/interface/TrackCountingTagInfo.h"
#include "RecoBTag/Analysis/interface/FlavourHistorgrams.h"
#include "RecoBTag/MCTools/interface/JetFlavour.h"


class TrackCountingTagPlotter : public BaseBTagPlotter {


 public:

  TrackCountingTagPlotter (JetTagPlotter *jetTagPlotter, bool update = false);

  ~TrackCountingTagPlotter () ;

  void analyzeTag (const reco::TrackCountingTagInfo & tagInfo,
  	const reco::JetTag & jetTag, const JetFlavour & jetFlavour);

  virtual void finalize ();

  // get "2d" histograms for misid. vs. b-eff
  virtual EffPurFromHistos * getEffPurFromHistos ()
	{return jetTagPlotter_->getEffPurFromHistos();}

  virtual void write ();

  void epsPlot(const TString & name);

  void psPlot(const TString & name);

 private:

  JetTagPlotter *jetTagPlotter_;

  FlavourHistorgrams<double> * hSignificance;
  FlavourHistorgrams<double> * tkcntHistosSig[4];
  EffPurFromHistos * effPurFromHistos2 ;
  EffPurFromHistos * effPurFromHistos3 ;
  bool finalized;
} ;

#endif
