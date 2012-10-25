#ifndef TauAnalysis_MCEmbeddingTools_MuonDetCleaner_h
#define TauAnalysis_MCEmbeddingTools_MuonDetCleaner_h

/** \class MuonDetCleaner
 *
 * Produce collections of recHits in muon detectors
 * from which the hits of the two muons from the Z -> mu+ mu- decay are subtracted
 * (the remaining hits may be due to pile-up/heavy quark decays in Z+jets events)
 *
 * \author Christian Veelken, LLR
 *
 * \version $Revision: 1.1 $
 *
 * $Id: MuonDetCleaner.h,v 1.1 2012/10/24 09:37:14 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <DataFormats/MuonReco/interface/Muon.h>

#include <map>

class MuonDetCleaner : public edm::EDProducer 
{
 public:
  explicit MuonDetCleaner(const edm::ParameterSet&);
  ~MuonDetCleaner();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  typedef std::map<uint32_t, int> detIdToIntMap;
  void fillHitMap(const reco::Muon*, detIdToIntMap&);
      
  edm::InputTag srcSelectedMuons_;

  int verbosity_;
};

#endif
