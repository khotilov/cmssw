#ifndef RecoMuon_L2MuonProducer_L2MuonProducer_H
#define RecoMuon_L2MuonProducer_L2MuonProducer_H

//-------------------------------------------------
//
/**  \class L2MuonProducer
 * 
 *   L2 muon reconstructor:
 *   reconstructs muons using DT, CSC and RPC
 *   information,<BR>
 *   starting from internal seeds (L2 muon track segments).
 *
 *
 *   $Date: 2006/05/19 15:23:20 $
 *   $Revision: 1.2 $
 *
 *   \author  R.Bellan - INFN TO
 */
//
//--------------------------------------------------

#include "FWCore/Framework/interface/EDProducer.h"

namespace edm {class ParameterSet; class Event; class EventSetup;}

class MuonTrackFinder;
class MuonServiceProxy;

class L2MuonProducer : public edm::EDProducer {

  public:

  /// constructor with config
  L2MuonProducer(const edm::ParameterSet&);
  
  /// destructor
  virtual ~L2MuonProducer(); 
  
  /// reconstruct muons
  virtual void produce(edm::Event&, const edm::EventSetup&);
  // ex virtual void reconstruct();
  
 private:

  // MuonSeed Collection Label
  std::string theSeedCollectionLabel;

  /// the track finder
  MuonTrackFinder* theTrackFinder; //It isn't the same as in ORCA

  /// the event setup proxy, it takes care the services update
  MuonServiceProxy *theService;
};

#endif
