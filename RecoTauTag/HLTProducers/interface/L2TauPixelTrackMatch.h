#ifndef L2TauPixelTrackMatch_h
#define L2TauPixelTrackMatch_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/Point3D.h"
#include <vector>


/** class L2TauPixelTrackMatch
 * this producer creates a new L2 tau jet collection with jets' vertices redefined 
 * from vertex z (relative to beamspot) of dr-matched pixel tracks 
 * that are above some pt threshold and beamline x & y.
 */
class L2TauPixelTrackMatch : public edm::EDProducer
{
public:

  explicit L2TauPixelTrackMatch(const edm::ParameterSet&);
  ~L2TauPixelTrackMatch();
  virtual void produce(edm::Event&, const edm::EventSetup&);

private:

  struct TinyTrack
  {
    float pt, eta, phi;
    math::XYZPoint vtx;
  };

  edm::InputTag m_jetSrc;
  float m_jetMinPt;
  float m_jetMaxEta;
  edm::InputTag m_trackSrc;
  float m_trackMinPt;
  float m_deltaR;
  edm::InputTag m_beamSpotTag;
};

#endif 
