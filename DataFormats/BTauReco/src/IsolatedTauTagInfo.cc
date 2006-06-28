#include "DataFormats/BTauReco/interface/IsolatedTauTagInfo.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/GenVector/VectorUtil.h>
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

using namespace edm;
using namespace reco;

const RefVector<TrackCollection> IsolatedTauTagInfo::tracksInCone( const math::XYZVector myVector, const float size,  const float pt_min) const { 
  
  RefVector<TrackCollection> tmp;
 
  RefVector<TrackCollection>::const_iterator myTrack = selectedTracks_.begin();
  for(;myTrack != selectedTracks_.end(); myTrack++)
    {
      const math::XYZVector trackMomentum = (*myTrack)->momentum() ;
      float pt_tk = (*myTrack)->pt();
      float deltaR = ROOT::Math::VectorUtil::DeltaR(myVector, trackMomentum);
      if ( deltaR < size && pt_tk > pt_min) tmp.push_back( *myTrack);
      }

  //  sort(tmp.begin(), tmp.end(), SortByDescendingTrackPt());
  return tmp;
}

const RefVector<TrackCollection> IsolatedTauTagInfo::tracksInCone( const math::XYZVector myVector, const float size,  const float pt_min, const float z_pv, const float dz_lt) const { 
  
  RefVector<TrackCollection> tmp;
 
  RefVector<TrackCollection>::const_iterator myTrack = selectedTracks_.begin();
  for(;myTrack != selectedTracks_.end(); myTrack++)
    {
      const math::XYZVector trackMomentum = (*myTrack)->momentum() ;
      float pt_tk = (*myTrack)->pt();
      float deltaR = ROOT::Math::VectorUtil::DeltaR(myVector, trackMomentum);
      if ( deltaR < size && pt_tk > pt_min && fabs((*myTrack)->dz() - z_pv) < dz_lt) tmp.push_back( *myTrack);
      }

  //  sort(tmp.begin(), tmp.end(), SortByDescendingTrackPt());
  return tmp;
}


const TrackRef IsolatedTauTagInfo::leadingSignalTrack(const float rm_cone, const float pt_min) const {

  const Jet & myjet = m_jetTag->jet(); 
  math::XYZVector jet3Vec   (myjet.px(),myjet.py(),myjet.pz()) ;

  const  RefVector<TrackCollection>  sTracks = tracksInCone(jet3Vec, rm_cone, pt_min);
  TrackRef leadTk;
  float pt_cut = pt_min;
  if (sTracks.size() >0) 
    {
      RefVector<TrackCollection>::const_iterator myTrack =sTracks.begin();
      for(;myTrack!=sTracks.end();myTrack++)
	{
	  if((*myTrack)->pt() > pt_cut) {
	    leadTk = *myTrack;
	    pt_cut = (*myTrack)->pt();
	  }
	}
    }
  return leadTk;
}


const TrackRef IsolatedTauTagInfo::leadingSignalTrack(const math::XYZVector myVector, const float rm_cone, const float pt_min) const {


const  RefVector<TrackCollection>  sTracks = tracksInCone(myVector, rm_cone, pt_min);
 TrackRef leadTk;
  float pt_cut = pt_min;
  if (sTracks.size() >0) 
    {
      RefVector<TrackCollection>::const_iterator myTrack =sTracks.begin();
      for(;myTrack!=sTracks.end();myTrack++)
	{
	  if((*myTrack)->pt() > pt_cut) {
	    leadTk = *myTrack;
	    pt_cut = (*myTrack)->pt();
	  }
	}
    }
  return leadTk;
}

double IsolatedTauTagInfo::discriminator(float m_cone, float sig_cone, float iso_cone, float pt_min_lt, float pt_min_tk, int nTracksIsoRing=0) const
{
  double myDiscriminator = 0;

const TrackRef leadTk = leadingSignalTrack(m_cone, pt_min_lt);

  if(!leadTk) {
    return myDiscriminator;
  }
  math::XYZVector trackMomentum = leadTk->momentum() ;
const  RefVector<TrackCollection> signalTracks = tracksInCone(trackMomentum, sig_cone , pt_min_tk);
const RefVector<TrackCollection> isolationTracks =tracksInCone(trackMomentum, iso_cone , pt_min_tk); 
  
  if(signalTracks.size() > 0 && (isolationTracks.size() - signalTracks.size())<nTracksIsoRing+1 )
    myDiscriminator=1;

  return myDiscriminator;
}
double IsolatedTauTagInfo::discriminator(math::XYZVector myVector, float m_cone, float sig_cone, float iso_cone, float pt_min_lt, float pt_min_tk, int nTracksIsoRing) const
{
  double myDiscriminator = 0;

const  TrackRef leadTk = leadingSignalTrack(myVector, m_cone, pt_min_lt);
  if(!leadTk) return myDiscriminator;

  math::XYZVector trackMomentum = leadTk->momentum() ;
const   RefVector<TrackCollection> signalTracks = tracksInCone(trackMomentum, sig_cone , pt_min_tk);
const  RefVector<TrackCollection> isolationTracks =tracksInCone(trackMomentum, iso_cone , pt_min_tk); 
  
  if(signalTracks.size() > 0 && (isolationTracks.size() - signalTracks.size())<nTracksIsoRing+1 )
    myDiscriminator=1;

  return myDiscriminator;
}

double IsolatedTauTagInfo::discriminator(float m_cone, float sig_cone, float iso_cone, float pt_min_lt, float pt_min_tk, int nTracksIsoRing, float dz_lt) const
{
  double myDiscriminator = 0;

const TrackRef leadTk = leadingSignalTrack(m_cone, pt_min_lt);

  if(!leadTk) {
    return myDiscriminator;
  }
  math::XYZVector trackMomentum = leadTk->momentum() ;
  float z_pv = leadTk->dz();
const  RefVector<TrackCollection> signalTracks = tracksInCone(trackMomentum, sig_cone , pt_min_tk, z_pv, dz_lt);
const RefVector<TrackCollection> isolationTracks =tracksInCone(trackMomentum, iso_cone , pt_min_tk, z_pv, dz_lt); 
  
  if(signalTracks.size() > 0 && (isolationTracks.size() - signalTracks.size())<nTracksIsoRing+1 )
    myDiscriminator=1;

  return myDiscriminator;
}
double IsolatedTauTagInfo::discriminator(math::XYZVector myVector, float m_cone, float sig_cone, float iso_cone, float pt_min_lt, float pt_min_tk, int nTracksIsoRing, float dz_lt) const
{
  double myDiscriminator = 0;

const  TrackRef leadTk = leadingSignalTrack(myVector, m_cone, pt_min_lt);
  if(!leadTk) return myDiscriminator;

  math::XYZVector trackMomentum = leadTk->momentum() ;
  float z_pv = leadTk->dz();
  const   RefVector<TrackCollection> signalTracks = tracksInCone(trackMomentum, sig_cone , pt_min_tk, z_pv, dz_lt);
  const  RefVector<TrackCollection> isolationTracks =tracksInCone(trackMomentum, iso_cone , pt_min_tk, z_pv, dz_lt); 
  
  if(signalTracks.size() > 0 && (isolationTracks.size() - signalTracks.size())<nTracksIsoRing+1 )
    myDiscriminator=1;

  return myDiscriminator;




}
