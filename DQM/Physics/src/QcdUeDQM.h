
/*
    This is the DQM code for UE physics plots
    11/12/2009 Sunil Bansal
*/

#ifndef QcdUeDQM_H
#define QcdUeDQM_H

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "AnalysisDataFormats/TrackInfo/interface/TrackInfo.h"
#include "AnalysisDataFormats/TrackInfo/interface/TrackInfoTrackAssociation.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h" 
#include <TMath.h>
#include <vector>
#define PI 3.141592654
class DQMStore;
class MonitorElement;
class TrackerGeometry;
class TH1F;
class TH2F;
class TH3F;
class TProfile;

class PtSorter {
public:
  template <class T> bool operator() ( const T& a, const T& b ) {
    return ( a.pt() > b.pt() );
  }
};


class QcdUeDQM : public edm::EDAnalyzer 
{
  public:

    QcdUeDQM(const edm::ParameterSet &parameters);
    virtual ~QcdUeDQM();

    void                          analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup);
    void                          beginJob(void);
    void                          beginLuminosityBlock(const edm::LuminosityBlock &l, 
                                                       const edm::EventSetup &iSetup);
    void                          beginRun(const edm::Run &r, const edm::EventSetup &iSetup);
    void                          endJob(void);
    void                          endRun(const edm::Run &r, const edm::EventSetup &iSetup);
    void                          endLuminosityBlock(const edm::LuminosityBlock &l, 
                                                     const edm::EventSetup &iSetup);

  private:
    void                          book1D(std::vector<MonitorElement*> &mes, 
                                         const std::string &name, const std::string &title, 
                                         int nx, double x1, double x2, bool sumw2=1, bool sbox=1);
    void                          book2D(std::vector<MonitorElement*> &mes, 
                                         const std::string &name, const std::string &title, 
                                         int nx, double x1, double x2, int ny, double y1, double y2,
                                         bool sumw2=1, bool sbox=1);
    void                          bookProfile(std::vector<MonitorElement*> &mes, 
                                         const std::string &name, const std::string &title, 
                                         int nx, double x1, double x2,  double y1, double y2,
                                         bool sumw2=1, bool sbox=1); 
    void                          create1D(std::vector<TH1F*> &mes, 
                                           const std::string &name, const std::string &title, 
                                           int nx, double x1, double x2, bool sumw2=1, bool sbox=1);
    void                          create2D(std::vector<TH2F*> &mes, 
                                           const std::string &name, const std::string &title, 
                                           int nx, double x1, double x2, int ny, double y1, double y2,
                                           bool sumw2=1, bool sbox=1);
    void                          createProfile(std::vector<TProfile*> &mes, 
                                           const std::string &name, const std::string &title, 
                                           int nx, double x1, double x2,  double y1, double y2,
                                           bool sumw2=1, bool sbox=1);  
    void                          createHistos();
    void                          fill1D(std::vector<TH1F*> &hs, double val, double w=1.);
    void                          fill1D(std::vector<MonitorElement*> &mes, double val, double w=1.);
    void                          fill2D(std::vector<TH2F*> &hs, 
                                         double valx, double valy, double w=1.);
    void                          fill2D(std::vector<MonitorElement*> &mes, 
                                         double valx, double valy, double w=1.);
    void                          fillProfile(std::vector<TProfile*> &hs, 
                                         double valx, double valy, double w=1.);
    void                          fillProfile(std::vector<MonitorElement*> &mes, 
                                        double valx, double valy, double w=1.);
    void                          fill3D(std::vector<TH3F*> &hs, int gbin, double w=1.);
    void                          setLabel1D(std::vector<MonitorElement*> &mes);
   


    void                          fillHltBits(const edm::Event &iEvent);
    
    bool                          fillVtxPlots(const edm::Event &iEvent, const edm::Handle< reco::VertexCollection > vtxColl);
    void                          fillpTMaxRelated(const edm::Event &iEvent, const edm::Handle<reco::TrackCollection> &track, const edm::Handle< reco::VertexCollection > vtxColl);

    void                          fillChargedJetSpectra(const edm::Event &iEvent, const edm::Handle<reco::CandidateView> trackJets);
  
    void                          fillCaloJetSpectra(const edm::Event &iEvent, const edm::Handle<reco::CaloJetCollection> caloJets);

    void                          fillUE_with_ChargedJets(const edm::Event &iEvent, const reco::TrackCollection &track,  const edm::Handle<reco::CandidateView> &trackJets, const edm::Handle< reco::VertexCollection > vtxColl);
    void                          fillUE_with_CaloJets(const edm::Event &iEvent, const reco::TrackCollection &track, const edm::Handle<reco::CaloJetCollection> &caloJets, const edm::Handle< reco::VertexCollection > vtxColl);
    void                          fillUE_with_MaxpTtrack(const edm::Event &iEvent, const reco::TrackCollection &track, const edm::Handle< reco::VertexCollection > vtxColl); 

   
    template <typename TYPE>
    void                          getProduct(const std::string name, edm::Handle<TYPE> &prod,
                                             const edm::Event &event) const;    
    template <typename TYPE>
    bool                          getProductSafe(const std::string name, edm::Handle<TYPE> &prod,
                                                 const edm::Event &event) const;

    HLTConfigProvider hltConfig;   

    std::string                   hltResName_;         //HLT trigger results name
    std::vector<std::string>      hltProcNames_;       //HLT process name(s)
    std::vector<std::string>      hltTrgNames_;        //HLT trigger name(s)
    
   
    std::vector<int>              hltTrgBits_;         //HLT trigger bit(s)
    std::vector<bool>             hltTrgDeci_;         //HLT trigger descision(s)
    std::vector<std::string>      hltTrgUsedNames_;    //HLT used trigger name(s)
    std::string                   hltUsedResName_;     //used HLT trigger results name 
    int                           verbose_;            //verbosity (0=debug,1=warn,2=error,3=throw)
    const TrackerGeometry        *tgeo_;               //tracker geometry
    DQMStore                     *theDbe_;             //dqm store
    MonitorElement               *repSumMap_;          //report summary map
    MonitorElement               *repSummary_;         //report summary
    MonitorElement               *h2TrigCorr_;         //trigger correlation plot


    std::vector<MonitorElement*>  hEvtSel_pTMax_;        // Event selection efficiency
    std::vector<MonitorElement*>  hEvtSel_ChargedJet_;        // Event selection efficiency
    std::vector<MonitorElement*>  hEvtSel_CaloJet_;        // Event selection efficiency
    std::vector<MonitorElement*>  hNvertices_;           // Number of vertices
    std::vector<MonitorElement*>  hVertex_z_;            // z-position of vertex
    std::vector<MonitorElement*>  hNTrack500_;           //number of tracks with pT > 500 MeV 


    std::vector<MonitorElement*>  hTrack_pTSpectrum_;      //pt spectrum of leading track
    std::vector<MonitorElement*>  hTrack_etaSpectrum_;      //eta spectrum of leading track
    std::vector<MonitorElement*>  hTrack_phiSpectrum_;      //phi spectrum of leading track

    
   
   

    std::vector<MonitorElement*>   hChargedJetMulti_;         // Number of charged jets
    std::vector<MonitorElement*>   hChargedJetConstituent_;         // Number of constituent of charged jets 
    std::vector<MonitorElement*>   hChargedJet_pTSpectrum_;         // pT spectrum of charged jets 

    std::vector<MonitorElement*>   hCaloJetMulti_;         // Number of calo jets
    std::vector<MonitorElement*>   hCaloJetConstituent_;         // Number of constituent of calo jets
    std::vector<MonitorElement*>   hCaloJet_pTSpectrum_;         // pT spectrum of calo jets
    std::vector<MonitorElement*>  hCaloJet_etaSpectrum_;      //eta spectrum of leading calo jet
    std::vector<MonitorElement*>  hCaloJet_phiSpectrum_;      //phi spectrum of leading calo jet
   
    std::vector<MonitorElement*>  hdPhi_maxpTTrack_tracks_;  // delta phi between leading track and tracks
    std::vector<MonitorElement*>  hdPhi_caloJet_tracks_;  // delta phi between leading calo jet and tracks  
    std::vector<MonitorElement*>  hdPhi_chargedJet_tracks_;  // delta phi between leading charged jet and tracks
   

    std::vector<MonitorElement*>  hChargedJet_etaSpectrum_;      //eta spectrum of leading charged jet
    std::vector<MonitorElement*>  hChargedJet_phiSpectrum_;      //phi spectrum of leading charged jet
  
    std::vector<MonitorElement*>  hdNdEtadPhi_pTMax_Toward500_;    // number of tracks in toward region of leadin track (pT > 500 MeV)
    std::vector<MonitorElement*>  hdNdEtadPhi_pTMax_Transverse500_;    // number of tracks in transverse region of leadin track (pT > 500 MeV)  
    std::vector<MonitorElement*>  hdNdEtadPhi_pTMax_Away500_;    // number of tracks in away region of leadin track (pT > 500 MeV)
    std::vector<MonitorElement*>  hdNdEtadPhi_caloJet_Toward500_;    // number of tracks in toward region of leadin calo Jet (pT > 500 MeV)
    std::vector<MonitorElement*>  hdNdEtadPhi_caloJet_Transverse500_;    // number of tracks in transverse region of leadin calo Jet (pT > 500 MeV)  
    std::vector<MonitorElement*>  hdNdEtadPhi_caloJet_Away500_;    // number of tracks in away region of leadin calo Jet (pT > 500 MeV)
    std::vector<MonitorElement*>  hdNdEtadPhi_trackJet_Toward500_;    // number of tracks in toward region of leadin calo Jet (pT > 500 MeV)
    std::vector<MonitorElement*>  hdNdEtadPhi_trackJet_Transverse500_;    // number of tracks in transverse region of leadin calo Jet (pT > 500 MeV)  
    std::vector<MonitorElement*>  hdNdEtadPhi_trackJet_Away500_;    // number of tracks in away region of leadin calo Jet (pT > 500 MeV) 


    
    std::vector<MonitorElement*>  hpTSumdEtadPhi_pTMax_Toward500_;    // pT sum of tracks in toward region of leadin track (pT > 500 MeV)
    std::vector<MonitorElement*>  hpTSumdEtadPhi_pTMax_Transverse500_;    // pT sum of tracks in transverse region of leadin track (pT > 500 MeV)  
    std::vector<MonitorElement*>  hpTSumdEtadPhi_pTMax_Away500_;    // pT sum of tracks in away region of leadin track (pT > 500 MeV)
    std::vector<MonitorElement*>  hpTSumdEtadPhi_caloJet_Toward500_;    // pT sum of tracks in toward region of leadin calo Jet (pT > 500 MeV)
    std::vector<MonitorElement*>  hpTSumdEtadPhi_caloJet_Transverse500_;    // pT sum of tracks in transverse region of leadin calo Jet (pT > 500 MeV)  
    std::vector<MonitorElement*>  hpTSumdEtadPhi_caloJet_Away500_;    // pT sum of tracks in away region of leadin calo Jet (pT > 500 MeV)
    std::vector<MonitorElement*>  hpTSumdEtadPhi_trackJet_Toward500_;    // pT sum of tracks in toward region of leadin calo Jet (pT > 500 MeV)
    std::vector<MonitorElement*>  hpTSumdEtadPhi_trackJet_Transverse500_;    // pT sum of tracks in transverse region of leadin calo Jet (pT > 500 MeV)  
    std::vector<MonitorElement*>  hpTSumdEtadPhi_trackJet_Away500_;    // pT sum of tracks in away region of leadin calo Jet (pT > 500 MeV)
   
    

 
   

    edm::InputTag caloJetLabel_;
    edm::InputTag chargedJetLabel_;
    edm::InputTag trackLabel_;
    edm::InputTag vtxLabel_;
    
};

//--------------------------------------------------------------------------------------------------
template <typename TYPE>
inline void QcdUeDQM::getProduct(const std::string name, edm::Handle<TYPE> &prod,
                                    const edm::Event &event) const
{
  // Try to access data collection from EDM file. We check if we really get just one
  // product with the given name. If not we throw an exception.

  event.getByLabel(edm::InputTag(name),prod);
  if (!prod.isValid()) 
    throw edm::Exception(edm::errors::Configuration, "QcdUeDQM::GetProduct()\n")
      << "Collection with label " << name << " is not valid" <<  std::endl;
}

//--------------------------------------------------------------------------------------------------
template <typename TYPE>
inline bool QcdUeDQM::getProductSafe(const std::string name, edm::Handle<TYPE> &prod,
                                        const edm::Event &event) const
{
  // Try to safely access data collection from EDM file. We check if we really get just one
  // product with the given name. If not, we return false.

  if (name.size()==0)
    return false;

  try {
    event.getByLabel(edm::InputTag(name),prod);
    if (!prod.isValid()) 
      return false;
  } catch (...) {
    return false;
  }
  return true;
}

//--------------------------------------------------------------------------------------------------
#endif
