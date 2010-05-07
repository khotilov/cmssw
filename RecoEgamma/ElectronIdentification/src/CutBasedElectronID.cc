#include "RecoEgamma/ElectronIdentification/interface/CutBasedElectronID.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
//#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"
//#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <algorithm>

void CutBasedElectronID::setup(const edm::ParameterSet& conf) {
  
  // Get all the parameters
  baseSetup(conf);
  
  type_ = conf.getParameter<std::string>("electronIDType");
  quality_ = conf.getParameter<std::string>("electronQuality");
  version_ = conf.getParameter<std::string>("electronVersion");
  verticesCollection = conf.getParameter<edm::InputTag>("verticesCollection");
  
  if (type_ == "robust" || type_ == "classbased") {
    //if (quality_ == "loose" || quality_ == "tight" ||
    //    quality_ == "medium" || quality_ == "highenergy" ) {
    std::string stringCut = type_+quality_+"EleIDCuts"+version_;
    cuts_ = conf.getParameter<edm::ParameterSet>(stringCut);
    //}
    //else {
    //throw cms::Exception("Configuration")
    // << "Invalid electronQuality parameter in CutBasedElectronID: must be loose, tight or highenergy.\n";
    //}
  } 
  else {
    throw cms::Exception("Configuration")
      << "Invalid electronType parameter in CutBasedElectronID: must be robust or classbased\n";
  }
}

int CutBasedElectronID::classify(const reco::GsfElectron* electron) {
  
  double eta = fabs(electron->superCluster()->eta());
  double eOverP = electron->eSuperClusterOverP();
  //double pin  = electron->trackMomentumAtVtx().R(); 
  //double pout = electron->trackMomentumOut().R(); 
  //double fBrem = (pin-pout)/pin;
  double fBrem = electron->fbrem();

  int cat = -1;
  if (version_ == "V00" || version_ == "V01") {
    if((electron->isEB() && fBrem<0.06) || (electron->isEE() && fBrem<0.1)) 
      cat=1;
    else if (eOverP < 1.2 && eOverP > 0.8) 
      cat=0;
    else 
      cat=2;
    
    return cat;

  } else if (version_ == "V02") {
    if (electron->isEB()) {       // BARREL
      if(fBrem < 0.12)
        cat=1;
      else if (eOverP < 1.2 && eOverP > 0.9) 
        cat=0;
      else 
        cat=2;
    } else {                     // ENDCAP
      if(fBrem < 0.2)
        cat=1;
      else if (eOverP < 1.22 && eOverP > 0.82) 
        cat=0;
      else 
        cat=2;
    }
    
    return cat;
  } else if (version_ == "V03" || version_ == "V04" ||version_ == "") {
    if (electron->isEB()) {
      if ((fBrem >= 0.12) and (eOverP > 0.9) and (eOverP < 1.2))
        cat = 0;
      else if ((eta >  .445   and eta <  .45  ) ||
               (eta >  .79    and eta <  .81  ) ||
               (eta > 1.137   and eta < 1.157 ) ||
               (eta > 1.47285 and eta < 1.4744))
        cat = 6;
      else if ((electron->trackerDrivenSeed()) and (!electron->ecalDrivenSeed()))
        cat = 8;
      else if (fBrem < 0.12)
        cat = 1;
      else
        cat = 2;
    } else {
      if ((fBrem >= 0.2) and (eOverP > 0.82) and (eOverP < 1.22))
        cat = 3;
      else if (eta > 1.5 and eta <  1.58)
        cat = 7;
      else if ((electron->trackerDrivenSeed()) and (!electron->ecalDrivenSeed()))
        cat = 8;
      else if (fBrem < 0.2)
        cat = 4;
      else
        cat = 5;
    }

    return cat;
  }

  return -1;
}

double CutBasedElectronID::result(const reco::GsfElectron* electron ,
                                  const edm::Event& e ,
                                  const edm::EventSetup& es) { 

  double scTheta = (2*atan(exp(-electron->superCluster()->eta())));
  double scEt = electron->superCluster()->energy()*sin(scTheta);
 
  double eta = electron->p4().Eta();
  double eOverP = electron->eSuperClusterOverP();
  double eSeed = electron->superCluster()->seed()->energy();
  double pin  = electron->trackMomentumAtVtx().R();   
  double pout = electron->trackMomentumOut().R(); 
  double eSeedOverPin = eSeed/pin; 
  double fBrem = (pin-pout)/pin;
  double hOverE = electron->hadronicOverEm();
  EcalClusterLazyTools lazyTools = getClusterShape(e,es);
  std::vector<float> vLocCov = lazyTools.localCovariances(*(electron->superCluster()->seed())) ;
  double sigmaee = sqrt(vLocCov[0]);
  double e25Max = lazyTools.e2x5Max(*(electron->superCluster()->seed()))  ;
  double e15 = lazyTools.e1x5(*(electron->superCluster()->seed()))  ;
  double e55 = lazyTools.e5x5(*(electron->superCluster()->seed())) ;
  double e25Maxoe55 = e25Max/e55 ;
  double e15oe55 = e15/e55 ;
  double deltaPhiIn = electron->deltaPhiSuperClusterTrackAtVtx();
  double deltaEtaIn = electron->deltaEtaSuperClusterTrackAtVtx();

  double ip = 0;
  int mishits = electron->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
  double tkIso = electron->dr03TkSumPt();
  double ecalIso = electron->dr04EcalRecHitSumEt();
  double ecalIsoPed = (electron->isEB())?std::max(0.,ecalIso-1.):ecalIso;
  double hcalIso = electron->dr04HcalTowerSumEt();
  double hcalIso1 = electron->dr04HcalDepth1TowerSumEt();
  double hcalIso2 = electron->dr04HcalDepth2TowerSumEt();
  //
  // calculate the conversion track partner related criteria
  // calculate the reference point of the track
  const math::XYZPoint tpoint = electron->gsfTrack()->referencePoint();
  //const GlobalPoint gp(tpoint.x(), tpoint.y(), tpoint.z());
  // calculate the magnetic field for that point
  edm::ESHandle<MagneticField> magneticField;
  es.get<IdealMagneticFieldRecord>().get(magneticField);
  const  MagneticField *mField = magneticField.product();
  double bfield = mField->inTesla(GlobalPoint(0.,0.,0.)).z();//mField->inTesla(gp).mag();
  //
  edm::Handle<reco::TrackCollection> ctfTracks;
  e.getByLabel("generalTracks", ctfTracks); 
  ConversionFinder convFinder;

  if (version_ == "V00") {
     std::vector<float> vCov = lazyTools.covariances(*(electron->superCluster()->seed())) ;
     sigmaee = sqrt(vCov[0]);  
     if (electron->isEE())
       sigmaee = sigmaee - 0.02*(fabs(eta) - 2.3);   //correct sigmaetaeta dependence on eta in endcap
  }

  if (version_ != "V01" and version_ != "V00") {
    edm::Handle<reco::VertexCollection> vtxH;
    e.getByLabel(verticesCollection, vtxH);
    if (vtxH->size() != 0) {
      reco::VertexRef vtx(vtxH, 0);
      ip = fabs(electron->gsfTrack()->dxy(math::XYZPoint(vtx->x(),vtx->y(),vtx->z())));
    } else
      ip = fabs(electron->gsfTrack()->dxy());
    
    if (electron->isEB()) {
      std::vector<float> vCov = lazyTools.scLocalCovariances(*(electron->superCluster()));
      sigmaee = sqrt(vCov[0]); 
    } 
  }

  if (version_ == "V03" and type_ == "robust") {
    edm::Handle<reco::BeamSpot> pBeamSpot;
    // uses the same name for the vertex collection to avoid adding more new names
    e.getByLabel(verticesCollection, pBeamSpot); 
    if (pBeamSpot.isValid()) {
      const reco::BeamSpot *bspot = pBeamSpot.product();
      const math::XYZPoint bspotPosition = bspot->position();
      ip = fabs(electron->gsfTrack()->dxy(bspotPosition));
    } else
      ip = fabs(electron->gsfTrack()->dxy());
  }

  // .....................................................................................
  std::vector<double> cut;
  // ROBUST Selection
  if (type_ == "robust") {
    
    double result = 0;

    // hoe, sigmaEtaEta, dPhiIn, dEtaIn
    if (electron->isEB())
      cut = cuts_.getParameter<std::vector<double> >("barrel");
    else
      cut = cuts_.getParameter<std::vector<double> >("endcap");
    // check isolations: if only isolation passes result = 2   
    if (quality_ == "highenergy") {
      if ((tkIso > cut[6] || hcalIso2 > cut[11]) ||
          (electron->isEB() && ((ecalIso + hcalIso1) > cut[7]+cut[8]*scEt)) ||
          (electron->isEE() && (scEt >= 50.) && ((ecalIso + hcalIso1) > cut[7]+cut[8]*(scEt-50))) ||
          (electron->isEE() && (scEt < 50.) && ((ecalIso + hcalIso1) > cut[9]+cut[10]*(scEt-50))))
        result = 0;
      else
        result = 2;
    } else {

      //if (electron->isEB() && (version_ == "V03" || version_ == "V04" || version_ == "")) 
      //  ecalIso = std::max(0., ecalIso - 1.);
        
      if ((tkIso > cut[6]) || (ecalIso > cut[7]) || (hcalIso > cut[8]) || (hcalIso1 > cut[9]) || (hcalIso2 > cut[10]) || 
	  (tkIso/electron->p4().Pt() > cut[11]) || (ecalIso/electron->p4().Pt() > cut[12]) || (hcalIso/electron->p4().Pt() > cut[13]) ||
	  ((tkIso+ecalIso+hcalIso)>cut[14]) || (((tkIso+ecalIso+hcalIso)/ electron->p4().Pt()) > cut[15]) || 
	  ((tkIso+ecalIsoPed+hcalIso)>cut[16]) || (((tkIso+ecalIsoPed+hcalIso)/ electron->p4().Pt()) > cut[17])  )
        result = 0.;
      else
        result = 2.;
    }


    if (hOverE > cut[0]) 
      return result;    

    if (sigmaee > cut[1]) 
      return result;    

    if (fabs(deltaPhiIn) > cut[2]) 
      return result;    

    if (fabs(deltaEtaIn) > cut[3]) 
      return result;    
    
    if (e25Maxoe55 < cut[4] && e15oe55 < cut[5])
      return result;
    // some extra electron id cuts
    if (sigmaee < cut[18]) // inverted sigmaee cut - spike removal related
      return result;

    if (  eOverP < cut[19] ||  eOverP > cut[20]) // lower and upper cut in E/P
      return result;
    
    result = result + 1;

    if (ip > cut[21])
      return result;
    if (mishits > cut[22]) // expected missing hits
      return result;
    // positive cut[23] means to demand a valid hit in 1st layer PXB 
    if (cut[23] >0 && not (electron->gsfTrack()->hitPattern().hasValidHitInFirstPixelBarrel()))
      return result;
    // cut[24]: Dist cut[25]: dcot
    if (convFinder.isElFromConversion(*electron, ctfTracks, bfield, cut[24], cut[25]))
      return result;
      
    result += 4;
    
    return result;
  }
  
  int cat = classify(electron);
  int eb;

  if (electron->isEB()) 
    eb = 0;
  else 
    eb = 1; 

  // LOOSE and TIGHT Selections
  if (type_ == "classbased" && (version_ == "V01" || version_ == "V00")) {
    
    if ((eOverP < 0.8) && (fBrem < 0.2)) 
      return 0.;
    
    cut = cuts_.getParameter<std::vector<double> >("hOverE");
    if (hOverE > cut[cat+4*eb]) 
      return 0.;    
    
    cut = cuts_.getParameter<std::vector<double> >("sigmaEtaEta");
    if (sigmaee > cut[cat+4*eb]) 
      return 0.;    
    
    cut = cuts_.getParameter<std::vector<double> >("deltaPhiIn");
    if (eOverP < 1.5) {
      if (fabs(deltaPhiIn) > cut[cat+4*eb]) 
        return 0.;    
    } else {
      if (fabs(deltaPhiIn) > cut[3+4*eb])
        return 0.;
    }
    
    cut = cuts_.getParameter<std::vector<double> >("deltaEtaIn");
    if (fabs(deltaEtaIn) > cut[cat+4*eb]) 
      return 0.;    
    
    cut = cuts_.getParameter<std::vector<double> >("eSeedOverPin");
    if (eSeedOverPin < cut[cat+4*eb]) 
      return 0.;  
    
    if (quality_ == "tight")
      if (eOverP < 0.9*(1-fBrem))
        return 0.;
    
    return 1.;
  }
  
  if (type_ == "classbased" and version_ == "V02") {
    double result = 0.;

    int bin = 0;

    if (scEt < 20.)
      bin = 2;
    else if (scEt > 30.)
      bin = 0;
    else
      bin = 1;

    if (fBrem > 0)
      eSeedOverPin = eSeedOverPin + fBrem;
    
    if (bin != 2) {     
      tkIso = tkIso*pow(40./scEt, 2); 
      ecalIso = ecalIso*pow(40./scEt, 2); 
      hcalIso = hcalIso*pow(40./scEt, 2); 
    }

    std::vector<double> cutTk = cuts_.getParameter<std::vector<double> >("cutisotk");
    std::vector<double> cutEcal = cuts_.getParameter<std::vector<double> >("cutisoecal");
    std::vector<double> cutHcal = cuts_.getParameter<std::vector<double> >("cutisohcal");
    if ((tkIso > cutTk[cat+3*eb+bin*6]) ||
        (ecalIso > cutEcal[cat+3*eb+bin*6]) ||
        (hcalIso > cutHcal[cat+3*eb+bin*6]))
      result = 0.;
    else
      result = 2.;

    if (fBrem > -2) {

      //std::cout << "hoe" << hOverE << std::endl;
      std::vector<double> cuthoe = cuts_.getParameter<std::vector<double> >("cuthoe");
      //std::cout << "see" << sigmaee << std::endl;
      std::vector<double> cutsee = cuts_.getParameter<std::vector<double> >("cutsee");
      //std::cout << "dphiin" << fabs(deltaPhiIn) << std::endl;
      std::vector<double> cutdphi = cuts_.getParameter<std::vector<double> >("cutdphi");
      //std::cout << "detain" << fabs(deltaEtaIn) << std::endl;
      std::vector<double> cutdeta = cuts_.getParameter<std::vector<double> >("cutdeta");
      //std::cout << "eseedopin " << eSeedOverPin << std::endl;
      std::vector<double> cuteopin = cuts_.getParameter<std::vector<double> >("cuteopin");
      //std::cout << "ip" << ip << std::endl;
      std::vector<double> cutip = cuts_.getParameter<std::vector<double> >("cutip");
      std::vector<double> cutmishits = cuts_.getParameter<std::vector<double> >("cutmishits");

      if ((hOverE < cuthoe[cat+3*eb+bin*6]) and
          (sigmaee < cutsee[cat+3*eb+bin*6]) and
          (fabs(deltaPhiIn) < cutdphi[cat+3*eb+bin*6]) and
          (fabs(deltaEtaIn) < cutdeta[cat+3*eb+bin*6]) and
          (eSeedOverPin > cuteopin[cat+3*eb+bin*6]) and
          (ip < cutip[cat+3*eb+bin*6]) and
          (mishits < cutmishits[cat+3*eb+bin*6]))
        result = result + 1.;
    }
    return result;
  }

  if (type_ == "classbased" && (version_ == "V03" || version_ == "V04" || version_ == "")) {
    double result = 0.;
    
    int bin = 0;
    
    if (scEt < 20.)
      bin = 2;
    else if (scEt > 30.)
      bin = 0;
    else
      bin = 1;

    if (fBrem > 0)
      eSeedOverPin = eSeedOverPin + fBrem;
    
    float iso_sum = tkIso + ecalIso + hcalIso;
    float iso_sum_corrected = iso_sum*pow(40./scEt, 2);

    std::vector<double> cutIsoSum = cuts_.getParameter<std::vector<double> >("cutiso_sum");
    std::vector<double> cutIsoSumCorr = cuts_.getParameter<std::vector<double> >("cutiso_sumoet");
    if ((iso_sum < cutIsoSum[cat+bin*9]) and
        (iso_sum_corrected < cutIsoSumCorr[cat+bin*9]))
      result += 2.;

    if (fBrem > -2) {
      //std::cout << "hoe" << hOverE << std::endl;
      std::vector<double> cuthoe = cuts_.getParameter<std::vector<double> >("cuthoe");
      //std::cout << "see" << sigmaee << std::endl;
      std::vector<double> cutsee = cuts_.getParameter<std::vector<double> >("cutsee");
      //std::cout << "dphiin" << fabs(deltaPhiIn) << std::endl;
      std::vector<double> cutdphi = cuts_.getParameter<std::vector<double> >("cutdphiin");
      //std::cout << "detain" << fabs(deltaEtaIn) << std::endl;
      std::vector<double> cutdeta = cuts_.getParameter<std::vector<double> >("cutdetain");
      //std::cout << "eseedopin " << eSeedOverPin << std::endl;
      std::vector<double> cuteopin = cuts_.getParameter<std::vector<double> >("cuteseedopcor");
      std::vector<double> cutet = cuts_.getParameter<std::vector<double> >("cutet");

      if ((hOverE < cuthoe[cat+bin*9]) and
          (sigmaee < cutsee[cat+bin*9]) and
          (fabs(deltaPhiIn) < cutdphi[cat+bin*9]) and
          (fabs(deltaEtaIn) < cutdeta[cat+bin*9]) and
          (eSeedOverPin > cuteopin[cat+bin*9]) and
          (scEt > cutet[cat+bin*9]))
        result += 1.;
    }

    //std::cout << "ip" << ip << std::endl;
    std::vector<double> cutip = cuts_.getParameter<std::vector<double> >("cutip_gsf");
    if (ip < cutip[cat+bin*9])
      result += 8;
    
    std::vector<double> cutmishits = cuts_.getParameter<std::vector<double> >("cutfmishits");
    std::vector<double> cutdcotdist = cuts_.getParameter<std::vector<double> >("cutdcotdist");
    
    std::pair<double, double> convParam(9999., 9999.);
    reco::TrackRef convTk = convFinder.getConversionPartnerTrack(*electron, ctfTracks, bfield, 0.04, 0.04, 0.45);
    if (convTk.isNonnull()) {
      math::XYZTLorentzVector tklv(convTk->px(), convTk->py(),
                                   convTk->pz(), convTk->p());
      
      convParam = convFinder.getConversionInfo(electron->p4(), electron->charge(), 
                                                electron->gsfTrack()->d0(), 
                                                tklv, convTk->charge(), 
                                                convTk->d0(), bfield);  
    }
    
    Float_t dcotdistcomb = ((0.4 - std::max(convParam.first, convParam.second)) > 0?(0.4 - std::max(convParam.first, convParam.second)):0);
    
    if ((mishits < cutmishits[cat+bin*9]) and 
        (dcotdistcomb < cutdcotdist[cat+bin*9]))
      result += 4;
    
    return result;
  }

  return -1.;
}
