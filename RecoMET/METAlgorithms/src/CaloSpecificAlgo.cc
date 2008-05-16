#include "DataFormats/Math/interface/LorentzVector.h"
#include "RecoMET/METAlgorithms/interface/CaloSpecificAlgo.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

using namespace reco;
using namespace std;

//-------------------------------------------------------------------------
// This algorithm adds calorimeter specific global event information to 
// the MET object which may be useful/needed for MET Data Quality Monitoring
// and MET cleaning.  This list is not exhaustive and additional 
// information will be added in the future. 
//-------------------------------------
//reco::CaloMET CaloSpecificAlgo::addInfo(const CandidateCollection *towers, CommonMETData met)
reco::CaloMET CaloSpecificAlgo::addInfo(edm::Handle<edm::View<Candidate> > towers, CommonMETData met, bool noHF)
{ 
  // Instantiate the container to hold the calorimeter specific information
  SpecificCaloMETData specific;
  // Initialise the container 
  specific.MaxEtInEmTowers = 0.0;         // Maximum energy in EM towers
  specific.MaxEtInHadTowers = 0.0;        // Maximum energy in HCAL towers
  specific.HadEtInHO = 0.0;          // Hadronic energy fraction in HO
  specific.HadEtInHB = 0.0;          // Hadronic energy in HB
  specific.HadEtInHF = 0.0;          // Hadronic energy in HF
  specific.HadEtInHE = 0.0;          // Hadronic energy in HE
  specific.EmEtInEB = 0.0;           // Em energy in EB
  specific.EmEtInEE = 0.0;           // Em energy in EE
  specific.EmEtInHF = 0.0;           // Em energy in HF
  specific.EtFractionHadronic = 0.0; // Hadronic energy fraction
  specific.EtFractionEm = 0.0;       // Em energy fraction
  specific.METSignificance = -1.0;    // MET Significance
  specific.CaloSETInpHF = 0.0;        // CaloSET in HF+ 
  specific.CaloSETInmHF = 0.0;        // CaloSET in HF- 
  specific.CaloMETInpHF = 0.0;        // CaloMET in HF+ 
  specific.CaloMETInmHF = 0.0;        // CaloMET in HF- 
  specific.CaloMETPhiInpHF = 0.0;     // CaloMET-phi in HF+ 
  specific.CaloMETPhiInmHF = 0.0;     // CaloMET-phi in HF- 
  
  double totalEt = 0.0; 
  double totalEm     = 0.0;
  double totalHad    = 0.0;
  double MaxTowerEm  = 0.0;
  double MaxTowerHad = 0.0;
  double sumEtInpHF = 0.0;
  double sumEtInmHF = 0.0;
  double MExInpHF = 0.0;
  double MEyInpHF = 0.0;
  double MExInmHF = 0.0;
  double MEyInmHF = 0.0;

  if( towers->size() == 0 )  // if there are no towers, return specific = 0
    {
   //   LogDebug("CaloMET") << "Number of Candidate CaloTowers is zero : Unable to calculate calo specific info. " ;
      const LorentzVector p4( met.mex, met.mey, 0.0, met.met );
      const Point vtx( 0.0, 0.0, 0.0 );
      CaloMET specificmet( specific, met.sumet, p4, vtx );
      return specificmet;
    }
  /*
  //retreive calo tower information from candidates
  //start with the first element of the candidate list
  CandidateCollection::const_iterator tower = towers->begin();
  //get the EDM references to the CaloTowers from the candidate list
  edm::Ref<CaloTowerCollection> towerRef = tower->get<CaloTowerRef>();
  */
  
  /*
  edm::Ref<CaloTowerCollection> towerRef = towers->begin()->get<CaloTowerRef>();
  //finally instantiate now, a list of pointers to the CaloTowers
  const CaloTowerCollection *towerCollection = towerRef.product();
  //iterate over all CaloTowers and record information
  CaloTowerCollection::const_iterator calotower = towerCollection->begin();
  */

  //edm::Ref<CaloTowerCollection> towerRef = towers->begin();
  //finally instantiate now, a list of pointers to the CaloTowers
  //const CaloTowerCollection *towerCollection = towerRef.product();
  //iterate over all CaloTowers and record information
  //CaloTowerCollection::const_iterator calotower = towerCollection->begin();

  edm::View<Candidate>::const_iterator towerCand = towers->begin();


  for( ; towerCand != towers->end(); towerCand++ ) 
       {
	 const Candidate* candidate = &(*towerCand);
	 if (candidate) {
	   const CaloTower* calotower = dynamic_cast<const CaloTower*> (candidate);
	   if (calotower)
	     {
	       totalEt  += calotower->et();
	       totalEm  += calotower->emEt();
	       totalHad += calotower->hadEt();
	       
	       bool hadIsDone = false;
	       bool emIsDone = false;
	       int cell = calotower->constituentsSize();
	       while ( --cell >= 0 && (!hadIsDone || !emIsDone) ) 
		 {
		   DetId id = calotower->constituent( cell );
		   if( !hadIsDone && id.det() == DetId::Hcal ) 
		     {
		       HcalSubdetector subdet = HcalDetId(id).subdet();
		       if( subdet == HcalBarrel || subdet == HcalOuter )
			 {
			   if( calotower->hadEt() > MaxTowerHad ) MaxTowerHad = calotower->hadEt();
			   specific.HadEtInHB   += calotower->hadEt();
			   specific.HadEtInHO   += calotower->outerEt();
			 }
		       else if( subdet == HcalEndcap )
			 {
			   if( calotower->hadEt() > MaxTowerHad ) MaxTowerHad = calotower->hadEt();
			   specific.HadEtInHE   += calotower->hadEt();
			 }
		       else if( subdet == HcalForward )
			 {
			   if (!noHF)
			     {
			       if( calotower->hadEt() > MaxTowerHad ) MaxTowerHad = calotower->hadEt();
			       if( calotower->emEt()  > MaxTowerEm  ) MaxTowerEm  = calotower->emEt();
			     }
			   
			   specific.HadEtInHF   += calotower->hadEt();
			   specific.EmEtInHF    += calotower->emEt(); 
			   
			   if (calotower->eta()>=0)
			     {
			       sumEtInpHF            += calotower->et();
			       MExInpHF -= (calotower->et() * cos(calotower->phi()));
			       MEyInpHF -= (calotower->et() * sin(calotower->phi()));
			     }
			   else
			     {
			       sumEtInmHF            += calotower->et();
			       MExInmHF -= (calotower->et() * cos(calotower->phi()));
			       MEyInmHF -= (calotower->et() * sin(calotower->phi()));
			     }
			 }
		       hadIsDone = true;
		     }
		   else if( !emIsDone && id.det() == DetId::Ecal )
		     {
		       EcalSubdetector subdet = EcalSubdetector( id.subdetId() );
		       if( calotower->emEt()  > MaxTowerEm  ) MaxTowerEm  = calotower->emEt();
		       if( subdet == EcalBarrel )
			 {
			   specific.EmEtInEB    += calotower->emEt(); 
			 }
		       else if( subdet == EcalEndcap ) 
			 {
			   specific.EmEtInEE    += calotower->emEt();
			 }
		       emIsDone = true;
		     }
		 }
	     }
	 }
       }
  
  // Form sub-det specific MET-vectors
  LorentzVector METpHF(MExInpHF, MEyInpHF, 0, sqrt(MExInpHF*MExInpHF + MEyInpHF*MEyInpHF));
  LorentzVector METmHF(MExInmHF, MEyInmHF, 0, sqrt(MExInmHF*MExInmHF + MEyInmHF*MEyInmHF));

  specific.CaloMETInpHF = METpHF.pt();
  specific.CaloMETInmHF = METmHF.pt();

  specific.CaloMETPhiInpHF = METpHF.Phi();
  specific.CaloMETPhiInmHF = METmHF.Phi();

  specific.CaloSETInpHF = sumEtInpHF;
  specific.CaloSETInmHF = sumEtInmHF;
 
  // Instantiate containers for the MET candidate and initialise them with
  // the MET information in "met" (of type CommonMETData)

  // remove HF from MET calculation if specified, but keep HF specific calculations
  if (noHF)
    {
      met.mex   -= (MExInmHF + MExInpHF);
      met.mey   -= (MEyInmHF + MEyInpHF);
      met.sumet -= (sumEtInpHF + sumEtInmHF);
      totalEm   -= specific.EmEtInHF;
      totalHad  -= specific.HadEtInHF;
      totalEt   -= (sumEtInpHF + sumEtInmHF);
      met.met    = sqrt(met.mex*met.mex + met.mey*met.mey);   
 } 

  specific.MaxEtInEmTowers         = MaxTowerEm;  
  specific.MaxEtInHadTowers        = MaxTowerHad;         
  specific.EtFractionHadronic = totalHad / totalEt; 
  specific.EtFractionEm       =  totalEm / totalEt;       

  const LorentzVector p4( met.mex, met.mey, 0.0, met.met );
  const Point vtx( 0.0, 0.0, 0.0 );
  // Create and return an object of type CaloMET, which is a MET object with 
  // the extra calorimeter specfic information added
  CaloMET specificmet( specific, met.sumet, p4, vtx );
  return specificmet;
}
//-------------------------------------------------------------------------
