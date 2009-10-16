///////////////////////////////////////////////////////////////////////////////
// File: CastorSD.cc
// Date: 02.04
// UpDate: 07.04 - C3TF & C4TF semi-trapezoid added
// Description: Sensitive Detector class for Castor
///////////////////////////////////////////////////////////////////////////////

//  Added by WC 
#include "SimG4Core/Notification/interface/TrackInformation.h"

#include "SimG4CMS/Forward/interface/CastorSD.h"
//#include "SimDataFormats/CaloHit/interface/CastorShowerEvent.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

#include "G4ios.hh"
#include "G4Cerenkov.hh"
#include "G4LogicalVolumeStore.hh"

#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "Randomize.hh"
#include "G4Poisson.hh"

#define debugLog

CastorSD::CastorSD(G4String name, const DDCompactView & cpv,
		   SensitiveDetectorCatalog & clg, 
		   edm::ParameterSet const & p, 
		   const SimTrackManager* manager) : 
  CaloSD(name, cpv, clg, p, manager), numberingScheme(0), lvC3EF(0),
  lvC3HF(0), lvC4EF(0), lvC4HF(0) {
  
  edm::ParameterSet m_CastorSD = p.getParameter<edm::ParameterSet>("CastorSD");
  useShowerLibrary  = m_CastorSD.getParameter<bool>("useShowerLibrary");                  // Add useSh... to PSet (WC)
  energyThresholdSL = m_CastorSD.getParameter<double>("minEnergyInGeVforUsingSLibrary");  // Add useSh... to PSet (WC)
  energyThresholdSL = energyThresholdSL*GeV;   //  Convert GeV => MeV 
  
  if (useShowerLibrary) showerLibrary = new CastorShowerLibrary(name, p);
  
  setNumberingScheme(new CastorNumberingScheme());
  
  edm::LogInfo("ForwardSim") 
    << "***************************************************\n"
    << "*                                                 *\n" 
    << "* Constructing a CastorSD  with name " << GetName() << "\n"
    << "*                                                 *\n"
    << "***************************************************";
    
  const G4LogicalVolumeStore * lvs = G4LogicalVolumeStore::GetInstance();
  std::vector<G4LogicalVolume*>::const_iterator lvcite;
  for (lvcite = lvs->begin(); lvcite != lvs->end(); lvcite++) {
    if (strcmp(((*lvcite)->GetName()).c_str(),"C3EF") == 0) lvC3EF = (*lvcite);
    if (strcmp(((*lvcite)->GetName()).c_str(),"C3HF") == 0) lvC3HF = (*lvcite);
    if (strcmp(((*lvcite)->GetName()).c_str(),"C4EF") == 0) lvC4EF = (*lvcite);
    if (strcmp(((*lvcite)->GetName()).c_str(),"C4HF") == 0) lvC4HF = (*lvcite);
    if (lvC3EF != 0 && lvC3HF != 0 && lvC4EF != 0 && lvC4HF != 0) break;
  }
  edm::LogInfo("ForwardSim") << "CastorSD:: LogicalVolume pointers\n"
			     << lvC3EF << " for C3EF; " << lvC3HF 
			     << " for C3HF; " << lvC4EF << " for C4EF; " 
			     << lvC4HF << " for C4HF.";

  //  if(useShowerLibrary) edm::LogInfo("ForwardSim") << "\n Using Castor Shower Library \n";

}

//=============================================================================================

CastorSD::~CastorSD() {
  if (useShowerLibrary) delete showerLibrary;
}

//=============================================================================================

void CastorSD::initRun(){
  if (useShowerLibrary) {
    // showerLibrary = new CastorShowerLibrary(name, cpv, p);
    G4ParticleTable * theParticleTable = G4ParticleTable::GetParticleTable();
    showerLibrary->initParticleTable(theParticleTable);
    edm::LogInfo("ForwardSim") << "CastorSD::initRun: Using Castor Shower Library \n";
  }
}

//=============================================================================================

double CastorSD::getEnergyDeposit(G4Step * aStep) {
  
  float NCherPhot = 0.;

  if (aStep == NULL) {
    return 0;
  } else {
    // preStepPoint information *********************************************
  
    G4SteppingControl  stepControlFlag = aStep->GetControlFlag();
    G4StepPoint*       preStepPoint = aStep->GetPreStepPoint();
    G4VPhysicalVolume* currentPV    = preStepPoint->GetPhysicalVolume();
    G4LogicalVolume*   currentLV    = currentPV->GetLogicalVolume();
    G4String           name   = currentPV->GetName();
    std::string        nameVolume;
    nameVolume.assign(name,0,4);
    
#ifdef DebugLog
    if (aStep->IsFirstStepInVolume()) 
      LogDebug("ForwardSim") << "CastorSD::getEnergyDeposit:"
			     << "\n IsFirstStepInVolume " ; 
#endif
    
    // Get theTrack 
    G4Track*        theTrack = aStep->GetTrack();
    
#ifdef DebugLog
    if (useShowerLibrary) {
      edm::LogInfo("ForwardSim") << "CastorSD::getEnergyDeposit:"
				 << "\n TrackID , ParentID , ParticleName ,"
				 << " eta , phi , z , time ,"
				 << " K , E , Mom " 
				 << "\n  TRACKINFO: " 
				 << theTrack->GetTrackID() 
				 << " , " 
				 << theTrack->GetParentID() 
				 << " , "
				 << theTrack->GetDefinition()->GetParticleName() 
				 << " , "
				 << theTrack->GetPosition().eta() 
				 << " , "
				 << theTrack->GetPosition().phi() 
				 << " , "
				 << theTrack->GetPosition().z() 
				 << " , "
				 << theTrack->GetGlobalTime() 
				 << " , "
				 << theTrack->GetKineticEnergy() 
				 << " , "
				 << theTrack->GetTotalEnergy() 
				 << " , "
                                 << theTrack->GetMomentum().mag() ;
      if(theTrack->GetTrackID() != 1) 
	edm::LogInfo("ForwardSim") << "CastorSD::getEnergyDeposit:"
				   << "\n CurrentStepNumber , VertexPosition ,"
				   << " LogicalVolumeAtVertex , CreatorProcess"
				   << "\n  TRACKINFO2: " 
				   << theTrack->GetCurrentStepNumber() 
				   << " , " 
				   << theTrack->GetVertexPosition() 
				   << " , "
				   << theTrack->GetLogicalVolumeAtVertex()->GetName() 
				   << " , " 
				   << theTrack->GetCreatorProcess()->GetProcessName() ;
    } // end of if(useShowerLibrary)
#endif
    
    // if particle moves from interaction point or "backwards (halo)
    int backward = 0;
    G4ThreeVector  hitPoint = preStepPoint->GetPosition();	
    G4ThreeVector  hit_mom  = preStepPoint->GetMomentumDirection();
    double zint = hitPoint.z();
    double pz   = hit_mom.z();
    // Check if theTrack moves backward 
    if (pz * zint < 0.) backward = 1;
    
    // Check that theTrack is above the energy threshold to use Shower Library 
    bool aboveThreshold = false;
    if(theTrack->GetKineticEnergy() > energyThresholdSL) aboveThreshold = true;
    
    // Check if theTrack is a muon (if so, DO NOT use Shower Library) 
    bool notaMuon = true;
    G4int mumPDG  =  13;
    G4int mupPDG  = -13;
    G4int parCode = theTrack->GetDefinition()->GetPDGEncoding();
    if (parCode == mupPDG || parCode == mumPDG ) notaMuon = false;
    
    if(useShowerLibrary && aboveThreshold && notaMuon && (!backward)) {
      // Use Castor shower library if energy is above threshold, is not a muon 
      // and is not moving backward 
      getFromLibrary(aStep);
    } else {
      // Usual calculations
      // G4ThreeVector      hitPoint = preStepPoint->GetPosition();	
      // G4ThreeVector      hit_mom = preStepPoint->GetMomentumDirection();
      G4double           stepl    = aStep->GetStepLength()/cm;
      G4double           beta     = preStepPoint->GetBeta();
      G4double           charge   = preStepPoint->GetCharge();
      //        G4VProcess*        curprocess   = preStepPoint->GetProcessDefinedStep();
      //        G4String           namePr   = preStepPoint->GetProcessDefinedStep()->GetProcessName();
      //        std::string nameProcess;
      //        nameProcess.assign(namePr,0,4);

      //        G4LogicalVolume*   lv    = currentPV->GetLogicalVolume();
      //        G4Material*        mat   = lv->GetMaterial();
      //        G4double           rad   = mat->GetRadlen();


      // postStepPoint information *********************************************
      G4StepPoint* postStepPoint= aStep->GetPostStepPoint();   
      G4VPhysicalVolume* postPV    = postStepPoint->GetPhysicalVolume();
      G4String           postname   = postPV->GetName();
      std::string        postnameVolume;
      postnameVolume.assign(postname,0,4);

      // theTrack information  *************************************************
      // G4Track*        theTrack = aStep->GetTrack();   
      G4double        entot    = theTrack->GetTotalEnergy();
      G4ThreeVector   vert_mom = theTrack->GetVertexMomentumDirection();

      G4ThreeVector  localPoint = theTrack->GetTouchable()->GetHistory()->
	GetTopTransform().TransformPoint(hitPoint);

      // calculations...       *************************************************

      float costheta =vert_mom.z()/sqrt(vert_mom.x()*vert_mom.x()+
					vert_mom.y()*vert_mom.y()+
					vert_mom.z()*vert_mom.z());
      float theta = acos(std::min(std::max(costheta,float(-1.)),float(1.)));
      float eta = -log(tan(theta/2));
      float phi = -100.;
      if (vert_mom.x() != 0) phi = atan2(vert_mom.y(),vert_mom.x()); 
      if (phi < 0.) phi += twopi;
      G4String       particleType = theTrack->GetDefinition()->GetParticleName();
      G4int          primaryID    = theTrack->GetTrackID();
      // *************************************************


      // *************************************************
      double edep   = aStep->GetTotalEnergyDeposit();

      // *************************************************


      // *************************************************
      // take into account light collection curve for plate
      //      double weight = curve_Castor(nameVolume, preStepPoint);
      //      double edep   = aStep->GetTotalEnergyDeposit() * weight;
      // *************************************************


      // *************************************************
      /*    comments for sensitive volumes:      
	    C001 ...-... CP06,CBU1 ...-...CALM --- > fibres and bundle 
	    for first release of CASTOR
	    CASF  --- > quartz plate  for first and second releases of CASTOR  
	    GF2Q, GFNQ, GR2Q, GRNQ 
	    for tests with my own test geometry of HF (on ask of Gavrilov)
	    C3TF, C4TF - for third release of CASTOR
      */
      double meanNCherPhot=0;

      if (currentLV == lvC3EF || currentLV == lvC4EF || currentLV == lvC3HF ||
	  currentLV == lvC4HF) {
	//      if(nameVolume == "C3EF" || nameVolume == "C4EF" || nameVolume == "C3HF" || nameVolume == "C4HF") {

	float bThreshold = 0.67;
	float nMedium = 1.4925;
	//     float photEnSpectrDL = (1./400.nm-1./700.nm)*10000000.cm/nm; /* cm-1  */
	//     float photEnSpectrDL = 10714.285714;

	float photEnSpectrDE = 1.24;    /* see below   */
	/*     E = 2pi*(1./137.)*(eV*cm/370.)/lambda  =     */
	/*       = 12.389184*(eV*cm)/lambda                 */
	/*     Emax = 12.389184*(eV*cm)/400nm*10-7cm/nm  = 3.01 eV     */
	/*     Emin = 12.389184*(eV*cm)/700nm*10-7cm/nm  = 1.77 eV     */
	/*     delE = Emax - Emin = 1.24 eV  --> */
	/*   */
	/* default for Castor nameVolume  == "CASF" or (C3TF & C4TF)  */
	float effPMTandTransport = 0.15;
	/* for test HF geometry volumes:   
	   if(nameVolume == "GF2Q" || nameVolume == "GFNQ" || 
	   nameVolume == "GR2Q" || nameVolume == "GRNQ")
	   effPMTandTransport = 0.15;
	*/

	float thFullRefl = 23.;  /* 23.dergee */
	float thFullReflRad = thFullRefl*pi/180.;

	/* default for Castor nameVolume  == "CASF" or (C3TF & C4TF)  */
	float thFibDir = 45.;  /* .dergee */
	/* for test HF geometry volumes:   
	   if(nameVolume == "GF2Q" || nameVolume == "GFNQ" ||
	   nameVolume == "GR2Q" || nameVolume == "GRNQ")
	   thFibDir = 0.0; // .dergee
	*/
	float thFibDirRad = thFibDir*pi/180.;
	/*   */
	/*   */

	// at which theta the point is located:
	//     float th1    = hitPoint.theta();

	// theta of charged particle in LabRF(hit momentum direction):
	float costh =hit_mom.z()/sqrt(hit_mom.x()*hit_mom.x()+
				      hit_mom.y()*hit_mom.y()+
				      hit_mom.z()*hit_mom.z());
	if (zint < 0) costh = -costh;
	float th = acos(std::min(std::max(costh,float(-1.)),float(1.)));

	// just in case (can do bot use):
	if (th < 0.) th += twopi;

	float thgrad = th*180./pi;


	// theta of cone with Cherenkov photons w.r.t.direction of charged part.:
	float costhcher =1./(nMedium*beta);
	float thcher = acos(std::min(std::max(costhcher,float(-1.)),float(1.)));
	float thchergrad = thcher*180./pi;

	// diff thetas of charged part. and quartz direction in LabRF:
	float DelFibPart = fabs(th - thFibDirRad);
	float DelFibPartgrad = DelFibPart*180./pi;



	// define real distances:
	float d = fabs(tan(th)-tan(thFibDirRad));   

	//       float a = fabs(tan(thFibDirRad)-tan(thFibDirRad+thFullReflRad));   
	//       float r = fabs(tan(th)-tan(th+thcher));   
	
	float a = tan(thFibDirRad)+tan(fabs(thFibDirRad-thFullReflRad));   
	float r = tan(th)+tan(fabs(th-thcher));   


	// define losses d_qz in cone of full reflection inside quartz direction
	float d_qz;
	float variant;

	if(DelFibPart > (thFullReflRad + thcher) ) {d_qz = 0.; variant=0.;}
	// if(d > (r+a) ) {d_qz = 0.; variant=0.;}
	else {
	  if((th + thcher) < (thFibDirRad+thFullReflRad) && 
	     (th - thcher) > (thFibDirRad-thFullReflRad) 
	     ) {d_qz = 1.; variant=1.;}
	  // if((DelFibPart + thcher) < thFullReflRad ) {d_qz = 1.; variant=1.;}
	  // if((d+r) < a ) {d_qz = 1.; variant=1.;}
	  else {
	    if((thFibDirRad + thFullReflRad) < (th + thcher) && 
	       (thFibDirRad - thFullReflRad) > (th - thcher) ) {
	      // if((thcher - DelFibPart ) > thFullReflRad ) 
	      // if((r-d ) > a ) 
	      d_qz = 0.; variant=2.;

	    } else {
	      // if((thcher + DelFibPart ) > thFullReflRad && 
	      //    thcher < (DelFibPart+thFullReflRad) ) 
	      //      {
	      d_qz = 0.; variant=3.;


	      // use crossed length of circles(cone projection)
	      // dC1/dC2 : 
	      float arg_arcos = 0.;
	      float tan_arcos = 2.*a*d;
	      if(tan_arcos != 0.) arg_arcos =(r*r-a*a-d*d)/tan_arcos; 
	      arg_arcos = fabs(arg_arcos);
	      float th_arcos = acos(std::min(std::max(arg_arcos,float(-1.)),float(1.)));
	      d_qz = th_arcos/pi/2.;
	      d_qz = fabs(d_qz);



	      //	    }
	      //             else
	      //  	     {
	      //                d_qz = 0.; variant=4.;
	      //#ifdef debugLog
	      // std::cout <<" ===============>variant 4 information: <===== " <<std::endl;
	      // std::cout <<" !!!!!!!!!!!!!!!!!!!!!!  variant = " << variant  <<std::endl;
	      //#endif 
	      //
	      // 	     }
	    }
	  }
	}


	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if(charge != 0. && beta > bThreshold && d_qz != 0. ) {

	  meanNCherPhot = 370.*charge*charge* 
	    ( 1. - 1./(nMedium*nMedium*beta*beta) )*
	    photEnSpectrDE*stepl;

	  // dLamdX:
	  //   meanNCherPhot = (2.*pi/137.)*charge*charge* 
	  //                     ( 1. - 1./(nMedium*nMedium*beta*beta) )*
	  //                     photEnSpectrDL*stepl;


	  //     NCherPhot = meanNCherPhot;
	  // Poisson:
	  G4int poissNCherPhot = (G4int) G4Poisson(meanNCherPhot);

	  if(poissNCherPhot < 0) poissNCherPhot = 0; 
	  
	  NCherPhot = poissNCherPhot * effPMTandTransport * d_qz;



#ifdef debugLog
	  LogDebug("ForwardSim") << " ==============================> start all "
				 << "information:<========= \n" << " =====> for "
				 << "test:<===  \n" << " variant = " << variant  
				 << "\n thgrad = " << thgrad  << "\n thchergrad "
				 << "= " << thchergrad  << "\n DelFibPartgrad = "
				 << DelFibPartgrad << "\n d_qz = " << d_qz  
				 << "\n =====> Start Step Information <===  \n"
				 << " ===> calo preStepPoint info <===  \n" 
				 << " hitPoint = " << hitPoint  << "\n"
				 << " hitMom = " << hit_mom  << "\n"
				 << " stepControlFlag = " << stepControlFlag 
	    // << "\n curprocess = " << curprocess << "\n"
	    // << " nameProcess = " << nameProcess 
				 << "\n charge = " << charge << "\n"
				 << " beta = " << beta << "\n"
				 << " bThreshold = " << bThreshold << "\n"
				 << " thgrad =" << thgrad << "\n"
				 << " effPMTandTransport=" << effPMTandTransport 
	    // << "\n volume = " << name 
				 << "\n nameVolume = " << nameVolume << "\n"
				 << " nMedium = " << nMedium << "\n"
	    //  << " rad length = " << rad << "\n"
	    //  << " material = " << mat << "\n"
				 << " stepl = " << stepl << "\n"
				 << " photEnSpectrDE = " << photEnSpectrDE <<"\n"
				 << " edep = " << edep << "\n"
				 << " ===> calo theTrack info <=== " << "\n"
				 << " particleType = " << particleType << "\n"
				 << " primaryID = " << primaryID << "\n"
				 << " entot= " << entot << "\n"
				 << " vert_eta= " << eta  << "\n"
				 << " vert_phi= " << phi << "\n"
				 << " vert_mom= " << vert_mom  << "\n"
				 << " ===> calo hit preStepPointinfo <=== "<<"\n"
				 << " local point = " << localPoint << "\n"
				 << " ==============================> final info"
				 << ":  <=== \n" 
				 << " meanNCherPhot = " << meanNCherPhot << "\n"
				 << " poissNCherPhot = " << poissNCherPhot <<"\n"
				 << " NCherPhot = " << NCherPhot;
#endif 
	
	  // Included by WC
	  //	     std::cout << "\n volume = "         << name 
	  //	          << "\n nameVolume = "     << nameVolume << "\n"
	  //	          << "\n postvolume = "     << postname 
	  //	          << "\n postnameVolume = " << postnameVolume << "\n"
	  //	          << "\n particleType = "   << particleType 
	  //	          << "\n primaryID = "      << primaryID << "\n";

	}
      }


#ifdef debugLog
      LogDebug("ForwardSim") << "CastorSD:: " << nameVolume 
	//      << " Light Collection Efficiency " << weight
			     << " Weighted Energy Deposit " << edep/MeV 
			     << " MeV\n";
#endif
      // Temporary member for testing purpose only...
      // unit_id = setDetUnitId(aStep);
      // if(NCherPhot != 0) std::cout << "\n  UnitID = " << unit_id << "  ;  NCherPhot = " << NCherPhot ;
       
      return NCherPhot;
    }
  } 
  return 0;
}

//=======================================================================================

uint32_t CastorSD::setDetUnitId(G4Step* aStep) {
  return (numberingScheme == 0 ? 0 : numberingScheme->getUnitID(aStep));
}

//=======================================================================================

void CastorSD::setNumberingScheme(CastorNumberingScheme* scheme) {

  if (scheme != 0) {
    edm::LogInfo("ForwardSim") << "CastorSD: updates numbering scheme for " 
			       << GetName();
    if (numberingScheme) delete numberingScheme;
    numberingScheme = scheme;
  }
}

//=======================================================================================

int CastorSD::setTrackID (G4Step* aStep) {

  theTrack     = aStep->GetTrack();

  double etrack = preStepPoint->GetKineticEnergy();
  TrackInformation * trkInfo = (TrackInformation *)(theTrack->GetUserInformation());
  int      primaryID = trkInfo->getIDonCaloSurface();
  if (primaryID == 0) {
#ifdef DebugLog
    LogDebug("ForwardSim") << "CastorSD: Problem with primaryID **** set by force "
			   << "to TkID **** " << theTrack->GetTrackID();
#endif
    primaryID = theTrack->GetTrackID();
  }

  if (primaryID != previousID.trackID())
    resetForNewPrimary(preStepPoint->GetPosition(), etrack);

  return primaryID;
}

//=======================================================================================

uint32_t CastorSD::rotateUnitID(uint32_t unitID, G4Track* track, CastorShowerEvent shower) {
// ==============================================================
//
//   o   Exploit Castor phi symmetry to return newUnitID for  
//       shower hits based on track phi coordinate 
//
// ==============================================================
  
  const float pi = 3.141592654;
  
  // Get 'track' phi:
  float   trackPhi = theTrack->GetPosition().phi(); 
  // Get phi from primary that gave rise to SL 'shower':
  float  showerPhi = shower.getPrimPhi(); 
  // Delta phi:
  
  //  Find the HexSector for which 'track' and 'shower' belong
  int  trackHexSector = (int) (  trackPhi / (pi/8) ) ;
  int showerHexSector = (int) ( showerPhi / (pi/8) ) ;
  
  uint32_t  newUnitID;
  uint32_t         sec = ( ( unitID>>4 ) & 0xF ) ;
  uint32_t  complement = ( unitID & 0xFFFFFF0F ) ;
  
  // edm::LogInfo("ForwardSim") << "\n CastorSD::rotateUnitID:  " 
  //                            << "\n      unitID = " << unitID 
  //                            << "\n         sec = " << sec 
  //                            << "\n  complement = " << complement ; 
  
  int aux ;
  int dSec = trackHexSector - showerHexSector ;
  if( dSec<0 ) {
     sec += 16 ;
     sec += dSec ;
     aux  = (int) (sec/16) ;
     sec -= aux*16 ;
     sec  = sec<<4 ;
     newUnitID = complement | sec ;
  } else {
     sec += dSec ;
     aux  = (int) (sec/16) ;
     sec -= aux*16 ;
     sec  = sec<<4 ;
     newUnitID = complement | sec ;
  }
  
#ifdef DebugLog
  if(newUnitID != unitID) {
    LogDebug("ForwardSim") << "\n CastorSD::rotateUnitID:  " 
			   << "\n     unitID = " << unitID 
			   << "\n  newUnitID = " << newUnitID ; 
  }
#endif
  
  return newUnitID ;

}

//=======================================================================================

void CastorSD::getFromLibrary (G4Step* aStep) {

/////////////////////////////////////////////////////////////////////
//
//   Method to get hits from the Shower Library
//
//   CastorShowerEvent hits returned by getShowerHits are used to  
//   replace the full simulation of the shower from theTrack
//    
//   "updateHit" save the Hits to a CaloG4Hit container
//
/////////////////////////////////////////////////////////////////////

  preStepPoint  = aStep->GetPreStepPoint(); 
  theTrack      = aStep->GetTrack();   
  bool ok;
  
  // ****    Call method to retrieve hits from the ShowerLibrary   ****
  CastorShowerEvent hits = showerLibrary->getShowerHits(aStep, ok);

  double etrack    = preStepPoint->GetKineticEnergy();
  int    primaryID = setTrackID(aStep);
  // int    primaryID = theTrack->GetTrackID();

  // Reset entry point for new primary
  posGlobal = preStepPoint->GetPosition();
  resetForNewPrimary(posGlobal, etrack);

  // Check whether track is EM or HAD
  G4int particleCode = theTrack->GetDefinition()->GetPDGEncoding();
  bool isEM , isHAD ;
  if (particleCode==emPDG || particleCode==epPDG || particleCode==gammaPDG) {
    isEM = true ; isHAD = false;
  } else {
    isEM = false; isHAD = true ;
  }

#ifdef DebugLog
  //  edm::LogInfo("ForwardSim") << "\n CastorSD::getFromLibrary:  " 
  LogDebug("ForwardSim") << "\n CastorSD::getFromLibrary:  " 
			 << hits.getNhit() << " hits for " << GetName() << " from " 
			 << theTrack->GetDefinition()->GetParticleName() << " of "
			 << preStepPoint->GetKineticEnergy()/GeV << " GeV and trackID " 
			 << theTrack->GetTrackID()  ;
#endif

  //  Loop over hits retrieved from the library
  for (unsigned int i=0; i<hits.getNhit(); i++) {
    
    // Get nPhotoElectrons and set edepositEM / edepositHAD accordingly
    double nPhotoElectrons    = hits.getNphotons(i);  
    if(isEM)  {
       // edepositEM  = nPhotoElectrons*GeV; 
       edepositEM  = nPhotoElectrons; 
       edepositHAD = 0.;                 
    } else if(isHAD) {
       edepositEM  = 0.;                  
       edepositHAD = nPhotoElectrons;
       // edepositHAD = nPhotoElectrons*GeV;
    }
    
    // Get hit position and time
    double                time = hits.getTime(i);
    math::XYZPoint    position = hits.getHitPosition(i);
    
    // Get hit detID
    unsigned int        unitID = hits.getDetID(i);
    
    // Make the detID "rotation" from one sector to another taking into account the 
    // sectors of the impinging particle (theTrack) and of the particle that produced 
    // the 'hits' retrieved from shower library   
    unsigned int rotatedUnitID = rotateUnitID(unitID , theTrack , hits);
    currentID.setID(rotatedUnitID, time, primaryID, 0);
    // currentID.setID(unitID, time, primaryID, 0);
   
    // check if it is in the same unit and timeslice as the previous one
    if (currentID == previousID) {
      updateHit(currentHit);
    } else {
      if (!checkHit()) currentHit = createNewHit();
    }
  }  //  End of loop over hits

  //Now kill the current track
  if (ok) {
    theTrack->SetTrackStatus(fStopAndKill);
#ifdef DebugLog
    LogDebug("ForwardSim") << "CastorSD::getFromLibrary:"
			   << "\n \"theTrack\" with TrackID() = " 
			   << theTrack->GetTrackID() 
			   << " and with energy " 
			   << theTrack->GetTotalEnergy()
			   << " has been set to be killed" ;
#endif
    G4TrackVector tv = *(aStep->GetSecondary());
    for (unsigned int kk=0; kk<tv.size(); kk++) {
      if (tv[kk]->GetVolume() == preStepPoint->GetPhysicalVolume()) {
	tv[kk]->SetTrackStatus(fStopAndKill);
#ifdef DebugLog
	LogDebug("ForwardSim") << "CastorSD::getFromLibrary:"
			       << "\n tv[" << kk << "]->GetTrackID() = " 
			       << tv[kk]->GetTrackID() 
			       << " with energy " 
			       << tv[kk]->GetTotalEnergy()
			       << " has been set to be killed" ;
#endif
      }
    }
  }
}
