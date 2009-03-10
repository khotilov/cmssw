///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// modified by P. Biallass 29.03.2006 to implement new cosmic generator (CMSCGEN.cc) and new normalization of flux (CMSCGENnorm.cc)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 04.12.2008 sonne: replaced Min/MaxE by Min/MaxP to get cos_sf/ug scripts working again

#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/JamesRandom.h>

#include "GeneratorInterface/CosmicMuonGenerator/interface/CosmicMuonGenerator.h"

CosmicMuonGenerator::CosmicMuonGenerator() : delRanGen(false)
{
  //initialize class which normalizes flux (added by P.Biallass 29.3.2006)
  Norm = new CMSCGENnorm();
  //initialize class which produces the cosmic muons  (modified by P.Biallass 29.3.2006)
  Cosmics = new CMSCGEN();
  // set default control parameters
  NumberOfEvents = 100;
  RanSeed = 135799468;
  MinP =     3.;
  MinP_CMS =     MinP;
  MaxP =   3000.;
  MinTheta =  0.*Deg2Rad;
  MaxTheta = 84.26*Deg2Rad;
  MinPhi =    0.*Deg2Rad;
  MaxPhi =  360.*Deg2Rad;
  MinT0  = -12.5;
  MaxT0  =  12.5;
  ElossScaleFactor = 1.0;
  RadiusOfTarget = 8000.;
  ZDistOfTarget = 15000.;
  ZCentrOfTarget = 0.;
  TrackerOnly = false;
  TIFOnly_constant = false;
  TIFOnly_linear = false;
  MTCCHalf = false;
  EventRate = 0.;
  rateErr_stat = 0.;
  rateErr_syst = 0.;

  SumIntegrals = 0.;
  Ngen = 0.;
  Nsel = 0.;
  Ndiced = 0.;
  NotInitialized = true;
  Target3dRadius = 0.;
  SurfaceRadius = 0.;
  //set plug as default onto PX56 shaft
  PlugVx = PlugOnShaftVx;
  PlugVz = PlugOnShaftVz;

  std::cout << std::endl;
  std::cout << "*********************************************************" << std::endl;
  std::cout << "*********************************************************" << std::endl;
  std::cout << "***                                                   ***" << std::endl;
  std::cout << "***  C O S M I C  M U O N  G E N E R A T O R  (vC++)  ***" << std::endl;
  std::cout << "***                                                   ***" << std::endl;
  std::cout << "*********************************************************" << std::endl;
  std::cout << "*********************************************************" << std::endl;
  std::cout << std::endl;
}

CosmicMuonGenerator::~CosmicMuonGenerator()
{
  if (delRanGen)
    delete RanGen;
  delete Norm; 
  delete Cosmics;
}

void CosmicMuonGenerator::runCMG(){
  initialize();
  for (unsigned int iGen=0; iGen<NumberOfEvents; ++iGen){ nextEvent(); }
  terminate();
}

void CosmicMuonGenerator::initialize(CLHEP::HepRandomEngine *rng){
  if (delRanGen)
    delete RanGen;
  if (!rng) {
    RanGen = new CLHEP::HepJamesRandom;
    RanGen->setSeed(RanSeed, 0); //set seed for Random Generator (seed can be controled by config-file)
    delRanGen = true;
  } else {
    RanGen = rng;
    delRanGen = false;
  }
  checkIn();
  if (NumberOfEvents > 0){
    // set up "surface geometry" dimensions
    double RadiusTargetEff = RadiusOfTarget; //get this from cfg-file
    double Z_DistTargetEff = ZDistOfTarget;  //get this from cfg-file
    double Z_CentrTargetEff = ZCentrOfTarget;  //get this from cfg-file
    if(TrackerOnly==true){
    RadiusTargetEff = RadiusTracker;
    Z_DistTargetEff = Z_DistTracker;
    }
    Target3dRadius = sqrt(RadiusTargetEff*RadiusTargetEff + Z_DistTargetEff*Z_DistTargetEff) + MinStepSize;
    if (Debug) std::cout << "  radius of sphere  around  target = " << Target3dRadius << " mm" << std::endl;

    if (MinTheta > 90.*Deg2Rad) //upgoing muons from neutrinos
      SurfaceRadius = (RadiusCMS)*(-tan(MinTheta)) + MinStepSize;
    else
      SurfaceRadius = (SurfaceOfEarth+PlugWidth+RadiusTargetEff)*tan(MaxTheta) + Target3dRadius;  
    if (Debug) std::cout << "  starting point radius at Surface + PlugWidth = " << SurfaceRadius << " mm" << std::endl;

    OneMuoEvt.PlugVx = PlugVx;
    OneMuoEvt.PlugVz = PlugVz;
    //set energy and angle limits for CMSCGEN, give same seed as above 
    if (MinTheta >= 90.*Deg2Rad) //upgoing muons from neutrinos
      Cosmics->initializeNuMu(MinP, MaxP, MinTheta, MaxTheta, MinEnu, MaxEnu, 
			      MinPhi, MaxPhi, RanSeed);
    else 
      Cosmics->initialize(MinP, MaxP, MinTheta, MaxTheta, RanSeed, TIFOnly_constant, TIFOnly_linear);
   
#if ROOT_INTERACTIVE
  // book histos
  TH1D* ene = new TH1D("ene","generated energy",210,0.,1050.);
  TH1D* the = new TH1D("the","generated theta",90,0.,90.);
  TH1D* phi = new TH1D("phi","generated phi",120,0.,360.);
  TH3F* ver = new TH3F("ver","Z-X-Y coordinates",50,-25.,25.,20,-10.,10.,20,-10.,10.);
#endif
    if (EventDisplay) initEvDis();
    std::cout << std::endl;
    std::cout << "  generating " << NumberOfEvents << " events with random seed " << RanSeed << std::endl;
    NotInitialized = false;
  }
}

void CosmicMuonGenerator::nextEvent(){

  double E = 0.; double Theta = 0.; double Phi = 0.; double RxzV = 0.; double PhiV = 0.;
  if (int(Nsel)%100 == 0) std::cout << "    generated " << int(Nsel) << " events" << std::endl;
  // generate cosmic (E,theta,phi)
  bool   notSelected = true;
  while (notSelected){
	bool   badMomentumGenerated = true;
	while (badMomentumGenerated){

	  if (MinTheta > 90.*Deg2Rad) //upgoing muons from neutrinos
	    Cosmics->generateNuMu();
	  else
	    Cosmics->generate(); //dice one event now
	  
	  E = sqrt(Cosmics->momentum_times_charge()*Cosmics->momentum_times_charge() + MuonMass*MuonMass);
	  Theta = TMath::ACos( Cosmics->cos_theta() ) ; //angle has to be in RAD here
	  Ngen+=1.;   //count number of initial cosmic events (in surface area), vertices will be added later
	  badMomentumGenerated = false;
	  Phi = RanGen->flat()*(MaxPhi-MinPhi) + MinPhi;
	}
	Norm->events_n100cos(E, Theta); //test if this muon is in normalization range
	Ndiced += 1; //one more cosmic is diced
	//std::cout << "diced: E=" << E << " Theta=" << Theta << " Phi=" << Phi << std::endl;
    // generate vertex
    double Nver = 0.;
    bool   badVertexGenerated = true;
    while (badVertexGenerated){
      RxzV = sqrt(RanGen->flat())*SurfaceRadius;
      PhiV = RanGen->flat()*TwoPi;
      // check phi range (for a sphere with Target3dRadius around the target)
      double dPhi = Pi; if (RxzV > Target3dRadius) dPhi = asin(Target3dRadius/RxzV);
      double rotPhi = PhiV + Pi; if (rotPhi > TwoPi) rotPhi -= TwoPi;
      double disPhi = fabs(rotPhi - Phi); if (disPhi > Pi) disPhi = TwoPi - disPhi;
      if (disPhi < dPhi) badVertexGenerated = false;
      Nver+=1.;
    }
    Ngen += (Nver-1.); //add number of generated vertices to initial cosmic events
    
    // complete event at surface
    int                             id =  13; // mu-
    if (Cosmics->momentum_times_charge() >0.) id = -13; // mu+
    double absMom = sqrt(E*E - MuonMass*MuonMass);
    double verMom = absMom*cos(Theta);
    double horMom = absMom*sin(Theta);
    double Px = horMom*sin(Phi); // [GeV/c]
    double Py = -verMom;         // [GeV/c]
    double Pz = horMom*cos(Phi); // [GeV/c]
    double Vx = RxzV*sin(PhiV);  // [mm]

    double Vy;
    if (MinTheta > 90.*Deg2Rad) //upgoing muons from neutrinos
      Vy = -RadiusCMS;
    else
      Vy = SurfaceOfEarth + PlugWidth;  // [mm]

    double Vz = RxzV*cos(PhiV);  // [mm]
    double T0 = (RanGen->flat()*(MaxT0-MinT0) + MinT0)*SpeedOfLight; // [mm/c];
    //std::cout << "Vx=" << Vx << " Vy=" << Vy << " Vz=" << Vz << std::endl;
    OneMuoEvt.create(id, Px, Py, Pz, E, MuonMass, Vx, Vy, Vz, T0); 
    // if angles are ok, propagate to target
    if (goodOrientation()) { 
      if (MinTheta > 90.*Deg2Rad) //upgoing muons from neutrinos
	OneMuoEvt.propagate(0., RadiusOfTarget, ZDistOfTarget, ZCentrOfTarget, TrackerOnly, MTCCHalf);
      else
	OneMuoEvt.propagate(ElossScaleFactor, RadiusOfTarget, ZDistOfTarget, ZCentrOfTarget, TrackerOnly, MTCCHalf);
    }
      // if cosmic hits target test also if P>Pmin_CMS; the default is MinP_surface=MinP_CMS, thus no bias from access shaft

    //std::cout << "OneMuoEvt.hitTarget()=";
    //if (OneMuoEvt.hitTarget()) std::cout << "true";
    //else std::cout << "false";
    //std::cout << " P(Mu)=" << sqrt(OneMuoEvt.e()*OneMuoEvt.e() - MuonMass*MuonMass)
    //      << " MinP_CMS=" << MinP_CMS << std::endl;

    if (OneMuoEvt.hitTarget() && sqrt(OneMuoEvt.e()*OneMuoEvt.e() - MuonMass*MuonMass) > MinP_CMS){
      Nsel+=1.; //count number of generated and accepted events  
      notSelected = false;
    }
  }
  // plot variables of selected events
#if ROOT_INTERACTIVE
  ene->Fill(OneMuoEvt.e());
  the->Fill((OneMuoEvt.theta()*Rad2Deg));
  phi->Fill((OneMuoEvt.phi()*Rad2Deg));
  ver->Fill((OneMuoEvt.vz()/1000.),(OneMuoEvt.vx()/1000.),(OneMuoEvt.vy()/1000.));
#endif
  if (Debug){
    std::cout << "new event" << std::endl;
    std::cout << "  Px,Py,Pz,E,m = " << OneMuoEvt.px() << ", " << OneMuoEvt.py() << ", "
         << OneMuoEvt.pz() << ", " << OneMuoEvt.e() << ", " << OneMuoEvt.m() << " GeV" << std::endl;
    std::cout << "  Vx,Vy,Vz,t0  = " << OneMuoEvt.vx() << ", " << OneMuoEvt.vy() << ", " 
         << OneMuoEvt.vz() << ", " << OneMuoEvt.t0() << " mm" << std::endl;
  }
  if (EventDisplay) displayEv();
  
}

void CosmicMuonGenerator::terminate(){
  if (NumberOfEvents > 0){
    std::cout << std::endl;
    std::cout << "*********************************************************" << std::endl;
    std::cout << "*********************************************************" << std::endl;
    std::cout << "***                                                   ***" << std::endl;
    std::cout << "***    C O S M I C   M U O N   S T A T I S T I C S    ***" << std::endl;
    std::cout << "***                                                   ***" << std::endl;
    std::cout << "*********************************************************" << std::endl;
    std::cout << "*********************************************************" << std::endl;
    std::cout << std::endl;  
    std::cout << "       number of initial cosmic events:  " << int(Ngen) << std::endl;
    std::cout << "       number of actually diced events:  " << int(Ndiced) << std::endl;
    std::cout << "       number of generated and accepted events:  " << int(Nsel) << std::endl;
    double selEff = Nsel/Ngen; // selection efficiency
    std::cout << "       event selection efficiency:  " << selEff*100. << "%" << std::endl;
    int n100cos =  Norm->events_n100cos(0., 0.); //get final amount of cosmics in defined range for normalisation of flux
    std::cout << "       events with ~100 GeV and 1 - cos(theta) < 1/2pi: " << n100cos << std::endl;
    std::cout << std::endl;
    std::cout << "       momentum range: " << MinP             << " ... " << MaxP << " GeV" << std::endl;
    std::cout << "       theta  range:   " << MinTheta*Rad2Deg << " ... " << MaxTheta*Rad2Deg << " deg" << std::endl; 
    std::cout << "       phi    range:   " << MinPhi*Rad2Deg   << " ... " << MaxPhi*Rad2Deg << " deg" << std::endl;
    std::cout << "       time   range:   " << MinT0            << " ... " << MaxT0 << " ns" << std::endl;
    std::cout << "       energy  loss:   " << ElossScaleFactor*100. << "%" << std::endl;
    std::cout << std::endl;
    double area = 1.e-6*Pi*SurfaceRadius*SurfaceRadius; // area on surface [m^2] 
    if (MinTheta > 90.*Deg2Rad) //upgoing muons from neutrinos)
      std::cout << "       area of initial cosmics at CMS detector bottom surface:   " << area << " m^2" << std::endl;
    else
      std::cout << "       area of initial cosmics on Surface + PlugWidth:   " << area << " m^2" << std::endl;
    std::cout << "       depth of CMS detector (from Surface):   " << SurfaceOfEarth/1000 << " m" << std::endl;
       
    if(n100cos>0){
      // rate: corrected for area and selection-Eff. and normalized to known flux, integration over solid angle (dOmega) is implicit
      // flux is normalised with respect to known flux of vertical 100GeV muons in area at suface level 
      // rate seen by detector is lower than rate at surface area, so has to be corrected for selection-Eff.
      // normalisation factor has unit [1/s/m^2] 
      // rate = N/time --> normalization factor gives 1/runtime/area 
      // normalization with respect to number of actually diced events (Ndiced)

      if (MinTheta > 90.*Deg2Rad) {//upgoing muons from neutrinos) 
	double Omega = (cos(MinTheta)-cos(MaxTheta)) * (MaxPhi-MinPhi);
	EventRate = (Ndiced * 3.e-13) * Omega * area*1.e4 * selEff;//area in cm, flux=3.e-13cm^-2s^-1sr^-1
	rateErr_stat = EventRate/sqrt( (double) n100cos);  // stat. rate error 
	rateErr_syst = EventRate/3.e-13 * 1.0e-13;  // syst. rate error, from error of known flux 
      }
      else {
	EventRate= (Ndiced * Norm->norm(n100cos)) * area * selEff;
	rateErr_stat = EventRate/sqrt( (double) n100cos);  // stat. rate error 
	rateErr_syst = EventRate/2.63e-3 * 0.06e-3;  // syst. rate error, from error of known flux 
      }

      // normalisation in region 1.-cos(theta) < 1./(2.*Pi), if MaxTheta even lower correct for this
      if(MaxTheta<0.572){
	double spacean = 2.*Pi*(1.-cos(MaxTheta));
	EventRate= (Ndiced * Norm->norm(n100cos)) * area * selEff * spacean;
	rateErr_stat = EventRate/sqrt( (double) n100cos);  // rate error 
	rateErr_syst = EventRate/2.63e-3 * 0.06e-3;  // syst. rate error, from error of known flux 
      }

    }else{
      EventRate=Nsel; //no info as no muons at 100 GeV
      rateErr_stat =Nsel;
      rateErr_syst =Nsel;
      std::cout << std::endl;
      if (MinP > 100.)
	std::cout << " !!! MinP > 100 GeV. Cannot apply normalisation!" << std::endl;
      else
	std::cout << " !!! Not enough statistics to apply normalisation (rate=1 +- 1) !!!" << std::endl;
    } 
    
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "       rate is " << EventRate << " +-" << rateErr_stat <<" (stat) " << "+-" << 
      rateErr_syst << " (syst) " <<" muons per second" << std::endl;
    if(EventRate!=0) std::cout << "       number of events corresponds to " << Nsel/EventRate << " s" << std::endl;  //runtime at CMS = Nsel/rate
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << std::endl;
    std::cout << "*********************************************************" << std::endl;
    std::cout << "*********************************************************" << std::endl;
  }
}

void CosmicMuonGenerator::checkIn(){
  if (MinP < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: min.energy is out of range (0 GeV ... inf]" << std::endl << std::endl; }
  if (MaxP < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.energy is out of range (0 GeV ... inf]" << std::endl << std::endl; }
  if (MaxP <= MinP){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.energy is not greater than min.energy" << std::endl << std::endl; }
  if (MinTheta < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: min.theta is out of range [0 deg ... 90 deg)" << std::endl << std::endl; }
  if (MaxTheta < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.theta is out of range [0 deg ... 90 deg)" << std::endl << std::endl; }
  if (MaxTheta <= MinTheta){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.theta is not greater than min.theta" << std::endl << std::endl; }
  if (MinPhi < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: min.phi is out of range [0 deg ... 360 deg]" << std::endl << std::endl; }
  if (MaxPhi < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.phi is out of range [0 deg ... 360 deg]" << std::endl << std::endl; }
  if (MaxPhi <= MinPhi){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.phi is not greater than min.phi" << std::endl << std::endl; }
  if (MaxT0 <= MinT0){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.t0 is not greater than min.t0" << std::endl << std::endl; }
  if (ElossScaleFactor < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: E-loss scale factor is out of range [0 ... inf)" << std::endl << std::endl; }
  if (MinEnu < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: min.Enu is out of range [0 GeV ... inf]" << std::endl << std::endl; }
  if (MaxEnu < 0.){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.Enu is out of range [0 GeV ... inf]" << std::endl << std::endl; }
  if (MaxEnu <= MinEnu){ NumberOfEvents = 0;
    std::cout << "  CMG-ERR: max.Enu is not greater than min.Enu" << std::endl << std::endl; }

}

bool CosmicMuonGenerator::goodOrientation(){
  // check angular range (for a sphere with Target3dRadius around the target)
  bool goodAngles = false;
  bool phiaccepted = false;
  bool thetaaccepted = false;
  double RxzV = sqrt(OneMuoEvt.vx()*OneMuoEvt.vx() + OneMuoEvt.vz()*OneMuoEvt.vz());

  double rVY;
  if (MinTheta > 90.*Deg2Rad) //upgoing muons from neutrinos
    rVY = -sqrt(RxzV*RxzV + RadiusCMS*RadiusCMS);
  else 
    rVY = sqrt(RxzV*RxzV + (SurfaceOfEarth+PlugWidth)*(SurfaceOfEarth+PlugWidth));

  double Phi = OneMuoEvt.phi();
  double PhiV = atan2(OneMuoEvt.vx(),OneMuoEvt.vz()) + Pi; if (PhiV > TwoPi) PhiV -= TwoPi;
  double disPhi = fabs(PhiV - Phi); if (disPhi > Pi) disPhi = TwoPi - disPhi;
  double dPhi = Pi; if (RxzV > Target3dRadius) dPhi = asin(Target3dRadius/RxzV);
  if (disPhi < dPhi) phiaccepted = true;
  double Theta = OneMuoEvt.theta();
  double ThetaV = asin(RxzV/rVY);
  double dTheta = Pi; if (fabs(rVY) > Target3dRadius) dTheta = asin(Target3dRadius/fabs(rVY));
  //std::cout << "    dPhi = " <<   dPhi << "  (" <<   Phi << " <p|V> " <<   PhiV << ")" << std::endl;
  //std::cout << "  dTheta = " << dTheta << "  (" << Theta << " <p|V> " << ThetaV << ")" << std::endl;

  if (!phiaccepted && RxzV < Target3dRadius)
  //if (RxzV < Target3dRadius)
    std::cout << "Rejected phi=" << Phi << "  PhiV="  << PhiV 
	 << "  dPhi=" << dPhi << "  disPhi=" << disPhi
	 << "  RxzV=" << RxzV << "  Target3dRadius=" << Target3dRadius 
	 << "  Theta=" << Theta << std::endl;

  if (fabs(Theta-ThetaV) < dTheta) thetaaccepted = true;
  if (phiaccepted && thetaaccepted) goodAngles = true;
  return goodAngles;
}

void CosmicMuonGenerator::initEvDis(){
#if ROOT_INTERACTIVE
  float rCMS = RadiusCMS/1000.;
  float zCMS = Z_DistCMS/1000.;
  if(TrackerOnly==true){
    rCMS = RadiusTracker/1000.;
    zCMS = Z_DistTracker/1000.;
}
  TH2F* disXY = new TH2F("disXY","X-Y view",160,-rCMS,rCMS,160,-rCMS,rCMS);
  TH2F* disZY = new TH2F("disZY","Z-Y view",150,-zCMS,zCMS,160,-rCMS,rCMS);
  gStyle->SetPalette(1,0);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1.5);
  TCanvas *disC = new TCanvas("disC","Cosmic Muon Event Display",0,0,800,410);
  disC->Divide(2,1);
  disC->cd(1);
  gPad->SetTicks(1,1);
  disXY->SetMinimum(log10(MinP));
  disXY->SetMaximum(log10(MaxP));
  disXY->GetXaxis()->SetLabelSize(0.05);
  disXY->GetXaxis()->SetTitleSize(0.05);
  disXY->GetXaxis()->SetTitleOffset(1.0);
  disXY->GetXaxis()->SetTitle("X [m]");
  disXY->GetYaxis()->SetLabelSize(0.05);
  disXY->GetYaxis()->SetTitleSize(0.05);
  disXY->GetYaxis()->SetTitleOffset(0.8);
  disXY->GetYaxis()->SetTitle("Y [m]");
  disC->cd(2);
  gPad->SetGrid(1,1);
  gPad->SetTicks(1,1);
  disZY->SetMinimum(log10(MinP));
  disZY->SetMaximum(log10(MaxP));
  disZY->GetXaxis()->SetLabelSize(0.05);
  disZY->GetXaxis()->SetTitleSize(0.05);
  disZY->GetXaxis()->SetTitleOffset(1.0);
  disZY->GetXaxis()->SetTitle("Z [m]");
  disZY->GetYaxis()->SetLabelSize(0.05);
  disZY->GetYaxis()->SetTitleSize(0.05);
  disZY->GetYaxis()->SetTitleOffset(0.8);
  disZY->GetYaxis()->SetTitle("Y [m]");
#endif
}

void CosmicMuonGenerator::displayEv(){
#if ROOT_INTERACTIVE
  double RadiusDet=RadiusCMS;
  double Z_DistDet=Z_DistCMS;
  if(TrackerOnly==true){
    RadiusDet = RadiusTracker;
    Z_DistDet = Z_DistTracker;
  }
  disXY->Reset();
  disZY->Reset();
  TMarker* InteractionPoint = new TMarker(0.,0.,2);
  TArc* r8m = new TArc(0.,0.,(RadiusDet/1000.));
  TLatex* logEaxis = new TLatex(); logEaxis->SetTextSize(0.05);
  float energy = float(OneMuoEvt.e());
  float verX = float(OneMuoEvt.vx()/1000.); // [m]
  float verY = float(OneMuoEvt.vy()/1000.); // [m]
  float verZ = float(OneMuoEvt.vz()/1000.); // [m]
  float dirX = float(OneMuoEvt.px())/fabs(OneMuoEvt.py());
  float dirY = float(OneMuoEvt.py())/fabs(OneMuoEvt.py());
  float dirZ = float(OneMuoEvt.pz())/fabs(OneMuoEvt.py());
  float yStep = disXY->GetYaxis()->GetBinWidth(1);
  int   NbinY = disXY->GetYaxis()->GetNbins();
  for (int iy=0; iy<NbinY; ++iy){
    verX += dirX*yStep;
    verY += dirY*yStep;
    verZ += dirZ*yStep;
    float rXY = sqrt(verX*verX + verY*verY)*1000.; // [mm]
    float absZ = fabs(verZ)*1000.;                 // [mm]
    if (rXY < RadiusDet && absZ < Z_DistDet){
      disXY->Fill(verX,verY,log10(energy));
      disZY->Fill(verZ,verY,log10(energy));
      disC->cd(1); disXY->Draw("COLZ"); InteractionPoint->Draw("SAME"); r8m->Draw("SAME");
      logEaxis->DrawLatex((0.65*RadiusDet/1000.),(1.08*RadiusDet/1000.),"log_{10}E(#mu^{#pm})");
      disC->cd(2); disZY->Draw("COL"); InteractionPoint->Draw("SAME");
      gPad->Update();
    }
  }
#endif
}

void CosmicMuonGenerator::setNumberOfEvents(unsigned int N){ if (NotInitialized) NumberOfEvents = N; }

void CosmicMuonGenerator::setRanSeed(int N){ if (NotInitialized) RanSeed = N; }

void CosmicMuonGenerator::setMinP(double P){ if (NotInitialized) MinP = P; }

void CosmicMuonGenerator::setMinP_CMS(double P){ if (NotInitialized) MinP_CMS = P; }

void CosmicMuonGenerator::setMaxP(double P){ if (NotInitialized) MaxP = P; }

void CosmicMuonGenerator::setMinTheta(double Theta){ if (NotInitialized) MinTheta = Theta*Deg2Rad; }

void CosmicMuonGenerator::setMaxTheta(double Theta){ if (NotInitialized) MaxTheta = Theta*Deg2Rad; } 

void CosmicMuonGenerator::setMinPhi(double Phi){ if (NotInitialized) MinPhi = Phi*Deg2Rad; }

void CosmicMuonGenerator::setMaxPhi(double Phi){ if (NotInitialized) MaxPhi = Phi*Deg2Rad; }

void CosmicMuonGenerator::setMinT0(double T0){ if (NotInitialized) MinT0 = T0; }

void CosmicMuonGenerator::setMaxT0(double T0){ if (NotInitialized) MaxT0 = T0; }

void CosmicMuonGenerator::setElossScaleFactor(double ElossScaleFact){ if (NotInitialized) ElossScaleFactor = ElossScaleFact; }

void CosmicMuonGenerator::setRadiusOfTarget(double R){ if (NotInitialized) RadiusOfTarget = R; }

void CosmicMuonGenerator::setZDistOfTarget(double Z){ if (NotInitialized) ZDistOfTarget = Z; }

void CosmicMuonGenerator::setZCentrOfTarget(double Z){ if (NotInitialized) ZCentrOfTarget = Z; }

void CosmicMuonGenerator::setTrackerOnly(bool Tracker){ if (NotInitialized) TrackerOnly = Tracker; }

void CosmicMuonGenerator::setTIFOnly_constant(bool TIF){ if (NotInitialized) TIFOnly_constant = TIF; }

void CosmicMuonGenerator::setTIFOnly_linear(bool TIF){ if (NotInitialized) TIFOnly_linear = TIF; }

void CosmicMuonGenerator::setMTCCHalf(bool MTCC){ if (NotInitialized) MTCCHalf = MTCC; }

void CosmicMuonGenerator::setPlugVx(double PlugVtx){ if (NotInitialized) PlugVx = PlugVtx; }
void CosmicMuonGenerator::setPlugVz(double PlugVtz){ if (NotInitialized) PlugVz = PlugVtz; }

void CosmicMuonGenerator::setMinEnu(double MinEn){ if (NotInitialized) MinEnu = MinEn; }
void CosmicMuonGenerator::setMaxEnu(double MaxEn){ if (NotInitialized) MaxEnu = MaxEn; }

double CosmicMuonGenerator::getRate(){ return EventRate; }
