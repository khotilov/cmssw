#include "Validation/Geometry/interface/MaterialBudgetHcalHistos.h"

#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDLogicalPart.h"
#include "DetectorDescription/Core/interface/DDSplit.h"
#include "DetectorDescription/Core/interface/DDValue.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"

MaterialBudgetHcalHistos::MaterialBudgetHcalHistos(const edm::ParameterSet &p,
						   TestHistoMgr* mgr) : hmgr(mgr) {
  theFileName = p.getUntrackedParameter<std::string>("HistoFile","MBHC.root");
  binEta      = p.getUntrackedParameter<int>("NBinEta", 260);
  binPhi      = p.getUntrackedParameter<int>("NBinPhi", 180);
  maxEta      = p.getUntrackedParameter<double>("MaxEta", 5.2);
  edm::LogInfo("MaterialBudget") << "MaterialBudgetHcalHistos: Output file "
				 << theFileName;
  book();

}

void MaterialBudgetHcalHistos::fillBeginJob(const DDCompactView & cpv) {

  std::string attribute = "ReadOutName";
  std::string value     = "HcalHits";
  DDSpecificsFilter filter1;
  DDValue           ddv1(attribute,value,0);
  filter1.setCriteria(ddv1,DDSpecificsFilter::equals);
  DDFilteredView fv1(cpv);
  fv1.addFilter(filter1);
  sensitives = getNames(fv1);
  edm::LogInfo("MaterialBudget") << "MaterialBudgetHcalHistos: Names to be "
				 << "tested for " << attribute << " = " 
				 << value << " has " << sensitives.size()
				 << " elements";
  for (unsigned int i=0; i<sensitives.size(); i++) 
    edm::LogInfo("MaterialBudget") << "MaterialBudgetHcalHistos:  Sensitive[" 
				   << i << "] = " << sensitives[i];

  attribute = "Volume";
  value     = "HF";
  DDSpecificsFilter filter2;
  DDValue           ddv2(attribute,value,0);
  filter2.setCriteria(ddv2,DDSpecificsFilter::equals);
  DDFilteredView fv2(cpv);
  fv2.addFilter(filter2);
  hfNames = getNames(fv2);
  fv2.firstChild();
  DDsvalues_type sv(fv2.mergedSpecifics());
  std::vector<double> temp =  getDDDArray("Levels",sv);
  edm::LogInfo("MaterialBudget") << "MaterialBudgetHcalHistos: Names to be "
				 << "tested for " << attribute << " = " 
				 << value << " has " << hfNames.size()
				 << " elements";
  for (unsigned int i=0; i<hfNames.size(); i++) {
    int level = static_cast<int>(temp[i]);
    hfLevels.push_back(level);
    edm::LogInfo("MaterialBudget") << "MaterialBudgetHcalHistos:  HF[" << i 
				   << "] = " << hfNames[i] << " at level " 
				   << hfLevels[i];
  }

}

void MaterialBudgetHcalHistos::fillStartTrack(const G4Track* aTrack) {

  id     = layer  = steps   = 0;
  radLen = intLen = stepLen = 0;

  const G4ThreeVector& dir = aTrack->GetMomentum() ;
  if (dir.theta() != 0 ) {
    eta = dir.eta();
  } else {
    eta = -99;
  }
  phi = dir.phi();
  double theEnergy = aTrack->GetTotalEnergy();
  int    theID     = (int)(aTrack->GetDefinition()->GetPDGEncoding());

  LogDebug("MaterialBudget") << "MaterialBudgetHcalHistos: Track " 
			     << aTrack->GetTrackID() << " Code " << theID
			     << " Energy " << theEnergy/GeV << " GeV; Eta "
			     << eta << " Phi " << phi/deg << " PT "
			     << dir.perp()/GeV << " GeV";
}


void MaterialBudgetHcalHistos::fillPerStep(const G4Step* aStep) {

  G4Material * material = aStep->GetPreStepPoint()->GetMaterial();
  double step    = aStep->GetStepLength();
  double radl    = material->GetRadlen();
  double intl    = material->GetNuclearInterLength();
  double density = material->GetDensity() / (g/cm3);

  int    idOld   = id;
  const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();
  std::string         name  = touch->GetVolume(0)->GetName();
  LogDebug("MaterialBudget") << "MaterialBudgetHcalHistos: Step at " << name
			     << " Length " << step << " in " 
			     << material->GetName() << " of density "
			     << density << " g/cc; Radiation Length "
			     << radl << " mm; Interaction Length " << intl
			     << " mm";

  int det=0, lay=0;
  if (isSensitive(name)) {
    if (isItHF(touch)) {
      det = 5;
      lay = 20;
    } else {
      det   = (touch->GetReplicaNumber(1))/1000;
      lay   = (touch->GetReplicaNumber(0)/10)%100 + 1;
    }
    LogDebug("MaterialBudget") << "MaterialBudgetHcalHistos: Det " << det
			       << " Layer " << lay << " Eta " << eta << " Phi "
			       << phi/deg;
    if (det > 0) {
      if (layer != lay) {
	id++;
	layer = lay;
      }
    }
  }

  if (id != idOld) {
    fillHisto(idOld);
    idOld = id;
  }

  stepLen += step;
  radLen  += step/radl;
  intLen  += step/intl;
  if (layer == 20) {
    if (!isItHF(aStep->GetPostStepPoint()->GetTouchable())) {
      LogDebug("MaterialBudget") << "MaterialBudgetHcalHistos: After HF in " 
				 << aStep->GetPostStepPoint()->GetTouchable()->GetVolume(0)->GetName();
      fillHisto(id);
      id++;
    }
  }
}


void MaterialBudgetHcalHistos::fillEndTrack() {

  fillHisto(maxSet-1);
}

void MaterialBudgetHcalHistos::book() {
  
  double maxPhi=pi;
  edm::LogInfo("MaterialBudget") << "MaterialBudgetHcalHistos: Booking user "
				 << "histos === with " << binEta << " bins "
				 << "in eta from " << -maxEta << " to "
				 << maxEta << " and " << binPhi << " bins "
				 << "in phi from " << -maxPhi << " to " 
				 << maxPhi;
  
  char  name[10], title[40];
  // total X0
  for (int i=0; i<maxSet; i++) {
    sprintf(name, "%d", i+100);
    sprintf(title, "MB(X0) prof Eta in region %d", i);
    hmgr->addHistoProf1( new TProfile(name, title, binEta, -maxEta, maxEta) );
    sprintf(name, "%d", i+200);
    sprintf(title, "MB(L0) prof Eta in region %d", i);
    hmgr->addHistoProf1( new TProfile(name, title, binEta, -maxEta, maxEta) );
    sprintf(name, "%d", i+300);
    sprintf(title, "MB(Step) prof Eta in region %d", i);
    hmgr->addHistoProf1( new TProfile(name, title, binEta, -maxEta, maxEta) );
    sprintf(name, "%d", i+400);
    sprintf(title, "Eta in region %d", i);
    hmgr->addHisto1( new TH1F(name, title, binEta, -maxEta, maxEta) );
    sprintf(name, "%d", i+500);
    sprintf(title, "MB(X0) prof Ph in region %d", i);
    hmgr->addHistoProf1( new TProfile(name, title, binPhi, -maxPhi, maxPhi) );
    sprintf(name, "%d", i+600);
    sprintf(title, "MB(L0) prof Ph in region %d", i);
    hmgr->addHistoProf1( new TProfile(name, title, binPhi, -maxPhi, maxPhi) );
    sprintf(name, "%d", i+700);
    sprintf(title, "MB(Step) prof Ph in region %d", i);
    hmgr->addHistoProf1( new TProfile(name, title, binPhi, -maxPhi, maxPhi) );
    sprintf(name, "%d", i+800);
    sprintf(title, "Phi in region %d", i);
    hmgr->addHisto1( new TH1F(name, title, binPhi, -maxPhi, maxPhi) );
    sprintf(name, "%d", i+900);
    sprintf(title, "MB(X0) prof Eta Phi in region %d", i);
    hmgr->addHistoProf2( new TProfile2D(name, title, binEta, -maxEta, maxEta, 
					binPhi, -maxPhi, maxPhi) );
    sprintf(name, "%d", i+1000);
    sprintf(title, "MB(L0) prof Eta Phi in region %d", i);
    hmgr->addHistoProf2( new TProfile2D(name, title, binEta, -maxEta, maxEta, 
					binPhi, -maxPhi, maxPhi) );
    sprintf(name, "%d", i+1100);
    sprintf(title, "MB(Step) prof Eta Phi in region %d", i);
    hmgr->addHistoProf2( new TProfile2D(name, title, binEta, -maxEta, maxEta, 
					binPhi, -maxPhi, maxPhi) );
    sprintf(name, "%d", i+1200);
    sprintf(title, "Eta vs Phi in region %d", i);
    hmgr->addHisto2( new TH2F(name, title, binEta, -maxEta, maxEta, 
			      binPhi, -maxPhi, maxPhi) );
  }

  edm::LogInfo("MaterialBudget") << "MaterialBudgetHcalHistos: Booking user "
				 << "histos done ===";

}

void MaterialBudgetHcalHistos::fillHisto(int ii) {

  LogDebug("MaterialBudget") << "MaterialBudgetHcalHistos:FillHisto called "
			     << "with index " << ii << " integrated  step "
			     << stepLen << " X0 " << radLen << " Lamda " 
			     << intLen;
  
  if (ii >=0 && ii < maxSet) {
    hmgr->getHisto1(400+ii)->Fill(eta);
    hmgr->getHisto1(800+ii)->Fill(phi);
    hmgr->getHisto2(1200+ii)->Fill(eta, phi);
  
    hmgr->getHistoProf1(100+ii)->Fill(eta, radLen);
    hmgr->getHistoProf1(500+ii)->Fill(phi, radLen);
    hmgr->getHistoProf2(900+ii)->Fill(eta, phi, radLen);
  
    hmgr->getHistoProf1(200+ii)->Fill(eta, intLen);
    hmgr->getHistoProf1(600+ii)->Fill(phi, intLen);
    hmgr->getHistoProf2(1000+ii)->Fill(eta, phi, intLen);
  
    hmgr->getHistoProf1(300+ii)->Fill(eta, stepLen);
    hmgr->getHistoProf1(700+ii)->Fill(phi, stepLen);
    hmgr->getHistoProf2(1100+ii)->Fill(eta, phi, stepLen);
  }
}

void MaterialBudgetHcalHistos::hend() {
  edm::LogInfo("MaterialBudget") << "MaterialBudgetHcalHistos: Save user "
				 << "histos ===";
  hmgr->save( theFileName );

}

std::vector<std::string> MaterialBudgetHcalHistos::getNames(DDFilteredView& fv) {

  std::vector<std::string> tmp;
  bool dodet = fv.firstChild();
  while (dodet) {
    const DDLogicalPart & log = fv.logicalPart();
    std::string namx = DDSplit(log.name()).first;
    bool ok = true;
    for (unsigned int i=0; i<tmp.size(); i++)
      if (namx == tmp[i]) ok = false;
    if (ok) tmp.push_back(namx);
    dodet = fv.next();
  }
  return tmp;
}

std::vector<double> MaterialBudgetHcalHistos::getDDDArray(const std::string & str,
							  const DDsvalues_type & sv) {

  LogDebug("MaterialBudget") << "MaterialBudgetHcalHistos:getDDDArray called "
			     << "for " << str;
  DDValue value(str);
  if (DDfetch(&sv,value)) {
    LogDebug("MaterialBudget") << value;
    const std::vector<double> & fvec = value.doubles();
    int nval = fvec.size();
    if (nval < 1) {
      edm::LogError("MaterialBudget") << "MaterialBudgetHcalHistos : # of " 
				      << str << " bins " << nval
				      << " < 1 ==> illegal";
      throw cms::Exception("Unknown", "MaterialBudgetHcalHistos")
        << "nval < 1 for array " << str << "\n";
    }

    return fvec;
  } else {
    edm::LogError("MaterialBudget") << "MaterialBudgetHcalHistos : cannot get "
				    << "array " << str;
    throw cms::Exception("Unknown", "MaterialBudgetHcalHistos")
      << "cannot get array " << str << "\n";
  }
}

bool MaterialBudgetHcalHistos::isSensitive (std::string name) {

  std::vector<std::string>::const_iterator it = sensitives.begin();
  for (; it != sensitives.end(); it++)
    if (name == *it) return true;
  return false;
}

bool MaterialBudgetHcalHistos::isItHF (const G4VTouchable* touch) {

  std::vector<std::string>::const_iterator it = hfNames.begin();
  int levels = ((touch->GetHistoryDepth())+1);
  for (unsigned int it=0; it < hfNames.size(); it++) {
    if (levels >= hfLevels[it]) {
      std::string name = touch->GetVolume(levels-hfLevels[it])->GetName();
      if (name == hfNames[it]) return true;
    }
  }
  return false;
}

