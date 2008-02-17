#include "TopQuarkAnalysis/TopEventSelection/interface/TtHadLRSignalSelObservables.h"

using namespace reco;
using namespace std;
using namespace math;

TtHadLRSignalSelObservables::TtHadLRSignalSelObservables()
{
}

TtHadLRSignalSelObservables::~TtHadLRSignalSelObservables()
{
}

void TtHadLRSignalSelObservables::operator() (TtHadEvtSolution &TS)
{
  evtselectVarVal.clear();
  vector<pat::Jet> topJets;
  topJets.clear();
  topJets.push_back(TS.getHadp());
  topJets.push_back(TS.getHadq());
  topJets.push_back(TS.getHadj());
  topJets.push_back(TS.getHadk());
  topJets.push_back(TS.getHadb());
  topJets.push_back(TS.getHadbbar());
  
  //Assume the highest Et jets are the b-jets, the others the jets from the W - but how to work out which two are from which W? FIXME!!! for now don't sort jets in Et before after this...
  TLorentzVector *Hadp = new TLorentzVector();	
  TLorentzVector *Hadq = new TLorentzVector();
  TLorentzVector *Hadj = new TLorentzVector();	
  TLorentzVector *Hadk = new TLorentzVector();
  TLorentzVector *Hadb = new TLorentzVector();	
  TLorentzVector *Hadbbar = new TLorentzVector();
  Hadp->SetPxPyPzE(topJets[0].recJet().px(),topJets[0].recJet().py(),topJets[0].recJet().pz(),topJets[3].recJet().energy());
  Hadq->SetPxPyPzE(topJets[1].recJet().px(),topJets[1].recJet().py(),topJets[1].recJet().pz(),topJets[2].recJet().energy());	
  Hadj->SetPxPyPzE(topJets[2].recJet().px(),topJets[2].recJet().py(),topJets[2].recJet().pz(),topJets[2].recJet().energy());
  Hadk->SetPxPyPzE(topJets[3].recJet().px(),topJets[3].recJet().py(),topJets[3].recJet().pz(),topJets[3].recJet().energy());
  Hadb->SetPxPyPzE(topJets[4].recJet().px(),topJets[4].recJet().py(),topJets[4].recJet().pz(),topJets[1].recJet().energy());
  Hadbbar->SetPxPyPzE(topJets[4].recJet().px(),topJets[4].recJet().py(),topJets[4].recJet().pz(),topJets[1].recJet().energy());
  
  //sort the topJets in Et
  std::sort(topJets.begin(),topJets.end(),EtComparator);
  
  //Et-Sum of the lightest jets
  double EtSum = topJets[5].recJet().et()+topJets[5].recJet().et();
  double Obs1 = (EtSum>0 ? EtSum : -1);
  evtselectVarVal.push_back(pair<unsigned int,double>(1,Obs1));
  
  //Log of the difference in Bdisc between the 2nd and the 3rd jets (ordered in Bdisc) - does that still
  //make sense for the fully hadronic channel as well?FIXME!!!
  
  //sort the jets in Bdiscriminant
  std::sort(topJets.begin(),topJets.end(),BdiscComparator);
  
  double BGap = topJets[1].bDiscriminator("trackCountingJetTags") - topJets[2].bDiscriminator("trackCountingJetTags");
  double Obs2 = (BGap>0 ? log(BGap) : -1);
  evtselectVarVal.push_back(pair<unsigned int,double>(2,Obs2));
  
  //Circularity of the event
  double N=0,D=0,C_tmp=0,C=1000;
  double nx,ny,nz;
  
  // C = 2min(E(pt.n)^2/E(pt)^2) = 2*N/D but it is theorically preferable to use C'=PI/2*min(E|pt.n|/E|pt|), sum over all jets+lepton+MET (cf PhysRevD 48 R3953(Nov 1993))
  
  for(unsigned int i=0;i<6;i++){
    D += fabs(topJets[i].recJet().pt());
  }
  
  if((D>0)){
    // Loop over all the unit vectors in the transverse plane in order to find the miminum : 
    for(unsigned int i=0; i<360; i++){
      nx = cos((2*PI/360)*i);
      ny = sin((2*PI/360)*i);
      nz = 0;
      N=0;
      for(unsigned int i=0;i<4;i++){
	N += fabs(topJets[i].recJet().px()*nx+topJets[i].recJet().py()*ny+topJets[i].recJet().pz()*nz);
      }
      C_tmp = 2*N/D;
      if(C_tmp<C) C = C_tmp;
    }
  }
  
  double Obs3 = ( C!=1000 ? C : -1);
  evtselectVarVal.push_back(pair<unsigned int,double>(3,Obs3));
  
  
  //HT variable (Et-sum of the six jets)
  double HT=0;
  for(unsigned int i=0;i<6;i++){
    HT += topJets[i].recJet().et();
  }
  
  double Obs4 = ( HT!=0 ? HT : -1);
  evtselectVarVal.push_back(pair<unsigned int,double>(4,Obs4));
  
  //Transverse Mass of the system
  XYZTLorentzVector pjets;
  // for the six jets 
  for(unsigned int i=0;i<6;i++){
    pjets += topJets[i].recJet().p4();
  }
  double MT = sqrt(pow(pjets.mass(),2));
  double Obs5 = ( MT>0 ? MT : -1);
  evtselectVarVal.push_back(pair<unsigned int,double>(5,Obs5));
  
  //CosTheta(Hadp,Hadq) and CosTheta(Hadj,Hadq) 
  //sort the lightJets in Et
  std::sort(topJets.begin(),topJets.end(),EtComparator);
  
  double px1 = topJets[2].recJet().px();     double px2 = topJets[3].recJet().px();
  double py1 = topJets[2].recJet().py();     double py2 = topJets[3].recJet().py();
  double pz1 = topJets[2].recJet().pz();     double pz2 = topJets[3].recJet().pz();
  double E1  = topJets[2].recJet().energy(); double E2  = topJets[3].recJet().energy();
  double px3 = topJets[4].recJet().px();     double px4 = topJets[5].recJet().px();
  double py3 = topJets[4].recJet().py();     double py4 = topJets[5].recJet().py();
  double pz3 = topJets[4].recJet().pz();     double pz4 = topJets[5].recJet().pz();
  double E3  = topJets[4].recJet().energy(); double E4  = topJets[5].recJet().energy();
  
  TLorentzVector *LightJet1 = new TLorentzVector();
  TLorentzVector *LightJet2 = new TLorentzVector();
  TLorentzVector *LightJet3 = new TLorentzVector();
  TLorentzVector *LightJet4 = new TLorentzVector();
  
  LightJet1->SetPxPyPzE(px1,py1,pz1,E1);
  //LightJet1->Boost(BoostBackToCM);
  LightJet2->SetPxPyPzE(px2,py2,pz2,E2);
  //LightJet2->Boost(BoostBackToCM);
  LightJet3->SetPxPyPzE(px3,py3,pz3,E3);
  //LightJet1->Boost(BoostBackToCM);
  LightJet4->SetPxPyPzE(px4,py4,pz4,E4);
  //LightJet2->Boost(BoostBackToCM);
  
  double CosTheta1 = cos(LightJet2->Angle(LightJet1->Vect()));
  double CosTheta2 = cos(LightJet4->Angle(LightJet3->Vect()));
  
  double Obs6 = ( -1<CosTheta1 ? CosTheta1 : -2);
  evtselectVarVal.push_back(pair<unsigned int,double>(6,Obs6));
  double Obs7 = ( -1<CosTheta2 ? CosTheta2 : -2);
  evtselectVarVal.push_back(pair<unsigned int,double>(7,Obs7));
  
  delete LightJet1;
  delete LightJet2;
  delete LightJet3;
  delete LightJet4;
  
  // try to find out more powerful observables related to the b-disc
  //sort the jets in Bdiscriminant
  std::sort(topJets.begin(),topJets.end(),BdiscComparator);
  
  double BjetsBdiscSum = topJets[0].bDiscriminator("trackCountingJetTags") + topJets[1].bDiscriminator("trackCountingJetTags");
  double Ljets1BdiscSum = topJets[2].bDiscriminator("trackCountingJetTags") + topJets[3].bDiscriminator("trackCountingJetTags");
  double Ljets2BdiscSum = topJets[4].bDiscriminator("trackCountingJetTags") + topJets[5].bDiscriminator("trackCountingJetTags");
  //cout<<"BjetsBdiscSum = "<<BjetsBdiscSum<<endl;
  //cout<<"LjetsBdiscSum = "<<LjetsBdiscSum<<endl;
  
  double Obs8 = (Ljets1BdiscSum !=0 ? log(BjetsBdiscSum/Ljets1BdiscSum) : -1);
  evtselectVarVal.push_back(pair<unsigned int,double>(8,Obs8));	
  double Obs9 = (Ljets2BdiscSum !=0 ? log(BjetsBdiscSum/Ljets2BdiscSum) : -1);
  evtselectVarVal.push_back(pair<unsigned int,double>(9,Obs9));
  double Obs10 = (BGap>0 ? log(BjetsBdiscSum*BGap) : -1);
  evtselectVarVal.push_back(pair<unsigned int,double>(10,Obs10));
  
  
  // Et-Ratio between the two highest in Et jets and four highest jets	
  double Obs11 = (HT!=0 ? (HT-EtSum)/HT : -1);
  evtselectVarVal.push_back(pair<unsigned int,double>(11,Obs11));
  
  
  //Sphericity and Aplanarity without boosting back the system to CM frame
  
  TMatrixDSym Matrix(3);
  double PX2 = pow(Hadp->Px(),2)+pow(Hadq->Px(),2)+pow(Hadb->Px(),2)+
    pow(Hadj->Px(),2)+pow(Hadk->Px(),2)+pow(Hadbbar->Px(),2);
  double PY2 = pow(Hadp->Py(),2)+pow(Hadq->Py(),2)+pow(Hadb->Py(),2)+
    pow(Hadj->Py(),2)+pow(Hadk->Py(),2)+pow(Hadbbar->Py(),2);
  double PZ2 = pow(Hadp->Pz(),2)+pow(Hadq->Pz(),2)+pow(Hadb->Pz(),2)+ 
    pow(Hadj->Pz(),2)+pow(Hadk->Pz(),2)+pow(Hadbbar->Pz(),2);
  double P2  = PX2+PY2+PZ2;
  
  double PXY = Hadp->Px()*Hadp->Py()+Hadq->Px()*Hadq->Py()+Hadb->Px()*Hadb->Py()+
    Hadj->Px()*Hadj->Py()+Hadk->Px()*Hadk->Py()+Hadbbar->Px()*Hadbbar->Py();
  double PXZ = Hadp->Px()*Hadp->Pz()+Hadq->Px()*Hadq->Pz()+Hadb->Px()*Hadb->Pz()+
    Hadj->Px()*Hadj->Pz()+Hadk->Px()*Hadk->Pz()+Hadbbar->Px()*Hadbbar->Pz();
  double PYZ = Hadp->Py()*Hadp->Pz()+Hadq->Py()*Hadq->Pz()+Hadb->Py()*Hadb->Pz()+
    Hadj->Py()*Hadj->Pz()+Hadk->Py()*Hadk->Pz()+Hadbbar->Py()*Hadbbar->Pz();
  
  Matrix(0,0) = PX2/P2; Matrix(0,1) = PXY/P2; Matrix(0,2) = PXZ/P2;
  Matrix(1,0) = PXY/P2; Matrix(1,1) = PY2/P2; Matrix(1,2) = PYZ/P2;
  Matrix(2,0) = PXZ/P2; Matrix(2,1) = PYZ/P2; Matrix(2,2) = PZ2/P2;
  
  TMatrixDSymEigen pTensor(Matrix);
  
  vector<double> EigValues;
  EigValues.clear();
  for(int i=0;i<3;i++){
    EigValues.push_back(pTensor.GetEigenValues()[i]);
  }
  
  std::sort(EigValues.begin(),EigValues.end(),dComparator);
  
  double Sphericity = 1.5*(EigValues[1]+EigValues[2]);
  double Aplanarity = 1.5*EigValues[2];
  
  double Obs12 = (isnan(Sphericity) ? -1 : Sphericity);
  evtselectVarVal.push_back(pair<unsigned int,double>(12,Obs12));
  
  double Obs13 = (isnan(Aplanarity) ? -1 : Aplanarity);
  evtselectVarVal.push_back(pair<unsigned int,double>(13,Obs13));
  
  
  //Sphericity and Aplanarity with boosting back the system to CM frame	
  TLorentzVector *TtbarSystem = new TLorentzVector();
  TtbarSystem->SetPx(Hadp->Px()+Hadq->Px()+Hadb->Px()+Hadj->Px()+Hadk->Px()+Hadbbar->Px());
  TtbarSystem->SetPy(Hadp->Py()+Hadq->Py()+Hadb->Py()+Hadj->Py()+Hadk->Py()+Hadbbar->Py());
  TtbarSystem->SetPz(Hadp->Pz()+Hadq->Pz()+Hadb->Pz()+Hadj->Pz()+Hadk->Pz()+Hadbbar->Pz());
  TtbarSystem->SetE(Hadp->E()+Hadq->E()+Hadb->E()+Hadj->E()+Hadk->E()+Hadbbar->E());
  
  TVector3 BoostBackToCM = -(TtbarSystem->BoostVector());
  Hadp->Boost(BoostBackToCM);
  Hadq->Boost(BoostBackToCM);
  Hadb->Boost(BoostBackToCM);
  Hadj->Boost(BoostBackToCM);
  Hadk->Boost(BoostBackToCM);
  Hadbbar->Boost(BoostBackToCM);
  
  double BOOST_PX2 = pow(Hadp->Px(),2)+pow(Hadq->Px(),2)+pow(Hadb->Px(),2)+
    pow(Hadj->Px(),2)+pow(Hadk->Px(),2)+pow(Hadbbar->Px(),2);
  double BOOST_PY2 = pow(Hadp->Py(),2)+pow(Hadq->Py(),2)+pow(Hadb->Py(),2)+
    pow(Hadj->Py(),2)+pow(Hadk->Py(),2)+pow(Hadbbar->Py(),2);
  double BOOST_PZ2 = pow(Hadp->Pz(),2)+pow(Hadq->Pz(),2)+pow(Hadb->Pz(),2)+
    pow(Hadj->Pz(),2)+pow(Hadk->Pz(),2)+pow(Hadbbar->Pz(),2);
  
  double BOOST_P2  = BOOST_PX2+BOOST_PY2+BOOST_PZ2;
  
  double BOOST_PXY = Hadp->Px()*Hadp->Py()+Hadq->Px()*Hadq->Py()+Hadb->Px()*Hadb->Py()+
    Hadj->Px()*Hadj->Py()+Hadk->Px()*Hadk->Py()+Hadbbar->Px()*Hadbbar->Py();
  double BOOST_PXZ = Hadp->Px()*Hadp->Pz()+Hadq->Px()*Hadq->Pz()+Hadb->Px()*Hadb->Pz()+
    Hadj->Px()*Hadj->Pz()+Hadk->Px()*Hadk->Pz()+Hadbbar->Px()*Hadbbar->Pz();
  double BOOST_PYZ = Hadp->Py()*Hadp->Pz()+Hadq->Py()*Hadq->Pz()+Hadb->Py()*Hadb->Pz()+
    Hadj->Py()*Hadj->Pz()+Hadk->Py()*Hadk->Pz()+Hadbbar->Py()*Hadbbar->Pz();
  
  TMatrixDSym BOOST_Matrix(3);
  
  BOOST_Matrix(0,0) = BOOST_PX2/BOOST_P2; BOOST_Matrix(0,1) = BOOST_PXY/BOOST_P2; BOOST_Matrix(0,2) = BOOST_PXZ/BOOST_P2;
  BOOST_Matrix(1,0) = BOOST_PXY/BOOST_P2; BOOST_Matrix(1,1) = BOOST_PY2/BOOST_P2; BOOST_Matrix(1,2) = BOOST_PYZ/BOOST_P2;
  BOOST_Matrix(2,0) = BOOST_PXZ/BOOST_P2; BOOST_Matrix(2,1) = BOOST_PYZ/BOOST_P2; BOOST_Matrix(2,2) = BOOST_PZ2/BOOST_P2;
  
  TMatrixDSymEigen BOOST_pTensor(BOOST_Matrix);
  
  vector<double> BOOST_EigValues;
  BOOST_EigValues.clear();
  for(int i=0;i<3;i++){
    BOOST_EigValues.push_back(BOOST_pTensor.GetEigenValues()[i]);
  }
  
  std::sort(BOOST_EigValues.begin(),BOOST_EigValues.end(),dComparator);
  
  double BOOST_Sphericity = 1.5*(BOOST_EigValues[1]+BOOST_EigValues[2]);
  double BOOST_Aplanarity = 1.5*BOOST_EigValues[2];
  
  double Obs14 = ( isnan(BOOST_Sphericity) ? -1 : BOOST_Sphericity );
  evtselectVarVal.push_back(pair<unsigned int,double>(14,Obs14));
  
  double Obs15 = ( isnan(BOOST_Aplanarity) ? -1 : BOOST_Aplanarity );
  evtselectVarVal.push_back(pair<unsigned int,double>(15,Obs15));
  
  
  // Centrality of the six jets
  double H=0;
  for(unsigned int i=0;i<6;i++){
    H += topJets[i].recJet().energy();
  }
  double Obs16 = ( H != 0 ? HT/H : -1 );
  evtselectVarVal.push_back(pair<unsigned int,double>(16,Obs16));
  
  
  // distance from the origin in the (BjetsBdiscSum, Ljets1BdiscSum) and  (BjetsBdiscSum, Ljets2BdiscSum)
  double Obs17 = ( BjetsBdiscSum != 0 && Ljets1BdiscSum != 0 ? 0.707*(BjetsBdiscSum-Ljets1BdiscSum) : -1 );
  evtselectVarVal.push_back(pair<unsigned int,double>(17,Obs17));
  double Obs18 = ( BjetsBdiscSum != 0 && Ljets2BdiscSum != 0 ? 0.707*(BjetsBdiscSum-Ljets2BdiscSum) : -1 );
  evtselectVarVal.push_back(pair<unsigned int,double>(18,Obs18));
  
  // Ratio of the Et Sum of the two lightest jets over HT=Sum of the Et of the six highest jets in Et
  double Obs19 = ( HT != 0 ? EtSum/HT : -1 );
  evtselectVarVal.push_back(pair<unsigned int,double>(19,Obs19));
  
  
  // Put the vector in the TtHadEvtSolution
  TS.setLRSignalEvtObservables(evtselectVarVal);
}
