#include "OHltMenu.h"

using namespace std;

void OHltMenu::AddTrigger(TString trign, int presc, float eventS) {
  names.push_back(trign);
  prescales[trign] 	       	= presc;
  eventSizes[trign] 	       	= eventS;
}

void OHltMenu::AddL1forPreLoop(TString trign, int presc) {
  L1names.push_back(trign);
  L1prescales[trign] 	       	= presc;
}

void OHltMenu::print() {
  cout << "Menu - isL1Menu="<<isL1Menu << " - doL1preloop="<<doL1preloop<<  endl;
  cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" <<  endl;
  for (unsigned int i=0;i<names.size();i++) {
    cout<<names[i]<<" "<<prescales[names[i]]<<" "<<eventSizes[names[i]]<<" "<<endl;
  }
  cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" <<  endl;

  if (doL1preloop) {
    cout << endl << "L1 Menu - for L1preloop"<<  endl;
    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" <<  endl;
    for (unsigned int i=0;i<L1names.size();i++) {
      cout<<L1names[i]<<" "<<L1prescales[L1names[i]]<<endl;
    }
    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" <<  endl;
  }
  cout << endl;
}

void OHltMenu::SetMapL1SeedsOfStandardHLTPath() {
  typedef vector<TString> myvec;
  typedef pair< TString, vector<TString> > mypair;
  typedef map< TString, vector<TString> >  mymap;

  myvec vtmp;  

  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_L1Jet15", myvec(1,"L1_SingleJet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Jet30", myvec(1,"L1_SingleJet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Jet50", myvec(1,"L1_SingleJet10_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Jet80", myvec(1,"L1_SingleJet20_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Jet110", myvec(1,"L1_SingleJet40_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Jet180", myvec(1,"L1_SingleJet40_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Jet250", myvec(1,"L1_SingleJet40_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_FwdJet20", myvec(1,"L1_IsoEG10_Jet15_ForJet10")));

  vtmp.clear(); vtmp.push_back("L1_SingleJet60_00001"); vtmp.push_back("L1_DoubleJet40_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleJet150", vtmp));
  vtmp.clear(); vtmp.push_back("L1_SingleJet60_00001"); vtmp.push_back("L1_DoubleJet40_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleJet125_Aco", vtmp));

  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleFwdJet50", myvec(1,"L1_SingleJet10_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DiJetAve15", myvec(1,"L1_SingleJet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DiJetAve30", myvec(1,"L1_SingleJet10_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DiJetAve50", myvec(1,"L1_SingleJet20_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DiJetAve70", myvec(1,"L1_SingleJet40_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DiJetAve130", myvec(1,"L1_SingleJet40_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DiJetAve220", myvec(1,"L1_SingleJet40_00001")));

  vtmp.clear();
  vtmp.push_back("L1_SingleJet60_00001"); vtmp.push_back("L1_DoubleJet40_00001"); vtmp.push_back("L1_TripleJet30_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_TripleJet85", vtmp));
    
  vtmp.clear();
  vtmp.push_back("L1_SingleJet60_00001"); vtmp.push_back("L1_DoubleJet40_00001");
  vtmp.push_back("L1_TripleJet30_00001"); vtmp.push_back("L1_QuadJet20_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_QuadJet60", vtmp));
    
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_QuadJet30", myvec(1,"L1_QuadJet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_SumET120", myvec(1,"L1_ETT60_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_L1MET20", myvec(1,"L1_ETM20_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_MET25", myvec(1,"L1_ETM20_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_MET35", myvec(1,"L1_ETM30_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_MET50", myvec(1,"L1_ETM40_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_MET65", myvec(1,"L1_ETM50_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_MET75", myvec(1,"L1_ETM50_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_MET35_HT350", myvec(1,"L1_HTT300_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Jet180_MET60", myvec(1,"L1_SingleJet60_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Jet60_MET70_Aco", myvec(1,"L1_SingleJet60_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Jet100_MET60_Aco", myvec(1,"L1_SingleJet60_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleJet125_MET60", myvec(1,"L1_SingleJet60_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleFwdJet40_MET60", myvec(1,"L1_ETM40_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleJet60_MET60_Aco", myvec(1,"L1_SingleJet60_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleJet50_MET70_Aco", myvec(1,"L1_SingleJet60_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleJet40_MET70_Aco", myvec(1,"L1_SingleJet60_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_TripleJet60_MET60", myvec(1,"L1_SingleJet60_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_QuadJet35_MET60", myvec(1,"L1_SingleJet60_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoEle15_L1I", myvec(1,"L1_SingleIsoEG12_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoEle18_L1R", myvec(1,"L1_SingleEG15_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoEle15_LW_L1I", myvec(1,"L1_SingleIsoEG12_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_LooseIsoEle15_LW_L1R", myvec(1,"L1_SingleEG12_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Ele10_SW_L1R", myvec(1,"L1_SingleEG8_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Ele15_SW_L1R", myvec(1,"L1_SingleEG12_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Ele15_LW_L1R", myvec(1,"L1_SingleEG10_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_EM80", myvec(1,"L1_SingleEG15_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_EM200", myvec(1,"L1_SingleEG15_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleIsoEle10_L1I", myvec(1,"L1_DoubleIsoEG8_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleIsoEle12_L1R", myvec(1,"L1_DoubleEG10_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleIsoEle10_LW_L1I", myvec(1,"L1_DoubleIsoEG8_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleIsoEle12_LW_L1R", myvec(1,"L1_DoubleEG10_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleEle5_SW_L1R", myvec(1,"L1_DoubleEG5_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleEle10_LW_OnlyPixelM_L1R", myvec(1,"L1_DoubleEG5_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleEle10_Z", myvec(1,"L1_DoubleIsoEG8_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoPhoton30_L1I", myvec(1,"L1_SingleIsoEG12_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoPhoton10_L1R", myvec(1,"L1_SingleEG8_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoPhoton15_L1R", myvec(1,"L1_SingleEG12_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoPhoton20_L1R", myvec(1,"L1_SingleEG15_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoPhoton25_L1R", myvec(1,"L1_SingleEG15_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoPhoton40_L1R", myvec(1,"L1_SingleEG15_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Photon15_L1R", myvec(1,"L1_SingleEG8_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Photon15_LooseEcalIso_L1R", myvec(1,"L1_SingleEG10_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Photon25_L1R", myvec(1,"L1_SingleEG15_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleIsoPhoton20_L1I", myvec(1,"L1_DoubleIsoEG8_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleIsoPhoton20_L1R", myvec(1,"L1_DoubleEG10_00001")));

  vtmp.clear(); vtmp.push_back("L1_SingleMu7"); vtmp.push_back("L1_DoubleMu3"); 
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_L1Mu", vtmp));
  
  vtmp.clear(); vtmp.push_back("L1_SingleMuOpen");
  vtmp.push_back("L1_SingleMu3"); vtmp.push_back("L1_SingleMu5");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_L1MuOpen", vtmp));
  
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_L2Mu9", myvec(1,"L1_SingleMu7")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoMu9", myvec(1,"L1_SingleMu7")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoMu11", myvec(1,"L1_SingleMu7")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoMu13", myvec(1,"L1_SingleMu10")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoMu15", myvec(1,"L1_SingleMu10")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_NoTrackerIsoMu15", myvec(1,"L1_SingleMu7")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Mu3", myvec(1,"L1_SingleMu3")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Mu5", myvec(1,"L1_SingleMu5")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Mu7", myvec(1,"L1_SingleMu7")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Mu9", myvec(1,"L1_SingleMu7")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Mu11", myvec(1,"L1_SingleMu7")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Mu13", myvec(1,"L1_SingleMu10")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Mu15", myvec(1,"L1_SingleMu10")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Mu15_L1Mu7", myvec(1,"L1_SingleMu7")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Mu15_Vtx2cm", myvec(1,"L1_SingleMu7")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Mu15_Vtx2mm", myvec(1,"L1_SingleMu7")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleIsoMu3", myvec(1,"L1_DoubleMu3")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleMu3", myvec(1,"L1_DoubleMu3")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleMu3_Vtx2cm", myvec(1,"L1_DoubleMu3")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleMu3_Vtx2mm", myvec(1,"L1_DoubleMu3")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleMu3_JPsi", myvec(1,"L1_DoubleMu3")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleMu3_Upsilon", myvec(1,"L1_DoubleMu3")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleMu7_Z", myvec(1,"L1_DoubleMu3")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleMu3_SameSign", myvec(1,"L1_DoubleMu3")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleMu3_Psi2S", myvec(1,"L1_DoubleMu3")));

  vtmp.clear();
  vtmp.push_back("L1_SingleJet60_00001"); vtmp.push_back("L1_DoubleJet60_00001");
  vtmp.push_back("L1_TripleJet30_00001"); vtmp.push_back("L1_QuadJet20_00001"); vtmp.push_back("L1_HTT300_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagIP_Jet180", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleJet60_00001"); vtmp.push_back("L1_DoubleJet40_00001");
  vtmp.push_back("L1_TripleJet30_00001"); vtmp.push_back("L1_QuadJet20_00001"); vtmp.push_back("L1_HTT300_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagIP_Jet120_Relaxed", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleJet60_00001"); vtmp.push_back("L1_DoubleJet60_00001");
  vtmp.push_back("L1_TripleJet30_00001"); vtmp.push_back("L1_QuadJet20_00001"); vtmp.push_back("L1_HTT300_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagIP_DoubleJet120", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleJet60_00001"); vtmp.push_back("L1_DoubleJet40_00001");
  vtmp.push_back("L1_TripleJet30_00001"); vtmp.push_back("L1_QuadJet20_00001"); vtmp.push_back("L1_HTT300_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagIP_DoubleJet60_Relaxed", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleJet60_00001"); vtmp.push_back("L1_DoubleJet60_00001");
  vtmp.push_back("L1_TripleJet30_00001"); vtmp.push_back("L1_QuadJet20_00001"); vtmp.push_back("L1_HTT300_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagIP_TripleJet70", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleJet60_00001"); vtmp.push_back("L1_DoubleJet40_00001");
  vtmp.push_back("L1_TripleJet30_00001"); vtmp.push_back("L1_QuadJet20_00001"); vtmp.push_back("L1_HTT300_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagIP_TripleJet40_Relaxed", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleJet60_00001"); vtmp.push_back("L1_DoubleJet60_00001");
  vtmp.push_back("L1_TripleJet30_00001"); vtmp.push_back("L1_QuadJet20_00001"); vtmp.push_back("L1_HTT300_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagIP_QuadJet40", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleJet60_00001"); vtmp.push_back("L1_DoubleJet40_00001");
  vtmp.push_back("L1_TripleJet30_00001"); vtmp.push_back("L1_QuadJet20_00001"); vtmp.push_back("L1_HTT300_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagIP_QuadJet30_Relaxed", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleJet60_00001"); vtmp.push_back("L1_DoubleJet60_00001");
  vtmp.push_back("L1_TripleJet30_00001"); vtmp.push_back("L1_QuadJet20_00001"); vtmp.push_back("L1_HTT300_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagIP_HT470", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleJet60_00001"); vtmp.push_back("L1_DoubleJet40_00001");
  vtmp.push_back("L1_TripleJet30_00001"); vtmp.push_back("L1_QuadJet20_00001"); vtmp.push_back("L1_HTT300_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagIP_HT320_Relaxed", vtmp));

  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagMu_DoubleJet120", myvec(1,"L1_Mu5_Jet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagMu_DoubleJet60_Relaxed", myvec(1,"L1_Mu5_Jet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagMu_TripleJet70", myvec(1,"L1_Mu5_Jet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagMu_TripleJet40_Relaxed", myvec(1,"L1_Mu5_Jet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagMu_QuadJet40", myvec(1,"L1_Mu5_Jet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagMu_QuadJet30_Relaxed", myvec(1,"L1_Mu5_Jet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagMu_HT370", myvec(1,"L1_HTT300_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagMu_HT250_Relaxed", myvec(1,"L1_HTT200_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleMu3_BJPsi", myvec(1,"L1_DoubleMu3")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleMu4_BJPsi", myvec(1,"L1_DoubleMu3")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_TripleMu3_TauTo3Mu", myvec(1,"L1_DoubleMu3")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoTau_MET65_Trk20", myvec(1,"L1_SingleTauJet50_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoTau_MET35_Trk15_L1MET", myvec(1,"L1_TauJet10_ETM30_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_LooseIsoTau_MET30", myvec(1,"L1_SingleTauJet50_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_LooseIsoTau_MET30_L1MET", myvec(1,"L1_TauJet10_ETM30_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleIsoTau_Trk3", myvec(1,"L1_DoubleTauJet20_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_DoubleLooseIsoTau", myvec(1,"L1_DoubleTauJet8_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoEle8_IsoMu7", myvec(1,"L1_Mu3_IsoEG5_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoEle10_Mu10_L1R", myvec(1,"L1_Mu3_EG12_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoEle12_IsoTau_Trk3", myvec(1,"L1_IsoEG10_TauJet8_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoEle10_BTagIP_Jet35", myvec(1,"L1_IsoEG10_Jet8_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoEle12_Jet40", myvec(1,"L1_IsoEG10_Jet12_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoEle12_DoubleJet80", myvec(1,"L1_IsoEG10_Jet12_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoElec5_TripleJet30", myvec(1,"L1_EG5_TripleJet15")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoEle12_TripleJet60", myvec(1,"L1_IsoEG10_Jet12_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoEle12_QuadJet35", myvec(1,"L1_IsoEG10_Jet12_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoMu14_IsoTau_Trk3", myvec(1,"L1_Mu5_TauJet8_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoMu7_BTagIP_Jet35", myvec(1,"L1_Mu5_Jet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoMu7_BTagMu_Jet20", myvec(1,"L1_Mu5_Jet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_IsoMu7_Jet40", myvec(1,"L1_Mu5_Jet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_NoL2IsoMu8_Jet40", myvec(1,"L1_Mu5_Jet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Mu14_Jet50", myvec(1,"L1_Mu5_Jet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_Mu5_TripleJet30", myvec(1,"L1_Mu3_TripleJet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_BTagMu_Jet20_Calib", myvec(1,"L1_Mu5_Jet6_00001")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_ZeroBias", myvec(1,"OpenL1_ZeroBias")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_MinBias", myvec(1,"L1_MinBias_HTT10_00001")));

  vtmp.clear();
  vtmp.push_back("L1_SingleHfBitCountsRing1_1"); vtmp.push_back("L1_DoubleHfBitCountsRing1_P1N1");
  vtmp.push_back("L1_SingleHfRingEtSumsRing1_4"); vtmp.push_back("L1_DoubleHfRingEtSumsRing1_P4N4");
  vtmp.push_back("L1_SingleHfRingEtSumsRing2_4"); vtmp.push_back("L1_DoubleHfRingEtSumsRing2_P4N4");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_MinBiasHcal", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleEG1"); vtmp.push_back("L1_DoubleEG1_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_MinBiasEcal", vtmp));

  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_MinBiasPixel", myvec(1,"OpenL1_ZeroBias")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_MinBiasPixel_Trk5", myvec(1,"OpenL1_ZeroBias")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_CSCBeamHalo", myvec(1,"L1_SingleMuBeamHalo")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_CSCBeamHaloOverlapRing1", myvec(1,"L1_SingleMuBeamHalo")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_CSCBeamHaloOverlapRing2", myvec(1,"L1_SingleMuBeamHalo")));
  map_L1SeedsOfStandardHLTPath.insert(mypair("HLT_CSCBeamHaloRing2or3", myvec(1,"L1_SingleMuBeamHalo")));

  vtmp.clear();
  vtmp.push_back("L1_SingleJet10_00001");
  vtmp.push_back("L1_SingleJet20_00001"); vtmp.push_back("L1_SingleJet40_00001");
  vtmp.push_back("L1_SingleTauJet10_00001"); vtmp.push_back("L1_SingleTauJet20_00001");
  vtmp.push_back("L1_SingleTauJet30_00001"); vtmp.push_back("L1_SingleTauJet50_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("AlCa_IsoTrack", vtmp));
  
  vtmp.clear();
  vtmp.push_back("OpenL1_ZeroBias");
  vtmp.push_back("L1_SingleEG1"); vtmp.push_back("L1_DoubleEG1_00001");
  vtmp.push_back("L1_SingleHfBitCountsRing1_1"); vtmp.push_back("L1_DoubleHfBitCountsRing1_P1N1");
  vtmp.push_back("L1_SingleHfRingEtSumsRing1_4"); vtmp.push_back("L1_DoubleHfRingEtSumsRing1_P4N4");
  vtmp.push_back("L1_SingleHfRingEtSumsRing2_4"); vtmp.push_back("L1_DoubleHfRingEtSumsRing2_P4N4");
  map_L1SeedsOfStandardHLTPath.insert(mypair("AlCa_EcalPhiSym", vtmp));
  
  vtmp.clear();
  vtmp.push_back("OpenL1_ZeroBias");
  vtmp.push_back("L1_SingleEG1"); vtmp.push_back("L1_DoubleEG1_00001");
  vtmp.push_back("L1_SingleJet20_00001");vtmp.push_back("L1_SingleJet30_00001");
  vtmp.push_back("L1_SingleJet40_00001"); vtmp.push_back("L1_SingleJet50_00001");
  vtmp.push_back("L1_SingleJet60_00001"); vtmp.push_back("L1_SingleJet10_00001");
  vtmp.push_back("L1_SingleIsoEG5_00001"); vtmp.push_back("L1_SingleIsoEG10_00001"); vtmp.push_back("L1_SingleIsoEG12_00001");
  vtmp.push_back("L1_SingleIsoEG15_00001"); vtmp.push_back("L1_SingleIsoEG15"); 
  vtmp.push_back("L1_SingleEG10_00001"); vtmp.push_back("L1_SingleEG12_00001");
  vtmp.push_back("L1_SingleEG15_00001"); vtmp.push_back("L1_SingleEG20_00001"); 
  vtmp.push_back("L1_DoubleEG10_00001"); vtmp.push_back("L1_DoubleEG15_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("AlCa_EcalPi0", vtmp));

    // For measuring L1 rates, also add L1 bits to the map!
  /* New L1s from the L1Menu_startup2_v2 */
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleEG10_00001", myvec(1,"L1_DoubleEG10_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleEG1_00001", myvec(1,"L1_DoubleEG1_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleEG5_00001", myvec(1,"L1_DoubleEG5_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleForJet20", myvec(1,"L1_DoubleForJet20")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleHfBitCountsRing1_P1N1", myvec(1,"L1_DoubleHfBitCountsRing1_P1N1")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleHfBitCountsRing2_P1N1", myvec(1,"L1_DoubleHfBitCountsRing2_P1N1")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleHfRingEtSumsRing1_P200N200", myvec(1,"L1_DoubleHfRingEtSumsRing1_P200N200")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleHfRingEtSumsRing1_P4N4", myvec(1,"L1_DoubleHfRingEtSumsRing1_P4N4")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleHfRingEtSumsRing2_P200N200", myvec(1,"L1_DoubleHfRingEtSumsRing2_P200N200")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleHfRingEtSumsRing2_P4N4", myvec(1,"L1_DoubleHfRingEtSumsRing2_P4N4")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleIsoEG05_TopBottom", myvec(1,"L1_DoubleIsoEG05_TopBottom")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleIsoEG05_TopBottomCen", myvec(1,"L1_DoubleIsoEG05_TopBottomCen")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleIsoEG10_00001", myvec(1,"L1_DoubleIsoEG10_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleIsoEG8_00001", myvec(1,"L1_DoubleIsoEG8_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleJet40_00001", myvec(1,"L1_DoubleJet40_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleJet60_00001", myvec(1,"L1_DoubleJet60_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleMu3", myvec(1,"L1_DoubleMu3")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleMuOpen", myvec(1,"L1_DoubleMuOpen")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleMuTopBottom", myvec(1,"L1_DoubleMuTopBottom")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleNoIsoEG05_TopBottom", myvec(1,"L1_DoubleNoIsoEG05_TopBottom")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleNoIsoEG05_TopBottomCen", myvec(1,"L1_DoubleNoIsoEG05_TopBottomCen")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleTauJet20_00001", myvec(1,"L1_DoubleTauJet20_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_DoubleTauJet8_00001", myvec(1,"L1_DoubleTauJet8_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_EG12_Jet40_00001", myvec(1,"L1_EG12_Jet40_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_EG5_TripleJet6_00001", myvec(1,"L1_EG5_TripleJet6_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_ETM20_00001", myvec(1,"L1_ETM20_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_ETM30_00001", myvec(1,"L1_ETM30_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_ETM40_00001", myvec(1,"L1_ETM40_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_ETM50_00001", myvec(1,"L1_ETM50_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_ETT60_00001", myvec(1,"L1_ETT60_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_HTT100_00001", myvec(1,"L1_HTT100_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_HTT200_00001", myvec(1,"L1_HTT200_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_HTT300_00001", myvec(1,"L1_HTT300_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_IsoEG10_Jet12_00001", myvec(1,"L1_IsoEG10_Jet12_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_IsoEG10_Jet6_00001", myvec(1,"L1_IsoEG10_Jet6_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_IsoEG10_Jet6_ForJet6_00001", myvec(1,"L1_IsoEG10_Jet6_ForJet6_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_IsoEG10_Jet8_00001", myvec(1,"L1_IsoEG10_Jet8_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_IsoEG10_TauJet8_00001", myvec(1,"L1_IsoEG10_TauJet8_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_MinBias_ETT10_00001", myvec(1,"L1_MinBias_ETT10_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_MinBias_HTT10_00001", myvec(1,"L1_MinBias_HTT10_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_Mu3_EG12_00001", myvec(1,"L1_Mu3_EG12_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_Mu3_IsoEG5_00001", myvec(1,"L1_Mu3_IsoEG5_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_Mu3_Jet6_00001", myvec(1,"L1_Mu3_Jet6_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_Mu3_TripleJet6_00001", myvec(1,"L1_Mu3_TripleJet6_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_Mu5_IsoEG10_00001", myvec(1,"L1_Mu5_IsoEG10_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_Mu5_Jet6_00001", myvec(1,"L1_Mu5_Jet6_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_Mu5_TauJet8_00001", myvec(1,"L1_Mu5_TauJet8_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_QuadJet20_00001", myvec(1,"L1_QuadJet20_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_QuadJet6_00001", myvec(1,"L1_QuadJet6_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleEG1", myvec(1,"L1_SingleEG1")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleEG10_00001", myvec(1,"L1_SingleEG10_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleEG12_00001", myvec(1,"L1_SingleEG12_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleEG15_00001", myvec(1,"L1_SingleEG15_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleEG20_00001", myvec(1,"L1_SingleEG20_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleEG5_00001", myvec(1,"L1_SingleEG5_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleEG5_Endcap_00001", myvec(1,"L1_SingleEG5_Endcap_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleEG8_00001", myvec(1,"L1_SingleEG8_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleForJet10", myvec(1,"L1_SingleForJet10")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleForJet6", myvec(1,"L1_SingleForJet6")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleHfBitCountsRing1_1", myvec(1,"L1_SingleHfBitCountsRing1_1")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleHfBitCountsRing2_1", myvec(1,"L1_SingleHfBitCountsRing2_1")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleHfRingEtSumsRing1_200", myvec(1,"L1_SingleHfRingEtSumsRing1_200")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleHfRingEtSumsRing1_4", myvec(1,"L1_SingleHfRingEtSumsRing1_4")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleHfRingEtSumsRing2_200", myvec(1,"L1_SingleHfRingEtSumsRing2_200")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleHfRingEtSumsRing2_4", myvec(1,"L1_SingleHfRingEtSumsRing2_4")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleIsoEG10_00001", myvec(1,"L1_SingleIsoEG10_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleIsoEG12_00001", myvec(1,"L1_SingleIsoEG12_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleIsoEG15_00001", myvec(1,"L1_SingleIsoEG15_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleIsoEG5_00001", myvec(1,"L1_SingleIsoEG5_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleIsoEG5_Endcap_00001", myvec(1,"L1_SingleIsoEG5_Endcap_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleIsoEG8_00001", myvec(1,"L1_SingleIsoEG8_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet10_00001", myvec(1,"L1_SingleJet10_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet10_Barrel_00001", myvec(1,"L1_SingleJet10_Barrel_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet10_Central", myvec(1,"L1_SingleJet10_Central")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet10_Endcap", myvec(1,"L1_SingleJet10_Endcap")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet20_00001", myvec(1,"L1_SingleJet20_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet20_Barrel_00001", myvec(1,"L1_SingleJet20_Barrel_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet30_00001", myvec(1,"L1_SingleJet30_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet30_Barrel_00001", myvec(1,"L1_SingleJet30_Barrel_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet40_00001", myvec(1,"L1_SingleJet40_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet40_Barrel_00001", myvec(1,"L1_SingleJet40_Barrel_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet50_00001", myvec(1,"L1_SingleJet50_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet60_00001", myvec(1,"L1_SingleJet60_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet6_00001", myvec(1,"L1_SingleJet6_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet6_Barrel_00001", myvec(1,"L1_SingleJet6_Barrel_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet6_Central", myvec(1,"L1_SingleJet6_Central")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleJet6_Endcap", myvec(1,"L1_SingleJet6_Endcap")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleMu0", myvec(1,"L1_SingleMu0")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleMu10", myvec(1,"L1_SingleMu10")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleMu14", myvec(1,"L1_SingleMu14")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleMu3", myvec(1,"L1_SingleMu3")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleMu5", myvec(1,"L1_SingleMu5")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleMu7", myvec(1,"L1_SingleMu7")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleMuBeamHalo", myvec(1,"L1_SingleMuBeamHalo")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleMuOpen", myvec(1,"L1_SingleMuOpen")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleTauJet10_00001", myvec(1,"L1_SingleTauJet10_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleTauJet10_Barrel_00001", myvec(1,"L1_SingleTauJet10_Barrel_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleTauJet20_00001", myvec(1,"L1_SingleTauJet20_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleTauJet20_Barrel_00001", myvec(1,"L1_SingleTauJet20_Barrel_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleTauJet30_00001", myvec(1,"L1_SingleTauJet30_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleTauJet30_Barrel_00001", myvec(1,"L1_SingleTauJet30_Barrel_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleTauJet50_00001", myvec(1,"L1_SingleTauJet50_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleTauJet8_00001", myvec(1,"L1_SingleTauJet8_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_SingleTauJet8_Barrel_00001", myvec(1,"L1_SingleTauJet8_Barrel_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_TauJet10_ETM30_00001", myvec(1,"L1_TauJet10_ETM30_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_TauJet10_ETM40_00001", myvec(1,"L1_TauJet10_ETM40_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_TripleJet30_00001", myvec(1,"L1_TripleJet30_00001")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("L1_TripleMu3", myvec(1,"L1_TripleMu3")));  


  /* New Taus */
  // This is openhlt and not standard hlt,
  // but the same mechanism can also be used here!
  // L1 prescales can be checked in the same way as 
  // for standard hlt in CheckOpenHlt(). 
  // Look for "New Taus" in OHltTreeOpen.cpp!
  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet20_00001"); vtmp.push_back("L1_SingleJet30_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau20", vtmp));

  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet30_00001"); vtmp.push_back("L1_SingleJet60_00001");
  //vtmp.push_back("L1_SingleTauJet50_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau30_Trk5", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet20_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau20_Trk5", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet30_00001"); vtmp.push_back("L1_SingleJet60_00001");
  //vtmp.push_back("L1_SingleTauJet50_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau30_Trk10", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet20_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau20_Trk10", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_DoubleTauJet20_00001"); vtmp.push_back("L1_DoubleJet40_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_DoubleLooseIsoTau15", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_DoubleTauJet20_00001"); vtmp.push_back("L1_DoubleJet40_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_DoubleLooseIsoTau15_Trk5", vtmp));

  
  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet20_00001"); vtmp.push_back("L1_SingleJet40_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau20_L2R", vtmp));

  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet30_00001"); vtmp.push_back("L1_SingleJet60_00001");
  //vtmp.push_back("L1_SingleTauJet50_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau30_Trk5_L2R", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet20_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau20_Trk5_L2R", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet30_00001"); vtmp.push_back("L1_SingleJet60_00001");
  //vtmp.push_back("L1_SingleTauJet50_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau30_Trk10_L2R", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet20_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau20_Trk10_L2R", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_DoubleTauJet20_00001"); vtmp.push_back("L1_DoubleJet40_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_DoubleLooseIsoTau15_L2R", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_DoubleTauJet20_00001"); vtmp.push_back("L1_DoubleJet40_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_DoubleLooseIsoTau15_Trk5_L2R", vtmp));

  // L3 iso
  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet30_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau30_Trk5_L3I", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet30_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau30_Trk5_L2R_L3I", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet30_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau30_Trk10_L3I", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet30_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau30_Trk10_L2R_L3I", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_DoubleTauJet20_00001"); vtmp.push_back("L1_DoubleJet40_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_DoubleLooseIsoTau15_Trk5_L3I", vtmp));

  vtmp.clear();
  vtmp.push_back("L1_DoubleTauJet20_00001"); vtmp.push_back("L1_DoubleJet40_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_DoubleLooseIsoTau15_Trk5_L2R_L3I", vtmp));

  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet20_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau20_Trk5_L3I", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet20_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau20_Trk5_L2R_L3I", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet20_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau20_Trk10_L3I", vtmp));
  
  vtmp.clear();
  vtmp.push_back("L1_SingleTauJet20_00001"); vtmp.push_back("L1_SingleJet60_00001");
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_SingleLooseIsoTau20_Trk10_L2R_L3I", vtmp));
  
  vtmp.clear(); 
  vtmp.push_back("L1_SingleTauJet20_00001"); vtmp.push_back("L1_SingleMu14"); 
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_L1Mu14", vtmp)); 


  // Muons
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_IsoMu3", myvec(1,"L1_SingleMu3"))); 
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_IsoMu9", myvec(1,"L1_SingleMu7")));  
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_Mu5", myvec(1,"L1_SingleMu3")));   
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_L2Mu11", myvec(1,"L1_SingleMu7")));   
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_L1Mu20", myvec(1,"L1_SingleMu20")));   
  map_L1SeedsOfStandardHLTPath.insert(mypair("OpenHLT_L2Mu9_1Jet30", myvec(1,"L1_Mu5_Jet6_00001")));   

}
