#include "FastSimulation/Validation/test/JetComparison.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"


JetComparison::JetComparison(edm::ParameterSet const& conf)
{

  fMinEnergy  = conf.getParameter<double>("MinEnergy");





  nEvent=0;
  //  output file
  outputFile_ = conf.getUntrackedParameter<string>("outputFile", "myfile.root");

  
  if ( outputFile_.size() != 0 ) {
    LogInfo("OutputInfo") << " jet histograms will be saved to '" << outputFile_.c_str() << "'";
  } else {
    LogInfo("OutputInfo") << " jet histograms will NOT be saved";
  }
  

  Char_t histo[20];

  sprintf (histo, "number_of_jets" ) ;
  meNumberJet = new TH1F(histo, histo, 500 , 0. , 500.);   
  
  sprintf (histo, "jet_et" ) ;
  meEtJet = new TH1F(histo, histo, 1000 , 0. , 1000.);   

  sprintf (histo, "gen_et" ) ;
  meEtGen = new TH1F(histo, histo, 1000 , 0. , 1000.);   

  sprintf (histo, "jet_et_two_matched" ) ;
  meEtJetMatched = new TH1F(histo, histo, 1000 , 0. , 1000.);   

  sprintf (histo, "gen_et_two_matched" ) ;
  meEtGenMatched = new TH1F(histo, histo, 1000 , 0. , 1000.);   

  sprintf (histo, "jet_eta" ) ;
  meEtaJet = new TH1F(histo, histo, 500 , 0. , 5.);   

  sprintf (histo, "gen_eta" ) ;
  meEtaGen = new TH1F(histo, histo, 500 , 0. , 5.);   
    
  sprintf (histo, "ratio_jet_gen_vs_eta" ) ;
  meRatio =  new TH2F(histo, histo, 50,  0., 5., 20000 , 0., 5.);   

  sprintf (histo, "hadronic_fraction_vs_eta" ) ;
  meHadronicFrac_vs_eta =  new TH2F(histo, histo, 50,  0., 5., 1000 , 0., 1.);   

  sprintf (histo, "number_of_towers_with_60_percents_energy_vs_eta" ) ;
  meNTowers60_vs_eta =  new TH2F(histo, histo, 50,  0., 5., 100 , -0.5, 99.5);   

  sprintf (histo, "number_of_towers_with_90_percents_energy_vs_eta" ) ;
  meNTowers90_vs_eta =  new TH2F(histo, histo, 50,  0., 5., 100 , -0.5, 99.5);  

  sprintf (histo, "distance_R_between_matched_jets") ;
  meDistR =  new TH1F(histo, histo, 300,0.,0.3);  

  sprintf (histo, "distance_R_between_matched_jets_vs_eta") ;
  meDistR_vs_eta =  new TH2F(histo, histo,50,  0., 5., 300,0.,0.3);  
    
  sprintf (histo, "number_of_calotowers_vs_eta") ;
  meNTowers_vs_eta =  new TH2F(histo, histo,50,  0., 5., 51, -0.5, 50.5);  

  fFile = new TFile(outputFile_.c_str(),"RECREATE");
}


JetComparison::~JetComparison() {
   
  cout << " outputFile_.size() =  " << outputFile_.size() << endl;
  
}

void JetComparison::endJob() {
  
  cout << " outputFile_.size() =  " << outputFile_.size() << endl;

  fFile->cd();

  meEtJet->Write();
  meEtGen ->Write();
  meEtJetMatched ->Write();
  meEtGenMatched ->Write();
  meEtaJet->Write();
  meEtaGen->Write();
  meRatio->Write();
    
  meNumberJet->Scale(1./nEvent);
  meNumberJet->Write();
  meHadronicFrac_vs_eta->Write();
  meNTowers60_vs_eta->Write();
  meNTowers90_vs_eta->Write();
  meDistR->Write();
  meDistR_vs_eta->Write();
  meNTowers_vs_eta->Write();
  Float_t   x[50];
  Float_t   y[50];
  Float_t y_er[50];
  Float_t x_er[50];
  Float_t mean;
  Float_t sigm;

  for (Int_t i=1;i<51;i++){
    mean= (meRatio->ProjectionY("",i,i,"")->GetMean());
    sigm= (meRatio->ProjectionY("",i,i,"")->GetRMS());
    cout << " Mean " << mean <<endl;
    x[i-1] = (i-0.5)*0.1;
    y[i-1] = mean;
    y_er[i-1]=sigm;
    x_er[i-1]=0.;
  }

  gr = new TGraphErrors(50,x,y,x_er,y_er);
  gr->SetMarkerStyle(21);
  gr->SetMarkerColor(44);
  gr->Write();


  fFile->Close();
}
void JetComparison::beginJob(const edm::EventSetup& c){

}
void JetComparison::analyze(edm::Event const& event, edm::EventSetup const& c) {

  nEvent++;

  Handle<CaloJetCollection> jets;
  event.getByLabel("iterativeCone5CaloJets", jets);
  meNumberJet->Fill(jets->size());

  for(unsigned int ijet=0; ijet< jets->size();ijet++)
    {
      //  cout << " JETS energy = " << (*jets)[ijet].et() <<endl;
      double etJet = (*jets)[ijet].et();
      double etaJet = (*jets)[ijet].eta();
      
      meEtJet->Fill(etJet);
      meEtaJet->Fill(etaJet);
    
    }


  Handle<GenJetCollection> jetsgen;
  event.getByLabel("iterativeCone5GenJets", jetsgen);
  for(unsigned int igen =0; igen < jetsgen->size();igen++)
    {
      //     cout << " GENS energy = " << (*jetsgen)[igen].et() <<endl;
      double etGen = (*jetsgen)[igen].et();
      double etaGen = (*jetsgen)[igen].eta();
      meEtGen->Fill(etGen);
      meEtaGen->Fill(etaGen);
    }

  for(unsigned int ijet=0; ijet< jets->size();ijet++)
    {
     
      double etaJet = (*jets)[ijet].eta();
      double phiJet = (*jets)[ijet].phi();
      double etJet = (*jets)[ijet].et();
      double fracHadr = (*jets)[ijet].energyFractionHadronic();
      int n60 = (*jets)[ijet].n60();
      int n90 = (*jets)[ijet].n90();
      int nCaloTower = (*jets)[ijet].nConstituents ();
      double etaGen;
      double phiGen;
      double etGen;
    
     
      unsigned int isize = jetsgen->size();
      if (isize>2) isize = 2;
      for(unsigned int igen =0; igen < isize; igen++)
	{
	  etaGen = (*jetsgen)[igen].eta();
	  phiGen = (*jetsgen)[igen].phi();
	  etGen = (*jetsgen)[igen].et();
	  
	  double rDist = sqrt(deltaR2(etaJet, phiJet, etaGen, phiGen));

	  if (rDist<0.3 && etGen > fMinEnergy ) {
	    meNTowers_vs_eta->Fill(fabs(etaGen),nCaloTower);
	    meDistR->Fill(rDist);
	    meDistR_vs_eta->Fill(fabs(etaGen), rDist); 
	    cout << "  JETS energy = " << (*jets)[ijet].et() << " ETA = " << etaJet << " PHI = "<<phiJet <<endl;
	    cout << "  GENS energy = " << (*jetsgen)[igen].et() << " ETA = " << etaGen << " PHI = "<<phiGen <<endl;
	    cout << "Jets matched " << endl;
	    meRatio->Fill(fabs(etaGen),etJet/etGen);
	    meEtJetMatched->Fill(etJet);
	    meEtGenMatched->Fill(etGen);
	    meNTowers60_vs_eta->Fill(fabs(etaJet),  n60); 
	    meNTowers90_vs_eta->Fill(fabs(etaJet),  n90); 
	    meHadronicFrac_vs_eta -> Fill(fabs(etaJet), fracHadr);
	  }
	}
      


    }
    

}

double JetComparison::deltaR2(double eta0, double phi0, double eta, double phi){
  double dphi=phi-phi0;
  if(dphi>M_PI) dphi-=2*M_PI;
  else if(dphi<=-M_PI) dphi+=2*M_PI;
  return dphi*dphi+(eta-eta0)*(eta-eta0);
}



#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(JetComparison);

 
