#include "Validation/HcalRecHits/interface/HcalRecHitsValidation.h"


HcalRecHitsValidation::HcalRecHitsValidation(edm::ParameterSet const& conf) {
  // DQM ROOT output
  outputFile_ = conf.getUntrackedParameter<std::string>("outputFile", "myfile.root");

  if ( outputFile_.size() != 0 ) {
    edm::LogInfo("OutputInfo") << " Hcal RecHit Task histograms will be saved to '" << outputFile_.c_str() << "'";
  } else {
    edm::LogInfo("OutputInfo") << " Hcal RecHit Task histograms will NOT be saved";
  }

  nevtot = 0;

  dbe_ = 0;
  // get hold of back-end interface
  dbe_ = edm::Service<DQMStore>().operator->();
   
  Char_t histo[20];

  hcalselector_ = conf.getUntrackedParameter<std::string>("hcalselector", "all");
  ecalselector_ = conf.getUntrackedParameter<std::string>("ecalselector", "yes");
  eventype_     = conf.getUntrackedParameter<std::string>("eventype", "single");
  sign_         = conf.getUntrackedParameter<std::string>("sign", "*");
  mc_           = conf.getUntrackedParameter<std::string>("mc", "yes");
  famos_        = conf.getUntrackedParameter<bool>("Famos", false);
  //  std::cout << "*** famos_ = " << famos_ << std::endl; 

  subdet_ = 5;
  if (hcalselector_ == "noise") subdet_ = 0;
  if (hcalselector_ == "HB"   ) subdet_ = 1;
  if (hcalselector_ == "HE"   ) subdet_ = 2;
  if (hcalselector_ == "HO"   ) subdet_ = 3;
  if (hcalselector_ == "HF"   ) subdet_ = 4;
  if (hcalselector_ == "all"  ) subdet_ = 5;
  if (hcalselector_ == "ZS"   ) subdet_ = 6;

  etype_ = 1;
  if (eventype_ == "multi") etype_ = 2;

  iz = 1;
  if(sign_ == "-") iz = -1;
  if(sign_ == "*") iz = 0;

  imc = 1;
  if(mc_ == "no") imc = 0;

  if ( dbe_ ) {
    std::cout << " dbe_->setCurrentFolder" << std::endl; 
    dbe_->setCurrentFolder("HcalRecHitsV/HcalRecHitTask");


    // General counters
    sprintf  (histo, "N_HB" );
    Nhb = dbe_->book1D(histo, histo, 2600,0.,2600.);
    sprintf  (histo, "N_HE" );
    Nhe = dbe_->book1D(histo, histo, 2600,0.,2600.);
    sprintf  (histo, "N_HO" );
    Nho = dbe_->book1D(histo, histo, 2200,0.,2200.);
    sprintf  (histo, "N_HF" );
    Nhf = dbe_->book1D(histo, histo, 1800,0., 1800.);

    // ZS
    if(subdet_ == 6) {

      for (unsigned int i1 = 0;  i1 < 82; i1++) {
	for (unsigned int i2 = 0;  i2 < 72; i2++) {
	  for (unsigned int i3 = 0;  i3 < 4;  i3++) {
	    for (unsigned int i4 = 0;  i4 < 4;  i4++) {
	      emap_min [i1][i2][i3][i4] = 99999.;     
	    }
	  }
	}
      }
      
      sprintf  (histo, "ZSmin_map_depth1" );
      map_depth1 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      sprintf  (histo, "ZSmin_map_depth2" );
      map_depth2 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      sprintf  (histo, "ZSmin_map_depth3" );
      map_depth3 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      sprintf  (histo, "ZSmin_map_depth4" );
      map_depth4 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);

      
      sprintf  (histo, "ZS_Nreco_HB1" );
      ZS_nHB1 = dbe_->book1D(histo, histo, 2500, 0., 2500.);
      sprintf  (histo, "ZS_Nreco_HB2" );
      ZS_nHB2 = dbe_->book1D(histo, histo,  500, 0.,  500.);
      sprintf  (histo, "ZS_Nreco_HE1" );
      ZS_nHE1 = dbe_->book1D(histo, histo, 2000, 0., 2000.);
      sprintf  (histo, "ZS_Nreco_HE2" );
      ZS_nHE2 = dbe_->book1D(histo, histo, 2000, 0., 2000.);
      sprintf  (histo, "ZS_Nreco_HE3" );
      ZS_nHE3 = dbe_->book1D(histo, histo,  500, 0.,  500.);
      sprintf  (histo, "ZS_Nreco_HO" );
      ZS_nHO  = dbe_->book1D(histo, histo, 2500, 0., 2500.);
      sprintf  (histo, "ZS_Nreco_HF1" );
      ZS_nHF1 = dbe_->book1D(histo, histo, 1000, 0., 1000.);
      sprintf  (histo, "ZS_Nreco_HF2" );
      ZS_nHF2 = dbe_->book1D(histo, histo, 1000, 0., 1000.);
      
      sprintf  (histo, "ZSmin_simple1D_HB1" );
      ZS_HB1 = dbe_->book1D(histo, histo,120, -2., 10.);
      sprintf  (histo, "ZSmin_simple1D_HB2" );
      ZS_HB2 = dbe_->book1D(histo, histo,120, -2., 10.);
      sprintf  (histo, "ZSmin_simple1D_HE1" );
      ZS_HE1 = dbe_->book1D(histo, histo,120, -2., 10.);
      sprintf  (histo, "ZSmin_simple1D_HE2" );
      ZS_HE2 = dbe_->book1D(histo, histo,120, -2., 10.);
      sprintf  (histo, "ZSmin_simple1D_HE3" );
      ZS_HE3 = dbe_->book1D(histo, histo,120, -2., 10.);
      sprintf  (histo, "ZSmin_simple1D_HO" );
      ZS_HO = dbe_->book1D(histo, histo,120, -2., 10.);
      sprintf  (histo, "ZSmin_simple1D_HF1" );
      ZS_HF1 = dbe_->book1D(histo, histo,200, -10., 10.);
      sprintf  (histo, "ZSmin_simple1D_HF2" );
      ZS_HF2 = dbe_->book1D(histo, histo,200, -10., 10.);

      sprintf  (histo, "ZSmin_sequential1D_HB1" );
      ZS_seqHB1 = dbe_->book1D(histo, histo,2400, -1200., 1200.);
      sprintf  (histo, "ZSmin_sequential1D_HB2" );
      ZS_seqHB2 = dbe_->book1D(histo, histo,2400, -1200., 1200.);
      sprintf  (histo, "ZSmin_sequential1D_HE1" );
      ZS_seqHE1 = dbe_->book1D(histo, histo,4400, -2200., 2200.);
      sprintf  (histo, "ZSmin_sequential1D_HE2" );
      ZS_seqHE2 = dbe_->book1D(histo, histo,4400, -2200., 2200.);
      sprintf  (histo, "ZSmin_sequential1D_HE3" );
      ZS_seqHE3 = dbe_->book1D(histo, histo,4400, -2200., 2200.);
      sprintf  (histo, "ZSmin_sequential1D_HO" );
      ZS_seqHO  = dbe_->book1D(histo, histo,2400, -1200., 1200.);
      sprintf  (histo, "ZSmin_sequential1D_HF1" );
      ZS_seqHF1 = dbe_->book1D(histo, histo,6000, -3000., 3000.);
      sprintf  (histo, "ZSmin_sequential1D_HF2" );
      ZS_seqHF2 = dbe_->book1D(histo, histo,6000, -3000., 3000.);
    }

    // ALL others, except ZS
    else {
      if (ecalselector_ == "yes") {
	sprintf  (histo, "map_ecal" );
	map_ecal = dbe_->book2D(histo, histo, 70, -3.045, 3.045, 72, -3.1415926536, 3.1415926536);
      }
      
      sprintf  (histo, "emap_HBdepth1" );
      emap_HBdepth1 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      sprintf  (histo, "emap_HBdepth2" );
      emap_HBdepth2 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      sprintf  (histo, "emap_HEdepth1" );
      emap_HEdepth1 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      sprintf  (histo, "emap_HEdepth2" );
      emap_HEdepth2 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      sprintf  (histo, "emap_HEdepth3" );
      emap_HEdepth3 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      sprintf  (histo, "emap_HO" );
      emap_HO = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      sprintf  (histo, "emap_HFdepth1" );
      emap_HFdepth1 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      sprintf  (histo, "emap_HFdepth2" );
      emap_HFdepth2 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      
      sprintf  (histo, "emean_vs_ieta_HBdepth1" );
      emean_vs_ieta_HBdepth1 = dbe_->bookProfile(histo, histo, 82, -41., 41., 2010, -10.,2000. );
      sprintf  (histo, "emean_vs_ieta_HBdepth2" );
      emean_vs_ieta_HBdepth2 = dbe_->bookProfile(histo, histo, 82, -41., 41., 2010, -10.,2000. );
      sprintf  (histo, "emean_vs_ieta_HEdepth1" );
      emean_vs_ieta_HEdepth1 = dbe_->bookProfile(histo, histo, 82, -41., 41., 2010, -10.,2000. );
      sprintf  (histo, "emean_vs_ieta_HEdepth2" );
      emean_vs_ieta_HEdepth2 = dbe_->bookProfile(histo, histo, 82, -41., 41., 2010, -10.,2000. );
      sprintf  (histo, "emean_vs_ieta_HEdepth3" );
      emean_vs_ieta_HEdepth3 = dbe_->bookProfile(histo, histo, 82, -41., 41., 2010, -10.,2000. );
      sprintf  (histo, "emean_vs_ieta_HO" );
      emean_vs_ieta_HO = dbe_->bookProfile(histo, histo, 82, -41., 41., 2010, -10.,2000. );
      sprintf  (histo, "emean_vs_ieta_HFdepth1" );
      emean_vs_ieta_HFdepth1 = dbe_->bookProfile(histo, histo, 82, -41., 41., 2010, -10.,2000. );
      sprintf  (histo, "emean_vs_ieta_HFdepth2" );
      emean_vs_ieta_HFdepth2 = dbe_->bookProfile(histo, histo, 82, -41., 41., 2010, -10.,2000. );



      sprintf  (histo, "occupancy_map_HBdepth1" );
      occupancy_map_HBdepth1 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      sprintf  (histo, "occupancy_map_HBdepth2" );
      occupancy_map_HBdepth2 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      sprintf  (histo, "occupancy_map_HEdepth1" );
      occupancy_map_HEdepth1 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      sprintf  (histo, "occupancy_map_HEdepth2" );
      occupancy_map_HEdepth2 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);      
      sprintf  (histo, "occupancy_map_HEdepth3" );
      occupancy_map_HEdepth3 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      sprintf  (histo, "occupancy_map_HO" );
      occupancy_map_HO = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);      
      sprintf  (histo, "occupancy_map_HFdepth1" );
      occupancy_map_HFdepth1 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      sprintf  (histo, "occupancy_map_HFdepth2" );
      occupancy_map_HFdepth2 = dbe_->book2D(histo, histo, 82, -41., 41., 72, 0., 72.);
      
      sprintf  (histo, "occupancy_vs_ieta_HBdepth1" );
      occupancy_vs_ieta_HBdepth1 = dbe_->book1D(histo, histo, 82, -41., 41.);
      sprintf  (histo, "occupancy_vs_ieta_HBdepth2" );
      occupancy_vs_ieta_HBdepth2 = dbe_->book1D(histo, histo, 82, -41., 41.);
      sprintf  (histo, "occupancy_vs_ieta_HEdepth1" );
      occupancy_vs_ieta_HEdepth1 = dbe_->book1D(histo, histo, 82, -41., 41.);
      sprintf  (histo, "occupancy_vs_ieta_HEdepth2" );
      occupancy_vs_ieta_HEdepth2 = dbe_->book1D(histo, histo, 82, -41., 41.);
      sprintf  (histo, "occupancy_vs_ieta_HEdepth3" );
      occupancy_vs_ieta_HEdepth3 = dbe_->book1D(histo, histo, 82, -41., 41.);
      sprintf  (histo, "occupancy_vs_ieta_HO" );
      occupancy_vs_ieta_HO = dbe_->book1D(histo, histo, 82, -41., 41.);
      sprintf  (histo, "occupancy_vs_ieta_HFdepth1" );
      occupancy_vs_ieta_HFdepth1 = dbe_->book1D(histo, histo, 82, -41., 41.);
      sprintf  (histo, "occupancy_vs_ieta_HFdepth2" );
      occupancy_vs_ieta_HFdepth2 = dbe_->book1D(histo, histo, 82, -41., 41.);
      
      
      sprintf  (histo, "profile_z1" );
      profile_z1 = dbe_->bookProfile(histo, histo, 82, -41., 41., 2000, -2000., 2000.);
      sprintf  (histo, "profile_z2" );
      profile_z2 = dbe_->bookProfile(histo, histo, 82, -41., 41., 2000, -2000., 2000.);
      sprintf  (histo, "profile_z3" );
      profile_z3 = dbe_->bookProfile(histo, histo, 82, -41., 41., 2000, -2000., 2000.);
      sprintf  (histo, "profile_z4" );
      profile_z4 = dbe_->bookProfile(histo, histo, 82, -41., 41., 2000, -2000., 2000.);

      if(imc !=0) { 
	sprintf  (histo, "map_econe_depth1" );
	map_econe_depth1 =
	  dbe_->book2D(histo, histo, 520, -5.2, 5.2, 72, -3.1415926536, 3.1415926536);
	sprintf  (histo, "map_econe_depth2" );
	map_econe_depth2 =
	  dbe_->book2D(histo, histo, 520, -5.2, 5.2, 72, -3.1415926536, 3.1415926536);
	sprintf  (histo, "map_econe_depth3" );
	map_econe_depth3 =
	  dbe_->book2D(histo, histo, 520, -5.2, 5.2, 72, -3.1415926536, 3.1415926536);
	sprintf  (histo, "map_econe_depth4" );
	map_econe_depth4 =
	  dbe_->book2D(histo, histo, 520, -5.2, 5.2, 72, -3.1415926536, 3.1415926536);
      }
    }  // end-of (subdet_ =! 6)

    //======================= Now various cases one by one ===================

    if(subdet_ != 0) { // just not for noise  
      
      //    meEnConeEtaProfiel_depth1->Fill(eta_RecHit, HcalCone_d1);
      
      sprintf (histo, "HcalRecHitTask_En_rechits_cone_profile_vs_ieta_depth1");
      meEnConeEtaProfile_depth1 = dbe_->bookProfile(histo, histo, 82, -41., 41., 210, -10., 200.);   
      
      sprintf (histo, "HcalRecHitTask_En_rechits_cone_profile_vs_ieta_depth2");
      meEnConeEtaProfile_depth2 = dbe_->bookProfile(histo, histo, 82, -41., 41., 210, -10., 200.);  
      
      sprintf (histo, "HcalRecHitTask_En_rechits_cone_profile_vs_ieta_depth3");
      meEnConeEtaProfile_depth3 = dbe_->bookProfile(histo, histo, 82, -41., 41., 210, -10., 200.);  
      
      sprintf (histo, "HcalRecHitTask_En_rechits_cone_profile_vs_ieta_depth4");
      meEnConeEtaProfile_depth4 = dbe_->bookProfile(histo, histo, 82, -41., 41., 210, -10., 200.);  
      
    }
    
    if(etype_ == 1 && subdet_ != 0) { // single part., not for noise
      
      sprintf  (histo, "Delta_phi_cluster-MC");
      meDeltaPhi =  dbe_->book2D(histo, histo, 520, -5.2, 5.2, 60, -0.6, 0.6);
      
      sprintf  (histo, "Delta_eta_cluster-MC");
      meDeltaEta =  dbe_->book2D(histo, histo, 520, -5.2, 5.2, 60, -0.6, 0.6);
      
      sprintf  (histo, "Delta_phi_simcluster-MC");
      meDeltaPhiS =  dbe_->book2D(histo, histo, 520, -5.2, 5.2, 60, -0.6, 0.6);
      
      sprintf  (histo, "Delta_eta_simcluster-MC");
      meDeltaEtaS =  dbe_->book2D(histo, histo, 520, -5.2, 5.2, 60, -0.6, 0.6);
      
      
    }
    // NOISE-specific
    
    if (hcalselector_ == "noise" ){
      
      sprintf  (histo, "e_hb" ) ;
      e_hb = dbe_->book1D(histo, histo,1000, -5., 5.);
      sprintf  (histo, "e_he" ) ;
      e_he = dbe_->book1D(histo, histo,1000, -5., 5.);
      sprintf  (histo, "e_ho" ) ;
      e_ho = dbe_->book1D(histo, histo,1000, -5., 5.);
      sprintf  (histo, "e_hfl" ) ;
      e_hfl = dbe_->book1D(histo, histo,2000, -10., 10.);
      sprintf  (histo, "e_hfs" ) ;
      e_hfs = dbe_->book1D(histo, histo,2000, -10., 10.);
      
    }
    // ************** HB **********************************
    if (subdet_ == 1 || subdet_ == 5 ){
      
      if(etype_ == 1 && subdet_ == 1 ) { 
	sprintf (histo, "HcalRecHitTask_number_of_rechits_in_cone_HB" ) ;
	meNumRecHitsConeHB    = dbe_->book1D(histo, histo, 100, 0., 100.);
	
	sprintf (histo, "HcalRecHitTask_number_of_rechits_above_1GeV_HB");
	meNumRecHitsThreshHB = dbe_->book1D(histo, histo,  30, 0., 30.); 
	
	sprintf (histo, "HcalRecHitTask_sum_of_rechits_energy_HB" ) ;
	meSumRecHitsEnergyHB = dbe_->book1D(histo,histo, 60 , -20., 280.);
	
	sprintf (histo, "HcalRecHitTask_sum_of_rechits_energy_in_cone_HB" ) ;
	meSumRecHitsEnergyConeHB = dbe_->book1D(histo,histo, 60 , -20., 280.); 

	if (ecalselector_ == "yes") {  
	  sprintf (histo, "HcalRecHitTask_energy_hcal_vs_ecal_HB");
	  meEnergyHcalVsEcalHB = dbe_->book2D(histo, histo, 300, 0., 150., 300, 0., 150.);  	
	  sprintf (histo, "HcalRecHitTask_number_of_ecalrechits_in_cone_HB");
	  meNumEcalRecHitsConeHB = dbe_->book1D(histo, histo,   300, 0., 300.);
	  
	  sprintf (histo, "HcalRecHitTask_energy_ecal_plus_hcal_HB" ) ;
	  meEcalHcalEnergyHB = dbe_->book1D(histo,histo, 60 , -20., 280.);
	  
	  sprintf (histo, "HcalRecHitTask_energy_ecal_plus_hcal_in_cone_HB" ) ;
	  meEcalHcalEnergyConeHB =  dbe_->book1D(histo,histo, 60 , -20., 280.);
	}
      }
      
      sprintf (histo, "HcalRecHitTask_energy_of_rechits_HB" ) ;
      meRecHitsEnergyHB = dbe_->book1D(histo, histo, 210 , -10. , 200.); 
      
      sprintf (histo, "HcalRecHitTask_timing_HB" ) ;
      meTimeHB = dbe_->book1D(histo, histo, 2000 , -100. , 100.); 
      
      sprintf (histo, "HcalRecHitTask_timing_vs_energy_HB" ) ;
      meTE_HB = dbe_->book2D(histo, histo, 200, -5., 95.,  220, -10., 100.);
      
      sprintf (histo, "HcalRecHitTask_timing_vs_energy_HB_depth1" ) ;
      meTE_HB1 = dbe_->book2D(histo, histo, 200, -5., 95.,  220, -10., 100.);
      
      sprintf (histo, "HcalRecHitTask_timing_vs_energy_HB_depth2" ) ;
      meTE_HB2 = dbe_->book2D(histo, histo, 200, -5., 95.,  220, -10., 100.);
      
      sprintf (histo, "HcalRecHitTask_timing_vs_energy_profile_HB" ) ;
      meTEprofileHB = dbe_->bookProfile(histo, histo, 100, -5., 95., 110, -10., 100.); 
      
      sprintf (histo, "HcalRecHitTask_energy_rechits_vs_simhits_HB");
      meRecHitSimHitHB = dbe_->book2D(histo, histo, 120, 0., 1.2,  300, 0., 150.);
      
      sprintf (histo, "HcalRecHitTask_energy_rechits_vs_simhits_profile_HB");
      meRecHitSimHitProfileHB = dbe_->bookProfile(histo, histo, 120, 0., 1.2, 500, 0., 500.);  
      
    }
    
    // ********************** HE ************************************
    if ( subdet_ == 2 || subdet_ == 5 ){
      
      if(etype_ == 1 && subdet_ == 2 ) { 
	
	sprintf (histo, "HcalRecHitTask_number_of_rechits_in_cone_HE" ) ;
	meNumRecHitsConeHE    = dbe_->book1D(histo, histo, 100, 0., 100.);
	
	sprintf (histo, "HcalRecHitTask_number_of_rechits_above_1GeV_HE");
	meNumRecHitsThreshHE = dbe_->book1D(histo, histo,  30, 0., 30.);  
	
	sprintf (histo, "HcalRecHitTask_sum_of_rechits_energy_HE" ) ;
	meSumRecHitsEnergyHE = dbe_->book1D(histo,histo, 60 , -20., 280.);
	
	sprintf (histo, "HcalRecHitTask_sum_of_rechits_energy_in_cone_HE" ) ;
	meSumRecHitsEnergyConeHE = dbe_->book1D(histo,histo, 60 , -20., 280.);
	
	if (ecalselector_ == "yes") {  	
	  sprintf (histo, "HcalRecHitTask_energy_ecal_plus_hcal_HE" ) ;
	  meEcalHcalEnergyHE = dbe_->book1D(histo,histo, 80, -20., 380.);
	  
	  sprintf (histo, "HcalRecHitTask_energy_hcal_vs_ecal_HE");
	  meEnergyHcalVsEcalHE = dbe_->book2D(histo, histo, 300, 0., 150., 300, 0., 150.);
	  sprintf (histo, "HcalRecHitTask_number_of_ecalrechits_in_cone_HE");
	  meNumEcalRecHitsConeHE = dbe_->book1D(histo, histo, 300, 0., 300.);   
	  sprintf (histo, "HcalRecHitTask_energy_ecal_plus_hcal_in_cone_HE" ) ;
	  meEcalHcalEnergyConeHE =  dbe_->book1D(histo,histo, 60, -20., 280.);
	}	      
      }
      
      sprintf (histo, "HcalRecHitTask_energy_of_rechits_HE" ) ;
      meRecHitsEnergyHE = dbe_->book1D(histo, histo, 210, -10., 200.); 
      
      sprintf (histo, "HcalRecHitTask_timing_HE" ) ;
      meTimeHE = dbe_->book1D(histo, histo, 200 , -100. , 100.); 
      
      sprintf (histo, "HcalRecHitTask_timing_vs_energy_HE" ) ;
      meTE_HE = dbe_->book2D(histo, histo, 200, -5., 95.,  220, -10., 100.);
      
      sprintf (histo, "HcalRecHitTask_timing_vs_energy_HE_depth1" ) ;
      meTE_HE1 = dbe_->book2D(histo, histo, 200, -5., 95.,  220, -10., 100.);
      
      sprintf (histo, "HcalRecHitTask_timing_vs_energy_HE_depth2" ) ;
      meTE_HE2 = dbe_->book2D(histo, histo, 200, -5., 95.,  220, -10., 100.);
      
      sprintf (histo, "HcalRecHitTask_timing_vs_energy_profile_HE" ) ;
      meTEprofileHE = dbe_->bookProfile(histo, histo, 100, -5., 95., 110, -10., 100.); 
      
      sprintf (histo, "HcalRecHitTask_energy_rechits_vs_simhits_HE");
      meRecHitSimHitHE = dbe_->book2D(histo, histo, 100, 0., 1.0,  300, 0., 150.);
      
      sprintf (histo, "HcalRecHitTask_energy_rechits_vs_simhits_profile_HE");
      meRecHitSimHitProfileHE = dbe_->bookProfile(histo, histo, 100, 0., 1.0, 500, 0., 500.);  
   
      
    }

    // ************** HO ****************************************
    if ( subdet_ == 3 || subdet_ == 5  ){
      
      if(etype_ == 1 && subdet_ == 3) { 
      
	sprintf (histo, "HcalRecHitTask_number_of_rechits_in_cone_HO" ) ;
	meNumRecHitsConeHO    = dbe_->book1D(histo, histo, 100, 0 , 100.);
	
	sprintf (histo, "HcalRecHitTask_number_of_rechits_above_1GeV_HO");
	meNumRecHitsThreshHO = dbe_->book1D(histo, histo,   30, 0., 30.);   
	
	sprintf (histo, "HcalRecHitTask_sum_of_rechits_energy_HO" ) ;
	meSumRecHitsEnergyHO = dbe_->book1D(histo,histo, 60 , -20., 280.);
	
	sprintf (histo, "HcalRecHitTask_sum_of_rechits_energy_in_cone_HO" ) ;
	meSumRecHitsEnergyConeHO = dbe_->book1D(histo,histo, 60 , -20., 280.);   
      }      

      sprintf (histo, "HcalRecHitTask_energy_of_rechits_HO" ) ;
      meRecHitsEnergyHO = dbe_->book1D(histo, histo, 210 , -10. , 200.); 
      
      sprintf (histo, "HcalRecHitTask_timing_HO" ) ;
      meTimeHO = dbe_->book1D(histo, histo, 2000 , -100., 100.); 
      
      sprintf (histo, "HcalRecHitTask_timing_vs_energy_HO" ) ;
      meTE_HO= dbe_->book2D(histo, histo, 200, -5., 95., 220, -10., 100.);
      
      sprintf (histo, "HcalRecHitTask_timing_vs_energy_profile_HO" ) ;
      meTEprofileHO = dbe_->bookProfile(histo, histo, 100, -5., 95.,  110, -10., 100.); 
      
      
      sprintf (histo, "HcalRecHitTask_energy_rechits_vs_simhits_HO");
      meRecHitSimHitHO = dbe_->book2D(histo, histo, 120, 0., 1.2,  300, 0., 300.);
      
      sprintf (histo, "HcalRecHitTask_energy_rechits_vs_simhits_profile_HO");
      meRecHitSimHitProfileHO = dbe_->bookProfile(histo, histo, 120, 0., 1.2, 500, 0., 500.);  
        
    }   
  
    // ********************** HF ************************************
    if ( subdet_ == 4 || subdet_ == 5 ){

      if(etype_ == 1 &&  subdet_ == 4) { 

	sprintf (histo, "HcalRecHitTask_number_of_rechits_in_cone_HF" ) ;
	meNumRecHitsConeHF    = dbe_->book1D(histo, histo, 100, 0 , 100.);
	
	sprintf (histo, "HcalRecHitTask_sum_of_rechits_energy_HF" ) ;
	meSumRecHitsEnergyHF = dbe_->book1D(histo,histo, 60 , -20., 280.);
	
	sprintf (histo, "HcalRecHitTask_sum_of_rechits_energy_in_cone_HF" ) ;
	meSumRecHitsEnergyConeHF = dbe_->book1D(histo,histo, 60 , -20., 280.);   
	sprintf (histo, "HcalRecHitTask_sum_of_rechits_energy_in_cone_HFL" ) ;
	meSumRecHitsEnergyConeHFL = dbe_->book1D(histo,histo, 60 , -20., 280.);   
	sprintf (histo, "HcalRecHitTask_sum_of_rechits_energy_in_cone_HFS" ) ;
	meSumRecHitsEnergyConeHFS = dbe_->book1D(histo,histo, 60 , -20., 280.);   

      }

      sprintf (histo, "HcalRecHitTask_energy_of_rechits_HF" ) ;
      meRecHitsEnergyHF = dbe_->book1D(histo, histo, 210 , -10. , 200.); 

      sprintf (histo, "HcalRecHitTask_timing_HF" ) ;
      meTimeHF = dbe_->book1D(histo, histo, 2000 , -100. , 100.); 
      
      sprintf (histo, "HcalRecHitTask_timing_vs_energy_HF" ) ;
      meTE_HF = dbe_->book2D(histo, histo, 200, -5., 95., 220, -10., 100.);
      
      sprintf (histo, "HcalRecHitTask_timing_vs_energy_HFL" ) ;
      meTE_HFL = dbe_->book2D(histo, histo, 200, -5., 95., 220, -10., 100.);
      
      sprintf (histo, "HcalRecHitTask_timing_vs_energy_HFS" ) ;
      meTE_HFS = dbe_->book2D(histo, histo, 200, -5., 95., 220, -10., 100.);
      
      sprintf (histo, "HcalRecHitTask_timing_vs_energy_profile_HF" ) ;
      meTEprofileHF = dbe_->bookProfile(histo, histo, 100, -5., 95., 110, -10., 100.); 
            
      sprintf (histo, "HcalRecHitTask_energy_rechits_vs_simhits_HF");
      meRecHitSimHitHF  = dbe_->book2D(histo, histo, 120, 0., 60., 200, 0., 200.);      
      sprintf (histo, "HcalRecHitTask_energy_rechits_vs_simhits_HFL");
      meRecHitSimHitHFL = dbe_->book2D(histo, histo, 120, 0., 60., 200, 0., 200.);      
      sprintf (histo, "HcalRecHitTask_energy_rechits_vs_simhits_HFS");
      meRecHitSimHitHFS = dbe_->book2D(histo, histo, 120, 0., 60., 200, 0., 200.);

      sprintf (histo, "HcalRecHitTask_energy_rechits_vs_simhits_profile_HF");
      meRecHitSimHitProfileHF  = dbe_->bookProfile(histo, histo, 60, 0., 60., 500, 0., 500.);  
      sprintf (histo, "HcalRecHitTask_energy_rechits_vs_simhits_profile_HFL");
      meRecHitSimHitProfileHFL = dbe_->bookProfile(histo, histo, 60, 0., 60., 500, 0., 500.);  
      sprintf (histo, "HcalRecHitTask_energy_rechits_vs_simhits_profile_HFS");
      meRecHitSimHitProfileHFS = dbe_->bookProfile(histo, histo, 60, 0., 60., 500, 0., 500.);  
   
    }
  }  //end-of if(_dbe) 

}


HcalRecHitsValidation::~HcalRecHitsValidation() {

  if (hcalselector_ == "ZS"  ) {

    for (unsigned int i1 = 0;  i1 < 82; i1++) {
      for (unsigned int i2 = 0;  i2 < 72; i2++) {
	
	int index = (i1-41) * 72 + i2;
	
	double e = emap_min[i1][i2][0][0];
	if( e < 10000.) {
	  ZS_HB1->Fill(e);
	  ZS_seqHB1->Fill(double(index),e);
	}
	e = emap_min[i1][i2][1][0];
	if( e < 10000.) {
	  ZS_HB2->Fill(e);
	  ZS_seqHB2->Fill(double(index),e);
	}
	
	e = emap_min[i1][i2][0][1];
	if( e < 10000.) {
	  ZS_HE1->Fill(e);
	  ZS_seqHE1->Fill(double(index),e);
	}
	e = emap_min[i1][i2][1][1];
	if( e < 10000.) {
	  ZS_HE2->Fill(e);
	  ZS_seqHE2->Fill(double(index),e);
	}
	e = emap_min[i1][i2][2][1];
	if( e < 10000.) {
	  ZS_HE3->Fill(e);
	  ZS_seqHE3->Fill(double(index),e);
	}
	
	
	e = emap_min[i1][i2][3][2];
	if( e < 10000.) {
	  ZS_HO->Fill(e);
	  ZS_seqHO->Fill(double(index),e);
	}
	
	e = emap_min[i1][i2][0][3];
	if( e < 10000.) {
	  ZS_HF1->Fill(e);
	  ZS_seqHF1->Fill(double(index),e);
	}
	
	e = emap_min[i1][i2][1][3];
	if( e < 10000.) {
	  ZS_HF2->Fill(e);
	  ZS_seqHF2->Fill(double(index),e);
	}
	
	for (unsigned int i3 = 0;  i3 < 4;  i3++) {  // depth
	  double emin = 100000.;
	  for (unsigned int i4 = 0;  i4 < 4;  i4++) {  // subdet
	    if ( emin > emap_min [i1][i2][i3][i4]) 
	      emin = emap_min [i1][i2][i3][i4];
	  }
	  if( i3 == 0 && emin < 10000.)
	    map_depth1->Fill(double(i1-41),double(i2),emin);
	  if( i3 == 1 && emin < 10000.)
	    map_depth2->Fill(double(i1-41),double(i2),emin);
	  if( i3 == 2 && emin < 10000.) 
	    map_depth3->Fill(double(i1-41),double(i2),emin);
	  if( i3 == 3 && emin < 10000.) 
	    map_depth4->Fill(double(i1-41),double(i2),emin);
	}
      }
    } 
  }
  
  // mean energies and occupancies evaluation
  else {

    int nx = emap_HBdepth1->getNbinsX();    
    int ny = emap_HBdepth1->getNbinsY();
    float cnorm;
    float fev = float (nevtot);
    std::cout << "*** nevtot " <<  nevtot << std::endl; 

    float sumphi_hb1, sumphi_hb2, sumphi_he1, sumphi_he2, sumphi_he3,
      sumphi_ho, sumphi_hf1, sumphi_hf2;

    if(nx != 82)  std::cout << "*** problem with binning " << std::endl;
    float phi_factor;  

    for (int i = 1; i <= nx; i++) {
      sumphi_hb1 = 0.;
      sumphi_hb2 = 0.;
      sumphi_he1 = 0.;
      sumphi_he2 = 0.;
      sumphi_he3 = 0.;
      sumphi_ho  = 0.; 
      sumphi_hf1 = 0.;
      sumphi_hf2 = 0.;
      for (int j = 1; j <= ny; j++) {
	cnorm = emap_HBdepth1->getBinContent(i,j) / fev;
        emap_HBdepth1->setBinContent(i,j,cnorm);
	cnorm = emap_HBdepth2->getBinContent(i,j) / fev; 
        emap_HBdepth2->setBinContent(i,j,cnorm);
	cnorm = emap_HEdepth1->getBinContent(i,j) / fev;
        emap_HEdepth1->setBinContent(i,j,cnorm);
	cnorm = emap_HEdepth2->getBinContent(i,j) / fev;
        emap_HEdepth2->setBinContent(i,j,cnorm);
	cnorm = emap_HEdepth3->getBinContent(i,j) / fev;
        emap_HEdepth3->setBinContent(i,j,cnorm);
	cnorm = emap_HO->getBinContent(i,j) / fev;
        emap_HO->setBinContent(i,j,cnorm);
	cnorm = emap_HFdepth1->getBinContent(i,j) / fev;
        emap_HFdepth1->setBinContent(i,j,cnorm);
	cnorm = emap_HFdepth2->getBinContent(i,j) / fev;
        emap_HFdepth2->setBinContent(i,j,cnorm);

	cnorm = occupancy_map_HBdepth1->getBinContent(i,j) / fev;   
	occupancy_map_HBdepth1->setBinContent(i,j,cnorm);

	/*
	std::cout << "*** occupancy depth1  i j = " << i << " " << j 
		  << " bin cont. " 
		  << occupancy_map_depth1->getBinContent(i,j) 
		  << "  cnorm " <<  cnorm << std::endl;
	*/
	cnorm = occupancy_map_HBdepth2->getBinContent(i,j) / fev;   
	occupancy_map_HBdepth2->setBinContent(i,j,cnorm);
	cnorm = occupancy_map_HEdepth1->getBinContent(i,j) / fev;   
	occupancy_map_HEdepth1->setBinContent(i,j,cnorm);
	cnorm = occupancy_map_HEdepth2->getBinContent(i,j) / fev;   
	occupancy_map_HEdepth2->setBinContent(i,j,cnorm);
	cnorm = occupancy_map_HEdepth3->getBinContent(i,j) / fev;   
	occupancy_map_HEdepth3->setBinContent(i,j,cnorm);
	cnorm = occupancy_map_HO->getBinContent(i,j) / fev;   
	occupancy_map_HO->setBinContent(i,j,cnorm);
	cnorm = occupancy_map_HFdepth1->getBinContent(i,j) / fev;   
	occupancy_map_HFdepth1->setBinContent(i,j,cnorm);
	cnorm = occupancy_map_HFdepth2->getBinContent(i,j) / fev;   
	occupancy_map_HFdepth2->setBinContent(i,j,cnorm);
     
        sumphi_hb1 += occupancy_map_HBdepth1->getBinContent(i,j);
        sumphi_hb2 += occupancy_map_HBdepth2->getBinContent(i,j);
        sumphi_he1 += occupancy_map_HEdepth1->getBinContent(i,j);
        sumphi_he2 += occupancy_map_HEdepth2->getBinContent(i,j);
        sumphi_he3 += occupancy_map_HEdepth3->getBinContent(i,j);
        sumphi_ho  += occupancy_map_HO->getBinContent(i,j);
        sumphi_hf1 += occupancy_map_HFdepth1->getBinContent(i,j);
        sumphi_hf2 += occupancy_map_HFdepth2->getBinContent(i,j);
      }

      int ieta = i - 42;        // -41 -1, 0 40 
      if(ieta >=0 ) ieta +=1;   // -41 -1, 1 41 

      if(ieta >= -20 && ieta <= 20 )
	{phi_factor = 72.;}
      else {
	if(ieta >= 40 || ieta <= -40 ) {phi_factor = 18.;}
        else 
	  phi_factor = 36.;
      }  
      if(ieta >= 0) ieta -= 1; // -41 -1, 0 40
	      
      /*
      std::cout << "*** ieta = " << ieta << "  sumphi_hb1, sumphi_hb2, sumphi_he1, sumphi_he2, simphi_he3, sumphi_ho, simphi_hf1, sumphi_hf2" << std::endl 
		<< sumphi_hb1 << " " << sumphi_hb2 << " " << sumphi_he1 << " "
		<< sumphi_he2 << " " << simphi_he3 << " " << sumphi_ho  << " " 
		<< simphi_hf1 << " " << sumphi_hf2 << std::endl << std::endl;
      */

      cnorm = sumphi_hb1 / phi_factor;
      occupancy_vs_ieta_HBdepth1->Fill(float(ieta), cnorm);
      cnorm = sumphi_hb2 / phi_factor;
      occupancy_vs_ieta_HBdepth2->Fill(float(ieta), cnorm);
      cnorm = sumphi_he1 / phi_factor;
      occupancy_vs_ieta_HEdepth1->Fill(float(ieta), cnorm);
      cnorm = sumphi_he2 / phi_factor;
      occupancy_vs_ieta_HEdepth2->Fill(float(ieta), cnorm);
      cnorm = sumphi_he3 / phi_factor;
      occupancy_vs_ieta_HEdepth3->Fill(float(ieta), cnorm);
      cnorm = sumphi_ho / phi_factor;
      occupancy_vs_ieta_HO->Fill(float(ieta), cnorm);
      cnorm = sumphi_hf1 / phi_factor;
      occupancy_vs_ieta_HFdepth1->Fill(float(ieta), cnorm);
      cnorm = sumphi_hf2 / phi_factor;
      occupancy_vs_ieta_HFdepth2->Fill(float(ieta), cnorm);

    }

  }


   
  std::cout << " outputFile_.size() =  " << outputFile_.size() << std::endl;
  std::cout << " dbe_ = " << dbe_ << std::endl;
  if ( outputFile_.size() != 0 && dbe_ ) dbe_->save(outputFile_);
  
}

void HcalRecHitsValidation::endJob() { }

void HcalRecHitsValidation::beginJob(const edm::EventSetup& c){ }

void HcalRecHitsValidation::analyze(edm::Event const& ev, edm::EventSetup const& c) {

  using namespace edm;

  // cuts for each subdet_ector mimiking  "Scheme B"
  //  double cutHB = 0.9, cutHE = 1.4, cutHO = 1.1, cutHFL = 1.2, cutHFS = 1.8; 

  // energy in HCAL
  double eHcal        = 0.;
  double eHcalCone    = 0.;  
  double eHcalConeHB  = 0.;  
  double eHcalConeHE  = 0.;  
  double eHcalConeHO  = 0.;  
  double eHcalConeHF  = 0.;  
  double eHcalConeHFL = 0.;  
  double eHcalConeHFS = 0.;  
  // Total numbet of RecHits in HCAL, in the cone, above 1 GeV theshold
  int nrechits       = 0;
  int nrechitsCone   = 0;
  int nrechitsThresh = 0;

  // energy in ECAL
  double eEcal       = 0.;
  double eEcalB      = 0.;
  double eEcalE      = 0.;
  double eEcalCone   = 0.;
  int numrechitsEcal = 0;

  // MC info 
  double phi_MC = -999999.;  // phi of initial particle from HepMC
  double eta_MC = -999999.;  // eta of initial particle from HepMC
  bool MC = false;

  // HCAL energy around MC eta-phi at all depths;
  double partR = 0.3;
  double ehcal_coneMC_1 = 0.;
  double ehcal_coneMC_2 = 0.;
  double ehcal_coneMC_3 = 0.;
  double ehcal_coneMC_4 = 0.;

  // Cone size for serach of the hottest HCAL cell around MC
  double searchR = 1.0; 
  double eps     = 0.001;

  // Single particle samples: actual eta-phi position of cluster around
  // hottest cell
  double etaHot  = 99999.; 
  double phiHot  = 99999.; 
  int    ietahot = 1000;
  int    iphihot = 1000;

  // MC information

  //  std::cout << "*** 1" << std::endl; 


  if(imc != 0) { 

  edm::Handle<edm::HepMCProduct> evtMC;
  //  ev.getByLabel("VtxSmeared",evtMC);
  ev.getByLabel("source",evtMC);
  if (!evtMC.isValid()) {
    std::cout << "no HepMCProduct found" << std::endl;    
  } else {
    MC=true;
    //    std::cout << "*** source HepMCProduct found"<< std::endl;
  }  

  // MC particle with highest pt is taken as a direction reference  
  double maxPt = -99999.;
  int npart    = 0;
  HepMC::GenEvent * myGenEvent = new  HepMC::GenEvent(*(evtMC->GetEvent()));
  for ( HepMC::GenEvent::particle_iterator p = myGenEvent->particles_begin();
	p != myGenEvent->particles_end(); ++p ) {
    double phip = (*p)->momentum().phi();
    double etap = (*p)->momentum().eta();
    //    phi_MC = phip;
    //    eta_MC = etap;
    double pt  = (*p)->momentum().perp();
    if(pt > maxPt) { npart++; maxPt = pt; phi_MC = phip; eta_MC = etap; }
  }
  //  std::cout << "*** Max pT = " << maxPt <<  std::endl;  

  }

  //   std::cout << "*** 2" << std::endl; 

  c.get<IdealGeometryRecord>().get (geometry);


  // Fill working vectors of HCAL RecHits quantities 
  fillRecHitsTmp(subdet_, ev); 

  //   std::cout << "*** 3" << std::endl; 


  //===========================================================================
  // IN ALL other CASES : ieta-iphi maps 
  //===========================================================================

  // ECAL 
  if(ecalselector_ == "yes" && (subdet_ == 1 || subdet_ == 2 || subdet_ == 5)) {
    Handle<EBRecHitCollection> rhitEB;

    if(famos_)
      ev.getByLabel("caloRecHits","EcalRecHitsEB", rhitEB);
    else
      ev.getByLabel("ecalRecHit","EcalRecHitsEB", rhitEB);

    EcalRecHitCollection::const_iterator RecHit = rhitEB.product()->begin();  
    EcalRecHitCollection::const_iterator RecHitEnd = rhitEB.product()->end();  
    
    for (; RecHit != RecHitEnd ; ++RecHit) {
      EBDetId EBid = EBDetId(RecHit->id());
       
      const CaloCellGeometry* cellGeometry =
	geometry->getSubdetectorGeometry (EBid)->getGeometry (EBid) ;
      double eta = cellGeometry->getPosition ().eta () ;
      double phi = cellGeometry->getPosition ().phi () ;
      double en  = RecHit->energy();
      eEcal  += en;
      eEcalB += en;

      map_ecal->Fill(eta, phi, en);

      double r   = dR(eta_MC, phi_MC, eta, phi);
      if( r < partR)  {
	eEcalCone += en;
	numrechitsEcal++; 
      }
    }

    
    Handle<EERecHitCollection> rhitEE;
 
    if(famos_)
      ev.getByLabel("caloRecHits", "EcalRecHitsEE", rhitEE );
    else
      ev.getByLabel("ecalRecHit","EcalRecHitsEE", rhitEE);

    RecHit = rhitEE.product()->begin();  
    RecHitEnd = rhitEE.product()->end();  
    
    for (; RecHit != RecHitEnd ; ++RecHit) {
      EEDetId EEid = EEDetId(RecHit->id());
      
      const CaloCellGeometry* cellGeometry =
	geometry->getSubdetectorGeometry (EEid)->getGeometry (EEid) ;
      double eta = cellGeometry->getPosition ().eta () ;
      double phi = cellGeometry->getPosition ().phi () ;	
      double en   = RecHit->energy();
      eEcal  += en;
      eEcalE += en;

      map_ecal->Fill(eta, phi, en);

      double r   = dR(eta_MC, phi_MC, eta, phi);
      if( r < partR)  {
	eEcalCone += en;
	numrechitsEcal++; 
      }
    }
  }     // end of ECAL selection 


  //    std::cout << "*** 4" << std::endl; 


  // Counting, including ZS items
  // Filling HCAL maps  ----------------------------------------------------
  double maxE = -99999.;
  
  int nhb1 = 0;
  int nhb2 = 0;
  int nhe1 = 0;
  int nhe2 = 0;
  int nhe3 = 0;
  int nho  = 0;
  int nhf1 = 0;
  int nhf2 = 0;  
  
  for (unsigned int i = 0; i < cen.size(); i++) {
    
    int sub    = csub[i];
    int depth  = cdepth[i];
    int ieta   = cieta[i]; 
    int iphi   = ciphi[i]; 
    double en  = cen[i]; 
    double eta = ceta[i]; 
    double phi = cphi[i]; 
    double z   = cz[i];
    
    /*   
	 std::cout << "*** point 4-1" << " ieta, iphi, depth, sub = "
	 << ieta << ", " << iphi << ", " << depth << ", " << sub  
	 << std::endl;
    */
    
    
    if( sub == 1 && depth == 1)  nhb1++;
    if( sub == 1 && depth == 2)  nhb2++;
    if( sub == 2 && depth == 1)  nhe1++;
    if( sub == 2 && depth == 2)  nhe2++;
    if( sub == 2 && depth == 3)  nhe3++;
    if( sub == 3 && depth == 4)  nho++;
    if( sub == 4 && depth == 1)  nhf1++;
    if( sub == 4 && depth == 2)  nhf2++;
    
    if( subdet_ == 6) {                                    // ZS specific
      if( en < emap_min[ieta+41][iphi][depth-1][sub-1] )
	emap_min[ieta+41][iphi][depth-1][sub-1] = en;
    }
    
    double emin = 1.;
    if(fabs(eta) > 3.) emin = 5.; 
    
    double r  = dR(eta_MC, phi_MC, eta, phi);
    if( r < searchR ) { // search for hottest cell in a big cone
      if(maxE < en && en > emin) {
	maxE    = en;
	etaHot  = eta;
	phiHot  = phi;
	ietahot = ieta;
	iphihot = iphi;
      }
    }

    /*   
    if(ieta == 27 ) { 
      std::cout << "*** ieta=28, iphi = " << iphi << "  det = " 
		<< sub << "  depth = " << depth << std::endl;
    }
    */

    if (depth == 1 ) profile_z1->Fill(double(ieta), z); 
    if (depth == 2 ) profile_z2->Fill(double(ieta), z); 
    if (depth == 3 ) profile_z3->Fill(double(ieta), z); 
    if (depth == 4 ) profile_z4->Fill(double(ieta), z); 


    if (depth == 1 && sub == 1 ) {
      emap_HBdepth1->Fill(double(ieta), double(iphi), en);
      emean_vs_ieta_HBdepth1->Fill(double(ieta), en);
      occupancy_map_HBdepth1->Fill(double(ieta), double(iphi));          
    }
    if (depth == 2  && sub == 1) {
      emap_HBdepth2->Fill(double(ieta), double(iphi), en);
      emean_vs_ieta_HBdepth2->Fill(double(ieta), en);
      occupancy_map_HBdepth2->Fill(double(ieta), double(iphi));          
    }
    if (depth == 1 && sub == 2) {
      emap_HEdepth1->Fill(double(ieta), double(iphi), en);
      emean_vs_ieta_HEdepth1->Fill(double(ieta), en);
      occupancy_map_HEdepth1->Fill(double(ieta), double(iphi));          
    }
    if (depth == 2 && sub == 2) {
      emap_HEdepth2->Fill(double(ieta), double(iphi), en);
      emean_vs_ieta_HEdepth2->Fill(double(ieta), en);
      occupancy_map_HEdepth2->Fill(double(ieta), double(iphi));          
    }
    if (depth == 3 && sub == 2) {
      emap_HEdepth3->Fill(double(ieta), double(iphi), en);
      emean_vs_ieta_HEdepth3->Fill(double(ieta), en);
      occupancy_map_HEdepth3->Fill(double(ieta), double(iphi));          
    }
    if (depth == 4 ) {
      emap_HO->Fill(double(ieta), double(iphi), en);
      emean_vs_ieta_HO->Fill(double(ieta), en);
      occupancy_map_HO->Fill(double(ieta), double(iphi));          
    }
    if (depth == 1 && sub == 4) {
      emap_HFdepth1->Fill(double(ieta), double(iphi), en);
      emean_vs_ieta_HFdepth1->Fill(double(ieta), en);
      occupancy_map_HFdepth1->Fill(double(ieta), double(iphi));          
    }
    if (depth == 2 && sub == 4) {
      emap_HFdepth2->Fill(double(ieta), double(iphi), en);
      emean_vs_ieta_HFdepth2->Fill(double(ieta), en);
      occupancy_map_HFdepth2->Fill(double(ieta), double(iphi));          
    }
    
    // maps with "Scheme B" cuts 
    /*
    double cut = cutHB;  
    if (sub == 2) cut = cutHE;
    if (sub == 3) cut = cutHO;
    if (sub == 4){
      if(depth == 1) cut = cutHFL;
      else  cut = cutHFS;
    }
    
    if(en > cut ) {
      if (depth == 1) {
	emap_depth1_cuts->Fill(double(ieta), double(iphi), en);
	emean_vs_ieta_depth1_cuts->Fill(double(ieta), en);
      }
      if (depth == 2) {
	emap_depth2_cuts->Fill(double(ieta), double(iphi), en);
	emean_vs_ieta_depth2_cuts->Fill(double(ieta), en);
      }
      if (depth == 3) {
	emap_depth3_cuts->Fill(double(ieta), double(iphi), en);
	emean_vs_ieta_depth3_cuts->Fill(double(ieta), en);
      }
      if (depth == 4) {
	emap_depth4_cuts->Fill(double(ieta), double(iphi), en);
	emean_vs_ieta_depth4_cuts->Fill(double(ieta), en);
      }
    }
    */    

    if( r < partR ) {
      if (depth == 1) ehcal_coneMC_1 += en; 
      if (depth == 2) ehcal_coneMC_2 += en; 
      if (depth == 3) ehcal_coneMC_3 += en; 
      if (depth == 4) ehcal_coneMC_4 += en; 
    }
    
  } 

  Nhb->Fill(double(nhb1 + nhb2));
  Nhe->Fill(double(nhe1 + nhe2 + nhe3));
  Nho->Fill(double(nho));
  Nhf->Fill(double(nhf1 + nhf2));
  
  if( subdet_ == 6) {               // ZS plots
    ZS_nHB1->Fill(double(nhb1));  
    ZS_nHB2->Fill(double(nhb2));  
    ZS_nHE1->Fill(double(nhe1));  
    ZS_nHE2->Fill(double(nhe2));  
    ZS_nHE3->Fill(double(nhe3));  
    ZS_nHO ->Fill(double(nho));  
    ZS_nHF1->Fill(double(nhf1));  
    ZS_nHF2->Fill(double(nhf2));  
  }

  //    std::cout << "*** 5" << std::endl; 
  
  
  // IN ALL THE CASES : sum of rechits around MC particle = f(eta,phi) 
  
  if(imc != 0) {
    map_econe_depth1->Fill(eta_MC, phi_MC, ehcal_coneMC_1);
    map_econe_depth2->Fill(eta_MC, phi_MC, ehcal_coneMC_2);
    map_econe_depth3->Fill(eta_MC, phi_MC, ehcal_coneMC_3);
    map_econe_depth4->Fill(eta_MC, phi_MC, ehcal_coneMC_4);
  }

  //  NOISE ================================================================= 
  
  if (hcalselector_ == "noise") {
    for (unsigned int i = 0; i < cen.size(); i++) {
      
      int sub   = csub[i];
      int depth = cdepth[i];
      double en = cen[i]; 
      
      if (sub == 1) e_hb->Fill(en);
      if (sub == 2) e_he->Fill(en);  
      if (sub == 3) e_ho->Fill(en);  
      if (sub == 4) {
	if(depth == 1)  
	  e_hfl->Fill(en);  
	else 
	  e_hfs->Fill(en);  
      }
    }
  }

  //===========================================================================
  // SUBSYSTEMS,  
  //===========================================================================
  
  else if ((subdet_ != 6) && (subdet_ != 0)) {

    //       std::cout << "*** 6" << std::endl; 
    
    
    double clusEta = 999.;
    double clusPhi = 999.; 
    double clusEn  = 0.;
    
    double HcalCone_d1 = 0.;
    double HcalCone_d2 = 0.;
    double HcalCone_d3 = 0.;
    double HcalCone_d4 = 0.;

    int ietaMax1  =  9999;
    int ietaMax2  =  9999;
    int ietaMax3  =  9999;
    int ietaMax4  =  9999;
    double enMax1 = -9999.;
    double enMax2 = -9999.;
    double enMax3 = -9999.;
    double enMax4 = -9999.;

    /*
    std::cout << "*** point 5-1" << "  eta_MC, phi_MC    etaHot,  phiHot = "
	      << eta_MC  << ", " << phi_MC << "   "
	      << etaHot  << ", " << phiHot  
	      << std::endl;
    */

    //   CYCLE over cells ====================================================

    for (unsigned int i = 0; i < cen.size(); i++) {
      int sub    = csub[i];
      int depth  = cdepth[i];
      double eta = ceta[i]; 
      double phi = cphi[i]; 
      double en  = cen[i]; 
      double t   = ctime[i];
      int   ieta = cieta[i];

      double rhot = dR(etaHot, phiHot, eta, phi); 
      if(rhot < partR && en > 1.) { 
	clusEta = (clusEta * clusEn + eta * en)/(clusEn + en);
    	clusPhi = phi12(clusPhi, clusEn, phi, en); 
        clusEn += en;
      }

      nrechits++;	    
      eHcal += en;
      if(en > 1. ) nrechitsThresh++;

      double r    = dR(eta_MC, phi_MC, eta, phi);
      if( r < partR ){
        if(sub == 1)   eHcalConeHB += en;
        if(sub == 2)   eHcalConeHE += en;
        if(sub == 3)   eHcalConeHO += en;
        if(sub == 4) {
	  eHcalConeHF += en;
	  if (depth == 1) eHcalConeHFL += en;
	  else            eHcalConeHFS += en;
	}
	eHcalCone += en;
	nrechitsCone++;

	// search for most energetic cell at the given depth in the cone
        if(depth == 1) {
	  HcalCone_d1 += en;
	  if(enMax1 < en) {
	    enMax1   = en;
	    ietaMax1 = ieta;
	  }
	}
        if(depth == 2) {
	  HcalCone_d2 += en;
	  if(enMax2 < en) {
	    enMax2   = en;
	    ietaMax2 = ieta;
	  }
	}
        if(depth == 3) {
	  HcalCone_d3 += en;
	  if(enMax3 < en) {
	    enMax3   = en;
	    ietaMax3 = ieta;
	  }
	}
        if(depth == 4) {
	  HcalCone_d4 += en;
	  if(enMax4 < en) {
	    enMax4   = en;
	    ietaMax4 = ieta;
	  }
	}


      }
      
      if(sub == 1 && (subdet_ == 1 || subdet_ == 5)) {  
	meTimeHB->Fill(t);
	meRecHitsEnergyHB->Fill(en);
	meTE_HB->Fill( en, t);
	if(depth == 1)      meTE_HB1->Fill( en, t);
	else if(depth == 2) meTE_HB2->Fill( en, t);
	meTEprofileHB->Fill(en, t);
      }
      
      if(sub == 2 && (subdet_ == 2 || subdet_ == 5)) {  
	meTimeHE->Fill(t);
	meRecHitsEnergyHE->Fill(en);
	meTE_HE->Fill( en, t);
	if(depth == 1)      meTE_HE1->Fill( en, t);
	else if(depth == 2) meTE_HE2->Fill( en, t);
	meTEprofileHE->Fill(en, t);
      }

      if(sub == 4 && (subdet_ == 4 || subdet_ == 5)) {  
	  meTimeHF->Fill(t);
	  meTE_HF->Fill(en, t);
	  if(depth == 1) meTE_HFL->Fill( en, t);
	  else           meTE_HFS->Fill( en, t);
	  meTEprofileHF->Fill(en, t);
	  meRecHitsEnergyHF->Fill(en);	  
      }

      if(sub == 3 && (subdet_ == 3 || subdet_ == 5)) {  
	meTimeHO->Fill(t);
	meRecHitsEnergyHO->Fill(en);
	meTE_HO->Fill( en, t);
	meTEprofileHO->Fill(en, t);
      }
    }

    meEnConeEtaProfile_depth1->Fill(double(ietaMax1), HcalCone_d1);
    meEnConeEtaProfile_depth2->Fill(double(ietaMax2), HcalCone_d2);
    meEnConeEtaProfile_depth3->Fill(double(ietaMax3), HcalCone_d3);
    meEnConeEtaProfile_depth4->Fill(double(ietaMax4), HcalCone_d4);


    //     std::cout << "*** 7" << std::endl; 

    
    // Single particle samples ONLY !  ======================================
    // Fill up some histos for "integrated" subsustems. 
    
    if(etype_ == 1) {

      /*
      std::cout << "*** point 7-1" << "  eta_MC, phi_MC   clusEta, clusPhi = "
                << eta_MC  << ", " << phi_MC << "   "
		<< clusEta << ", " << clusPhi 
		<< std::endl;
      */    

      double phidev = dPhiWsign(clusPhi, phi_MC);
      meDeltaPhi->Fill(eta_MC, phidev);
      double etadev = clusEta - eta_MC;
      meDeltaEta->Fill(eta_MC, etadev);

      if(subdet_ == 1) {
	meSumRecHitsEnergyHB->Fill(eHcal);
	meSumRecHitsEnergyConeHB->Fill(eHcalConeHB);    
	meNumRecHitsConeHB->Fill(double(nrechitsCone));
	meNumRecHitsThreshHB->Fill(double(nrechitsThresh));
      }

      if(subdet_ == 2) {
        meSumRecHitsEnergyHE->Fill(eHcal);
	meSumRecHitsEnergyConeHE->Fill(eHcalConeHE);    
	meNumRecHitsConeHE->Fill(double(nrechitsCone));
	meNumRecHitsThreshHE->Fill(double(nrechitsThresh));
      }

      if(subdet_ == 3) {
	meSumRecHitsEnergyHO->Fill(eHcal);
	meSumRecHitsEnergyConeHO->Fill(eHcalConeHO);    
	meNumRecHitsConeHO->Fill(double(nrechitsCone));
	meNumRecHitsThreshHO->Fill(double(nrechitsThresh));
      }

      if(subdet_ == 4) {
        if(eHcalConeHF > eps ) {
	meSumRecHitsEnergyHF ->Fill(eHcal);
	meSumRecHitsEnergyConeHF ->Fill(eHcalConeHF);    
	meSumRecHitsEnergyConeHFL ->Fill(eHcalConeHFL);    
	meSumRecHitsEnergyConeHFS ->Fill(eHcalConeHFS);    
	meNumRecHitsConeHF->Fill(double(nrechitsCone));
	}
      }

      //         std::cout << "*** 8" << std::endl; 


      // Also combine with ECAL if needed 
      if(subdet_ == 1  && ecalselector_ == "yes") {
		
	/*
	std::cout << "*** point 8-1" 
		  << "  eEcalB " << eEcalB << "  eHcal " << eHcal
		  << "  eEcalCone " <<  eEcalCone << "  eHcalCone " 
		  << eHcalCone
		  << "  numrechitsEcal " <<  numrechitsEcal
		  << std::endl;
	
	*/

       	meEcalHcalEnergyHB->Fill(eEcalB+eHcal);
      	meEcalHcalEnergyConeHB->Fill(eEcalCone+eHcalCone);
      	meNumEcalRecHitsConeHB->Fill(double(numrechitsEcal));

      }

      if(subdet_ == 2  && ecalselector_ == "yes"){
	
	/*
	std::cout << "*** point 8-2a" 
		  << "  eEcalE " << eEcalE << "  eHcal " << eHcal
		  << "  eEcalCone " <<  eEcalCone << "  eHcalCone " 
		  << eHcalCone
		  << "  numrechitsEcal " <<  numrechitsEcal
		  << std::endl;
	*/

	meEcalHcalEnergyHE->Fill(eEcalE+eHcal);
	meEcalHcalEnergyConeHE->Fill(eEcalCone+eHcalCone);
	meNumEcalRecHitsConeHE->Fill(double(numrechitsEcal));


      } 

      // Banana plots finally
      if(subdet_ == 1 && ecalselector_ == "yes")
	meEnergyHcalVsEcalHB -> Fill(eEcalCone,eHcalCone);
      if(subdet_ == 2 && ecalselector_ == "yes") 
	meEnergyHcalVsEcalHE -> Fill(eEcalCone,eHcalCone);
    }
  }

  //  std::cout << "*** 9" << std::endl; 


  //===========================================================================
  // Getting SimHits
  //===========================================================================

  if(subdet_ > 0 && subdet_ < 6 &&imc !=0 && !famos_  ) {  // not noise 

    double maxES = -9999.;
    double etaHotS = 1000.;
    double phiHotS = 1000.;
    
    edm::Handle<PCaloHitContainer> hcalHits;
    ev.getByLabel("g4SimHits","HcalHits",hcalHits);
    const PCaloHitContainer * SimHitResult = hcalHits.product () ;
    
    double enSimHits    = 0.;
    double enSimHitsHB  = 0.;
    double enSimHitsHE  = 0.;
    double enSimHitsHO  = 0.;
    double enSimHitsHF  = 0.;
    double enSimHitsHFL = 0.;
    double enSimHitsHFS = 0.;
    // sum of SimHits in the cone 
    
    for (std::vector<PCaloHit>::const_iterator SimHits = SimHitResult->begin () ; SimHits != SimHitResult->end(); ++SimHits) {
      HcalDetId cell(SimHits->id());
      int sub =  cell.subdet();
      const CaloCellGeometry* cellGeometry =
	geometry->getSubdetectorGeometry (cell)->getGeometry (cell) ;
      double etaS = cellGeometry->getPosition().eta () ;
      double phiS = cellGeometry->getPosition().phi () ;
      double en   = SimHits->energy();    

      double emin = 0.01;
      if(fabs(etaS) > 3.) emin = 1.;   

      double r  = dR(eta_MC, phi_MC, etaS, phiS);
      if( r < searchR ) { // search for hottest cell in a big cone
	if(maxES < en && en > emin ) {
	  maxES    = en;
	  etaHotS  = etaS;
	  phiHotS  = phiS;
	}
      }
       
      if ( r < partR ){ // just energy in the small cone
	enSimHits += en;
	if(sub == 1) enSimHitsHB += en; 
	if(sub == 2) enSimHitsHE += en; 
	if(sub == 3) enSimHitsHO += en; 
	if(sub == 4) {
	  enSimHitsHF += en;
	  int depth = cell.depth();
	  if(depth == 1) enSimHitsHFL += en;
	  else           enSimHitsHFS += en;
	} 
      }
    }


    // Second look over SimHits: cluster finding

    double clusEta = 999.;
    double clusPhi = 999.; 
    double clusEn  = 0.;

    for (std::vector<PCaloHit>::const_iterator SimHits = SimHitResult->begin () ; SimHits != SimHitResult->end(); ++SimHits) {
     HcalDetId cell(SimHits->id());

      const CaloCellGeometry* cellGeometry =
	geometry->getSubdetectorGeometry (cell)->getGeometry (cell) ;
      double etaS = cellGeometry->getPosition().eta () ;
      double phiS = cellGeometry->getPosition().phi () ;
      double en   =  SimHits->energy();    

      double emin = 0.01;
      if(fabs(etaS) > 3.) emin = 1.; 

      double rhot = dR(etaHotS, phiHotS, etaS, phiS); 
      if(rhot < partR && en > emin) { 
	clusEta = (clusEta * clusEn + etaS * en)/(clusEn + en);
    	clusPhi = phi12(clusPhi, clusEn, phiS, en); 
        clusEn += en;
      }
    }

    // SimHits cluster deviation from MC (eta, phi)
    if(etype_ == 1) {
      double phidev = dPhiWsign(clusPhi, phi_MC);
      meDeltaPhiS->Fill(eta_MC, phidev);
      double etadev = clusEta - eta_MC;
      meDeltaEtaS->Fill(eta_MC, etadev);
    }


    // Now some histos with SimHits
    
    if(subdet_ == 4 || subdet_ == 5) {
      if(eHcalConeHF > eps) {
      meRecHitSimHitHF->Fill( enSimHitsHF, eHcalConeHF );
      meRecHitSimHitProfileHF->Fill( enSimHitsHF, eHcalConeHF);
      
      meRecHitSimHitHFL->Fill( enSimHitsHFL, eHcalConeHFL );
      meRecHitSimHitProfileHFL->Fill( enSimHitsHFL, eHcalConeHFL);
      meRecHitSimHitHFS->Fill( enSimHitsHFS, eHcalConeHFS );
      meRecHitSimHitProfileHFS->Fill( enSimHitsHFS, eHcalConeHFS);       
      }
    }
    if(subdet_ == 1  || subdet_ == 5) { 
      meRecHitSimHitHB->Fill( enSimHitsHB,eHcalConeHB );
      meRecHitSimHitProfileHB->Fill( enSimHitsHB,eHcalConeHB);
    }
    if(subdet_ == 2  || subdet_ == 5) { 
      meRecHitSimHitHE->Fill( enSimHitsHE,eHcalConeHE );
      meRecHitSimHitProfileHE->Fill( enSimHitsHE,eHcalConeHE);
    }
    if(subdet_ == 3  || subdet_ == 5) { 
      meRecHitSimHitHO->Fill( enSimHitsHO,eHcalConeHO );
      meRecHitSimHitProfileHO->Fill( enSimHitsHO,eHcalConeHO);
    }

  }

  nevtot++;
}


///////////////////////////////////////////////////////////////////////////////
void HcalRecHitsValidation::fillRecHitsTmp(int subdet_, edm::Event const& ev){
  
  using namespace edm;
  
  
  // initialize data vectors
  csub.clear();
  cen.clear();
  ceta.clear();
  cphi.clear();
  ctime.clear();
  cieta.clear();
  ciphi.clear();
  cdepth.clear();


  if( subdet_ == 1 || subdet_ == 2  || subdet_ == 5 || subdet_ == 6 || subdet_ == 0) {
    
    //HBHE
    std::vector<edm::Handle<HBHERecHitCollection> > hbhecoll;
    ev.getManyByType(hbhecoll);

    std::vector<edm::Handle<HBHERecHitCollection> >::iterator i;
    
    int count = 0;
    for (i=hbhecoll.begin(); i!=hbhecoll.end(); i++) {
      
      count ++;  
    //      std::cout << "*** HBHE collection No. " <<  count << std::endl;     
      if ( count == 1) {
      for (HBHERecHitCollection::const_iterator j=(*i)->begin(); j!=(*i)->end(); j++) {
	HcalDetId cell(j->id());
	const CaloCellGeometry* cellGeometry =
	  geometry->getSubdetectorGeometry (cell)->getGeometry (cell) ;
	double eta  = cellGeometry->getPosition().eta () ;
	double phi  = cellGeometry->getPosition().phi () ;
        double zc   = cellGeometry->getPosition().z ();
	int sub     = cell.subdet();
	int depth   = cell.depth();
	int inteta  = cell.ieta();
	if(inteta > 0) inteta -= 1;
	int intphi  = cell.iphi()-1;
	double en   = j->energy();
	double t    = j->time();

	if((iz > 0 && eta > 0.) || (iz < 0 && eta <0.) || iz == 0) { 
	
	  csub.push_back(sub);
	  cen.push_back(en);
	  ceta.push_back(eta);
	  cphi.push_back(phi);
	  ctime.push_back(t);
	  cieta.push_back(inteta);
	  ciphi.push_back(intphi);
	  cdepth.push_back(depth);
	  cz.push_back(zc);
	}
      }

      }
    }
  }

  if( subdet_ == 4 || subdet_ == 5 || subdet_ == 6 || subdet_ == 0) {

    //HF
    std::vector<edm::Handle<HFRecHitCollection> > hfcoll;
    ev.getManyByType(hfcoll);
    std::vector<edm::Handle<HFRecHitCollection> >::iterator ihf;

    int count = 0;
    for (ihf=hfcoll.begin(); ihf!=hfcoll.end(); ihf++) {      
      count++;
      if(count == 1) {
      for (HFRecHitCollection::const_iterator j=(*ihf)->begin(); j!=(*ihf)->end(); j++) {
	HcalDetId cell(j->id());
	const CaloCellGeometry* cellGeometry =
	  geometry->getSubdetectorGeometry (cell)->getGeometry (cell) ;
	double eta   = cellGeometry->getPosition().eta () ;
	double phi   = cellGeometry->getPosition().phi () ;
        double zc     = cellGeometry->getPosition().z ();
	int sub      = cell.subdet();
	int depth    = cell.depth();
	int inteta   = cell.ieta();
	if(inteta > 0) inteta -= 1;
	int intphi   = cell.iphi()-1;
	double en    = j->energy();
	double t     = j->time();

	if((iz > 0 && eta > 0.) || (iz < 0 && eta <0.) || iz == 0) { 
	
	  csub.push_back(sub);
	  cen.push_back(en);
	  ceta.push_back(eta);
	  cphi.push_back(phi);
	  ctime.push_back(t);
	  cieta.push_back(inteta);
	  ciphi.push_back(intphi);
	  cdepth.push_back(depth);
	  cz.push_back(zc);
       
	}
      }
      }
    }
  }

  //HO

  if( subdet_ == 3 || subdet_ == 5 || subdet_ == 6 || subdet_ == 0) {
  
    std::vector<edm::Handle<HORecHitCollection> > hocoll;
    ev.getManyByType(hocoll);
    std::vector<edm::Handle<HORecHitCollection> >::iterator iho;
    
    int count = 0;
    for (iho=hocoll.begin(); iho!=hocoll.end(); iho++) {
      count++;
      if (count == 1) {
      for (HORecHitCollection::const_iterator j=(*iho)->begin(); j!=(*iho)->end(); j++) {
	HcalDetId cell(j->id());
	const CaloCellGeometry* cellGeometry =
	  geometry->getSubdetectorGeometry (cell)->getGeometry (cell) ;
	double eta   = cellGeometry->getPosition().eta () ;
	double phi   = cellGeometry->getPosition().phi () ;
        double zc    = cellGeometry->getPosition().z ();
	int sub      = cell.subdet();
	int depth    = cell.depth();
	int inteta   = cell.ieta();
	if(inteta > 0) inteta -= 1;
	int intphi   = cell.iphi()-1;
	double t     = j->time();
	double en    = j->energy();
	
	if((iz > 0 && eta > 0.) || (iz < 0 && eta <0.) || iz == 0) { 
	  csub.push_back(sub);
	  cen.push_back(en);
	  ceta.push_back(eta);
	  cphi.push_back(phi);
	  ctime.push_back(t);
	  cieta.push_back(inteta);
	  ciphi.push_back(intphi);
	  cdepth.push_back(depth);
	  cz.push_back(zc);
	}
      }
      }
    }
  }      
  
}
double HcalRecHitsValidation::dR(double eta1, double phi1, double eta2, double phi2) { 
  double PI = 3.1415926535898;
  double deltaphi= phi1 - phi2;
  if( phi2 > phi1 ) { deltaphi= phi2 - phi1;}
  if(deltaphi > PI) { deltaphi = 2.*PI - deltaphi;}
  double deltaeta = eta2 - eta1;
  double tmp = sqrt(deltaeta* deltaeta + deltaphi*deltaphi);
  return tmp;
}

double HcalRecHitsValidation::phi12(double phi1, double en1, double phi2, double en2) {
  // weighted mean value of phi1 and phi2
  
  double tmp;
  double PI = 3.1415926535898;
  double a1 = phi1; double a2  = phi2;

  if( a1 > 0.5*PI  && a2 < 0.) a2 += 2*PI; 
  if( a2 > 0.5*PI  && a1 < 0.) a1 += 2*PI; 
  tmp = (a1 * en1 + a2 * en2)/(en1 + en2);
  if(tmp > PI) tmp -= 2.*PI; 
 
  return tmp;

}

double HcalRecHitsValidation::dPhiWsign(double phi1, double phi2) {
  // clockwise      phi2 w.r.t phi1 means "+" phi distance
  // anti-clockwise phi2 w.r.t phi1 means "-" phi distance 

  double PI = 3.1415926535898;
  double a1 = phi1; double a2  = phi2;
  double tmp =  a2 - a1;
  if( a1*a2 < 0.) {
    if(a1 > 0.5 * PI)  tmp += 2.*PI;
    if(a2 > 0.5 * PI)  tmp -= 2.*PI;
  }
  return tmp;

}

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DQMServices/Core/interface/DQMStore.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(HcalRecHitsValidation);

