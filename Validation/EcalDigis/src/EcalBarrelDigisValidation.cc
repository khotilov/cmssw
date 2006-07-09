/*
 * \file EcalBarrelDigisValidation.cc
 *
 * $Date: 2006/06/23 17:29:58 $
 * $Revision: 1.6 $
 * \author F. Cossutti
 *
*/

#include <Validation/EcalDigis/interface/EcalBarrelDigisValidation.h>

EcalBarrelDigisValidation::EcalBarrelDigisValidation(const ParameterSet& ps)
  {
 
  // verbosity switch
  verbose_ = ps.getUntrackedParameter<bool>("verbose", false);
 
  if ( verbose_ ) {
    cout << " verbose switch is ON" << endl;
  } else {
    cout << " verbose switch is OFF" << endl;
  }
                                                                                                                                          
  dbe_ = 0;
                                                                                                                                          
  // get hold of back-end interface
  dbe_ = Service<DaqMonitorBEInterface>().operator->();
                                                                                                                                          
  if ( dbe_ ) {
    if ( verbose_ ) {
      dbe_->setVerbose(1);
    } else {
      dbe_->setVerbose(0);
    }
  }
                                                                                                                                          
  if ( dbe_ ) {
    if ( verbose_ ) dbe_->showDirStructure();
  }

  gainConv_[0] = 0.;
  gainConv_[1] = 1.;
  gainConv_[2] = 2.;
  gainConv_[3] = 12.;
  barrelADCtoGeV_ = 0.035;
  endcapADCtoGeV_ = 0.06;
 
  meEBDigiOccupancy_ = 0;

  meEBDigiADCGlobal_ = 0;

  for (int i = 0; i < 10 ; i++ ) {
    meEBDigiADCAnalog_[i] = 0;
    meEBDigiADCg1_[i] = 0;
    meEBDigiADCg6_[i] = 0;
    meEBDigiADCg12_[i] = 0;
    meEBDigiGain_[i] = 0;
  }

  meEBPedestal_ = 0;
                                 
  meEBMaximumgt100ADC_ = 0; 
                                 
  meEBMaximumgt10ADC_ = 0; 

  meEBnADCafterSwitch_ = 0;
 
  Char_t histo[20];
 
  
  if ( dbe_ ) {
    dbe_->setCurrentFolder("EcalDigiTask");
  
    sprintf (histo, "EcalDigiTask Barrel occupancy" ) ;
    meEBDigiOccupancy_ = dbe_->book2D(histo, histo, 360, 0., 360., 170, -85., 85.);
  
    sprintf (histo, "EcalDigiTask Barrel global pulse shape" ) ;
    meEBDigiADCGlobal_ = dbe_->bookProfile(histo, histo, 10, 0, 10, 10000, 0., 1000.) ;
    
    for (int i = 0; i < 10 ; i++ ) {

      sprintf (histo, "EcalDigiTask Barrel analog pulse %02d", i+1) ;
      meEBDigiADCAnalog_[i] = dbe_->book1D(histo, histo, 512, 0., 4096.);

      sprintf (histo, "EcalDigiTask Barrel ADC pulse %02d Gain 1", i+1) ;
      meEBDigiADCg1_[i] = dbe_->book1D(histo, histo, 512, 0., 4096);

      sprintf (histo, "EcalDigiTask Barrel ADC pulse %02d Gain 6", i+1) ;
      meEBDigiADCg6_[i] = dbe_->book1D(histo, histo, 512, 0., 4096);

      sprintf (histo, "EcalDigiTask Barrel ADC pulse %02d Gain 12", i+1) ;
      meEBDigiADCg12_[i] = dbe_->book1D(histo, histo, 512, 0., 4096);

      sprintf (histo, "EcalDigiTask Barrel gain pulse %02d", i+1) ;
      meEBDigiGain_[i] = dbe_->book1D(histo, histo, 4, 0, 4);

    }
    
    sprintf (histo, "EcalDigiTask Barrel pedestal for pre-sample" ) ;
    meEBPedestal_ = dbe_->book1D(histo, histo, 512, 0., 4096.) ;

    sprintf (histo, "EcalDigiTask Barrel maximum position gt 100 ADC" ) ;
    meEBMaximumgt100ADC_ = dbe_->book1D(histo, histo, 10, 0., 10.) ;

    sprintf (histo, "EcalDigiTask Barrel maximum position gt 10 ADC" ) ;
    meEBMaximumgt10ADC_ = dbe_->book1D(histo, histo, 10, 0., 10.) ;

    sprintf (histo, "EcalDigiTask Barrel ADC counts after gain switch" ) ;
    meEBnADCafterSwitch_ = dbe_->book1D(histo, histo, 10, 0., 10.) ;

  }
 
}

EcalBarrelDigisValidation::~EcalBarrelDigisValidation(){
 
}

void EcalBarrelDigisValidation::beginJob(const EventSetup& c){

}

void EcalBarrelDigisValidation::endJob(){

}

void EcalBarrelDigisValidation::analyze(const Event& e, const EventSetup& c){

  //LogInfo("EventInfo") << " Run = " << e.id().run() << " Event = " << e.id().event();

  Handle<EBDigiCollection> EcalDigiEB;

  e.getByType(EcalDigiEB);

  // BARREL

  // loop over Digis

  const EBDigiCollection * barrelDigi = EcalDigiEB.product () ;

  std::vector<double> ebAnalogSignal ;
  std::vector<double> ebADCCounts ;
  std::vector<double> ebADCGains ;
  ebAnalogSignal.reserve(EBDataFrame::MAXSAMPLES);
  ebADCCounts.reserve(EBDataFrame::MAXSAMPLES);
  ebADCGains.reserve(EBDataFrame::MAXSAMPLES);

  for (std::vector<EBDataFrame>::const_iterator digis = barrelDigi->begin () ;
       digis != barrelDigi->end () ;
       ++digis)
    {
    
      EBDetId ebid = digis->id () ;

      if (meEBDigiOccupancy_) meEBDigiOccupancy_->Fill( ebid.iphi(), ebid.ieta() );

      double Emax = 0. ;
      int Pmax = 0 ;
      double pedestalPreSample = 0.;
      double pedestalPreSampleAnalog = 0.;
      int countsAfterGainSwitch = -1;
      double higherGain = 1.;
      int higherGainSample = 0;

      for (int sample = 0 ; sample < digis->size () ; ++sample) {
        ebAnalogSignal[sample] = 0.;
        ebADCCounts[sample] = 0.;
        ebADCGains[sample] = 0.;
      }

      for (int sample = 0 ; sample < digis->size () ; ++sample)
        {
          ebADCCounts[sample] = (digis->sample (sample).adc ()) ;
          ebADCGains[sample] = (digis->sample (sample).gainId ()) ;
          ebAnalogSignal[sample] = (ebADCCounts[sample]*gainConv_[(int)ebADCGains[sample]]*barrelADCtoGeV_);
          if (Emax < ebAnalogSignal[sample] ) {
            Emax = ebAnalogSignal[sample] ;
            Pmax = sample ;
          }
          if ( sample < 3 ) {
            pedestalPreSample += ebADCCounts[sample] ;
            pedestalPreSampleAnalog += ebADCCounts[sample]*gainConv_[(int)ebADCGains[sample]]*barrelADCtoGeV_ ;
          }
          if ( sample > 0 && ebADCGains[sample] > ebADCGains[sample-1] ) {
            higherGain = ebADCGains[sample];
            higherGainSample = sample;
            countsAfterGainSwitch = 1;
          }
          if ( higherGain > 1 && higherGainSample != sample && ebADCGains[sample] == higherGain) countsAfterGainSwitch++ ;
        }
      pedestalPreSample /= 3. ; 
      pedestalPreSampleAnalog /= 3. ; 

      LogDebug("DigiInfo") << "Barrel Digi for EBDetId = " << ebid.rawId() << " eta,phi " << ebid.ieta() << " " << ebid.iphi() ;
      for ( int i = 0; i < 10 ; i++ ) {
        LogDebug("DigiInfo") << "sample " << i << " ADC = " << ebADCCounts[i] << " gain = " << ebADCGains[i] << " Analog = " << ebAnalogSignal[i];
      }
      LogDebug("DigiInfo") << "Maximum energy = " << Emax << " in sample " << Pmax << " Pedestal from pre-sample = " << pedestalPreSampleAnalog;
      if ( countsAfterGainSwitch > 0 ) LogDebug("DigiInfo") << "Counts after switch " << countsAfterGainSwitch;
        
      for ( int i = 0 ; i < 10 ; i++ ) {
        if (meEBDigiADCGlobal_) meEBDigiADCGlobal_->Fill( i , ebAnalogSignal[i] ) ;
        if (meEBDigiADCAnalog_[i]) meEBDigiADCAnalog_[i]->Fill( ebAnalogSignal[i]*100. ) ;
        if ( ebADCGains[i] == 3 ) {
          if (meEBDigiADCg1_[i]) meEBDigiADCg1_[i]->Fill( ebADCCounts[i] ) ;
        }
        else if ( ebADCGains[i] == 2 ) {
          if (meEBDigiADCg6_[i]) meEBDigiADCg6_[i]->Fill( ebADCCounts[i] ) ;
        }
        else if ( ebADCGains[i] == 1 ) {
          if (meEBDigiADCg12_[i]) meEBDigiADCg12_[i]->Fill( ebADCCounts[i] ) ;
        }
        if (meEBDigiGain_[i]) meEBDigiGain_[i]->Fill( ebADCGains[i] ) ;
      }

      if (meEBPedestal_) meEBPedestal_->Fill ( pedestalPreSample ) ;
      if (meEBMaximumgt10ADC_ && (Emax-pedestalPreSampleAnalog) > 10.*barrelADCtoGeV_) meEBMaximumgt10ADC_->Fill( Pmax ) ;
      if (meEBMaximumgt100ADC_ && (Emax-pedestalPreSampleAnalog) > 100.*barrelADCtoGeV_) meEBMaximumgt100ADC_->Fill( Pmax ) ;
      if (meEBnADCafterSwitch_) meEBnADCafterSwitch_->Fill( countsAfterGainSwitch ) ;
        
    } 

}

                                                                                                                                                             
