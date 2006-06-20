/*
 * \file EcalPreshowerDigisValidation.cc
 *
 * $Date: 2006/05/04 11:16:28 $
 * $Revision: 1.3 $
 * \author F. Cossutti
 *
*/

#include <Validation/EcalDigis/interface/EcalPreshowerDigisValidation.h>

EcalPreshowerDigisValidation::EcalPreshowerDigisValidation(const ParameterSet& ps)
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

  for (int i = 0; i < 3 ; i++ ) {
    meESDigiADC_[i] = 0;
  }

  Char_t histo[20];
 
  if ( dbe_ ) {
    dbe_->setCurrentFolder("EcalDigiTask");
  
    for ( int i = 0; i < 3 ; i++ ) {
      
      sprintf (histo, "EcalDigiTask Preshower ADC pulse %02d", i+1) ;
      meESDigiADC_[i] = dbe_->book1D(histo, histo, 512, 0., 4096.) ;
    }

  }
 
}

EcalPreshowerDigisValidation::~EcalPreshowerDigisValidation(){
 
}

void EcalPreshowerDigisValidation::beginJob(const EventSetup& c){

}

void EcalPreshowerDigisValidation::endJob(){

}

void EcalPreshowerDigisValidation::analyze(const Event& e, const EventSetup& c){

  //LogInfo("EventInfo") << " Run = " << e.id().run() << " Event = " << e.id().event();

  Handle<ESDigiCollection> EcalDigiES;

  e.getByType(EcalDigiES);

  // PRESHOWER
  
  // loop over Digis

  const ESDigiCollection * preshowerDigi = EcalDigiES.product () ;

  std::vector<double> esADCCounts ;
  esADCCounts.reserve(ESDataFrame::MAXSAMPLES);

  for (std::vector<ESDataFrame>::const_iterator digis = preshowerDigi->begin () ;
       digis != preshowerDigi->end () ;
       ++digis)
    {
       
      ESDetId esid = digis->id () ;

      for (int sample = 0 ; sample < digis->size () ; ++sample) {
        esADCCounts[sample] = 0.;
      }
       
      for (int sample = 0 ; sample < digis->size () ; ++sample)
        {
          esADCCounts[sample] = (digis->sample (sample).adc ()) ;
        }
      if (verbose_) {
        LogDebug("DigiInfo") << "Preshower Digi for ESDetId: z side " << esid.zside() << "  plane " << esid.plane() << esid.six() << ',' << esid.siy() << ':' << esid.strip();
        for ( int i = 0; i < 3 ; i++ ) {
          LogDebug("DigiInfo") << "sample " << i << " ADC = " << esADCCounts[i];
        }
      }
       
      for ( int i = 0 ; i < 3 ; i++ ) {
        if (meESDigiADC_[i]) meESDigiADC_[i]->Fill( esADCCounts[i] ) ;
      }
       
    } 

}

                                                                                                                                                             
