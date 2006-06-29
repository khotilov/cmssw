/*
 * \file EcalBarrelRecHitsValidation.cc
 *
 * $Date: 2006/05/22 $
 * \author C. Rovelli
 *
 */

#include <Validation/EcalRecHits/interface/EcalBarrelRecHitsValidation.h>
#include <DataFormats/EcalDetId/interface/EBDetId.h>

EcalBarrelRecHitsValidation::EcalBarrelRecHitsValidation(const ParameterSet& ps){
  
  // ---------------------- 
  digiProducer_              = ps.getParameter<std::string>("digiProducer");
  uncalibrecHitProducer_     = ps.getParameter<std::string>("uncalibrecHitProducer");
  EBdigiCollection_          = ps.getParameter<std::string>("EBdigiCollection");
  EBuncalibrechitCollection_ = ps.getParameter<std::string>("EBuncalibrechitCollection");
  
  // ---------------------- 
  // verbosity switch 
  verbose_ = ps.getUntrackedParameter<bool>("verbose", false);
  
  if ( verbose_ ) {
    cout << " verbose switch is ON" << endl;
  } else {
    cout << " verbose switch is OFF" << endl;
  }
  
  // ----------------------                 
  // get hold of back-end interface 
  dbe_ = 0;
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


  // ----------------------   
  meEBUncalibRecHitsOccupancy_             = 0;
  meEBUncalibRecHitsAmplitude_             = 0;
  meEBUncalibRecHitsPedestal_              = 0;
  meEBUncalibRecHitsJitter_                = 0;
  meEBUncalibRecHitsChi2_                  = 0;
  meEBUncalibRecHitMaxSampleRatio_         = 0;
  meEBUncalibRecHitsOccupancyGt100adc_     = 0;
  meEBUncalibRecHitsAmplitudeGt100adc_     = 0;
  meEBUncalibRecHitsPedestalGt100adc_      = 0;
  meEBUncalibRecHitsJitterGt100adc_        = 0;
  meEBUncalibRecHitsChi2Gt100adc_          = 0;
  meEBUncalibRecHitMaxSampleRatioGt100adc_ = 0;
  meEBUncalibRecHitsAmpFullMap_            = 0;
  meEBUncalibRecHitsPedFullMap_            = 0;
  for (int i=0; i<36 ; i++) 
    {
      meEBUncalibRecHitAmplMap_[i] = 0;
      meEBUncalibRecHitPedMap_[i]  = 0;
    }

  // ---------------------- 
  Char_t histo[200];
   
  if ( dbe_ ) 
    {
      dbe_->setCurrentFolder("EcalBarrelRecHitsTask");
      
      sprintf (histo, "EB Occupancy" );  
      meEBUncalibRecHitsOccupancy_ = dbe_->book2D(histo, histo, 170, -85., 85., 360, 0., 360.);
      
      sprintf (histo, "EB Amplitude" );
      meEBUncalibRecHitsAmplitude_ = dbe_->book1D(histo, histo, 200, 0., 4000.);
      
      sprintf (histo, "EB Pedestal" );
      meEBUncalibRecHitsPedestal_ = dbe_->book1D(histo, histo, 50, 194., 201.);
      
      sprintf (histo, "EB Jitter" );
      meEBUncalibRecHitsJitter_ = dbe_->book1D(histo, histo, 100, 0., 100.);
      
      sprintf (histo, "EB Chi2" );
      meEBUncalibRecHitsChi2_ = dbe_->book1D(histo, histo, 100, 0., 100.);

      sprintf (histo, "EB RecHit Max Sample Ratio"); 
      meEBUncalibRecHitMaxSampleRatio_ = dbe_->book1D(histo, histo, 80, 0.95, 1.05);

      sprintf (histo, "EB Occupancy gt 100 adc counts" );
      meEBUncalibRecHitsOccupancyGt100adc_ = dbe_->book2D(histo, histo, 170, -85., 85., 360, 0., 360.);
      
      sprintf (histo, "EB Amplitude gt 100 adc counts" );
      meEBUncalibRecHitsAmplitudeGt100adc_ = dbe_->book1D(histo, histo, 200, 0., 4000.);

      sprintf (histo, "EB Pedestal gt 100 adc counts" );
      meEBUncalibRecHitsPedestalGt100adc_ = dbe_->book1D(histo, histo, 50, 194., 201.);

      sprintf (histo, "EB Jitter gt 100 adc counts" );
      meEBUncalibRecHitsJitterGt100adc_ = dbe_->book1D(histo, histo, 100, 0., 100.);

      sprintf (histo, "EB Chi2 gt 100 adc counts" );
      meEBUncalibRecHitsChi2Gt100adc_ = dbe_->book1D(histo, histo, 100, 0., 100.);
    
      sprintf (histo, "EB RecHit Max Sample Ratio gt 100 adc counts"); 
      meEBUncalibRecHitMaxSampleRatioGt100adc_ = dbe_->book1D(histo, histo, 80, 0.95, 1.05);
      
      sprintf (histo, "EB Amplitude Full Map");
      meEBUncalibRecHitsAmpFullMap_ = dbe_->bookProfile2D(histo, histo, 170, -85., 85., 360, 0., 360., 200, 0., 4000.);

      sprintf (histo, "EB Pedestal Full Map");
      meEBUncalibRecHitsPedFullMap_ = dbe_->bookProfile2D(histo, histo, 170, -85., 85., 360, 0., 360., 50, 194., 201.);

      for (int i=0; i<36 ; i++) 
	{
	  sprintf(histo, "EB Amp SM%02d", i+1);
	  meEBUncalibRecHitAmplMap_[i] = dbe_->bookProfile2D(histo, histo, 85, 0., 85., 20, 0., 20., 200, 0., 4000.);
	  
	  sprintf(histo, "EB Ped SM%02d", i+1);
	  meEBUncalibRecHitPedMap_[i] = dbe_->bookProfile2D(histo, histo, 85, 0., 85., 20, 0., 20., 50, 194., 201.);
	}
    }
}

EcalBarrelRecHitsValidation::~EcalBarrelRecHitsValidation(){   

}

void EcalBarrelRecHitsValidation::beginJob(const EventSetup& c){  

}

void EcalBarrelRecHitsValidation::endJob(){

}

void EcalBarrelRecHitsValidation::analyze(const Event& e, const EventSetup& c){
  
  Handle< EBDigiCollection > EcalDigiEB;
  try {
    e.getByLabel( digiProducer_, EcalDigiEB);
  } catch ( std::exception& ex ) {
    edm::LogError("EcalRecHitsTaskError") << "Error! can't get the Digis " << std::endl;
  }
  
  Handle< EBUncalibratedRecHitCollection > EcalUncalibRecHitEB;
  try {
    e.getByLabel( uncalibrecHitProducer_, EBuncalibrechitCollection_, EcalUncalibRecHitEB);
  } catch ( std::exception& ex ) {
    edm::LogError("EcalRecHitsTaskError") << "Error! can't get the product " << EBuncalibrechitCollection_.c_str() << std::endl;
  }

  edm::ESHandle<EcalPedestals> ecalPeds; 
  try {
    c.get<EcalPedestalsRcd>().get(ecalPeds);
  } catch ( std::exception& ex ) {
    edm::LogError("EcalRecHitsTaskError") << "Error! can't get the Ecal pedestals" << std::endl;
  }


  // ---------------------- 
  // loop over UncalibRecHits
  const EBDigiCollection *EBDigi = EcalDigiEB.product();
  const EBUncalibratedRecHitCollection *EBUncalibRecHit = EcalUncalibRecHitEB.product();
  
  for (EcalUncalibratedRecHitCollection::const_iterator uncalibRecHit = EBUncalibRecHit->begin(); uncalibRecHit != EBUncalibRecHit->end() ; ++uncalibRecHit)
    {
      EBDetId EBid = EBDetId(uncalibRecHit->id());

      // general checks
      if (meEBUncalibRecHitsOccupancy_)  meEBUncalibRecHitsOccupancy_  -> Fill(EBid.ieta(), EBid.iphi());
      if (meEBUncalibRecHitsAmplitude_)  meEBUncalibRecHitsAmplitude_  -> Fill(uncalibRecHit->amplitude());
      if (meEBUncalibRecHitsPedestal_)   meEBUncalibRecHitsPedestal_   -> Fill(uncalibRecHit->pedestal());
      if (meEBUncalibRecHitsJitter_)     meEBUncalibRecHitsJitter_     -> Fill(uncalibRecHit->jitter());
      if (meEBUncalibRecHitsChi2_)       meEBUncalibRecHitsChi2_       -> Fill(uncalibRecHit->chi2());
      if (meEBUncalibRecHitsAmpFullMap_) meEBUncalibRecHitsAmpFullMap_ -> Fill(EBid.ieta(), EBid.iphi(), uncalibRecHit->amplitude()); 
      if (meEBUncalibRecHitsPedFullMap_) meEBUncalibRecHitsPedFullMap_ -> Fill(EBid.ieta(), EBid.iphi(), uncalibRecHit->pedestal()); 

      // general checks, with threshold at 3.5 GeV = 100 ADC counts
      if ( uncalibRecHit->amplitude() > 100 ) 
	{
	  if (meEBUncalibRecHitsOccupancyGt100adc_)  meEBUncalibRecHitsOccupancyGt100adc_  -> Fill(EBid.ieta(), EBid.iphi());
	  if (meEBUncalibRecHitsAmplitudeGt100adc_)  meEBUncalibRecHitsAmplitudeGt100adc_  -> Fill(uncalibRecHit->amplitude());
	  if (meEBUncalibRecHitsPedestalGt100adc_)   meEBUncalibRecHitsPedestalGt100adc_   -> Fill(uncalibRecHit->pedestal());
	  if (meEBUncalibRecHitsJitterGt100adc_)     meEBUncalibRecHitsJitterGt100adc_     -> Fill(uncalibRecHit->jitter());
	  if (meEBUncalibRecHitsChi2Gt100adc_)       meEBUncalibRecHitsChi2Gt100adc_       -> Fill(uncalibRecHit->chi2());
	}

      // supermodule maps
      int ic = EBid.ic();              
      int ie = (ic-1)/20 + 1;          
      int ip = (ic-1)%20 + 1;
      int ism = EBid.ism();
      float xie = ie - 0.5;
      float xip = ip - 0.5;      
      meEBUncalibRecHitPedMap_[ism-1] ->Fill(xie, xip, uncalibRecHit->pedestal());
      meEBUncalibRecHitAmplMap_[ism-1]->Fill(xie, xip, uncalibRecHit->amplitude());
      

      // find the rechit corresponding digi and the max sample
      EBDigiCollection::const_iterator myDigi = EBDigi->find(EBid);
      int sMax = -1;
      double eMax = 0.;
      if (myDigi != EBDigi->end())
	{
	  for (int sample = 0; sample < myDigi->size(); ++sample)
	    {
	      double analogSample = myDigi->sample(sample).adc();
	      if ( eMax < analogSample )
		{
		  eMax = analogSample;
		  sMax = sample;
		}
	    }
	}
      else
	continue;
      
      // ratio uncalibratedRecHit amplitude + ped / max energy digi  
      const EcalPedestals* myped = ecalPeds.product();
      std::map<const unsigned int,EcalPedestals::Item>::const_iterator it=myped->m_pedestals.find(EBid.rawId());
      if( it != myped->m_pedestals.end() )
	{
	  if (eMax > it->second.mean_x1 + 5 * it->second.rms_x1 ) //only real signal RecHit
	    {
	      if ( meEBUncalibRecHitMaxSampleRatio_ ) meEBUncalibRecHitMaxSampleRatio_->Fill( (uncalibRecHit->amplitude()+uncalibRecHit->pedestal())/eMax);
	      if ( meEBUncalibRecHitMaxSampleRatioGt100adc_ && (uncalibRecHit->amplitude()>100) ) meEBUncalibRecHitMaxSampleRatioGt100adc_->Fill( (uncalibRecHit->amplitude()+uncalibRecHit->pedestal())/eMax);
	      LogInfo("EcalRecHitsTaskInfo") << "barrel, eMax = " << eMax << " Amplitude = " << uncalibRecHit->amplitude()+uncalibRecHit->pedestal();  
	    }
	  else
	    continue;
	}
      else
	continue;
       
    }  // loop over the UncalibratedRecHitCollection

}
