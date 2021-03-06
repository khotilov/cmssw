// -*- C++ -*-
//
// Package:    HcalChannelMon
// Class:      HcalChannelMon
// 
/**\class HcalChannelMon HcalChannelMon.cc UserCode/HcalChannelMon/src/HcalChannelMon.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Kenneth Case Rossato
//         Created:  Fri Oct 17 16:17:21 CEST 2008
// $Id: HcalChannelMon.cc,v 1.1 2008/10/30 10:36:52 rossato Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "EventFilter/HcalRawToDigi/interface/HcalDCCHeader.h"
#include "EventFilter/HcalRawToDigi/interface/HcalHTRData.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"

#include "TH1D.h"
#include <vector>

//
// class decleration
//

class HcalChannelMon : public edm::EDAnalyzer {
   public:
      explicit HcalChannelMon(const edm::ParameterSet&);
      ~HcalChannelMon();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

  std::vector<int> ietaiphi_channels;
  const HcalElectronicsMap *readoutMap_;
  std::vector<HcalElectronicsId> channels;

  std::vector<TH1D*> hists;
  edm::Service<TFileService> fs;

  bool get_capid_info;
  bool get_sums;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HcalChannelMon::HcalChannelMon(const edm::ParameterSet& iConfig)
  : ietaiphi_channels(iConfig.getParameter<std::vector<int> >("ChannelList")),
    get_capid_info(iConfig.getParameter<bool>("GetCapIdInfo")),
    get_sums(iConfig.getParameter<bool>("GetSums"))
  //    get_tp_info(iConfig.getParameter<bool>("GetTPInfo"))
{
   //now do what ever initialization is needed


}


HcalChannelMon::~HcalChannelMon()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HcalChannelMon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<FEDRawDataCollection> rawraw;
   iEvent.getByType(rawraw);

   int sums   = 2;
   int caps   = sums + (get_sums ? 2 : 0);
   int nHists = caps + (get_capid_info ? 4 : 0);

   // iterate over Electronics Id
   std::vector<TH1D*>::iterator hist = hists.begin();
   for (std::vector<HcalElectronicsId>::const_iterator cit = channels.begin();
	cit != channels.end(); cit++, hist+=nHists) {

     // get dcc data
     const FEDRawData &raw = rawraw->FEDData(700 + cit->dccid());
     const HcalDCCHeader *dccHeader = (const HcalDCCHeader *)(raw.data());

     if (!dccHeader->getSpigotPresent(cit->spigot())) continue;
     
     // get htr data
     HcalHTRData htr;
     int retval=dccHeader->getSpigotData(cit->spigot(), htr, raw.size());
     if (retval!=0) {
       if (retval==-1)
	 //if (!silent)
	 edm::LogWarning("Invalid Data") << "Invalid HTR data (data beyond payload size) observed";
       continue;
     }

     if (dccHeader->getSpigotCRCError(cit->spigot())) {
       //       if (!silent) 
       edm::LogWarning("Invalid Data") << "CRC Error on HTR data observed";
       continue;
     }  
     if (!htr.check()) {
       //       if (!silent) 
       edm::LogWarning("Invalid Data") << "Invalid HTR data observed";
       continue;
     }  
     if (htr.isHistogramEvent()) {
       //       if (!silent) 
       edm::LogWarning("Invalid Data") << "Histogram data passed to non-histogram unpacker";
       continue;
     }

     // get fiber data

     //unsigned char daq_lengths[24];
     //unsigned short daq_samples[24*20];
     //unsigned char tp_lengths[24];
     //unsigned short tp_samples[24*20];
     //htr.unpack(daq_lengths, daq_samples, tp_lengths, tp_samples);

     const unsigned short *daq_first, *daq_last;
     const unsigned short *tp_first, *tp_last;
     htr.dataPointers(&daq_first, &daq_last, &tp_first, &tp_last);
     HcalQIESample *qie_begin = (HcalQIESample*)daq_first;
     HcalQIESample *qie_end   = (HcalQIESample*)(daq_last+1);

     //find *our* fiber data, be aware of possible ZS
     int timeslice = 0;
     double total = 0;
     double adc_last = -1;
     for (HcalQIESample *qie_work = qie_begin; qie_work != qie_end;
	  qie_work++) {
       
       if (qie_work->raw()==0xFFFF) {
         continue; // filler word
       }

       if (qie_work->fiber() != cit->fiberIndex() ||
	   qie_work->fiberChan() != cit->fiberChanId())
	 continue;

       //std::cout << "timeslice: " << timeslice << " adc: " << qie_work->adc() << std::endl;

       hist[0]->Fill(timeslice++, qie_work->adc());
       hist[1]->Fill(qie_work->adc());

       if (get_sums) {
	 if ( adc_last != -1 )
	   hist[sums]->Fill(qie_work->adc() + adc_last);
	 total += qie_work->adc();       
	 adc_last = qie_work->adc();
       }

       if (get_capid_info)
	 hist[caps+qie_work->capid()]->Fill(qie_work->adc());

     }

     if (get_sums)
       hist[sums+1]->Fill(total);

   }

}


// ------------ method called once each job just before starting event loop  ------------
void 
HcalChannelMon::beginJob(const edm::EventSetup&iSetup)
{
  edm::ESHandle<HcalDbService> pSetup;
  iSetup.get<HcalDbRecord>().get(pSetup);
  readoutMap_ = pSetup->getHcalMapping();

  for (std::vector<int>::const_iterator cit = ietaiphi_channels.begin();
       cit != ietaiphi_channels.end();) {
    int subdet = *cit; cit++;
    int ieta   = *cit; cit++;
    int iphi   = *cit; cit++;
    int depth  = *cit; cit++;

    HcalDetId here(HcalSubdetector(subdet), ieta, iphi, depth);

    HcalElectronicsId there(readoutMap_->lookup(here));

    channels.push_back(there);

    
    hists.push_back(fs->make<TH1D>(std::string("iEta")+(Long_t)ieta+"_iPhi"+(Long_t)iphi+"_depth"+(Long_t)depth+"_DigiShape", 
				   std::string("Digi Shape for iEta ")+(Long_t)ieta+" iPhi "+(Long_t)iphi, 10, -0.5, 9.5));
    hists.push_back(fs->make<TH1D>(std::string("iEta")+(Long_t)ieta+"_iPhi"+(Long_t)iphi+"_depth"+(Long_t)depth+"_Spectrum",
				   std::string("ADC Spectrum for iEta ")+(Long_t)ieta+" iPhi "+(Long_t)iphi, 128, -.5, 127.5));

    if (get_sums) {
      hists.push_back(fs->make<TH1D>(std::string("iEta")+(Long_t)ieta+"_iPhi"+(Long_t)iphi+"_depth"+(Long_t)depth+"_SumOfAdjacent",
				     std::string("Sum of adjacent ADC for iEta ")+(Long_t)ieta+" iPhi "+(Long_t)iphi, 128, -.5, 254.5));
      hists.push_back(fs->make<TH1D>(std::string("iEta")+(Long_t)ieta+"_iPhi"+(Long_t)iphi+"_depth"+(Long_t)depth+"_TotalSpectrum",
				     std::string("ADC Total Spectrum for iEta ")+(Long_t)ieta+" iPhi "+(Long_t)iphi, 128, -.5, 1270.5));
    }

    /*
    if (get_tp_info) {

    } */

    if (get_capid_info) {

      hists.push_back(fs->make<TH1D>(std::string("iEta")+(Long_t)ieta+"_iPhi"+(Long_t)iphi+"_depth"+(Long_t)depth+"_Cap0Spectrum",
				     std::string("CapId 0 ADC Spectrum for iEta ")+(Long_t)ieta+" iPhi "+(Long_t)iphi, 128, -.5, 127.5));
      hists.push_back(fs->make<TH1D>(std::string("iEta")+(Long_t)ieta+"_iPhi"+(Long_t)iphi+"_depth"+(Long_t)depth+"_Cap1Spectrum",
				     std::string("CapId 1 ADC Spectrum for iEta ")+(Long_t)ieta+" iPhi "+(Long_t)iphi, 128, -.5, 127.5));
      hists.push_back(fs->make<TH1D>(std::string("iEta")+(Long_t)ieta+"_iPhi"+(Long_t)iphi+"_depth"+(Long_t)depth+"_Cap2Spectrum",
				     std::string("CapId 2 ADC Spectrum for iEta ")+(Long_t)ieta+" iPhi "+(Long_t)iphi, 128, -.5, 127.5));
      hists.push_back(fs->make<TH1D>(std::string("iEta")+(Long_t)ieta+"_iPhi"+(Long_t)iphi+"_depth"+(Long_t)depth+"_Cap3Spectrum",
				     std::string("CapId 3 ADC Spectrum for iEta ")+(Long_t)ieta+" iPhi "+(Long_t)iphi, 128, -.5, 127.5));


    }
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HcalChannelMon::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalChannelMon);
