/**\class DumpADC

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// $Id: DumpADC.cc,v 1.7 2006/10/12 11:30:44 govoni Exp $
//
//

#include "Calibration/EcalTBProducers/interface/DumpADC.h"
#include "DataFormats/EcalRawData/interface/EcalRawDataCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "TBDataFormats/EcalTBObjects/interface/EcalTBHodoscopeRecInfo.h"
#include "TBDataFormats/EcalTBObjects/interface/EcalTBTDCRecInfo.h"
#include "TBDataFormats/EcalTBObjects/interface/EcalTBEventHeader.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "RecoEcal/EgammaClusterAlgos/interface/LogPositionCalc.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


//#include<fstream>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"

#include <iostream>
#include <string>
#include <stdexcept>
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//


//========================================================================
DumpADC::DumpADC( const edm::ParameterSet& iConfig ) :
  m_TB06Tree (iConfig.getUntrackedParameter<std::string>("rootfile","TB06Tree.root"),
               iConfig.getUntrackedParameter<std::string>("treename","Analysis"))
//========================================================================
{
   m_eventNum = 0 ; // FIXME
   m_minEvent = iConfig.getUntrackedParameter<int> ("minEvent",0) ; // FIXME
   m_beamEnergy = iConfig.getUntrackedParameter<double> ("beamEnergy",120) ;

   //now do what ever initialization is needed
//   rootfile_          = iConfig.getUntrackedParameter<std::string>("rootfile","ecalSimpleTBanalysis.root") ;
   ECALRawDataProducer_       = iConfig.getParameter<std::string> ("ECALRawDataProducer") ;
   ECALRawDataCollection_     = iConfig.getParameter<std::string> ("ECALRawDataCollection") ;
   EBuncalibRecHitCollection_ = iConfig.getParameter<std::string> ("EBuncalibRecHitCollection") ;
   uncalibRecHitProducer_     = iConfig.getParameter<std::string> ("uncalibRecHitProducer") ;
   hodoRecInfoCollection_     = iConfig.getParameter<std::string> ("hodoRecInfoCollection") ;
   hodoRecInfoProducer_       = iConfig.getParameter<std::string> ("hodoRecInfoProducer") ;
   tdcRecInfoCollection_      = iConfig.getParameter<std::string> ("tdcRecInfoCollection") ;
   tdcRecInfoProducer_        = iConfig.getParameter<std::string> ("tdcRecInfoProducer") ;
   eventHeaderCollection_     = iConfig.getParameter<std::string> ("eventHeaderCollection") ;
   eventHeaderProducer_       = iConfig.getParameter<std::string> ("eventHeaderProducer") ;
   // Use the maximum energy xtal instead of the xtal on beam information
   useMaxEnergyXtal_          = iConfig.getUntrackedParameter<bool> ("useMaxEnergyXtal") ;

   edm::LogInfo ("CTOR") << "Use Xtal in Beam information: " <<useMaxEnergyXtal_<< std::endl ;

   edm::LogInfo ("CTOR") << "DumpADC: fetching hitCollection: " << EBuncalibRecHitCollection_
                         << " produced by " << uncalibRecHitProducer_ << std::endl ;

}


//========================================================================
DumpADC::~DumpADC()
//========================================================================
{}

//========================================================================
void
DumpADC::beginJob (const edm::EventSetup & iSetup)
//========================================================================
{}

//========================================================================
void
DumpADC::endJob()
//========================================================================
{}

//
// member functions
//

//========================================================================
void
DumpADC::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
//========================================================================
{
//   std::cout << "PIETRO " << m_eventNum << "\t" << (m_eventNum < m_minEvent) << std::endl ;
   if (m_eventNum++ < m_minEvent) return ; // FIXME

   using namespace edm ;
   using namespace cms ;

   // fetch the uncalibrated digis
   // ----------------------------

   Handle< EBUncalibratedRecHitCollection > pEBUncalibRecHits ;
   const EBUncalibratedRecHitCollection*  EBuncalibRecHits = 0 ;

   try {
       iEvent.getByLabel (uncalibRecHitProducer_,
                          EBuncalibRecHitCollection_,
                          pEBUncalibRecHits) ;
       EBuncalibRecHits = pEBUncalibRecHits.product () ; // get a ptr to the product
     } catch ( std::exception& ex ) {
       edm::LogError ("missed") << "get collection " << EBuncalibRecHitCollection_ 
                                << " from producer " << uncalibRecHitProducer_
                                << " failed: " << std::endl ;
       edm::LogError ("missed") << ex.what () << std::endl ;
     }

   if (!EBuncalibRecHits) return ;
   if (EBuncalibRecHits->size () == 0) return ;

   Handle<EcalTBHodoscopeRecInfo> pHodo ;
   const EcalTBHodoscopeRecInfo* recHodo=0 ;

   try {
     iEvent.getByLabel(hodoRecInfoProducer_,
                       hodoRecInfoCollection_,
                       pHodo) ;
     recHodo = pHodo.product() ; // get a ptr to the product
   } catch ( std::exception& ex ) {
     edm::LogError ("missed") << ex.what () << std::endl ;
   }

   Handle<EcalTBTDCRecInfo> pTDC ;
   const EcalTBTDCRecInfo* recTDC=0 ;

   try {
     iEvent.getByLabel (tdcRecInfoProducer_, tdcRecInfoCollection_, pTDC) ;
     recTDC = pTDC.product () ; // get a ptr to the product
   } catch ( std::exception& ex ) {
     edm::LogError ("missed") << ex.what () << std::endl ;
   }

   Handle<EcalTBEventHeader> pEventHeader ;
   const EcalTBEventHeader* evtHeader=0 ;

   try {
     iEvent.getByLabel ( eventHeaderProducer_ , pEventHeader ) ;
     evtHeader = pEventHeader.product () ; // get a ptr to the product
   } catch ( std::exception& ex ) {
     edm::LogError ("missed") << ex.what () << std::endl ;
   }

   // FIXME if (!recHodo) return ;
   // FIXME if (!evtHeader) return ;

/*
   Handle<EcalRawDataCollection> dcchs ;
   try {
//     iEvent.getByLabel (ECALRawDataProducer_, ECALRawDataCollection_, dcchs) ;
     iEvent.getByLabel (ECALRawDataProducer_, dcchs) ;

     //FIXME questo e' un loop su che cosa?
     for ( EcalRawDataCollection::const_iterator dcchItr = dcchs->begin () ;
           dcchItr != dcchs->end () ;
           ++dcchItr )
       {
         EcalDCCHeaderBlock dcch = (*dcchItr) ;
         if ( dcch.getRunType () == EcalDCCHeaderBlock::BEAMH4 ||
              dcch.getRunType () == EcalDCCHeaderBlock::BEAMH2  )
           {
             std::cerr << "[DumpADC][analyze] " << dcch.getRunType () << std::endl ;
             int enable = true ;
           }
       }
   } catch ( std::exception& ex ) {
     std::cerr << ex.what () << std::endl ;
   }
*/

   // get the header info
   // -------------------

   int xtalNum = -1 ;
   int run = -1 ;
   int event = -1 ;
   int tableIsMoving = -1 ;
   int S6ADC = 0 ;
   std::string eventType = "0" ;
   if (evtHeader)
     {
     //D LogDebug("EvtHeaderDebug") << "Found evtHeader " << std::endl ;

       eventType = evtHeader->eventType () ;
       xtalNum = evtHeader->crystalInBeam () ;
       run = evtHeader->runNumber () ;
       event = evtHeader->eventNumber () ;
       tableIsMoving = evtHeader->tableIsMoving () ;
       S6ADC = evtHeader->S6ADC () ;
       /*
       LogDebug("EvtHeaderDebug") << "eventType "     << eventType     << "\n"
                                  << "xtalNum "       << xtalNum       << "\n"
                                  << "run "           << run           << "\n"
                                  << "event "         << event         << "\n"
                                  << "tableIsMoving " << tableIsMoving << "\n"
                                  << "S6ADC "         << S6ADC         << "\n" ;
        */         
     }

   if (not evtHeader or useMaxEnergyXtal_)
     {
     //D
       /*
       LogDebug("EvtHeaderDebug") << "evtHeader= " << evtHeader << " and "
                 << "useMaxEnergyXtal_= " << useMaxEnergyXtal_ << "\n" ;
       */          
       float maxAmplitude = -999999. ;
       DetId MExtal (0) ;
       // loop over uncalibrated rechits
       for (EBUncalibratedRecHitCollection::const_iterator EBUit
                                       = EBuncalibRecHits->begin () ;
            EBUit != EBuncalibRecHits->end () ;
            ++EBUit)
         {
           if (maxAmplitude < EBUit->amplitude ())
             {
               maxAmplitude = EBUit->amplitude () ;
               MExtal = EBUit->id () ;
             }
         } // loop over uncalibrated rechits
       if (maxAmplitude < -999990) return ;
       EBDetId EBMExtal (MExtal) ;
       xtalNum = EBMExtal.ic () ;
       //PG std::cout << "xtalNum " << xtalNum << std::endl ;
   }

   // get the hodoscope info
   // ----------------------

   float xhodo = -99. ;
   float yhodo = -99. ;
   float xslope = -99. ;
   float yslope = -99. ;
   float xqual = -99. ;
   float yqual = -99. ;

   if (recHodo)
     {
       xhodo=recHodo->posX () ;
       yhodo=recHodo->posY () ;
       xslope=recHodo->slopeX () ;
       yslope=recHodo->slopeY () ;
       xqual=recHodo->qualX () ;
       yqual=recHodo->qualY () ;
     }

//   edm::logDebug ("check") << "xhodo " << xhodo  << std::endl
//                           << "yhodo"  << yhodo  << std::endl
//                           << "xslope" << xslope << std::endl
//                           << "yslope" << yslope << std::endl
//                           << "xqual"  << xqual  << std::endl
//                           << "yqual"  << yqual  << std::endl ;

   // get the 7x7 matrix
   // ------------------

   double amplitude[49] ;
   //FIXME not sure about the "1"
   //EBDetId EBMExtal (1,xtalNum,EBDetId::SMCRYSTALMODE) ;

   EBDetId EBMExtal;
   try  {
         EBDetId dummy (1,xtalNum,EBDetId::SMCRYSTALMODE) ;
         EBMExtal = dummy;
        }
   catch ( cms::Exception &e )
        {
         edm::LogError ("missed") << "Problems in generating central crystal EBDetId object : "
                                  << "RunNb= " << run
                                  << " EvtNb= " << event
                                  << " XtalNb= " << xtalNum <<std::endl ;
         return;
        }

   // loop over the 7x7 matrix
   for (UInt_t icry=0 ; icry<49 ; ++icry)
     {
       UInt_t row = icry / 7 ;
       Int_t column = icry % 7 ;

       try {
             EBDetId tempo (EBMExtal.ieta ()+row-3,
                            EBMExtal.iphi ()+column-3, 
                            EBDetId::ETAPHIMODE) ;

             if (tempo.ism () == 1)
               {
                 amplitude [icry] = EBuncalibRecHits->find (tempo)->
                                                        amplitude () ;
//                 std::cout << "[CR][PG] Building element (" << tempo
//                           << " around " << EBMExtal
//                           << " ene " << amplitude [icry] << std::endl ;
               }
              else amplitude [icry] = 0. ;
           }
       catch ( cms::Exception &e )
	    {
             amplitude[icry] = 0. ;
             edm::LogError ("missed") << "Cannot get amplitude" << std::endl ;
           }
     } // loop over the 7x7 matrix

   /*
   std::cout << "Summary of TTree stored material: " << "\n"
             << "tableIsMoving " << tableIsMoving<<"\n"
             << "run" << run << "\n"
             << "event " << event << "\n"
             << "S6ADC " << S6ADC << "\n"
             << "xhodo " << xhodo << "\n"
             << "yhodo " << yhodo << "\n"
             << "xslope " << xslope << "\n"
             << "yslope " << yslope << "\n"
             << "xqual " << xqual << "\n"
             << "yqual " << yqual << "\n"
             << "xtalNum " << xtalNum << "\n"
             << "EBMExtal.ieta () " << EBMExtal.ieta () << "\n"
             << "EBMExtal.iphi () " << EBMExtal.iphi () << "\n"
             << "m_beamEnergy " << m_beamEnergy << "\n"
             << "amplitude " << amplitude << "\n" ;
   */ 

   m_TB06Tree.store (tableIsMoving,
                     run,event,S6ADC,
                     xhodo,yhodo,
                     xslope, yslope,
                     xqual, yqual,
                     xtalNum,
                     EBMExtal.ieta (), EBMExtal.iphi (),
                     m_beamEnergy,
                     amplitude) ;

//  m_TB06Tree.check () ;

}


