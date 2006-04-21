/**\class EcalSimpleTBAnalyzer

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// $Id: $
//
//

#include "RecoTBCalo/EcalSimpleTBAnalysis/interface/EcalSimpleTBAnalyzer.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "TBDataFormats/EcalTBObjects/interface/EcalTBHodoscopeRecInfo.h"
#include "TBDataFormats/EcalTBObjects/interface/EcalTBTDCRecInfo.h"
#include "TBDataFormats/EcalTBObjects/interface/EcalTBEventHeader.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

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
EcalSimpleTBAnalyzer::EcalSimpleTBAnalyzer( const edm::ParameterSet& iConfig )
//========================================================================
{
   //now do what ever initialization is needed
   rootfile_          = iConfig.getUntrackedParameter<std::string>("rootfile","ecalSimpleTBanalysis.root");
   hitCollection_     = iConfig.getParameter<std::string>("hitCollection");
   hitProducer_       = iConfig.getParameter<std::string>("hitProducer");
   hodoRecInfoCollection_     = iConfig.getParameter<std::string>("hodoRecInfoCollection");
   hodoRecInfoProducer_       = iConfig.getParameter<std::string>("hodoRecInfoProducer");
   tdcRecInfoCollection_     = iConfig.getParameter<std::string>("tdcRecInfoCollection");
   tdcRecInfoProducer_       = iConfig.getParameter<std::string>("tdcRecInfoProducer");
   eventHeaderCollection_     = iConfig.getParameter<std::string>("eventHeaderCollection");
   eventHeaderProducer_       = iConfig.getParameter<std::string>("eventHeaderProducer");

   std::cout << "EcalSimpleTBAnalyzer: fetching hitCollection: " << hitCollection_.c_str()
	<< " produced by " << hitProducer_.c_str() << std::endl;

}


//========================================================================
EcalSimpleTBAnalyzer::~EcalSimpleTBAnalyzer()
//========================================================================
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  // Amplitude vs TDC offset
//   if (h_ampltdc)
//   delete h_ampltdc;
  
//   // Reconstructed energies
//   delete h_e1x1;
//   delete h_e3x3; 
//   delete h_e5x5; 
  
//   delete h_bprofx; 
//   delete h_bprofy; 
  
//   delete h_qualx; 
//   delete h_qualy; 
  
//   delete h_slopex; 
//   delete h_slopey; 
  
//   delete h_mapx; 
//   delete h_mapy; 

}

//========================================================================
void
EcalSimpleTBAnalyzer::beginJob(edm::EventSetup const&) {
//========================================================================

  // Amplitude vs TDC offset
  h_ampltdc = new TH2F("h_ampltdc","Max Amplitude vs TDC offset", 100,0.,1.,2000, 0., 20000.);

  // Reconstructed energies
  h_e1x1 = new TH1F("h_e1x1","E1x1 energy", 1500, 0., 150.);
  h_e3x3 = new TH1F("h_e3x3","E3x3 energy", 1500, 0., 150.);
  h_e5x5 = new TH1F("h_e5x5","E5x5 energy", 1500, 0., 150.);

  h_e1e9 = new TH1F("h_e1e9","E1/E9 ratio", 600, 0., 1.2);
  h_e1e25 = new TH1F("h_e1e25","E1/E25 ratio", 600, 0., 1.2);
  h_e9e25 = new TH1F("h_e9e25","E9/E25 ratio", 600, 0., 1.2);

  h_bprofx = new TH1F("h_bprofx","Beam Profile X",100,-20.,20.);
  h_bprofy = new TH1F("h_bprofy","Beam Profile Y",100,-20.,20.);

  h_qualx = new TH1F("h_qualx","Beam Quality X",5000,0.,5.);
  h_qualy = new TH1F("h_qualy","Beam Quality X",5000,0.,5.);

  h_slopex = new TH1F("h_slopex","Beam Slope X",500, -5e-4 , 5e-4 );
  h_slopey = new TH1F("h_slopey","Beam Slope Y",500, -5e-4 , 5e-4 );

  h_mapx = new TH2F("h_mapx","Max Amplitude vs X",80,-20,20,1000,0.,200.);
  h_mapy = new TH2F("h_mapy","Max Amplitude vs Y",80,-20,20,1000,0.,200.);

  h_e1e9_mapx = new TH2F("h_e1e9_mapx","E1/E9 vs X",80,-20,20,600,0.,1.2);
  h_e1e9_mapy = new TH2F("h_e1e9_mapy","E1/E9 vs Y",80,-20,20,600,0.,1.2);

  h_e1e25_mapx = new TH2F("h_e1e25_mapx","E1/E25 vs X",80,-20,20,600,0.,1.2);
  h_e1e25_mapy = new TH2F("h_e1e25_mapy","E1/E25 vs Y",80,-20,20,600,0.,1.2);

  h_e9e25_mapx = new TH2F("h_e9e25_mapx","E9/E25 vs X",80,-20,20,600,0.,1.2);
  h_e9e25_mapy = new TH2F("h_e9e25_mapy","E9/E25 vs Y",80,-20,20,600,0.,1.2);
}

//========================================================================
void
EcalSimpleTBAnalyzer::endJob() {
//========================================================================

  TFile f(rootfile_.c_str(),"RECREATE");

  // Amplitude vs TDC offset
  h_ampltdc->Write(); 

  // Reconstructed energies
  h_e1x1->Write(); 
  h_e3x3->Write(); 
  h_e5x5->Write(); 

  h_e1e9->Write(); 
  h_e1e25->Write(); 
  h_e9e25->Write(); 

  h_bprofx->Write(); 
  h_bprofy->Write(); 

  h_qualx->Write(); 
  h_qualy->Write(); 

  h_slopex->Write(); 
  h_slopey->Write(); 

  h_mapx->Write(); 
  h_mapy->Write(); 

  h_e1e9_mapx->Write(); 
  h_e1e9_mapy->Write(); 

  h_e1e25_mapx->Write(); 
  h_e1e25_mapy->Write(); 

  h_e9e25_mapx->Write(); 
  h_e9e25_mapy->Write(); 

  f.Close();
}

//
// member functions
//

//========================================================================
void
EcalSimpleTBAnalyzer::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
//========================================================================

   using namespace edm;
   using namespace cms;



   // fetch the digis and compute signal amplitude
   Handle<EBRecHitCollection> phits;
   const EBRecHitCollection* hits=0;
   try {
     //std::cout << "EcalSimpleTBAnalyzer::analyze getting product with label: " << digiProducer_.c_str()<< " prodname: " << digiCollection_.c_str() << endl;
     iEvent.getByLabel( hitProducer_, hitCollection_,phits);
     hits = phits.product(); // get a ptr to the product
     //iEvent.getByLabel( hitProducer_, phits);
   } catch ( std::exception& ex ) {
     std::cerr << "Error! can't get the product " << hitCollection_.c_str() << std::endl;
   }

   Handle<EcalTBHodoscopeRecInfo> pHodo;
   const EcalTBHodoscopeRecInfo* recHodo=0;
   try {
     //std::cout << "EcalSimpleTBAnalyzer::analyze getting product with label: " << digiProducer_.c_str()<< " prodname: " << digiCollection_.c_str() << endl;
     iEvent.getByLabel( hodoRecInfoProducer_, hodoRecInfoCollection_, pHodo);
     recHodo = pHodo.product(); // get a ptr to the product
   } catch ( std::exception& ex ) {
     std::cerr << "Error! can't get the product " << hitCollection_.c_str() << std::endl;
   }

   Handle<EcalTBTDCRecInfo> pTDC;
   const EcalTBTDCRecInfo* recTDC=0;
   try {
     //std::cout << "EcalSimpleTBAnalyzer::analyze getting product with label: " << digiProducer_.c_str()<< " prodname: " << digiCollection_.c_str() << endl;
     iEvent.getByLabel( tdcRecInfoProducer_, tdcRecInfoCollection_, pTDC);
     recTDC = pTDC.product(); // get a ptr to the product
   } catch ( std::exception& ex ) {
     std::cerr << "Error! can't get the product " << hitCollection_.c_str() << std::endl;
   }

   Handle<EcalTBEventHeader> pEventHeader;
   const EcalTBEventHeader* evtHeader=0;
   try {
     //std::cout << "EcalSimpleTBAnalyzer::analyze getting product with label: " << digiProducer_.c_str()<< " prodname: " << digiCollection_.c_str() << endl;
     iEvent.getByLabel( eventHeaderProducer_ , pEventHeader );
     evtHeader = pEventHeader.product(); // get a ptr to the product
     std::cout << "Taken EventHeader " << std::endl;
   } catch ( std::exception& ex ) {
     std::cerr << "Error! can't get the product " << hitCollection_.c_str() << std::endl;
   }
   
   if (!hits)
     return;

   if (!recTDC)
     return;

   if (!recHodo)
     return;

   if (!evtHeader)
     return;

   if (hits->size() == 0)
     return;

   // Crystal hit by beam
   //   EBDetId maxHitId(1,evtHeader->crystalInBeam(),EBDetId::SMCRYSTALMODE);

   EBDetId maxHitId(0); 
   float maxHit= -999999.;


   for(EBRecHitCollection::const_iterator ithit = hits->begin(); ithit != hits->end(); ++ithit) 
     {
       if (ithit->energy()>=maxHit)
	 {
	   maxHit=ithit->energy();
	   maxHitId=ithit->id();
	 }
       
     }   

   if (maxHitId==EBDetId(0))
     return;
   //Find EBDetId in a 5x5 Matrix (to be substituted by the Selector code)
   // Something like 
   // EBFixedWindowSelector<EcalUncalibratedRecHit> Simple5x5Matrix(hits,maxHitId,5,5);
   // std::vector<EcalUncalibratedRecHit> Energies5x5 = Simple5x5Matrix.getHits();

   EBDetId Xtals5x5[25];
   for (UInt_t icry=0;icry<25;icry++)
     {
       UInt_t row = icry / 5;
       Int_t column= icry %5;
       try
	 {
	   Xtals5x5[icry]=EBDetId(maxHitId.ism(),maxHitId.ic()-85*(row-2)+column-2,EBDetId::SMCRYSTALMODE);
	   std::cout << "**** Xtal in the matrix **** row " << row  << ", column " << column << ", xtal " << Xtals5x5[icry].ic() << std::endl;
	 }
       catch ( std::runtime_error &e )
	 {
	   std::cout << "Cannot construct 5x5 matrix around EBDetId " << maxHitId << std::endl;
	   return;
	 }
     }
   
   double amplitude[25];
   
   double amplitude3x3=0;  
   double amplitude5x5=0;  

   for (Int_t icry=0;icry<25;icry++)
     {
       amplitude[icry]=(hits->find(Xtals5x5[icry]))->energy();
       amplitude5x5 += amplitude[icry];
       // Is in 3x3?
       if ( icry == 6  || icry == 7  || icry == 8 ||
	    icry == 11 || icry == 12 || icry ==13 ||
	    icry == 16 || icry == 17 || icry ==18   )
	 {
	   amplitude3x3+=amplitude[icry];
	 }
     }

   h_e1x1->Fill(amplitude[12]);
   h_e3x3->Fill(amplitude3x3);
   h_e5x5->Fill(amplitude5x5);

   h_e1e9->Fill(amplitude[12]/amplitude3x3);
   h_e1e25->Fill(amplitude[12]/amplitude5x5);
   h_e9e25->Fill(amplitude3x3/amplitude5x5);

   if (recTDC)
     h_ampltdc->Fill(recTDC->offset(),amplitude[12]);

   if (recHodo)
     {
       float x=recHodo->posX();
       float y=recHodo->posY();
       float xslope=recHodo->slopeX();
       float yslope=recHodo->slopeY();
       float xqual=recHodo->qualX();
       float yqual=recHodo->qualY();
       
       //Filling beam profiles
       h_bprofx->Fill(x);
       h_bprofy->Fill(y);
       h_qualx->Fill(xqual);
       h_qualy->Fill(yqual);
       h_slopex->Fill(xslope);
       h_slopey->Fill(yslope);
       
       h_mapx->Fill(x,amplitude[12]);
       h_mapy->Fill(y,amplitude[12]);

       h_e1e9_mapx->Fill(x,amplitude[12]/amplitude3x3);
       h_e1e9_mapy->Fill(y,amplitude[12]/amplitude3x3);

       h_e1e25_mapx->Fill(x,amplitude[12]/amplitude5x5);
       h_e1e25_mapy->Fill(y,amplitude[12]/amplitude5x5);

       h_e9e25_mapx->Fill(x,amplitude3x3/amplitude5x5);
       h_e9e25_mapy->Fill(y,amplitude3x3/amplitude5x5);
     }

}


