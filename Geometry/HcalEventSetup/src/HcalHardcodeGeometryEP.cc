// -*- C++ -*-
//
// Package:    HcalHardcodeGeometryEP
// Class:      HcalHardcodeGeometryEP
// 
/**\class HcalHardcodeGeometryEP HcalHardcodeGeometryEP.h tmp/HcalHardcodeGeometryEP/interface/HcalHardcodeGeometryEP.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jeremiah Mans
//         Created:  Mon Oct  3 11:35:27 CDT 2005
// $Id: HcalHardcodeGeometryEP.cc,v 1.5.2.1 2007/05/08 20:39:50 mansj Exp $
//
//

#include "Geometry/HcalEventSetup/src/HcalHardcodeGeometryEP.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HcalHardcodeGeometryEP::HcalHardcodeGeometryEP(const edm::ParameterSet& iConfig)
{
   //the following line is needed to tell the framework what
   // data is being produced
   setWhatProduced(this,"HCAL");
   loader_=0;
}


HcalHardcodeGeometryEP::~HcalHardcodeGeometryEP()
{ 
  if (loader_) delete loader_;
}


//
// member functions
//

// ------------ method called to produce the data  ------------
HcalHardcodeGeometryEP::ReturnType
HcalHardcodeGeometryEP::produce(const IdealGeometryRecord& iRecord)
{
   using namespace edm::es;

   //now do what ever other initialization is needed
   if (loader_==0) {
     edm::ESHandle<HcalTopology> topo;
     try {
       iRecord.get(topo);
       loader_=new HcalHardcodeGeometryLoader(*topo); 
     } catch (...) {
       loader_=new HcalHardcodeGeometryLoader();
       edm::LogInfo("HCAL") << "Using default HCAL topology";
     }
   }

   std::auto_ptr<CaloSubdetectorGeometry> pCaloSubdetectorGeometry(loader_->load()) ;

   return pCaloSubdetectorGeometry ;
}


