#include "RecoTBCalo/EcalTBHodoscopeReconstructor/interface/EcalTBHodoscopeRecInfoProducer.h"
#include "TBDataFormats/EcalTBObjects/interface/EcalTBHodoscopeRawInfo.h"
#include "TBDataFormats/EcalTBObjects/interface/EcalTBHodoscopeRecInfo.h"
#include "DataFormats/Common/interface/EDCollection.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/Selector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
    
EcalTBHodoscopeRecInfoProducer::EcalTBHodoscopeRecInfoProducer(edm::ParameterSet const& ps)
{
  rawInfoCollection_ = ps.getParameter<std::string>("rawInfoCollection");
  rawInfoProducer_   = ps.getParameter<std::string>("rawInfoProducer");
  recInfoCollection_        = ps.getParameter<std::string>("recInfoCollection");
  fitMethod_ = ps.getParameter<int>("fitMethod");

//   std::vector<double> planeShift_def;
//   planeShift_def.push_back( -0.333 );
//   planeShift_def.push_back( -0.333 );
//   planeShift_def.push_back( -0.333 );
//   planeShift_def.push_back( -0.333 );
  std::vector<double> planeShift = ps.getParameter< std::vector<double> >("planeShift");

//   std::vector<double> zPosition_def;
//   zPosition_def.push_back( -0.333 );
//   zPosition_def.push_back( -0.333 );
//   zPosition_def.push_back( -0.333 );
//   zPosition_def.push_back( -0.333 );
  std::vector<double> zPosition = ps.getParameter< std::vector<double> >("zPosition");
  
  produces<EcalTBHodoscopeRecInfo>(recInfoCollection_);
  
  algo_ = new EcalTBHodoscopeRecInfoAlgo(fitMethod_, planeShift, zPosition);
}

EcalTBHodoscopeRecInfoProducer::~EcalTBHodoscopeRecInfoProducer() {
  delete algo_;
}

void EcalTBHodoscopeRecInfoProducer::produce(edm::Event& e, const edm::EventSetup& es)
{
  // Get input
   edm::Handle<EcalTBHodoscopeRawInfo> ecalRawHodoscope;  
   try {
     //evt.getByLabel( digiProducer_, digiCollection_, pDigis);
     e.getByLabel( rawInfoProducer_, ecalRawHodoscope);
   } catch ( std::exception& ex ) {
     edm::LogError("EcalTBHodoscopeRecInfoError") << "Error! can't get the product " << rawInfoCollection_.c_str() ;
   }

  // Create empty output
  std::auto_ptr<EcalTBHodoscopeRecInfo> recInfo(new EcalTBHodoscopeRecInfo(algo_->reconstruct(*ecalRawHodoscope)));
  
  e.put(recInfo,recInfoCollection_);
} 


