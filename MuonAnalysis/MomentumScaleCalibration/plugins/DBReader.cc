// #include "CondFormats/PhysicsToolsObjects/interface/Histogram2D.h"
#include "CondFormats/MomentumScaleCalibrationObjects/interface/MuScleFitDBobject.h"
#include "CondFormats/DataRecord/interface/MuScleFitDBobjectRcd.h"

#include "DBReader.h"

#include <iostream>
#include <stdio.h>
#include <sys/time.h>
#include <string>

using namespace std;
using namespace cms;

DBReader::DBReader( const edm::ParameterSet& iConfig ) : type_(iConfig.getUntrackedParameter<string>("Type"))
{
}

void DBReader::beginJob ( const edm::EventSetup& iSetup ) {
  edm::ESHandle<MuScleFitDBobject> dbObject;
  iSetup.get<MuScleFitDBobjectRcd>().get(dbObject);
  edm::LogInfo("DBReader") << "[DBReader::analyze] End Reading MuScleFitDBobjectRcd" << endl;

  cout << "identifiers size from dbObject = " << dbObject->identifiers.size() << endl;
  cout << "parameters size from dbObject = " << dbObject->parameters.size() << endl;;

  // This string is one of: scale, resolution, background.
  // Create the corrector and set the parameters
  if( type_ == "scale" ) corrector_.reset(new MomentumScaleCorrector( dbObject.product() ) );
  else if( type_ == "resolution" ) resolution_.reset(new ResolutionFunction( dbObject.product() ) );
  else if( type_ == "background" ) background_.reset(new BackgroundFunction( dbObject.product() ) );
  else {
    cout << "Error: unrecognized type. Use one of those: 'scale', 'resolution', 'background'" << endl;
    exit(1);
  }

  cout << "pointer = " << corrector_.get() << endl;
}

//:  printdebug_(iConfig.getUntrackedParameter<uint32_t>("printDebug",1)){}

DBReader::~DBReader(){}

void DBReader::analyze( const edm::Event& e, const edm::EventSetup& iSetup){

//   cout << "checking size consistency" << endl;
//   if( corrector_->identifiers().size() != corrector_->parameters().size() ) {
//     cout << "Error: size of parameters("<<corrector_->parameters().size()<<") and identifiers("<<corrector_->identifiers().size()<<") don't match" << endl;
//     exit(1);
//   }

  // Looping directly on it does not work, because it is returned by value
  // and the iterator gets invalidated on the next line. Save it to a temporary object
  // and iterate on it.
  vector<double> parVecVec(corrector_->parameters());
  // vector<vector<double> >::const_iterator parVec = corrector_->parameters().begin();
  vector<double>::const_iterator parVec = parVecVec.begin();
  vector<int> functionId(corrector_->identifiers());
  vector<int>::const_iterator id = functionId.begin();
  cout << "total number of parameters read from database = parVecVec.size() = " << parVecVec.size() << endl;
  int iFunc = 0;
  for( ; id != functionId.end(); ++id, ++iFunc ) {
    int parNum = corrector_->function(iFunc)->parNum();
    cout << "For function id = " << *id << ", with "<<parNum<< " parameters: " << endl;
    for( int par=0; par<parNum; ++par ) {
      cout << "par["<<par<<"] = " << *parVec << endl;
      ++parVec;
    }
//     vector<double>::const_iterator par = parVec->begin();
//     int i=0;
//     for( ; par != parVec->end(); ++par, ++i ) {
//       cout << "par["<<i<<"] = " << *par << endl;
//     }
  }
}

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(DBReader);
