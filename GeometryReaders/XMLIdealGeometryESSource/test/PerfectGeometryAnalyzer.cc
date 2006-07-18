// -*- C++ -*-
//
// Package:    PerfectGeometryAnalyzer
// Class:      PerfectGeometryAnalyzer
// 
/**\class PerfectGeometryAnalyzer PerfectGeometryAnalyzer.cc test/PerfectGeometryAnalyzer/src/PerfectGeometryAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Tommaso Boccali
//         Created:  Tue Jul 26 08:47:57 CEST 2005
// $Id: PerfectGeometryAnalyzer.cc,v 1.5 2005/11/25 15:11:12 case Exp $
//
//


// system include files
#include <memory>
#include <iostream>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "DetectorDescription/Core/interface/DDExpandedView.h"
#include "DetectorDescription/Core/interface/DDSpecifics.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "DataSvc/RefException.h"
#include "CoralBase/Exception.h"

//
// class decleration
//

class PerfectGeometryAnalyzer : public edm::EDAnalyzer {
public:
  explicit PerfectGeometryAnalyzer( const edm::ParameterSet& );
  ~PerfectGeometryAnalyzer();

  
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
private:
  // ----------member data ---------------------------
  std::string label_;
  bool dumpHistory_;
  bool dumpSpecs_;
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
PerfectGeometryAnalyzer::PerfectGeometryAnalyzer( const edm::ParameterSet& iConfig ) :
   label_(iConfig.getUntrackedParameter<std::string>("label",""))
{
  dumpHistory_=iConfig.getUntrackedParameter<bool>("dumpGeoHistory");
  dumpSpecs_=iConfig.getUntrackedParameter<bool>("dumpSpecs");
}


PerfectGeometryAnalyzer::~PerfectGeometryAnalyzer()
{
 
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PerfectGeometryAnalyzer::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
   using namespace edm;

   std::cout << "Here I am " << std::endl;
   //
   // get the DDCompactView
   //
   edm::ESHandle<DDCompactView> pDD;
   iSetup.get<IdealGeometryRecord>().get(label_, pDD );

   try {
      DDExpandedView epv(*pDD);
      std::cout << " without firstchild or next... epv.logicalPart() =" << epv.logicalPart() << std::endl;
      if ( dumpHistory_ ) {
	typedef DDExpandedView::nav_type nav_type;
	typedef std::map<nav_type,int> id_type;
	id_type idMap;
	int id=0;
	std::ofstream dump("dumpGeoHistoryOnRead");
	do {
	  nav_type pos = epv.navPos();
	  idMap[pos]=id;
	  dump << id << " - " << epv.geoHistory() << std::endl;
	  ++id;
	} while (epv.next());
	dump.close();
      }
      if ( dumpSpecs_ ) {
	DDSpecifics::iterator<DDSpecifics> spit(DDSpecifics::begin()), spend(DDSpecifics::end());
	// ======= For each DDSpecific...
	std::ofstream dump("dumpSpecsOnRead");
	for (; spit != spend; ++spit) {
	  if ( !spit->isDefined().second ) continue;  
	  const DDSpecifics & sp = *spit;
	  dump << sp << std::endl;
	}
	dump.close();
      }
      std::cout << "finished" << std::endl;
     
   }catch(const DDLogicalPart& iException){
      throw cms::Exception("Geometry")
	<<"DDORAReader::readDB caught a DDLogicalPart exception: \""<<iException<<"\"";
   } catch (const coral::Exception& e) {
      throw cms::Exception("Geometry")
	<<"DDORAReader::readDB caught coral::Exception: \""<<e.what()<<"\"";
   } catch( const pool::RefException& er){
      throw cms::Exception("Geometry")
	<<"DDORAReader::readDB caught pool::RefException: \""<<er.what()<<"\"";
   } catch ( pool::Exception& e ) {
     throw cms::Exception("Geometry")
       << "DDORAReader::readDB caught pool::Exception: \"" << e.what() << "\"";
   } catch ( std::exception& e ) {
     throw cms::Exception("Geometry")
       <<  "DDORAReader::readDB caught std::exception: \"" << e.what() << "\"";
   } catch ( ... ) {
     throw cms::Exception("Geometry")
       <<  "DDORAReader::readDB caught UNKNOWN!!! exception." << std::endl;
   }
}


//define this as a plug-in
DEFINE_FWK_MODULE(PerfectGeometryAnalyzer)
