#include "Geometry/CommonTopologies/test/ValidateRadial.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "boost/lexical_cast.hpp"
#include "TProfile.h"

ValidateRadial::ValidateRadial(const edm::ParameterSet& cfg) 
  : epsilon_(cfg.getParameter<double>("Epsilon")),
    file_(new TFile(cfg.getParameter<std::string>("FileName").c_str(),"RECREATE")),
    printOut_(cfg.getParameter<bool>("PrintOut"))
{}

void ValidateRadial::
analyze(const edm::Event& e, const edm::EventSetup& es) {
  std::vector<const RadialStripTopology*> topologies = get_list_of_radial_topologies(e,es);
  for(unsigned i=0; i<topologies.size(); i++) {
    test_topology(topologies[i],i);
  }
  file_->Close();
}

std::vector<const RadialStripTopology*> ValidateRadial::
get_list_of_radial_topologies(const edm::Event&e, const edm::EventSetup& es) {
  std::vector<const RadialStripTopology*> topos;
  edm::ESHandle<TrackerGeometry> theTrackerGeometry;  es.get<TrackerDigiGeometryRecord>().get( theTrackerGeometry );  
  const uint32_t radial_detids[] = { 402666125,//TID r1
				     402668833,//TID r2
				     402673476,//TID r3
				     470066725,//TEC r1
				     470390853,//TEC r2
				     470114664,//TEC r3
				     470131344,//TEC r4
				     470079661,//TEC r5
				     470049476,//TEC r6
				     470045428}; //TEC r7
  for(unsigned i=0; i<10; i++) 
    topos.push_back(
	dynamic_cast<const RadialStripTopology*>(&dynamic_cast<const StripGeomDetUnit*>( theTrackerGeometry->idToDet( radial_detids[i] ) )->specificTopology())
      );
  return topos;
}

void ValidateRadial::
test_topology(const RadialStripTopology* t, unsigned i) {
  TProfile prof(("se2limit1"+boost::lexical_cast<std::string>(i)).c_str(),
		"Precision Limit of recoverable strip error (1st order);strip;strip error",
		t->nstrips()/8,0,t->nstrips());
  TProfile prof2(("se2limit2"+boost::lexical_cast<std::string>(i)).c_str(),
		 "Precision Limit of recoverable strip error (2nd order);strip;strip error",
		 t->nstrips()/8,0,t->nstrips());
  for(float strip = 0; strip<t->nstrips(); strip+=0.5) {
    for(float stripErr2 = 0.03; stripErr2>1e-10; stripErr2/=1.1)
      if(!pass_frame_change_test( t, strip, stripErr2, false ) ) {
	prof.Fill(strip,sqrt(stripErr2));
	break;
      }
    for(float stripErr2 = 0.03; stripErr2>1e-10; stripErr2/=1.1)
      if(!pass_frame_change_test( t, strip, stripErr2, true ) ) {
	prof2.Fill(strip,sqrt(stripErr2));
	break;
      }
}
  prof.Write();
  prof2.Write();
}


bool ValidateRadial::
pass_frame_change_test(const RadialStripTopology* t, const float strip, const float stripErr2, const bool secondOrder) {
  const LocalPoint lp = t->localPosition(strip);
  const LocalError le = t->localError(strip,stripErr2);
  const MeasurementPoint mp = t->measurementPosition(lp);
  const MeasurementError me = t->measurementError(lp, le);
  const float newstrip = t->strip(lp);
  const LocalPoint newlp = t->localPosition(mp);
  const LocalError newle = t->localError(mp,me);
  const MeasurementPoint newmp = t->measurementPosition(newlp);
  const MeasurementError newme = t->measurementError(newlp, newle);

  const bool pass1( fabs(strip - newstrip) < 0.001&&
		    fabs(strip - mp.x()) < 0.001 &&
		    fabs(mp.x() - newstrip) <0.001 && 
		    fabs(me.uu()-stripErr2)/stripErr2 < epsilon_ &&
		    fabs(me.uv()) < epsilon_ &&
		    me.uu()>0 &&  me.vv()>0 );
  const bool pass2( fabs(strip - newmp.x()) < 0.001 &&
		    fabs(newme.uu()-stripErr2)/stripErr2 < epsilon_ &&
		    fabs(newme.uv()) < epsilon_ &&
		    newme.uu()>0 &&  newme.vv()>0 );
  const bool pass = (secondOrder? pass2 : pass1);

  if(printOut_ && !pass)
    std::cout << "(" << strip << ", " << newstrip << ", " << mp.x() << ")\t"
	      << "(" << stripErr2 << ", " << me.uu() << ")\t\t"
	      << ( me.uv() ) << std::endl;
  return pass;
}
