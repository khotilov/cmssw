/// clhep
#include "CLHEP/Units/SystemOfUnits.h"

//#include "DetectorDescription/Core/interface/DDExpandedView.h"
//temporary
//#include "DetectorDescription/Core/interface/DDSolid.h"
//
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeomBuilderFromGeometricDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/GeomTopologyBuilder.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PlaneBuilderFromGeometricDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/PlaneBuilderForGluedDet.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/GeometrySurface/interface/MediumProperties.h"


#include <cfloat>
#include <vector>

using std::vector;
using std::string;


TrackerGeometry* TrackerGeomBuilderFromGeometricDet::build( const GeometricDet* gd){

  TrackerGeometry* tracker = new TrackerGeometry(gd);
  std::vector<const GeometricDet*> comp;
  gd->deepComponents(comp);

  std::vector<const GeometricDet*> dets[6];
  std::vector<const GeometricDet*> & pixB = dets[0]; pixB.reserve(comp.size());
  std::vector<const GeometricDet*> & pixF = dets[1]; pixF.reserve(comp.size());
  std::vector<const GeometricDet*> & tib  = dets[2];  tib.reserve(comp.size());
  std::vector<const GeometricDet*> & tid  = dets[3];  tid.reserve(comp.size());
  std::vector<const GeometricDet*> & tob  = dets[4];  tob.reserve(comp.size());
  std::vector<const GeometricDet*> & tec  = dets[5];  tec.reserve(comp.size());

  for(u_int32_t i = 0;i<comp.size();i++)
    dets[comp[i]->geographicalID().subdetId()-1].push_back(comp[i]);
  

  buildPixel(pixB,tracker,theDetIdToEnum.type(1), "barrel"); //"PixelBarrel" 
  buildPixel(pixF,tracker,theDetIdToEnum.type(2), "endcap"); //"PixelEndcap" 
  buildSilicon(tib,tracker,theDetIdToEnum.type(3), "barrel");// "TIB"	
  buildSilicon(tid,tracker,theDetIdToEnum.type(4), "endcap");//"TID" 
  buildSilicon(tob,tracker,theDetIdToEnum.type(5), "barrel");//"TOB"	
  buildSilicon(tec,tracker,theDetIdToEnum.type(6), "endcap");//"TEC"        
  buildGeomDet(tracker);//"GeomDet"
  return tracker;
}

void TrackerGeomBuilderFromGeometricDet::buildPixel(std::vector<const GeometricDet*>  const & gdv, 
						    TrackerGeometry* tracker,
						    GeomDetType::SubDetector& det,
						    const std::string& part){ 

  static std::map<std::string,PixelGeomDetType*> detTypeMap;

  for(u_int32_t i=0; i<gdv.size(); i++){

    std::string const & detName = gdv[i]->name();
    if (detTypeMap.find(detName) == detTypeMap.end()) {

      PixelTopology* t = 
	theTopologyBuilder->buildPixel(gdv[i]->bounds(),
				       gdv[i]->pixROCRows(),
				       gdv[i]->pixROCCols(),
				       gdv[i]->pixROCx(),
				       gdv[i]->pixROCy(),
				       part);
      
      detTypeMap[detName] = new PixelGeomDetType(t,detName,det);
      tracker->addType(detTypeMap[detName]);
    }

    PlaneBuilderFromGeometricDet::ResultType plane = buildPlaneWithMaterial(gdv[i]);
    GeomDetUnit* temp =  new PixelGeomDetUnit(&(*plane),detTypeMap[detName],gdv[i]);

    tracker->addDetUnit(temp);
    tracker->addDetUnitId(gdv[i]->geographicalID());
  }
}

void TrackerGeomBuilderFromGeometricDet::buildSilicon(std::vector<const GeometricDet*>  const & gdv, 
						      TrackerGeometry* tracker,
						      GeomDetType::SubDetector& det,
						      const std::string& part)
{ 
  static std::map<std::string,StripGeomDetType*> detTypeMap;
  
  for(u_int32_t i=0;i<gdv.size();i++){

    std::string const & detName = gdv[i]->name();
    if (detTypeMap.find(detName) == detTypeMap.end()) {

       StripTopology* t =
	theTopologyBuilder->buildStrip(gdv[i]->bounds(),
				       gdv[i]->siliconAPVNum(),
				       part);
      detTypeMap[detName] = new  StripGeomDetType( t,detName,det,
						   gdv[i]->stereo());
      tracker->addType(detTypeMap[detName]);
    }
    
    PlaneBuilderFromGeometricDet::ResultType plane = buildPlaneWithMaterial(gdv[i]);  
    GeomDetUnit* temp = new StripGeomDetUnit(&(*plane), detTypeMap[detName],gdv[i]);
    
    tracker->addDetUnit(temp);
    tracker->addDetUnitId(gdv[i]->geographicalID());
  }  
}


void TrackerGeomBuilderFromGeometricDet::buildGeomDet(TrackerGeometry* tracker){
  PlaneBuilderForGluedDet gluedplaneBuilder;
  std::vector<GeomDetUnit*> const & gdu= tracker->detUnits();
  std::vector<DetId> const & gduId = tracker->detUnitIds();

  for(u_int32_t i=0;i<gdu.size();i++){
    StripSubdetector sidet( gduId[i].rawId());
    tracker->addDet((GeomDet*) gdu[i]);
    tracker->addDetId(gduId[i]);      
    if(sidet.glued()!=0&&sidet.stereo()==1){
      int partner_pos=-1;
      for(u_int32_t jj=0;jj<gduId.size();jj++){
	if(sidet.partnerDetId()== gduId[jj]) {
	  partner_pos=jj;
	  break;
	}
      }
      const GeomDetUnit* dus = gdu[i];
      if(partner_pos==-1){
	throw cms::Exception("Configuration") <<"No partner detector found \n"
					<<"There is a problem on Tracker geometry configuration\n";
      }
      const GeomDetUnit* dum = gdu[partner_pos];
      std::vector<const GeomDetUnit *> glued(2);
      glued[0]=dum;
      glued[1]=dus;
      PlaneBuilderForGluedDet::ResultType plane = gluedplaneBuilder.plane(glued);
      GluedGeomDet* gluedDet = new GluedGeomDet(&(*plane),dum,dus);
      tracker->addDet((GeomDet*) gluedDet);
      tracker->addDetId(DetId(sidet.glued()));
    }
  }
}


// std::string TrackerGeomBuilderFromGeometricDet::getString(const std::string & s, DDExpandedView* ev) const
// {
//     DDValue val(s);
//     vector<const DDsvalues_type *> result;
//     ev->specificsV(result);
//     vector<const DDsvalues_type *>::iterator it = result.begin();
//     bool foundIt = false;
//     for (; it != result.end(); ++it)
//     {
// 	foundIt = DDfetch(*it,val);
// 	if (foundIt) break;

//     }    
//     if (foundIt)
//     { 
// 	const std::vector<std::string> & temp = val.strings(); 
// 	if (temp.size() != 1)
// 	{
// 	  throw cms::Exception("Configuration") << "I need 1 "<< s << " tags";
// 	}
// 	return temp[0]; 
//     }
//     return "NotFound";
// }

// double TrackerGeomBuilderFromGeometricDet::getDouble(const std::string & s,  DDExpandedView* ev) const 
// {
//   DDValue val(s);
//   vector<const DDsvalues_type *> result;
//   ev->specificsV(result);
//   vector<const DDsvalues_type *>::iterator it = result.begin();
//   bool foundIt = false;
//   for (; it != result.end(); ++it)
//     {
//       foundIt = DDfetch(*it,val);
//       if (foundIt) break;
//     }    
//   if (foundIt)
//     { 
//       const std::vector<std::string> & temp = val.strings(); 
//       if (temp.size() != 1)
// 	{
// 	  throw cms::Exception("Configuration") << "I need 1 "<< s << " tags";
// 	}
//       return double(atof(temp[0].c_str())); 
//     }
//   return 0;
// }

PlaneBuilderFromGeometricDet::ResultType
TrackerGeomBuilderFromGeometricDet::buildPlaneWithMaterial(const GeometricDet* gd) const
{
  PlaneBuilderFromGeometricDet planeBuilder;
  PlaneBuilderFromGeometricDet::ResultType plane = planeBuilder.plane(gd);  
  //
  // set medium properties (if defined)
  //
  plane->setMediumProperties(MediumProperties(gd->radLength(),gd->xi()));

  return plane;
}
