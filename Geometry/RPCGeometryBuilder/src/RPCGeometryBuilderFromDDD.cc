/** Implementation of the RPC Geometry Builder from DDD
 *
 *  \author Port of: MuDDDRPCBuilder (ORCA)
 *  \author M. Maggi - INFN Bari
 */
#include "Geometry/RPCGeometryBuilder/src/RPCGeometryBuilderFromDDD.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCRollSpecs.h"

#include <DetectorDescription/Core/interface/DDFilter.h>
#include <DetectorDescription/Core/interface/DDFilteredView.h>
#include <DetectorDescription/Core/interface/DDSolid.h>

#include "Geometry/MuonNumbering/interface/MuonDDDNumbering.h"
#include "Geometry/MuonNumbering/interface/MuonBaseNumber.h"
#include "Geometry/MuonNumbering/interface/RPCNumberingScheme.h"

#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"
#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"

#include "DataFormats/GeometryVector/interface/Basic3DVector.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include <iostream>
#include <algorithm>

RPCGeometryBuilderFromDDD::RPCGeometryBuilderFromDDD(bool comp11) : theComp11Flag(comp11)
{ }

RPCGeometryBuilderFromDDD::~RPCGeometryBuilderFromDDD() 
{ }

RPCGeometry* RPCGeometryBuilderFromDDD::build(const DDCompactView* cview, const MuonDDDConstants& muonConstants)
{
  std::string attribute = "ReadOutName"; // could come from .orcarc
  std::string value     = "MuonRPCHits";    // could come from .orcarc
  DDValue val(attribute, value, 0.0);

  // Asking only for the MuonRPC's
  DDSpecificsFilter filter;
  filter.setCriteria(val, // name & value of a variable 
		     DDSpecificsFilter::matches,
		     DDSpecificsFilter::AND, 
		     true, // compare strings otherwise doubles
		     true // use merged-specifics or simple-specifics
		     );
  DDFilteredView fview(*cview);
  fview.addFilter(filter);

  return this->buildGeometry(fview, muonConstants);
}

RPCGeometry* RPCGeometryBuilderFromDDD::buildGeometry(DDFilteredView& fview, const MuonDDDConstants& muonConstants)
{
#ifdef LOCAL_DEBUG  
  std::cout <<"Building the geometry service"<<std::endl;
#endif
  RPCGeometry* geometry = new RPCGeometry();

#ifdef LOCAL_DEBUG  
  std::cout << "About to run through the RPC structure" << std::endl;
  std::cout <<" First logical part "
  	    <<fview.logicalPart().name().name()<<std::endl;
#endif
  bool doSubDets = fview.firstChild();

#ifdef LOCAL_DEBUG  
  std::cout << "doSubDets = " << doSubDets << std::endl;
#endif
  while (doSubDets){

#ifdef LOCAL_DEBUG  
    std::cout <<"start the loop"<<std::endl; 
#endif

    // Get the Base Muon Number
    MuonDDDNumbering mdddnum(muonConstants);
#ifdef LOCAL_DEBUG  
    std::cout <<"Getting the Muon base Number"<<std::endl;
#endif
    MuonBaseNumber   mbn=mdddnum.geoHistoryToBaseNumber(fview.geoHistory());

#ifdef LOCAL_DEBUG  
    std::cout <<"Start the Rpc Numbering Schema"<<std::endl;
#endif
    // Get the The Rpc det Id 
    RPCNumberingScheme rpcnum(muonConstants);
    int detid = 0;

#ifdef LOCAL_DEBUG  
    std::cout <<"Getting the Unit Number"<<std::endl;
#endif
    detid = rpcnum.baseNumberToUnitNumber(mbn);
#ifdef LOCAL_DEBUG  
    std::cout <<"Getting the RPC det Id "<<detid <<std::endl;
#endif
    RPCDetId rpcid(detid);
    RPCDetId chid(rpcid.region(),rpcid.ring(),rpcid.station(),rpcid.sector(),rpcid.layer(),rpcid.subsector(),0);

#ifdef LOCAL_DEBUG  
    std::cout <<"The RPCDetid is "<<rpcid<<std::endl;
#endif

    DDValue numbOfStrips("nStrips");

    std::vector<const DDsvalues_type* > specs(fview.specifics());
    std::vector<const DDsvalues_type* >::iterator is=specs.begin();
    int nStrips=0;
    for (;is!=specs.end(); is++){
      if (DDfetch( *is, numbOfStrips)){
	nStrips=int(numbOfStrips.doubles()[0]);	
      }
    }
#ifdef LOCAL_DEBUG  
    if (nStrips == 0 )
      std::cout <<"No strip found!!"<<std::endl;
#endif
    
    std::vector<double> dpar=fview.logicalPart().solid().parameters();
    std::string name=fview.logicalPart().name().name();
    DDTranslation tran    = fview.translation();
    //removed .Inverse after comparing to DT...
    DDRotationMatrix rota = fview.rotation();//.Inverse();
    Surface::PositionType pos(tran.x()/cm,tran.y()/cm, tran.z()/cm);
    // CLHEP way
//     Surface::RotationType rot(rota.xx(),rota.xy(),rota.xz(),
// 			      rota.yx(),rota.yy(),rota.yz(),
// 			      rota.zx(),rota.zy(),rota.zz());

//ROOT::Math way
    DD3Vector x, y, z;
    rota.GetComponents(x,y,z);
    // doesn't this just re-inverse???
    Surface::RotationType rot (float(x.X()),float(x.Y()),float(x.Z()),
			       float(y.X()),float(y.Y()),float(y.Z()),
			       float(z.X()),float(z.Y()),float(z.Z())); 
    
    std::vector<float> pars;
    RPCRollSpecs* rollspecs= 0;
    Bounds* bounds = 0;



    if (dpar.size()==3){
      float width     = dpar[0]/cm;
      float length    = dpar[1]/cm;
      float thickness = dpar[2]/cm;
      //RectangularPlaneBounds* 
      bounds = 
	new RectangularPlaneBounds(width,length,thickness);
      pars.push_back(width);
      pars.push_back(length);
      pars.push_back(numbOfStrips.doubles()[0]); //h/2;

      if (!theComp11Flag) {
	//Correction of the orientation to get the REAL geometry.
        //Change of axes for the +z part only.
        //Including the 0 whell
        if (tran.z() >-1500. ){
          Basic3DVector<float> newX(-1.,0.,0.);
          Basic3DVector<float> newY(0.,-1.,0.);
          Basic3DVector<float> newZ(0.,0.,1.);
          rot.rotateAxes (newX, newY,newZ);
        }
      }
      
      rollspecs = new RPCRollSpecs(GeomDetEnumerators::RPCBarrel,name,pars);
#ifdef LOCAL_DEBUG  
      std::cout <<"Barrel "<<name
		<<" par "<<width
		<<" "<<length<<" "<<thickness;
#endif
    }else{
      float be = dpar[4]/cm;
      float te = dpar[8]/cm;
      float ap = dpar[0]/cm;
      float ti = 0.4/cm;
      //  TrapezoidalPlaneBounds* 
      bounds = 
	new TrapezoidalPlaneBounds(be,te,ap,ti);
      pars.push_back(dpar[4]/cm); //b/2;
      pars.push_back(dpar[8]/cm); //B/2;
      pars.push_back(dpar[0]/cm); //h/2;
      pars.push_back(numbOfStrips.doubles()[0]); //h/2;
      
#ifdef LOCAL_DEBUG  
      std::cout <<"Forward "<<name
		<<" par "<<dpar[4]/cm
		<<" "<<dpar[8]/cm<<" "<<dpar[3]/cm<<" "
		<<dpar[0];
#endif      

      rollspecs = new RPCRollSpecs(GeomDetEnumerators::RPCEndcap,name,pars);

      //Change of axes for the forward
      Basic3DVector<float> newX(1.,0.,0.);
      Basic3DVector<float> newY(0.,0.,1.);
      //      if (tran.z() > 0. )
      newY *= -1;
      Basic3DVector<float> newZ(0.,1.,0.);
      rot.rotateAxes (newX, newY,newZ);
      
    }
#ifdef LOCAL_DEBUG  
    std::cout <<"   Number of strips "<<nStrips<<std::endl;
#endif  

    BoundPlane::BoundPlanePointer surf = BoundPlane::build(pos, rot, bounds); 
    delete bounds; // bounds cloned by BoundPlane, so we can delete it

    RPCRoll* r=new RPCRoll(rpcid,surf,rollspecs);
    geometry->add(r);
    

    std::list<RPCRoll *> rls;
    if (chids.find(chid)!=chids.end()){
      rls = chids[chid];
    }
    rls.push_back(r);
    chids[chid]=rls;

    doSubDets = fview.nextSibling(); // go to next layer
  }
  // Create the RPCChambers and store them on the Geometry 
  for( std::map<RPCDetId, std::list<RPCRoll *> >::iterator ich=chids.begin();
       ich != chids.end(); ich++){
    RPCDetId chid = ich->first;
    std::list<RPCRoll * > rls = ich->second;

    // compute the overall boundplane. Distinguish between Barrel and Endcap

    std::vector<GlobalPoint> allP;
    float maxMajor=0;
    float minMinor=99999;

    RPCRoll fR=*(*rls.begin());

    for(std::list<RPCRoll *>::iterator rl=rls.begin();
	rl!=rls.end(); rl++){

      if ((*rl)->id().region() == 0){
	float x=(*rl)->surface().bounds().width()/2.;
	float y=(*rl)->surface().bounds().length()/2.;
	float z=(*rl)->surface().bounds().thickness()/2.;
	GlobalPoint gp1=(*rl)->toGlobal(LocalPoint(x,y,z));
	allP.push_back(gp1);
	GlobalPoint gp2=(*rl)->toGlobal(LocalPoint(-x,-y,-z));
	allP.push_back(gp2);
      }else{
	const TrapezoidalPlaneBounds  bTrap=*(static_cast<const TrapezoidalPlaneBounds *>(&(*rl)->surface().bounds()));
	std::vector<float> parsT=bTrap.parameters();
	if (parsT[0] < minMinor) minMinor=parsT[0];
	if (parsT[1] > maxMajor) maxMajor=parsT[1];
	float y = parsT[3];
	float z = parsT[2];
	GlobalPoint gp1=(*rl)->toGlobal(LocalPoint(0,y,z));
	allP.push_back(gp1);
	GlobalPoint gp2=(*rl)->toGlobal(LocalPoint(0,-y,-z));
	allP.push_back(gp2);
      }
   
    }

    Surface::PositionType pos = fR.position();     
    const Surface::RotationType rot = fR.rotation();
    const Bounds* bounds=0;

    if (chid.region()==0){
      double minX = 0; //it is surely a negative number;
      double maxX = 0;// it is surely a positiva number
      double minY = 0; //it is surely a negative number;
      double maxY = 0;// it is surely a positiva number
      double minZ = 0; //it is surely a negative number;
      double maxZ = 0;// it is surely a positiva number
      for (std::vector<GlobalPoint>::iterator p=allP.begin(); p<allP.end();p++){
	LocalPoint a = fR.toLocal(*p);
	if (a.x() < minX) minX=a.x();
	if (a.x() > maxX) maxX=a.x();
	if (a.y() < minY) minY=a.y();
	if (a.y() > maxY) maxY=a.y();
	if (a.z() < minZ) minZ=a.z();
	if (a.z() > maxZ) maxZ=a.z();
      }
      GlobalPoint pc = fR.toGlobal(LocalPoint((maxX+minX)/2.,(maxY+minY)/2.,(maxZ+minZ)/2.));
      pos = Surface::PositionType(pc);
      bounds = new RectangularPlaneBounds((maxX-minX)/2.,(maxY-minY)/2.,(maxZ-minZ)/2.);
    }else{
      double minY = 9999999; //it is surely a negative number;
      double maxY = 0;// it is surely a positiva number
      double minZ = 9999999; //it is surely a negative number;
      double maxZ = 0;// it is surely a positiva number
      for (std::vector<GlobalPoint>::iterator p=allP.begin(); p<allP.end();p++){
	LocalPoint a = fR.toLocal(*p);
	if (a.y() < minY) minY=a.y();
	if (a.y() > maxY) maxY=a.y();
	if (a.z() < minZ) minZ=a.z();
	if (a.z() > maxZ) maxZ=a.z();
      }
      GlobalPoint pc = fR.toGlobal(LocalPoint(0,(maxY+minY)/2.,(maxZ+minZ)/2.));
      pos = Surface::PositionType(pc);
      bounds = new TrapezoidalPlaneBounds(minMinor,maxMajor,(maxY-minY)/2.,(maxZ-minZ)/2.);
    }
    
    
    BoundPlane::BoundPlanePointer surf = BoundPlane::build(pos, rot, bounds); 
    RPCChamber* ch = new RPCChamber (chid, surf);     
 
    // Add the rolls to rhe chamber
    for(std::list<RPCRoll *>::iterator rl=rls.begin();
    rl!=rls.end(); rl++){
      ch->add(*rl);
    }
    // Add the chamber to the geometry
    geometry->add(ch);
    delete bounds; 
 } 
  return geometry;
}
