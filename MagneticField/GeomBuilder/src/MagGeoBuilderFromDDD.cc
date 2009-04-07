// #include "Utilities/Configuration/interface/Architecture.h"

/*
 *  See header file for a description of this class.
 *
 *  $Date: 2009/03/25 16:44:20 $
 *  $Revision: 1.20 $
 *  \author N. Amapane - INFN Torino
 */

#include "MagneticField/GeomBuilder/src/MagGeoBuilderFromDDD.h"
#include "MagneticField/GeomBuilder/src/volumeHandle.h"
#include "MagneticField/GeomBuilder/src/bSlab.h"
#include "MagneticField/GeomBuilder/src/bRod.h"
#include "MagneticField/GeomBuilder/src/bSector.h"
#include "MagneticField/GeomBuilder/src/bLayer.h"
#include "MagneticField/GeomBuilder/src/eSector.h"
#include "MagneticField/GeomBuilder/src/eLayer.h"
#include "MagneticField/GeomBuilder/src/FakeInterpolator.h"

#include "MagneticField/Layers/interface/MagBLayer.h"
#include "MagneticField/Layers/interface/MagESector.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDFilter.h"

#include "Utilities/BinningTools/interface/ClusterizingHistogram.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "MagneticField/Interpolation/interface/MagProviderInterpol.h"
#include "MagneticField/Interpolation/interface/MFGridFactory.h"
#include "MagneticField/Interpolation/interface/MFGrid3D.h"

#include "MagneticField/VolumeGeometry/interface/MagVolume6Faces.h"
#include "MagneticField/VolumeGeometry/interface/MagExceptions.h"
#include "MagneticField/Layers/interface/MagVerbosity.h"

#include "DataFormats/GeometryVector/interface/Pi.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <map>
#include <set>
#include "Utilities/General/interface/precomputed_value_sort.h"


bool MagGeoBuilderFromDDD::debug;

using namespace std;

MagGeoBuilderFromDDD::MagGeoBuilderFromDDD(string version_, bool debug_) :
  version (version_)
{  
  debug = debug_;
  if (debug) cout << "Constructing a MagGeoBuilderFromDDD" <<endl;
}

MagGeoBuilderFromDDD::~MagGeoBuilderFromDDD(){
  for (handles::const_iterator i=bVolumes.begin();
       i!=bVolumes.end(); ++i){
    delete (*i);
  }
  
  for (handles::const_iterator i=eVolumes.begin();
       i!=eVolumes.end(); ++i){
    delete (*i);
  }
}


void MagGeoBuilderFromDDD::summary(handles & volumes){  
  // The final countdown.
  int ivolumes  = volumes.size();  // number of volumes
  int isurfaces = ivolumes*6;       // number of individual surfaces
  int iassigned = 0;                // How many have been assigned
  int iunique   = 0;                // number of unique surfaces
  int iref_ass  = 0;
  int iref_nass = 0;

  set<const void *> ptrs;

  handles::const_iterator first = volumes.begin();
  handles::const_iterator last = volumes.end();

  for (handles::const_iterator i=first; i!=last; ++i){
    if (int((*i)->shape())>4) continue; // FIXME: missing shapes...
    for (int side = 0; side < 6; ++side) {
      int references = 	(*i)->references(side);
      if ((*i)->isPlaneMatched(side)) {
	++iassigned;
	bool firstOcc = (ptrs.insert(&((*i)->surface(side)))).second;
	if (firstOcc) iref_ass+=references;
	if (references<2){  
	  cout << "*** Only 1 ref, vol: " << (*i)->name << " # "
	       << (*i)->copyno << " side: " << side << endl;
	}	
      } else {
	iref_nass+=references;
	if (references>1){
	  cout << "*** Ref_nass >1 " <<endl;
	}
      }
    }
  }
  iunique = ptrs.size();

  cout << "    volumes   " << ivolumes  << endl
       << "    surfaces  " << isurfaces << endl
       << "    assigned  " << iassigned << endl
       << "    unique    " << iunique << endl
       << "    iref_ass  " << iref_ass << endl
       << "    iref_nass " << iref_nass << endl;
}


void MagGeoBuilderFromDDD::build(const DDCompactView & cpva)
{
//    DDCompactView cpv;
  DDExpandedView fv(cpva);

  if (debug) cout << "**********************************************************" <<endl;

  // The actual field interpolators
  map<string, MagProviderInterpol*> bInterpolators;
  map<string, MagProviderInterpol*> eInterpolators;
  
  // Counter of different (FIXME can be removed)
  int bVolCount = 0;
  int eVolCount = 0;


  // Look for MAGF tree (any better way to find out???)
  //  fv.reset();

  if (fv.logicalPart().name().name()!="MAGF") {
     std::string topNodeName(fv.logicalPart().name().name());

     //see if one of the children is MAGF
     bool doSubDets = fv.firstChild();
     
     bool go=true;
     while(go&& doSubDets) {
	if (fv.logicalPart().name().name()=="MAGF")
	   break;
	else
	   go = fv.nextSibling();
     }
     if (!go) {
	throw cms::Exception("NoMAGFinDDD")<<" Neither he top node, nor any child node of the DDCompactView is \"MAGF\" but the top node is instead \""<<topNodeName<<"\"";
     }
  }
  // Loop over MAGF volumes and create volumeHandles. 
  if (debug) { cout << endl << "*** MAGF: " << fv.geoHistory() << endl
		    << "translation: " << fv.translation() << endl
		    << " rotation: " << fv.rotation() << endl;
  }
  
  bool doSubDets = fv.firstChild();
  while (doSubDets){
    
    string name = fv.logicalPart().name().name();
    if (debug) cout << endl << "Name: " << name << endl
			       << "      " << fv.geoHistory() <<endl;

    // Build only the z-negative volumes, assuming symmetry
    // FIXME: should not use name but center...
    // even better, we should fix the XML!
    if (name.substr(2,2)=="ZP") {
      doSubDets = fv.nextSibling();
      continue;
    }
    
    bool mergeCylinders=true;

    // In the barrel, cylinders sectors will be skipped to build full 
    // cylinders out of sector copyno #1.
    // (these should be just volumes 1,2,4)
    bool expand = false;
    if (mergeCylinders) {
      if (name == "V_ZN_1"
	  || name == "V_ZN_2") {
	if (debug && fv.logicalPart().solid().shape()!=ddtubs) {
	  cout << "ERROR: MagGeoBuilderFromDDD::build: volume " << name
	       << " should be a cylinder" << endl;
	}
	if(fv.copyno()==1) {
	  // FIXME expand = true;
	} else {
	  //cout << "... to be skipped: "
	  //     << name << " " << fv.copyno() << endl;
	  //FIXME continue;
	}
      }
    }

    volumeHandle* v = new volumeHandle(fv, expand);
    
    // Select volumes, build volume handles.
    float Z = v->center().z();
    float R = v->center().perp();

    // v 85l: Barrel is everything up to |Z| = 661.0, excluding 
    // volume #7, centered at 6477.5
    // v 1103l: same numbers work fine. #16 instead of #7, same coords;
    // see comment below for V6,7
    // FIXME: misalignment?
    if ((fabs(Z)<647. || (R>350. && fabs(Z)<662.)) &&
	!(fabs(Z)>480 && R<172) // in 1103l we place V_6 and V_7 in the 
	                        // endcaps to preserve nice layer structure
	                        // in the barrel. This does not hurt in v85l
	                        // where there is a single V1 
	) { // Barrel
      if (debug) cout << " (Barrel)" <<endl;
      bVolumes.push_back(v);
      // Build the interpolator of the "master" volume (the one which is
      // not replicated, i.e. copy number #1)
      if (v->copyno==1) {
	buildInterpolator(v, bInterpolators);
	++bVolCount;
      }
    } else {               // Endcaps
      if (debug) cout << " (Endcaps)" <<endl;
      eVolumes.push_back(v);
      if (v->copyno==1) { 
	buildInterpolator(v, eInterpolators);
	++eVolCount;
      }
    }

    doSubDets = fv.nextSibling(); // end of loop over MAGF
  }
    
  if (debug) {
    cout << "Number of volumes (barrel): " << bVolumes.size() <<endl
		  << "Number of volumes (endcap): " << eVolumes.size() <<endl;
    cout << "**********************************************************" <<endl;
  }

  // Now all volumeHandles are there, and parameters for each of the planes
  // are calculated.

  //----------------------------------------------------------------------
  // Print summary information

  if (debug) {
    cout << "-----------------------" << endl;
    cout << "SUMMARY: Barrel " << endl;
    summary(bVolumes);
    
    cout << endl << "SUMMARY: Endcaps " << endl;
    summary(eVolumes);
    cout << "-----------------------" << endl;
  }


  //----------------------------------------------------------------------
  // Find barrel layers.

  vector<bLayer> layers; // the barrel layers
  precomputed_value_sort(bVolumes.begin(), bVolumes.end(), ExtractRN());

  // Find the layers (in R)
  const float resolution = 1.; // cm
  float rmin = bVolumes.front()->RN()-resolution;
  float rmax = bVolumes.back()->RN()+resolution;
  ClusterizingHistogram  hisR( int((rmax-rmin)/resolution) + 1, rmin, rmax);

  if (debug) cout << " R layers: " << rmin << " " << rmax << endl;

  handles::const_iterator first = bVolumes.begin();
  handles::const_iterator last = bVolumes.end();  

  for (handles::const_iterator i=first; i!=last; ++i){
    hisR.fill((*i)->RN());
  }
  vector<float> rClust = hisR.clusterize(resolution);

  handles::const_iterator ringStart = first;
  handles::const_iterator separ = first;

  for (unsigned int i=0; i<rClust.size() - 1; ++i) {
    if (debug) cout << " Layer at RN = " << rClust[i];
    float rSepar = (rClust[i] + rClust[i+1])/2.f;
    while ((*separ)->RN() < rSepar) ++separ;

    bLayer thislayer(ringStart, separ);
    layers.push_back(thislayer);
    ringStart = separ;
  }
  {
    if (debug) cout << " Layer at RN = " << rClust.back();
    bLayer thislayer(separ, last);
    layers.push_back(thislayer);
  }

  if (debug) cout << "Barrel: Found " << rClust.size() << " clusters in R, "
		  << layers.size() << " layers " << endl << endl;


  //----------------------------------------------------------------------
  // Find endcap sectors

  vector<eSector> sectors; // the endcap sectors
  precomputed_value_sort(eVolumes.begin(), eVolumes.end(), ExtractPhi()); 
 

  // ASSUMPTION: There are 12 sectors and each sector is 30 deg wide.
  for (int i = 0; i<12; ++i) {
    int offset = eVolumes.size()/12;
    //    int isec = (i+binOffset)%12;
    if (debug) cout << " Sector at phi = "
		    << (*(eVolumes.begin()+((i)*offset)))->center().phi()
		    << endl;
    sectors.push_back(eSector(eVolumes.begin()+((i)*offset),
			      eVolumes.begin()+((i+1)*offset)));
  }
   
  if (debug) cout << "Endcap: Found " 
		  << sectors.size() << " sectors " << endl;


  //----------------------------------------------------------------------  
  // Match surfaces.

//  cout << "------------------" << endl << "Now associating planes..." << endl;

//   // Loop on layers
//   for (vector<bLayer>::const_iterator ilay = layers.begin();
//        ilay!= layers.end(); ++ilay) {
//     cout << "On Layer: " << ilay-layers.begin() << " RN: " << (*ilay).RN()
// 	 <<endl;     

//     // Loop on wheels
//     for (vector<bWheel>::const_iterator iwheel = (*ilay).wheels.begin();
// 	 iwheel != (*ilay).wheels.end(); ++iwheel) {
//       cout << "  On Wheel: " << iwheel- (*ilay).wheels.begin()<< " Z: "
// 	   << (*iwheel).minZ() << " " << (*iwheel).maxZ() << " " 
// 	   << ((*iwheel).minZ()+(*iwheel).maxZ())/2. <<endl;

//       // Loop on sectors.
//       for (int isector = 0; isector<12; ++isector) {
// 	// FIXME: create new constructor...
// 	bSectorNavigator navy(layers,
// 			      ilay-layers.begin(),
// 			      iwheel-(*ilay).wheels.begin(),isector);
	
// 	const bSector & isect = (*iwheel).sector(isector);
	
// 	isect.matchPlanes(navy); //FIXME refcount
//       }
//     }
//   }


  //----------------------------------------------------------------------
  // Build MagVolumes and the MagGeometry hierarchy.

  //--- Barrel

  // Build MagVolumes and associate interpolators to them
  buildMagVolumes(bVolumes, bInterpolators);

  // Build MagBLayers
  for (vector<bLayer>::const_iterator ilay = layers.begin();
       ilay!= layers.end(); ++ilay) {
    mBLayers.push_back((*ilay).buildMagBLayer());
  }

  if (debug) {  
    cout << "*** BARREL ********************************************" << endl
	 << "Number of different volumes   = " << bVolCount << endl
	 << "Number of interpolators built = " << bInterpolators.size() << endl
    	 << "Number of MagBLayers built    = " << mBLayers.size() << endl;

    testInside(bVolumes); // Fixme: all volumes should be checked in one go.
  }
  
  //--- Endcap
  // Build MagVolumes  and associate interpolators to them
  buildMagVolumes(eVolumes, eInterpolators);

  // Build the MagESectors
  for (vector<eSector>::const_iterator isec = sectors.begin();
       isec!= sectors.end(); ++isec) {
    mESectors.push_back((*isec).buildMagESector());
  }

  if (debug) {
    cout << "*** ENDCAP ********************************************" << endl
	 << "Number of different volumes   = " << eVolCount << endl
	 << "Number of interpolators built = " << eInterpolators.size() << endl
    	 << "Number of MagESector built    = " << mESectors.size() << endl;

    testInside(eVolumes); // Fixme: all volumes should be checked in one go.
  }
}


void MagGeoBuilderFromDDD::buildMagVolumes(const handles & volumes, map<string, MagProviderInterpol*> & interpolators) {
  // Build all MagVolumes setting the MagProviderInterpol
  for (handles::const_iterator vol=volumes.begin(); vol!=volumes.end();
       ++vol){
    const MagProviderInterpol* mp = 0;
    if (interpolators.find((*vol)->magFile)!=interpolators.end()) {
      mp = interpolators[(*vol)->magFile];
    } else {
      cout << "No interpolator found for file " << (*vol)->magFile
	   << " vol: " << (*vol)->name << endl;
      cout << interpolators.size() <<endl;
      //FIXME08      continue;
    }
      

    // Get the volume number from the volume name.
    int volNum;
    string name = (*vol)->name;
    name.erase(0,name.rfind('_')+1);
    stringstream str;    
    str << name;
    str >> volNum;


    int key = volNum*100+(*vol)->copyno; 
    // FIXME: is copyno == sector??? should be the case for V1103l; but not
    // for 85l

    map<int, double>::const_iterator isf = theScalingFactors.find(key);
    if (isf == theScalingFactors.end()) {
      key = volNum*100;
      isf = theScalingFactors.find(key);
    }
    
    double sf = 1.;
    if (isf != theScalingFactors.end()) {
      sf = (*isf).second;

      //FIXME
      edm::LogInfo("MagneticField|VolumeBasedMagneticFieldESProducer") << "Applying scaling factor " << sf << " to "<< (*vol)->name << " (key:" << key << ")" << endl;
    }

    const GloballyPositioned<float> * gpos = (*vol)->placement();
    // FIXME check pos, rot corrsponds
    (*vol)->magVolume = new MagVolume6Faces(gpos->position(),
					    gpos->rotation(),
					    (*vol)->shape(),
					    (*vol)->sides(),
					    mp, sf);

    // FIXME temporary hack 
    (*vol)->magVolume->setIsIron((*vol)->isIron());

    // FIXME: debug, to be removed
    (*vol)->magVolume->name = (*vol)->name;
    (*vol)->magVolume->copyno = (*vol)->copyno;
  }
}


void MagGeoBuilderFromDDD::buildInterpolator(const volumeHandle * vol, map<string, MagProviderInterpol*> & interpolators){


  // In version grid_85l_030919, interpolators should be built only 
  // for volumes on NEGATIVE z 
  // (Z symmetry in field tables)
  if (version=="grid_85l_030919" && vol->center().z()>0) return;

  if (debug) {
    cout << "Building interpolator from "
	 << vol->name << " copyno " << vol->copyno
	 << " at " << vol->center()
	 << " phi: " << vol->center().phi()
	 << " file: " << vol->magFile
	 << endl;

    // In ver. grid_85l_030919, the master sector is sector 4 (along Y axis).
    // In ver. grid_1103l_071212, it is sector 1 (along X axis)
    double masterSectorPhi=0.;
    if (version=="grid_85l_030919") {
      masterSectorPhi=Geom::pi()/2.;
    }
    if ( fabs(vol->center().phi() - masterSectorPhi) > Geom::pi()/9.) {
      cout << "***WARNING wrong sector? " << endl;
    }
  }

  if (version == "fake") {
    interpolators[vol->magFile] = new magneticfield::FakeInterpolator();
    return;
  }

  string fullPath;

  try {
    edm::FileInPath mydata("MagneticField/Interpolation/data/"+version+"/"+vol->magFile);
    fullPath = mydata.fullPath();
  } catch (edm::Exception& exc) {
    cerr << "MagGeoBuilderFromDDD: exception in reading table; " << exc.what() << endl;
    //FIXME08    
    return;
    //throw;
  }
  
  
  try{
    if (vol->toExpand()){
      //FIXME
//       interpolators[vol->magFile] =
// 	MFGridFactory::build( fullPath, *(vol->placement()), vol->minPhi(), vol->maxPhi());
    } else {
      interpolators[vol->magFile] =
	MFGridFactory::build( fullPath, *(vol->placement()));
    }
  } catch (MagException& exc) {
    cout << exc.what() << endl;
    interpolators.erase(vol->magFile);
    throw; // FIXME, comment this to debug geometry in absence of interpolators
    //return;
  }


    if (debug) {
    // Check that all grid points of the interpolator are inside the volume.
      const MagVolume6Faces tempVolume(vol->placement()->position(),
				 vol->placement()->rotation(),
				 vol->shape(),
				 vol->sides(), 
				 interpolators[vol->magFile]);

      const MFGrid3D* grid = dynamic_cast<const MFGrid3D*>(interpolators[vol->magFile]);
      if (grid!=0) {
	
	vector<int> sizes = grid->dimensions();
	cout << "Grid has " << sizes.size() << " dimensions " 
	     << " number of nodes is " << sizes[0] << " " << sizes[1]
	     << " " << sizes[2] << endl;
      
	const double tolerance = 0.03;


	int dumpCount = 0;
	for (int j=0; j < sizes[1]; j++) {
	  for (int k=0; k < sizes[2]; k++) {
	    for (int i=0; i < sizes[0]; i++) {
	      MFGrid::LocalPoint lp = grid->nodePosition( i, j, k);
	      if (! tempVolume.inside(lp, tolerance)) {
		if (++dumpCount < 2) {
		  MFGrid::GlobalPoint gp = tempVolume.toGlobal(lp);
		  cout << "GRID ERROR: " << i << " " << j << " " << k
		       << " local: " << lp
		       << " global: " << gp
		       << " R= " << gp.perp() << " phi=" << gp.phi() << endl;
		}
	      }
	    }
	  }
	}
    
	cout << vol->name << " : Number of grid points outside the MagVolume: " << dumpCount << "/" << sizes[0]*sizes[1]*sizes[2] << endl;
      }
    }
}



void MagGeoBuilderFromDDD::testInside(handles & volumes) {
  // test inside() for all volumes.
  cout << "--------------------------------------------------" << endl;
  cout << " inside(center) test" << endl;
  for (handles::const_iterator vol=volumes.begin(); vol!=volumes.end();
       ++vol){
    for (handles::const_iterator i=volumes.begin(); i!=volumes.end();
	 ++i){
      if ((*i)==(*vol)) continue;
      //if ((*i)->magVolume == 0) continue;
      if ((*i)->magVolume->inside((*vol)->center())) {
	cout << "*** ERROR: center of " << (*vol)->name << " is inside " 
	     << (*i)->name <<endl;
      }
    }
    
    if ((*vol)->magVolume->inside((*vol)->center())) {
      cout << (*vol)->name << " OK " << endl;
    } else {
      cout << "*** ERROR: center of volume is not inside it, "
	   << (*vol)->name << endl;
    }
  }
  cout << "--------------------------------------------------" << endl;
}


vector<MagBLayer*> MagGeoBuilderFromDDD::barrelLayers() const{
  return mBLayers;
}

vector<MagESector*> MagGeoBuilderFromDDD::endcapSectors() const{
  return mESectors;
}

vector<MagVolume6Faces*> MagGeoBuilderFromDDD::barrelVolumes() const{
  vector<MagVolume6Faces*> v;
  v.reserve(bVolumes.size());
  for (handles::const_iterator i=bVolumes.begin();
       i!=bVolumes.end(); ++i){
    v.push_back((*i)->magVolume);
  }
  return v;
}

vector<MagVolume6Faces*> MagGeoBuilderFromDDD::endcapVolumes() const{
  vector<MagVolume6Faces*> v;
  v.reserve(eVolumes.size());
  for (handles::const_iterator i=eVolumes.begin();
       i!=eVolumes.end(); ++i){
    v.push_back((*i)->magVolume);
  }
  return v;
}


float MagGeoBuilderFromDDD::maxR() const{
  //FIXME: get from the actual geometry!!!
  if (version=="grid_85l_030919") return 1000.;
  else return 900.;
}

float MagGeoBuilderFromDDD::maxZ() const{
  //FIXME: get from the actual geometry!!!
  return 1600.;
}


void MagGeoBuilderFromDDD::setScaling(std::vector<int> keys, 
				      std::vector<double> values)
{
  if (keys.size() != values.size()) {
    throw cms::Exception("InvalidParameter") << "Invalid field scaling parameters 'scalingVolumes' and 'scalingFactors' ";
  }
  for (unsigned int i=0; i<keys.size(); ++i) {
    theScalingFactors[keys[i]] = values[i];
  }
}



