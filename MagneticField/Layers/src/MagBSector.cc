// #include "Utilities/Configuration/interface/Architecture.h"

/*
 *  See header file for a description of this class.
 *
 *  $Date: 2004/06/22 17:05:14 $
 *  $Revision: 1.2 $
 *  \author N. Amapane - INFN Torino
 */

#include "MagneticField/Layers/interface/MagBSector.h"
#include "MagneticField/Layers/interface/MagBRod.h"

// #include "MagneticField/MagLayers/interface/MagVerbosity.h"

#include <iostream>

using namespace std;

MagBSector::MagBSector(vector<MagBRod*>& rods, Geom::Phi<float> phiMin) : 
  theRods(rods),
  thePhiMin(phiMin)
{};

MagBSector::~MagBSector(){
  for (vector<MagBRod *>::const_iterator irod = theRods.begin();
       irod != theRods.end(); irod++) {
    delete (*irod);
  }
}

MagVolume * MagBSector::findVolume(const GlobalPoint & gp, double tolerance) const {
  MagVolume * result = 0;
  //  Geom::Phi<float> phi = gp.phi();

  // FIXME : use a binfinder
  for ( vector<MagBRod*>::const_iterator irod = theRods.begin();
	irod != theRods.end(); irod++) {
    // TOFIX
//     if (verbose.debugOut) cout << "     Trying rod at phi " << (*irod)->minPhi()
// 		    << " " << phi << endl ;
    result = (*irod)->findVolume(gp, tolerance);
    if (result!=0) return result;
  }

  return 0;
}



