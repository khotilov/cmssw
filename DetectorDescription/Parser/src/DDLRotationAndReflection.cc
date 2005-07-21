/***************************************************************************
                          DDLRotationAndReflection.cc  -  description
                             -------------------
    begin                : Tue Aug 6 2002
    email                : case@ucdhep.ucdavis.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *           DDDParser sub-component of DDD                                *
 *                                                                         *
 ***************************************************************************/

namespace std{} using namespace std;

// -------------------------------------------------------------------------
// Includes
// -------------------------------------------------------------------------
#include "DetectorDescription/Parser/interface/DDLRotationAndReflection.h"
#include "DetectorDescription/Parser/interface/DDLElementRegistry.h"

// DDCore dependencies
#include "DetectorDescription/Core/interface/DDPosPart.h"
#include "DetectorDescription/Core/interface/DDName.h"
#include "DetectorDescription/Core/interface/DDTransform.h"
#include "DetectorDescription/Base/interface/DDdebug.h"
#include "DetectorDescription/Base/interface/DDException.h"

#include "DetectorDescription/ExprAlgo/interface/ExprEvalSingleton.h"

// CLHEP dependencies
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <string>
#include <cmath>

// Default constructor
DDLRotationAndReflection::DDLRotationAndReflection() 
{
}

// Default destructor
DDLRotationAndReflection::~DDLRotationAndReflection()
{
}

void DDLRotationAndReflection::processElement (const string& name, const string& nmspace)
{

  DCOUT_V('P', "DDLRotationAndReflection::processElement started " << name);

  Hep3Vector x = makeX(nmspace);
  Hep3Vector y = makeY(nmspace);
  Hep3Vector z = makeZ(nmspace);

  DDXMLAttribute atts = getAttributeSet();

  try {
   if ((name == "Rotation") && isLeftHanded(x, y, z, nmspace) == 0)
     {
       HepRotation R;
       R.rotateAxes(x, y, z);      
       DDRotationMatrix* ddr = new DDRotationMatrix(R);
       DDRotation ddrot = DDrot(getDDName(nmspace), ddr);
       DCOUT_V ('p', "Rotation created: " << ddrot << endl);
     }
   else if ((name == "Rotation")  && isLeftHanded(x, y, z, nmspace) == 1)
     {
       string msg("\nDDLRotationAndReflection attempted to make a");
       msg += " left-handed rotation with a Rotation element. If";
       msg += " you meant to make a reflection, use ReflectionRotation";
       msg += " elements, otherwise, please check your matrix.  Other";
       msg += " errors may follow.  Rotation  matrix not created.";
       cout << msg << endl; // this could become a throwWarning or something.
     }
   else if (name == "ReflectionRotation" && isLeftHanded(x, y, z, nmspace) == 1) 
     {
       ExprEvalInterface & ev = ExprEvalSingleton::instance();
       DDRotation ddrot = 
	 DDrotReflect(getDDName(nmspace)
		      , ev.eval(nmspace, atts.find("thetaX")->second)
		      , ev.eval(nmspace, atts.find("phiX")->second)
		      , ev.eval(nmspace, atts.find("thetaY")->second)
		      , ev.eval(nmspace, atts.find("phiY")->second)
		      , ev.eval(nmspace, atts.find("thetaZ")->second)
		      , ev.eval(nmspace, atts.find("phiZ")->second));
       DCOUT_V ('p', "Rotation created: " << ddrot << endl);
     }
   else if (name == "ReflectionRotation" && isLeftHanded(x, y, z, nmspace) == 0)
     {
       string msg("WARNING:  Attempted to make a right-handed");
       msg += " rotation using a ReflectionRotation element. ";
       msg += " If you meant to make a Rotation, use Rotation";
       msg += " elements, otherwise, please check your matrix.";
       msg += "  Other errors may follow.  ReflectionRotation";
       msg += " matrix not created.";
       cout << msg << endl; // this could be a throwWarning or something.
     }
   else
     {
       string msg = "\nDDLRotationAndReflection::processElement tried to process wrong element.";
       throwError(msg);
     }
  } catch (DDException & e) {
    string msg(e.what());
    msg += "\nDDLRotationAndReflection failed to created requested DD object.";
    throwError(msg);
  }
  // after a rotation or reflection rotation has been processed, clear it
    clear();

  DCOUT_V('P', "DDLRotationAndReflection::processElement completed");
}


// returns 1 if it is a left-handed CLHEP rotation matrix, 0 if not, but is okay, -1 if 
// it is not an orthonormal matrix.
//
// Upon encountering the end tag of a Rotation element, we've got to feed
// the appropriate rotation in to the DDCore.  This is an attempt to do so.
//
// Basically, I cannibalized code from g3tog4 (see http link below) and then
// provided the information from our DDL to the same calls.  Tim Cox showed me
// how to build the rotation matrix (mathematically) and the g3tog4 code basically
// did the rest.
//

int DDLRotationAndReflection::isLeftHanded (Hep3Vector x, Hep3Vector y, Hep3Vector z, const string & nmspace)
{
  DCOUT_V('P', "DDLRotation::isLeftHanded started");

  int ret = 0;

  /**************** copied and cannibalized code:
 
  from g3tog4
 
  http://atlassw1.phy.bnl.gov/lxr/source/external/geant4.3.1/source/g3tog4/src/G4gsrotm.cc

 48         // Construct unit vectors 
 49     
 50     G4ThreeVector x(sin(th1r)*cos(phi1r), sin(th1r)*sin(phi1r), cos(th1r)->second;
 51     G4ThreeVector y(sin(th2r)*cos(phi2r), sin(th2r)*sin(phi2r), cos(th2r));
 52     G4ThreeVector z(sin(th3r)*cos(phi3r), sin(th3r)*sin(phi3r), cos(th3r));
 53 
 54         // check for orthonormality and left-handedness
 55 
 56     G4double check = (x.cross(y))*z;
 57     G4double tol = 1.0e-3;
 58         
 59     if (1-abs(check)>tol) {
 60         G4cerr << "Coordinate axes forming rotation matrix "
 61                << irot << " are not orthonormal.(" << 1-abs(check) << ")" 
 62          << G4endl;
 63         G4cerr << " thetaX=" << theta1;
 64         G4cerr << " phiX=" << phi1;
 65         G4cerr << " thetaY=" << theta2;
 66         G4cerr << " phiY=" << phi2;
 67         G4cerr << " thetaZ=" << theta3;
 68         G4cerr << " phiZ=" << phi3;
 69         G4cerr << G4endl;
 70         G4Exception("G4gsrotm error");
 71     }
 72     else if (1+check<=tol) {
 73         G4cerr << "G4gsrotm warning: coordinate axes forming rotation "
 74                << "matrix " << irot << " are left-handed" << G4endl;
 75     }   
 76
 77     G3toG4RotationMatrix* rotp = new G3toG4RotationMatrix;
 78 
 79     rotp->SetRotationMatrixByRow(x, y, z);

  ****************/


  // check for orthonormality and left-handedness
  
  double check = (x.cross(y))*z;
  double tol = 1.0e-3;
  ExprEvalInterface & ev = ExprEvalSingleton::instance();
  DDXMLAttribute atts = getAttributeSet();
  
  if (1-abs(check)>tol) {
    cout << "DDLRotationAndReflection Coordinate axes forming rotation matrix "
	 << getDDName(nmspace) 
	 << " are not orthonormal.(tolerance=" << tol 
	 << " check=" << abs(check)  << ")" 
	 << endl
	 << " thetaX=" << (atts.find("thetaX")->second) 
	 << ' ' << ev.eval(nmspace, atts.find("thetaX")->second)/deg << endl
	 << " phiX=" << (atts.find("phiX")->second) 
	 << ' ' << ev.eval(nmspace, atts.find("phiX")->second)/deg << endl
	 << " thetaY=" << (atts.find("thetaY")->second) 
	 << ' ' << ev.eval(nmspace, atts.find("thetaY")->second)/deg << endl
	 << " phiY=" << (atts.find("phiY")->second)
	 << ' ' << ev.eval(nmspace, atts.find("phiY")->second)/deg << endl
	 << " thetaZ=" << (atts.find("thetaZ")->second)
	 << ' ' << ev.eval(nmspace, atts.find("thetaZ")->second)/deg << endl
	 << " phiZ=" << (atts.find("phiZ")->second)
	 << ' ' << ev.eval(nmspace, atts.find("phiZ")->second)/deg 
	 << endl
	 << "  WAS NOT CREATED!" << endl;
    ret = -1;
  }
  else if (1+check<=tol) {
    ret = 1;    
  }
  DCOUT_V('P', "DDLRotation::isLeftHanded completed");
  return ret;
}

Hep3Vector DDLRotationAndReflection::makeX(string nmspace)
{
  Hep3Vector x = 0;
  DDXMLAttribute atts = getAttributeSet();
  if (atts.find("thetaX") != atts.end())
    {
      ExprEvalInterface & ev = ExprEvalSingleton::instance(); 
      double thetaX = ev.eval(nmspace, atts.find("thetaX")->second.c_str());
      double phiX = ev.eval(nmspace, atts.find("phiX")->second.c_str());
      // colx
      x[0] = sin(thetaX) * cos(phiX);
      x[1] = sin(thetaX) * sin(phiX);
      x[2] = cos(thetaX);
    }
  return x;
}

Hep3Vector DDLRotationAndReflection::makeY(string nmspace)
{
  Hep3Vector y = 0;
  DDXMLAttribute atts = getAttributeSet();
  if (atts.find("thetaY") != atts.end())
    {
      ExprEvalInterface & ev = ExprEvalSingleton::instance(); 
      double thetaY = ev.eval(nmspace, atts.find("thetaY")->second.c_str());
      double phiY = ev.eval(nmspace, atts.find("phiY")->second.c_str());
      
      // coly
      y[0] = sin(thetaY) * cos(phiY);
      y[1] = sin(thetaY) * sin(phiY);
      y[2] = cos(thetaY);
    }
  return y;
}

Hep3Vector DDLRotationAndReflection::makeZ(string nmspace)
{
  Hep3Vector z = 0;
  DDXMLAttribute atts = getAttributeSet();
  if (atts.find("thetaZ") != atts.end())
    {
      ExprEvalInterface & ev = ExprEvalSingleton::instance(); 
      double thetaZ = ev.eval(nmspace, atts.find("thetaZ")->second.c_str());
      double phiZ = ev.eval(nmspace, atts.find("phiZ")->second.c_str());
      
      // colz
      z[0] = sin(thetaZ) * cos(phiZ);
      z[1] = sin(thetaZ) * sin(phiZ);
      z[2] = cos(thetaZ);
    }
  return z;
}
