#ifndef DDSolidShapes_h
#define DDSolidShapes_h
#include "DetectorDescription/DDBase/interface/DDException.h"

enum DDSolidShape { dd_not_init,
                    ddbox, ddtubs, ddtrap, ddcons,
                    ddpolycone_rz, ddpolyhedra_rz,
		    ddpolycone_rrz, ddpolyhedra_rrz,
                    ddunion, ddsubtraction, ddintersection,
		    ddreflected,
		    ddshapeless,
		    ddpseudotrap, ddtrunctubs
		   };
		   
struct DDSolidShapesName {

  static const char * name(DDSolidShape s) 
  {
    static const char* c[] = { 
      "Solid not initialized",
      "Box", "Tube(section)", "Trapezoid", "Cone(section)",
      "Polycone_rz", "Polyhedra_rz",
      "Polycone_rrz", "Polyhedra_rrz",
      "UnionSolid", "SubtractionSolid", "IntersectionSolid",
      "ReflectedSolid", 
      "ShapelessSolid",
      "PseudoTrapezoid","TruncatedTube(section)"
    };
    return c[s];   			  
  }
  
  static const DDSolidShape index( const int& ind ) {
    switch (ind) {
    case 0:
      return dd_not_init;
      break;
    case 1:
      return ddbox;
      break;
    case 2:
      return ddtubs;
      break;
    case 3:
      return ddtrap;
      break;
    case 4:
      return ddcons;
      break;
    case 5:
      return ddpolycone_rz;
      break;
    case 6:
      return ddpolyhedra_rz;
      break;
    case 7:
      return ddpolycone_rrz;
      break;
    case 8:
      return ddpolyhedra_rrz;
      break;
    case 9:
      return ddunion;
      break;
    case 10:
      return ddsubtraction;
      break;
    case 11:
      return ddintersection;
      break;
    case 12:
      return ddreflected;
      break;
    case 13:
      return ddshapeless;
      break;
    case 14:
      return ddpseudotrap;
      break;
    case 15:
      return ddtrunctubs;
      break;
    default:
      throw DDException("DDSolidShapes:index wrong shape");   
      break;
    }
  }

};

		   

#endif
