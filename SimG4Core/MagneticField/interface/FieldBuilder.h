#ifndef SimG4Core_FieldBuilder_H
#define SimG4Core_FieldBuilder_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <memory>

class DDLogicalPart;
class MagneticField;

class G4FieldManager;
class G4Mag_UsualEqRhs;
class G4PropagatorInField;
class G4LogicalVolume;

namespace sim {
   class Field;
   class FieldBuilder
   {
      public:
	 FieldBuilder(const MagneticField*, const edm::ParameterSet&);
	 //~FieldBuilder();
	 void readFieldParameters(DDLogicalPart theLogicalPart,
				  const std::string& keywordField);
	 void configure(const std::string& keywordField,
			G4FieldManager * fM = 0,
			G4PropagatorInField * fP = 0);
	 G4LogicalVolume * fieldTopVolume();
      private:
	 void configureFieldManager(G4FieldManager * fM);
	 void configurePropagatorInField(G4PropagatorInField * fP);  
      private:
	 std::auto_ptr<Field> theField;
	 G4Mag_UsualEqRhs * theFieldEquation;
	 G4LogicalVolume * theTopVolume;
	 std::string keywordField;
	 std::string fieldType;
	 double fieldValue;
	 std::string stepper;
	 double minStep;
	 double dChord;
	 double dOneStep;
	 double dIntersection;
	 double dIntersectionAndOneStep;
	 double maxLoopCount;
	 double minEpsilonStep;
	 double maxEpsilonStep;
   };
}

#endif
