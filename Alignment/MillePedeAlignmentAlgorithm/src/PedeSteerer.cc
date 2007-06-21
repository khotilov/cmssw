/**
 * \file PedeSteerer.cc
 *
 *  \author    : Gero Flucke
 *  date       : October 2006
 *  $Revision: 1.11.2.4 $
 *  $Date: 2007/06/21 12:43:03 $
 *  (last update by $Author: flucke $)
 */

#include "PedeSteerer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Alignment/CommonAlignment/interface/Alignable.h"
#include "Alignment/CommonAlignment/interface/Utilities.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParameterStore.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParameterSelector.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/SelectionUserVariables.h"
#include "Alignment/CommonAlignmentParametrization/interface/RigidBodyAlignmentParameters.h"

#include "Alignment/TrackerAlignment/interface/TrackerAlignableId.h"
// for 'type identification' as Alignable
#include "Alignment/TrackerAlignment/interface/AlignableTracker.h"
#include "Alignment/MuonAlignment/interface/AlignableMuon.h"


#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include <fstream>
#include <sstream>
#include <algorithm>

// from ROOT
#include <TSystem.h>
#include <TMath.h>

// NOTE: The '+4', is backward incompatible for results with older binary files...
const unsigned int PedeSteerer::theMaxNumParam = RigidBodyAlignmentParameters::N_PARAM + 4;
const unsigned int PedeSteerer::theMinLabel = 1; // must be > 0

//___________________________________________________________________________
PedeSteerer::PedeSteerer(AlignableTracker *aliTracker, AlignableMuon *aliMuon,
			 AlignmentParameterStore *store,
			 const edm::ParameterSet &config, const std::string &defaultDir) :
  myParameterStore(store), myConfig(config),
  myDirectory(myConfig.getUntrackedParameter<std::string>("fileDir")),
  myParameterSign(myConfig.getUntrackedParameter<int>("parameterSign")),
  theCoordMaster(0)
{
  if (myParameterSign != 1 && myParameterSign != -1) {
    cms::Exception("BadConfig") << "Expect PedeSteerer.parameterSign = +/-1, "
				<< "found " << myParameterSign << ".";
  }

  // Correct directory, needed before asking for fileName(..):
  if (myDirectory.empty()) myDirectory = defaultDir;
  if (!myDirectory.empty() && myDirectory.find_last_of('/') != myDirectory.size() - 1) {
    myDirectory += '/'; // directory may need '/'
  }

  this->buildMap(aliTracker, aliMuon); //has to be done first
  const std::vector<Alignable*> &alis = myParameterStore->alignables();

  // Coordinate system selection and correction before fixing, at least for fixing at truth:
  theCoordDefiners = this->selectCoordinateAlis(alis);
  if (!theCoordDefiners.empty()) { // Create steering with constraints to define coordinate system:
    const std::string nameCoordFile(this->fileName("Coord"));
    theCoordMaster = aliTracker; // Both tracker and muon have global coordinates, though...
    if (!theCoordMaster) theCoordMaster = aliMuon; // ...which coordinates should not really matter anyway.
    // FIXME: With both tracker and muon, only tracker will be shifted in correctToReferenceSystem!
    // Could be fixed by introducing AlignableComposite with both, but must prevent double deletion.
    // If this has 'global frame', correctToReferenceSystem can use globalParameters (as it does)
    this->defineCoordinates(theCoordDefiners, theCoordMaster, nameCoordFile);
    edm::LogInfo("Alignment") << "@SUB=PedeSteerer" 
                              << theCoordDefiners.size() << " highest level objects define the "
			      << "coordinate system, steering file " << nameCoordFile << ".";
    this->correctToReferenceSystem();
  } 

  // Fixing parameters also checks for invalid parameter choices (should that be done first?):
  const std::string nameFixFile(this->fileName("FixPara"));
  const std::pair<unsigned int, unsigned int> nFixFixCor(this->fixParameters(alis, nameFixFile));
  if (nFixFixCor.first != 0 || nFixFixCor.second != 0) {
    edm::LogInfo("Alignment") << "@SUB=PedeSteerer" 
                              << nFixFixCor.first << " parameters fixed at 0. and "
                              << nFixFixCor.second << " at 'original' position, "
                              << "steering file " << nameFixFile << ".";
  } 

  if (this->buildNoHierarchyCollection(alis)) { // before hierarchyConstraints(..)
    edm::LogInfo("Alignment")
      << "@SUB=PedeSteerer" << myNoHieraCollection.size() << " alignables taken out of hierarchy.";
  }

  const std::string nameHierarchyFile(this->fileName("Hierarchy"));
  unsigned int nConstraint = this->hierarchyConstraints(alis, nameHierarchyFile);
  if (nConstraint) {
    edm::LogInfo("Alignment") << "@SUB=PedeSteerer" 
                              << "Hierarchy constraints for " << nConstraint << " alignables, "
                              << "steering file " << nameHierarchyFile << ".";
  }

  // Delete all SelectionUserVariables now? They will anyway be overwritten by MillePedeVariables...
}

//___________________________________________________________________________

PedeSteerer::~PedeSteerer()
{
}

//___________________________________________________________________________
/// Return 32-bit unique label for alignable, 0 indicates failure.
unsigned int PedeSteerer::alignableLabel(Alignable *alignable) const
{
  if (!alignable) return 0;

  AlignableToIdMap::const_iterator position = myAlignableToIdMap.find(alignable);
  if (position != myAlignableToIdMap.end()) {
    return position->second;
  } else {
    //throw cms::Exception("LogicError") 
    edm::LogWarning("LogicError")
      << "@SUB=PedeSteerer::alignableLabel" << "Alignable "
      << typeid(*alignable).name() << " not in map";
    return 0;
  }

  /*
// following ansatz does not work since the maximum label allowed by pede is 99 999 999...
//   TrackerAlignableId idProducer;
//   const DetId detId(alignable->geomDetId()); // does not work: only AlignableDet(Unit) has DetId...
//  if (detId.det() != DetId::Tracker) {
  const unsigned int detOffset = 28; // would like to use definition from DetId
  const TrackerAlignableId::UniqueId uniqueId(idProducer.alignableUniqueId(alignable));
  const uint32_t detId = uniqueId.first; // uniqueId is a pair...
  const uint32_t det = detId >> detOffset; // most significant bits are detector part
  if (det != DetId::Tracker) {
    //throw cms::Exception("LogicError") 
    edm::LogWarning("LogicError") << "@SUB=PedeSteerer::alignableLabel "
      << "Expecting DetId::Tracker (=" << DetId::Tracker << "), but found "
      << det << " which would make the pede labels ambigous. "
      << typeid(*alignable).name() << " " << detId;
    return 0;
  }
  // FIXME: Want: const AlignableObjectId::AlignableObjectIdType type = 
  const unsigned int aType = static_cast<unsigned int>(uniqueId.second);// alignable->alignableObjectId();
  if (aType != ((aType << detOffset) >> detOffset)) {
    // i.e. too many bits (luckily we are  not the muon system...)
    throw cms::Exception("LogicError")  << "@SUB=PedeSteerer::alignableLabel "
      << "Expecting alignableTypeId with at most " << 32 - detOffset
      << " bits, but the number is " << aType
      << " which would make the pede labels ambigous.";
    return 0;
  }

  const uint32_t detIdWithoutDet = (detId - (det << detOffset));
  return detIdWithoutDet + (aType << detOffset);
*/
}

//_________________________________________________________________________
unsigned int PedeSteerer::parameterLabel(unsigned int aliLabel, unsigned int parNum) const
{
  if (parNum >= theMaxNumParam) {
    throw cms::Exception("Alignment") << "@SUB=PedeSteerer::parameterLabel" 
                                      << "Parameter number " << parNum 
                                      << " out of range 0 <= num < " << theMaxNumParam;
  }
  return aliLabel + parNum;
}

//___________________________________________________________________________
unsigned int PedeSteerer::paramNumFromLabel(unsigned int paramLabel) const
{
  if (paramLabel < theMinLabel) {
    edm::LogError("LogicError") << "@SUB=PedeSteerer::paramNumFromLabel"
                                << "Label " << paramLabel << " should be >= " << theMinLabel;
    return 0;
  }
  return (paramLabel - theMinLabel) % theMaxNumParam;
}

//___________________________________________________________________________
unsigned int PedeSteerer::alignableLabelFromLabel(unsigned int paramLabel) const
{
  return paramLabel - this->paramNumFromLabel(paramLabel);
}

//___________________________________________________________________________
Alignable* PedeSteerer::alignableFromLabel(unsigned int label) const
{
  const unsigned int aliLabel = this->alignableLabelFromLabel(label);
  if (aliLabel < theMinLabel) return 0; // error already given
  
  if (myIdToAlignableMap.empty()) const_cast<PedeSteerer*>(this)->buildReverseMap();
  IdToAlignableMap::const_iterator position = myIdToAlignableMap.find(aliLabel);
  if (position != myIdToAlignableMap.end()) {
    return position->second;
  } else {
    edm::LogError("LogicError") << "@SUB=PedeSteerer::alignableFromLabel"
                                << "Alignable label " << aliLabel << " not in map";
    return 0;
  }
}

//_________________________________________________________________________
bool PedeSteerer::isNoHiera(const Alignable* ali) const
{
  return (myNoHieraCollection.find(ali) != myNoHieraCollection.end());
}

//_________________________________________________________________________
double PedeSteerer::cmsToPedeFactor(unsigned int parNum) const
{
  return 1.; // mmh, otherwise would need to FIXME hierarchyConstraint...

  switch (parNum) {
  case RigidBodyAlignmentParameters::dx:
  case RigidBodyAlignmentParameters::dy:
    return 1000.; // cm to mum *1/10 to get smaller values
  case RigidBodyAlignmentParameters::dz:
    return 2500.;   // cm to mum *1/4 
  case RigidBodyAlignmentParameters::dalpha:
  case RigidBodyAlignmentParameters::dbeta:
    return 1000.; // rad to mrad (no first guess for sensitivity yet)
  case RigidBodyAlignmentParameters::dgamma:
    return 10000.; // rad to mrad *10 to get larger values
  default:
    return 1.;
  }
}

//_________________________________________________________________________
unsigned int PedeSteerer::buildMap(Alignable *highestLevelAli1, Alignable *highestLevelAli2)
{

  myAlignableToIdMap.clear(); // just in case of re-use...

  std::vector<Alignable*> allComps;

  if (highestLevelAli1) {
    allComps.push_back(highestLevelAli1);
    highestLevelAli1->recursiveComponents(allComps);
  }
  if (highestLevelAli2) {
    allComps.push_back(highestLevelAli2);
    highestLevelAli2->recursiveComponents(allComps);
  }

  unsigned int id = theMinLabel;
  for (std::vector<Alignable*>::const_iterator iter = allComps.begin();
       iter != allComps.end(); ++iter) {
    myAlignableToIdMap.insert(AlignableToIdPair(*iter, id));
    id += theMaxNumParam;
  }
  
  return allComps.size();
}


//_________________________________________________________________________
unsigned int PedeSteerer::buildReverseMap()
{

  myIdToAlignableMap.clear();  // just in case of re-use...

  for (AlignableToIdMap::iterator it = myAlignableToIdMap.begin();
       it != myAlignableToIdMap.end(); ++it) {
    const unsigned int key = (*it).second;
    Alignable *ali = (*it).first;
    myIdToAlignableMap[key] = ali;
  }

  return myIdToAlignableMap.size();
}

//_________________________________________________________________________
unsigned int PedeSteerer::buildNoHierarchyCollection(const std::vector<Alignable*> &alis)
{
  myNoHieraCollection.clear();  // just in case of re-use...

  for (std::vector<Alignable*>::const_iterator iAli = alis.begin() ; iAli != alis.end(); ++iAli) {
    AlignmentParameters *params = (*iAli)->alignmentParameters();
    SelectionUserVariables *selVar = dynamic_cast<SelectionUserVariables*>(params->userVariables());
    if (!selVar) continue;
    // Now check whether taking out of hierarchy is selected - must be consistent!
    unsigned int numNoHieraPar = 0;
    unsigned int numHieraPar = 0;
    for (unsigned int iParam = 0; static_cast<int>(iParam) < params->size(); ++iParam) {
      const char selector = selVar->fullSelection()[iParam];
      if (selector == 'C' || selector == 'F' || selector == 'H') {
	++numNoHieraPar;
      } else if (selector == 'c' || selector == 'f' || selector == '1' || selector == 'r') {
	++numHieraPar;
      } // else ... accept '0' as undetermined
    }
    if (numNoHieraPar) { // Selected to be taken out.
      if (numHieraPar) { // Inconsistent: Some parameters still in hierarchy ==> exception!
	throw cms::Exception("BadConfig") 
	  << "[PedeSteerer::buildNoHierarchyCollection] All active parameters of alignables to be "
	  << " taken out of the hierarchy must be marked with capital letters 'C', 'F' or 'H'!";
      }
      bool isInHiera = false; // Check whether Alignable is really part of hierarchy:
      Alignable *mother = *iAli;
      while ((mother = mother->mother())) {
	if (mother->alignmentParameters()) isInHiera = true; // could 'break;', but loop is short
      }
      // Complain, but keep collection short if not in hierarchy:
      if (isInHiera) myNoHieraCollection.insert(*iAli);
      else edm::LogWarning("Alignment") << "@SUB=PedeSteerer::buildNoHierarchyCollection"
					<< "Alignable not in hierarchy, no need to remove it!";
    }
  } // end loop on alignables

  return myNoHieraCollection.size();
}


//_________________________________________________________________________
std::pair<unsigned int, unsigned int>
PedeSteerer::fixParameters(const std::vector<Alignable*> &alis, const std::string &fileName)
{
  // return number of parameters fixed at 0. and fixed at original position 
  std::pair<unsigned int, unsigned int> numFixNumFixCor(0, 0);

  std::ofstream *filePtr = 0;

  for (std::vector<Alignable*>::const_iterator iAli = alis.begin() ; iAli != alis.end(); ++iAli) {
    AlignmentParameters *params = (*iAli)->alignmentParameters();
    SelectionUserVariables *selVar = dynamic_cast<SelectionUserVariables*>(params->userVariables());
    if (!selVar) continue;

    for (unsigned int iParam = 0; static_cast<int>(iParam) < params->size(); ++iParam) {
      int whichFix = this->fixParameter(*iAli, iParam, selVar->fullSelection()[iParam], filePtr,
					fileName);
      if (whichFix == 1) {
        ++(numFixNumFixCor.first);
      } else if (whichFix == -1) {
        ++(numFixNumFixCor.second);
      }
    }
  }

  delete filePtr; // automatically flushes, no problem if NULL ptr.   

  return numFixNumFixCor;
}

//_________________________________________________________________________
int PedeSteerer::fixParameter(Alignable *ali, unsigned int iParam, char selector,
                              std::ofstream* &filePtr, const std::string &fileName)
{
  int result = 0;
  float fixAt = 0.;
  if (selector == 'c' || selector == 'C') {
    fixAt = -this->parameterSign() * RigidBodyAlignmentParameters(ali, true).parameters()[iParam];
    result = -1;
  } else if (selector == 'f' || selector == 'F') {
    result = 1;
  } else if (selector != '1' && selector != '0' && selector != 'r' && selector != 'H') {
    throw cms::Exception("BadConfig")
      << "[PedeSteerer::fixParameter] " << "Unexpected parameter selector '" << selector
      << "', use \n'f/F' (fix),\n'c/C' (fix at correct pos.),\n'1/H' (free),\n"
      << "'r' (free, but defining reference system) or \n'0' (ignore).\n"
      << "Capital letters mean that the Alignable is taken out of a possible hierarchy,\n"
      << "but must be used consistently for all its parameters.";
  }

  if (result) {
    if (!filePtr) {
      filePtr = this->createSteerFile(fileName, true);
      (*filePtr) << "Parameter\n";
    }
    std::ofstream &file = *filePtr;

    const unsigned int aliLabel = this->alignableLabel(ali);
    file << this->parameterLabel(aliLabel, iParam) << "  " 
         << fixAt * this->cmsToPedeFactor(iParam) << " -1.0";
    if (0) { // debug
      const GlobalPoint position(ali->globalPosition());
      file << "* eta " << position.eta() << ", z " << position.z()
           << ", r " << position.perp() << ", phi " << position.phi();
    }
    file << "\n";
  }

  return result;
}

//_________________________________________________________________________
std::vector<Alignable*> PedeSteerer::selectCoordinateAlis(const std::vector<Alignable*> &alis) const
{
  std::vector<Alignable*> coordAlis;

  for (std::vector<Alignable*>::const_iterator iAli = alis.begin() ; iAli != alis.end(); ++iAli) {
    AlignmentParameters *params = (*iAli)->alignmentParameters();
    SelectionUserVariables *selVar = dynamic_cast<SelectionUserVariables*>(params->userVariables());
    if (!selVar) continue;
    unsigned int refParam = 0;
    unsigned int nonRefParam = 0;
    for (unsigned int iParam = 0; static_cast<int>(iParam) < params->size(); ++iParam) {
      const char selector = selVar->fullSelection()[iParam];
      if (selector == 'r') {
	++refParam;
      } else if (selector != '0') { // allow inactive parameters
	++nonRefParam;
      }
    }
    // Check whether some 'r' selection string. If yes and selection makes sense, add to result:
    if (refParam) {
      if (nonRefParam) {
	throw cms::Exception("BadConfig") << "[PedeSteerer::selectCoordinateAlis] "
					  << "All active parameters of alignables defining "
					  << "the coordinate system must be marked with 'r'!";
      } else {
	Alignable *mother = *iAli;
	while ((mother = mother->mother())) {
	  if (mother->alignmentParameters()) {
	    throw cms::Exception("BadConfig") << "[PedeSteerer::selectCoordinateAlis] "
					      << "Alignables defining the coordinate system must "
					      << "be highest level!";
	  }
	}
	coordAlis.push_back(*iAli);
      }
    }
  } // end loop on alignables

  return coordAlis;
}


//_________________________________________________________________________
void PedeSteerer::defineCoordinates(const std::vector<Alignable*> &alis, Alignable *aliMaster,
				    const std::string &fileName)
{
  std::ofstream *filePtr = this->createSteerFile(fileName, true);
  (*filePtr) << "* Contraints to define coordinate system:\n";
  if (!aliMaster || aliMaster->alignmentParameters()) {
    throw cms::Exception("BadConfig")
      << "[PedeSteerer::defineCoordinates] " << "No master alignable or it has parameters!";
  }
  AlignmentParameters *par = new RigidBodyAlignmentParameters(aliMaster, false);
  aliMaster->setAlignmentParameters(par); // hierarchyConstraint needs parameters
  this->hierarchyConstraint(aliMaster, alis, *filePtr);
  aliMaster->setAlignmentParameters(0); // erase dummy parameters

  delete filePtr; // automatically flushes, no problem if NULL ptr.   
}

//_________________________________________________________________________
void PedeSteerer::correctToReferenceSystem()
{
  typedef RigidBodyAlignmentParameters RbPars;
  if (!theCoordMaster || theCoordDefiners.empty()) return; // nothing was defined

  std::vector<Alignable*> definerDets; // or ...DetUnits
  for (std::vector<Alignable*>::iterator it = theCoordDefiners.begin(), iE = theCoordDefiners.end();
       it != iE; ++it) {// find lowest level objects of alignables that define the coordinate system
    (*it)->deepComponents(definerDets);
  }

  for (unsigned int iLoop = 0; ; ++iLoop) { // iterate: shifts and rotations are not independent
    AlgebraicVector meanPars(RbPars::N_PARAM);
    for (std::vector<Alignable*>::iterator it = definerDets.begin(), iE = definerDets.end();
	 it != iE; ++it) { // sum up mean displacements/misrotations:
      meanPars += RbPars(*it, true).globalParameters();// requires theCoordMaster has global frame
    }
    meanPars /= definerDets.size();

    align::Scalar squareSum = 0.;
    for (unsigned int i = 0; i < RbPars::N_PARAM; ++i) {
      squareSum += meanPars[i] * meanPars[i] * this->cmsToPedeFactor(i) * this->cmsToPedeFactor(i);
    }

    if (true || iLoop == 0) {
      edm::LogInfo("Alignment") << "@SUB=PedeSteerer::correctToReferenceSystem"
				<< "Loop " << iLoop << " "
				<< "Mean misalignment of dets of defined coordinate system "
				<< "(will be iteratively corrected to about 0.):" << meanPars;
    }
    if (iLoop >=5 || squareSum < 1.e-20) { // sqrt(1.e-20)=1.e-10: close enough to stop iterating
      if (iLoop >=5) {
	edm::LogError("Alignment") << "@SUB=PedeSteerer::correctToReferenceSystem"
				   << "No convergence in " << iLoop << " iterations, " 
				   << "remaining misalignment: " << meanPars;
      }
      break;
    }

    const GlobalVector globalShift(meanPars[RbPars::dx],meanPars[RbPars::dy],meanPars[RbPars::dz]);
    align::EulerAngles globalAngles(3);
    globalAngles[0] = meanPars[RbPars::dalpha];
    globalAngles[1] = meanPars[RbPars::dbeta];
    globalAngles[2] = meanPars[RbPars::dgamma];
    theCoordMaster->move(-globalShift); // sign to revert
    theCoordMaster->rotateInGlobalFrame(align::toMatrix(-globalAngles)); // sign to revert
  }
  
}

//_________________________________________________________________________
unsigned int PedeSteerer::hierarchyConstraints(const std::vector<Alignable*> &alis,
					       const std::string &fileName)
{
  std::ofstream *filePtr = 0;

  unsigned int nConstraints = 0;
  std::vector<Alignable*> aliDaughts;
  for (std::vector<Alignable*>::const_iterator iA = alis.begin(), iEnd = alis.end();
       iA != iEnd; ++iA) {
    aliDaughts.clear();
    if (!(*iA)->firstCompsWithParams(aliDaughts)) {
      edm::LogError("Alignment") << "@SUB=PedeSteerer::hierarchyConstraints"
                                 << "Some but not all daughters with params!";
    }
//     edm::LogInfo("Alignment") << "@SUB=PedeSteerer::hierarchyConstraints"
// 			      << aliDaughts.size() << " ali param components";
    if (aliDaughts.empty()) continue;
//     edm::LogInfo("Alignment") << "@SUB=PedeSteerer::hierarchyConstraints"
// 			      << aliDaughts.size() << " alignable components ("
// 			      << (*iA)->size() << " in total) for " 
// 			      << aliId.alignableTypeName(*iA) 
// 			      << ", layer " << aliId.typeAndLayerFromAlignable(*iA).second
// 			      << ", position " << (*iA)->globalPosition()
// 			      << ", r = " << (*iA)->globalPosition().perp();
    if (!filePtr) filePtr = this->createSteerFile(fileName, true);
    ++nConstraints;
    this->hierarchyConstraint(*iA, aliDaughts, *filePtr);
  }

  delete filePtr; // automatically flushes, no problem if NULL ptr.   

  return nConstraints;
}

//_________________________________________________________________________
void PedeSteerer::hierarchyConstraint(const Alignable *ali,
                                      const std::vector<Alignable*> &components,
				      std::ofstream &file) const
{
  typedef AlignmentParameterStore::ParameterId ParameterId;
  typedef std::vector<Alignable*>::size_type IndexType;

  std::vector<std::vector<ParameterId> > paramIdsVec;
  std::vector<std::vector<float> > factorsVec;
  if (!myParameterStore->hierarchyConstraints(ali, components, paramIdsVec, factorsVec)) {
    edm::LogWarning("Alignment") << "@SUB=PedeSteerer::hierarchyConstraint"
				 << "Problems from store.";
  }

  for (unsigned int iConstr = 0; iConstr < paramIdsVec.size(); ++iConstr) {
    std::ostringstream aConstr;

    const std::vector<ParameterId> &parIds = paramIdsVec[iConstr];
    const std::vector<float> &factors = factorsVec[iConstr];
    // parIds.size() == factors.size() granted by myParameterStore->hierarchyConstraints
    for (unsigned int iParam = 0; iParam < parIds.size(); ++iParam) {
      Alignable *aliSubComp = parIds[iParam].first;
      const unsigned int compParNum = parIds[iParam].second;
      if (this->isNoHiera(aliSubComp)) {
	if (0) aConstr << "* Taken out of hierarchy: "; // conflict with !aConstr.str().empty()
	else continue;
      }
      const unsigned int aliLabel = this->alignableLabel(aliSubComp);
      const unsigned int paramLabel = this->parameterLabel(aliLabel, compParNum);
      // FIXME: multiply by cmsToPedeFactor(subcomponent)/cmsToPedeFactor(mother) (or vice a versa?)
      aConstr << paramLabel << "    " << factors[iParam];
      if (true) { // debug
	const TrackerAlignableId aliId;
	aConstr << "   ! for param " << compParNum << " of a " 
		<< aliId.alignableTypeName(aliSubComp) << " (label " << aliLabel << ")";
      }
      aConstr << "\n";
    } // end loop on params

    if (!aConstr.str().empty()) {
      if (true) { //debug
	const TrackerAlignableId aliId;
	file << "\n* Nr. " << iConstr << " of a '" << aliId.alignableTypeName(ali) << "' (label "
	     << this->alignableLabel(const_cast<Alignable*>(ali)) // ugly cast: FIXME!
	     << "), layer " << aliId.typeAndLayerFromAlignable(ali).second
	     << ", position " << ali->globalPosition()
	     << ", r = " << ali->globalPosition().perp();
      }
      file << "\nConstraint   0.\n" << aConstr.str(); // in future 'Wconstraint'?
    }
  } // end loop on constraints
}

//_________________________________________________________________________
std::ofstream* PedeSteerer::createSteerFile(const std::string &name, bool addToList)
{
  std::ofstream *result = new std::ofstream(name.c_str(), std::ios::out);
  if (!result || !result->is_open()) {
    delete result; // needed before exception in case just open failed
    throw cms::Exception("FileOpenProblem") << "[PedeSteerer::createSteerFile]" 
					    << "Could not open " << name 
					    << " as output file.";
  } else if (addToList) {
    mySteeringFiles.push_back(name); // keep track
  }

  return result;
}


//_________________________________________________________________________
std::string PedeSteerer::fileName(const std::string &addendum) const
{

  std::string name(myDirectory);
  name += myConfig.getParameter<std::string>("steerFile");
  name += addendum;
  name += ".txt";

  return name;
}

//_________________________________________________________________________
std::string PedeSteerer::buildMasterSteer(const std::vector<std::string> &binaryFiles)
{
  const std::string nameMasterSteer(this->fileName("Master"));
  std::ofstream *mainSteerPtr = this->createSteerFile(nameMasterSteer, false);
  if (!mainSteerPtr) return "";

  // add steering files to master steering file
  std::ofstream &mainSteerRef = *mainSteerPtr;
  for (unsigned int iFile = 0; iFile < mySteeringFiles.size(); ++iFile) {
    mainSteerRef << mySteeringFiles[iFile] << "\n";
  }

  // add binary files to master steering file
  mainSteerRef << "\nCfiles\n";
  for (unsigned int iFile = 0; iFile < binaryFiles.size(); ++iFile) {
    mainSteerRef << binaryFiles[iFile] << "\n";
  }

  // add method
  mainSteerRef << "\nmethod  " << myConfig.getParameter<std::string>("method") << "\n";

  // add further options
  const std::vector<std::string> opt(myConfig.getParameter<std::vector<std::string> >("options"));
  mainSteerRef << "\n* Outlier treatment and other options \n";
  for (unsigned int i = 0; i < opt.size(); ++i) {
    mainSteerRef << opt[i] << "\n";
  }

  delete mainSteerPtr;  // close (and flush) again

  return nameMasterSteer;
}

//_________________________________________________________________________
bool PedeSteerer::runPede(const std::string &masterSteer) const
{
  if (masterSteer.empty()) {
    edm::LogError("Alignment") << "@SUB=PedeSteerer::runPede" << "Empty master steer file, stop";
    return false;
  }

  std::string command(myConfig.getUntrackedParameter<std::string>("pedeCommand"));
  (command += " ") += masterSteer;
  const std::string dump(myConfig.getUntrackedParameter<std::string>("pedeDump"));
  if (!dump.empty()) {
    command += " > ";
    (command += myDirectory) += dump;
  }

  edm::LogInfo("Alignment") << "@SUB=PedeSteerer::runPede" << "Start running " << command;
  // FIXME: Recommended interface to system commands?
  int shellReturn = gSystem->Exec(command.c_str());
  if (shellReturn) {
    edm::LogError("Alignment") << "@SUB=PedeSteerer::runPede" << "Command returns " << shellReturn;
  } else {
    edm::LogInfo("Alignment") << "@SUB=PedeSteerer::runPede" << "Command returns " << shellReturn;
  }

  return !shellReturn;
}
