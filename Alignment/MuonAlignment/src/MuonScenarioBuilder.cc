/** \file
 *
 *  $Date: 2007/01/12 09:47:42 $
 *  $Revision: 1.1 $
 *  \author Andre Sznajder - UERJ(Brazil)
 */
 

#include <string>
#include <iostream>
#include <sstream>

// Framework
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Alignment

#include "Alignment/MuonAlignment/interface/MuonScenarioBuilder.h"

//__________________________________________________________________________________________________
MuonScenarioBuilder::MuonScenarioBuilder( Alignable* alignable )
{

  theAlignableMuon = dynamic_cast<AlignableMuon*>( alignable );

  if ( !theAlignableMuon )
    throw cms::Exception("TypeMismatch") << "Argument is not an AlignableMuon";

}


//__________________________________________________________________________________________________
void MuonScenarioBuilder::applyScenario( const edm::ParameterSet& scenario )
{

  // Apply the scenario to all main components of Muon.
  theScenario = scenario;
  theModifierCounter = 0;

  // Seed is set at top-level, and is mandatory
  if ( this->hasParameter_( "seed", theScenario ) )
	theModifier.setSeed( static_cast<long>(theScenario.getParameter<int>("seed")) );
  else
	throw cms::Exception("BadConfig") << "No generator seed defined!";  



  // DT Barrel
  std::vector<Alignable*> dtBarrel = theAlignableMuon->DTBarrel();
  this->decodeMovements_( theScenario, dtBarrel, "DTBarrel" );
  // CSC Endcap
  std::vector<Alignable*> cscEndcaps = theAlignableMuon->CSCEndcaps();
  this->decodeMovements_( theScenario, cscEndcaps, "CSCEndcap" );

  this->moveDTSectors(theScenario);
  this->moveCSCSectors(theScenario);
  this->moveMuon(theScenario);
  
  edm::LogInfo("TrackerScenarioBuilder") 
	<< "Applied modifications to " << theModifierCounter << " alignables";

}



std::vector<float> MuonScenarioBuilder::extractParameters(const edm::ParameterSet& pSet, char *blockId) {
  
  double scale_ = 0, scaleError_ = 0, phiX_ = 0, phiY_ = 0, phiZ_ = 0;
  double dX_ = 0, dY_ = 0, dZ_ = 0;
  
  std::ostringstream error;
  edm::ParameterSet Parameters = this->getParameterSet_((std::string)blockId, pSet);
  std::vector<std::string> parameterNames = Parameters.getParameterNames();
  for ( std::vector<std::string>::iterator iParam = parameterNames.begin(); iParam != parameterNames.end(); iParam++ ) {
    if ( (*iParam) == "scale" )    scale_ = Parameters.getParameter<double>( *iParam );
    else if ( (*iParam) == "scaleError" ) scaleError_ = Parameters.getParameter<double>( *iParam );
    else if ( (*iParam) == "phiX" )     phiX_     = Parameters.getParameter<double>( *iParam );
    else if ( (*iParam) == "phiY" )     phiY_     = Parameters.getParameter<double>( *iParam );
    else if ( (*iParam) == "phiZ" )     phiZ_     = Parameters.getParameter<double>( *iParam );
    else if ( (*iParam) == "dX" )       dX_       = Parameters.getParameter<double>( *iParam );
    else if ( (*iParam) == "dY" )       dY_       = Parameters.getParameter<double>( *iParam );
    else if ( (*iParam) == "dZ" )       dZ_       = Parameters.getParameter<double>( *iParam );
    else if ( Parameters.retrieve( *iParam ).typeCode() != 'P' )
      { // Add unknown parameter to list
	if ( !error.str().length() ) error << "Unknown parameter name(s): ";
	error << " " << *iParam;
      }
  }
  std::vector<float> param;
  param.push_back(scale_); param.push_back(scaleError_);
  param.push_back(phiX_); param.push_back(phiY_);
  param.push_back(phiZ_); param.push_back(dX_);
  param.push_back(dY_); param.push_back(dZ_);
  return param;

}

//_____________________________________________________________________________________________________
void MuonScenarioBuilder::moveDTSectors(const edm::ParameterSet& pSet) {
  
  std::vector<Alignable *> DTchambers = theAlignableMuon->DTChambers();
  //Take parameters
  std::vector<float> param = this->extractParameters(pSet, "DTsectors");
  float scale_ = param[0]; float scaleError_ = param[1];
  float phiX_ = param[2]; float phiY_ = param[3]; float phiZ_ = param[4];
  float dX_ = param[5]; float dY_ = param[6]; float dZ_ = param[7];
  
  float dx = scale_*dX_; float dy = scale_*dY_; float dz = scale_*dZ_;
  float phix = scale_*phiX_; float phiy = scale_*phiY_; float phiz = scale_*phiZ_;
  float errorx = scaleError_*dX_; float errory = scaleError_*dY_; float errorz = scaleError_*dZ_;
  float errorphix = scaleError_*phiX_; float errorphiy = scaleError_*phiY_; float errorphiz = scaleError_*phiZ_;
  std::vector<float> errorDisp;
  errorDisp.push_back(errorx); errorDisp.push_back(errory); errorDisp.push_back(errorz);
  std::vector<float> errorRotation;
  errorRotation.push_back(errorphix); errorRotation.push_back(errorphiy); errorRotation.push_back(errorphiz);
 
  int index[5][4][14];
  int counter = 0;
  //Create and index for the chambers in the Alignable vector
  for(std::vector<Alignable *>::iterator iter = DTchambers.begin(); iter != DTchambers.end(); ++iter) {
    DTChamberId myId((*iter)->geomDetId().rawId());
    index[myId.wheel()+2][myId.station()-1][myId.sector()-1] = counter;
    counter++;
  }
  for(int wheel = 0; wheel < 5; wheel++) {
    for(int sector = 0; sector < 12; sector++) {
      const std::vector<float> disp = theMuonModifier.gaussianRandomVector(dx, dy, dz);
      const std::vector<float> rotation = theMuonModifier.gaussianRandomVector(phix, phiy, phiz);
      for(int station = 0; station < 4; station++) {
        Alignable *myAlign = DTchambers.at(index[wheel][station][sector]);
        this->moveChamberInSector(myAlign, disp, rotation, errorDisp, errorRotation);
        if(sector == 3 && station == 3) {
	  Alignable *myAlignD = DTchambers.at(index[wheel][station][12]);
          this->moveChamberInSector(myAlignD, disp, rotation, errorDisp, errorRotation);
	} else if(sector == 9 && station == 3) {
	  Alignable *myAlignD = DTchambers.at(index[wheel][station][13]);
          this->moveChamberInSector(myAlignD, disp, rotation, errorDisp, errorRotation);
        }
      }
    }
  } 
}



//______________________________________________________________________________________________________
void MuonScenarioBuilder::moveCSCSectors(const edm::ParameterSet& pSet) {
  
  std::vector<Alignable *> CSCchambers = theAlignableMuon->CSCChambers();
  //Take Parameters
  std::vector<float> param = this->extractParameters(pSet, "CSCsectors");
  float scale_ = param[0]; float scaleError_ = param[1];
  float phiX_ = param[2]; float phiY_ = param[3]; float phiZ_ = param[4];
  float dX_ = param[5]; float dY_ = param[6]; float dZ_ = param[7];
  
  float dx = scale_*dX_; float dy = scale_*dY_; float dz = scale_*dZ_;
  float phix = scale_*phiX_; float phiy = scale_*phiY_; float phiz = scale_*phiZ_;
  float errorx = scaleError_*dX_; float errory = scaleError_*dY_; float errorz = scaleError_*dZ_;
  float errorphix = scaleError_*phiX_; float errorphiy = scaleError_*phiY_; float errorphiz = scaleError_*phiZ_;
  std::vector<float> errorDisp;
  errorDisp.push_back(errorx); errorDisp.push_back(errory); errorDisp.push_back(errorz);
  std::vector<float> errorRotation;
  errorRotation.push_back(errorphix); errorRotation.push_back(errorphiy); errorRotation.push_back(errorphiz);
  
  int index[2][4][4][36];
  int sector_index[2][4][4][36];
  int counter = 0;
  //Create an index for the chambers in the alignable vector
  for(std::vector<Alignable *>::iterator iter = CSCchambers.begin(); iter != CSCchambers.end(); ++iter) {
    CSCDetId myId((*iter)->geomDetId().rawId());
    index[myId.endcap()-1][myId.station()-1][myId.ring()-1][myId.chamber()-1] = counter;
    sector_index[myId.endcap()-1][myId.station()-1][myId.ring()-1][myId.chamber()-1] = CSCTriggerNumbering::sectorFromTriggerLabels(CSCTriggerNumbering::triggerSectorFromLabels(myId),CSCTriggerNumbering::triggerSubSectorFromLabels(myId) , myId.station());
    counter++;
  }
  for(int endcap = 0; endcap < 2; endcap++) {
    for(int ring = 0; ring < 2; ring++) {
      for(int sector = 1; sector < 7; sector++) {
	const std::vector<float> disp = theMuonModifier.gaussianRandomVector(dx, dy, dz);
	const std::vector<float> rotation = theMuonModifier.gaussianRandomVector(phix, phiy, phiz);
      	//Different cases are considered in order to fit endcap geometry
	for(int station = 0; station < 4; station++) {
	  if(station == 0) {
	    int r_ring[2];
	    if(ring == 0) {
	      r_ring[0] = 0; r_ring[1] = 3;
	    } else {
	      r_ring[0] = 1; r_ring[1] = 2;
	    }
	    for(int r_counter = 0; r_counter < 2; r_counter++) {
	      for(int chamber = 0; chamber < 36; chamber++) {
		if(sector == (sector_index[endcap][station][r_ring[r_counter]][chamber]+1)/2) {
		  Alignable *myAlign = CSCchambers.at(index[endcap][station][r_ring[r_counter]][chamber]);
                  this->moveChamberInSector(myAlign, disp, rotation, errorDisp, errorRotation);
		}
	      }
	    }
	  } else if(station == 3 && ring == 1) {
	    continue;
	  } else {
	    for(int chamber = 0; chamber < 36; chamber++) {
	      if(ring == 0 && chamber > 17) continue;
	      if(sector == sector_index[endcap][station][ring][chamber]) {
		Alignable *myAlign = CSCchambers.at(index[endcap][station][ring][chamber]);
                this->moveChamberInSector(myAlign, disp, rotation, errorDisp, errorRotation);
	      }
	    }
	  }
	}
      }
    }
  }
}


//______________________________________________________________________________________________________
void MuonScenarioBuilder::moveMuon(const edm::ParameterSet& pSet) {

  std::vector<Alignable *> DTbarrel = theAlignableMuon->DTBarrel();	
  std::vector<Alignable *> CSCendcaps = theAlignableMuon->CSCEndcaps();  
  //Take Parameters
  std::vector<float> param = this->extractParameters(pSet, "Muons");
  float scale_ = param[0]; float scaleError_ = param[1];
  float phiX_ = param[2]; float phiY_ = param[3]; float phiZ_ = param[4];
  float dX_ = param[5]; float dY_ = param[6]; float dZ_ = param[7];
  float dx = scale_*dX_; float dy = scale_*dY_; float dz = scale_*dZ_;
  float phix = scale_*phiX_; float phiy = scale_*phiY_; float phiz = scale_*phiZ_;
  float errorx = scaleError_*dX_; float errory = scaleError_*dY_; float errorz = scaleError_*dZ_;
  float errorphix = scaleError_*phiX_; float errorphiy = scaleError_*phiY_; float errorphiz = scaleError_*phiZ_;
  //Create an index for the chambers in the alignable vector
  const std::vector<float> disp = theMuonModifier.gaussianRandomVector(dx, dy, dz);
  const std::vector<float> rotation = theMuonModifier.gaussianRandomVector(phix, phiy, phiz);
  for(std::vector<Alignable *>::iterator iter = DTbarrel.begin(); iter != DTbarrel.end(); ++iter) {
    theMuonModifier.moveAlignable( *iter, false, true, disp[0], disp[1], disp[2] );
    theMuonModifier.rotateAlignable( *iter, false, true, rotation[0],  rotation[1], rotation[2] );
    theMuonModifier.addAlignmentPositionError( *iter, errorx, errory, errorz );
    theMuonModifier.addAlignmentPositionErrorFromRotation( *iter,  errorphix, errorphiy, errorphiz ); 
  }
  for(std::vector<Alignable *>::iterator iter = CSCendcaps.begin(); iter != CSCendcaps.end(); ++iter) {
    theMuonModifier.moveAlignable( *iter, false, true, disp[0], disp[1], disp[2] );
    theMuonModifier.rotateAlignable( *iter, false, true, rotation[0],  rotation[1], rotation[2] );
    theMuonModifier.addAlignmentPositionError( *iter, errorx, errory, errorz );
    theMuonModifier.addAlignmentPositionErrorFromRotation( *iter,  errorphix, errorphiy, errorphiz ); 
  }
  
}


//______________________________________________________________________________________________________
void MuonScenarioBuilder::moveChamberInSector(Alignable *chamber, std::vector<float> disp, std::vector<float>rotation, std::vector<float> dispError, std::vector<float> rotationError) {
   
    RotationType rotx( Basic3DVector<float>(1.0, 0.0, 0.0), rotation[0] );
    RotationType roty( Basic3DVector<float>(0.0, 1.0, 0.0), rotation[1] );
    RotationType rotz( Basic3DVector<float>(0.0, 0.0, 1.0), rotation[2] );
    RotationType rot = rotz * roty * rotx;
    GlobalPoint pos = chamber->globalPosition();
    GlobalPoint dispRot(pos.basicVector()-rot*pos.basicVector());
    disp[0] += dispRot.x(); disp[1] += dispRot.y(); disp[2] += dispRot.z();
    theMuonModifier.moveAlignable( chamber, false, true, disp[0], disp[1], disp[2] );
    theMuonModifier.rotateAlignable( chamber, false, true, rotation[0],  rotation[1], rotation[2] );
    theMuonModifier.addAlignmentPositionError( chamber, dispError[0], dispError[1], dispError[2] );
    theMuonModifier.addAlignmentPositionErrorFromRotation( chamber,  rotationError[0], rotationError[1], rotationError[2] );

}  
