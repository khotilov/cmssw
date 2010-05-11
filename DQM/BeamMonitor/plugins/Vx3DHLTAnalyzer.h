#ifndef Vx3DHLTAnalyzer_H
#define Vx3DHLTAnalyzer_H

// -*- C++ -*-
//
// Package:    Vx3DHLTAnalyzer
// Class:      Vx3DHLTAnalyzer
// 
/**\class Vx3DHLTAnalyzer Vx3DHLTAnalyzer.cc interface/Vx3DHLTAnalyzer.h

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mauro Dinardo,28 S-020,+41227673777,
//         Created:  Tue Feb 23 13:15:31 CET 2010
// $Id: Vx3DHLTAnalyzer.h,v 1.23 2010/04/25 07:09:08 dinardo Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include <TF3.h>

#include <iostream>
#include <fstream>
#include <vector>


using namespace std;


// #################
// # Fit variables #
// #################
#define DIM 3
void Gauss3DFunc(int& /*npar*/, double* /*gin*/, double& fval, double* par, int /*iflag*/);
typedef struct
{
  double x;
  double y;
  double z;
  double Covariance[DIM][DIM];
} VertexType;
vector<VertexType> Vertices;
bool considerVxCovariance;
unsigned int counterVx; // Counts the number of vertices taken into account for the fit
double maxTransRadius;  // Max transverse radius in which the vertices must be [cm]
double maxLongLength;   // Max longitudinal length in which the vertices must be [cm]
double xPos,yPos,zPos;  // x,y,z approximate positions of the beam spot
double pi;
// ######################
// # cfg file parameter #
// ######################
double VxErrCorr;       // Coefficient to compensate the under-estimation of the vertex errors


class Vx3DHLTAnalyzer : public edm::EDAnalyzer {
   public:
      explicit Vx3DHLTAnalyzer(const edm::ParameterSet&);
      ~Vx3DHLTAnalyzer();


   private:
      virtual void beginJob();
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual unsigned int HitCounter(const edm::Event& iEvent);
      virtual char* formatTime(const time_t& t);
      virtual int MyFit(vector<double>* vals);
      virtual void reset(string ResetType);
      virtual void writeToFile(vector<double>* vals,
			       edm::TimeValue_t BeginTimeOfFit,
			       edm::TimeValue_t EndTimeOfFit,
			       unsigned int BeginLumiOfFit,
			       unsigned int EndLumiOfFit,
			       int dataType);
      virtual void beginLuminosityBlock(const edm::LuminosityBlock& lumiBlock, 
					const edm::EventSetup& iSetup);
      virtual void endLuminosityBlock(const edm::LuminosityBlock& lumiBlock,
				      const edm::EventSetup& iSetup);
      virtual void endJob();


      // #######################
      // # cfg file parameters #
      // #######################
      edm::InputTag vertexCollection;
      bool debugMode;
      unsigned int nLumiReset;
      bool dataFromFit;
      unsigned int minNentries;
      double xRange;
      double xStep;
      double yRange;
      double yStep;
      double zRange;
      double zStep;
      string fileName;


      // ##############
      // # Histograms #
      // ##############
      MonitorElement* mXlumi;
      MonitorElement* mYlumi;
      MonitorElement* mZlumi;
      
      MonitorElement* sXlumi;
      MonitorElement* sYlumi;
      MonitorElement* sZlumi;
      
      MonitorElement* dxdzlumi;
      MonitorElement* dydzlumi;

      MonitorElement* Vx_X;
      MonitorElement* Vx_Y;
      MonitorElement* Vx_Z;
      
      MonitorElement* Vx_ZX;
      MonitorElement* Vx_ZY;
      MonitorElement* Vx_XY;
      
      MonitorElement* Vx_ZX_profile;
      MonitorElement* Vx_ZY_profile;

      MonitorElement* goodVxCounter;

      MonitorElement* fitResults;

      MonitorElement* hitCounter;

      MonitorElement* reportSummary;
      MonitorElement* reportSummaryMap;
      

      // ######################
      // # Internal variables #
      // ######################
      ofstream outputFile;
      ofstream outputDebugFile;
      edm::TimeValue_t beginTimeOfFit;
      edm::TimeValue_t endTimeOfFit;
      int nBinsHistoricalPlot;
      unsigned int runNumber;
      unsigned int lumiCounter;
      unsigned int lumiCounterHisto;
      unsigned int totalHits;
      unsigned int maxLumiIntegration;
      unsigned int numberGoodFits;
      unsigned int numberFits;
      unsigned int beginLumiOfFit;
      unsigned int endLumiOfFit;
      unsigned int lastLumiOfFit;
      double minVxDoF;
      double minVxWgt;
      bool internalDebug;
};

#endif
