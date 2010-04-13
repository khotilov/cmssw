// -*- C++ -*-
//
// Package:    Vx3DHLTAnalyzer
// Class:      Vx3DHLTAnalyzer
// 
/**\class Vx3DHLTAnalyzer Vx3DHLTAnalyzer.cc plugins/Vx3DHLTAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mauro Dinardo,28 S-020,+41227673777,
//         Created:  Tue Feb 23 13:15:31 CET 2010
// $Id: Vx3DHLTAnalyzer.cc,v 1.62 2010/04/07 09:46:18 dinardo Exp $
//
//


#include "DQM/BeamMonitor/interface/Vx3DHLTAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <TFitterMinuit.h>


using namespace std;
using namespace reco;
using namespace edm;


Vx3DHLTAnalyzer::Vx3DHLTAnalyzer(const ParameterSet& iConfig)
{
  vertexCollection = edm::InputTag("pixelVertices");
  debugMode        = false;
  nLumiReset       = 1;
  dataFromFit      = true;
  minNentries      = 35;
  xRange           = 4.;
  xStep            = 0.001;
  yRange           = 4.;
  yStep            = 0.001;
  zRange           = 40.;
  zStep            = 0.05;
  fileName         = "BeamPixelResults.txt";

  vertexCollection = iConfig.getParameter<InputTag>("vertexCollection");
  debugMode        = iConfig.getParameter<bool>("debugMode");
  nLumiReset       = iConfig.getParameter<unsigned int>("nLumiReset");
  dataFromFit      = iConfig.getParameter<bool>("dataFromFit");
  minNentries      = iConfig.getParameter<int>("minNentries");
  xRange           = iConfig.getParameter<double>("xRange");
  xStep            = iConfig.getParameter<double>("xStep");
  yRange           = iConfig.getParameter<double>("yRange");
  yStep            = iConfig.getParameter<double>("yStep");
  zRange           = iConfig.getParameter<double>("zRange");
  zStep            = iConfig.getParameter<double>("zStep");
  fileName         = iConfig.getParameter<string>("fileName");
}


Vx3DHLTAnalyzer::~Vx3DHLTAnalyzer()
{
}


void Vx3DHLTAnalyzer::analyze(const Event& iEvent, const EventSetup& iSetup)
{
  Handle<VertexCollection> Vx3DCollection;
  iEvent.getByLabel(vertexCollection,Vx3DCollection);

  unsigned int i,j;
  double det;
  VertexType MyVertex;

  if (runNumber != iEvent.id().run())
    {
      reset("scratch");
      runNumber = iEvent.id().run();

      if (debugMode == true)
	{
	  stringstream debugFile;
	  string tmp(fileName);

	  if (outputDebugFile.is_open() == true) outputDebugFile.close();
	  tmp.erase(strlen(fileName.c_str())-4,4);
	  debugFile << tmp.c_str() << "_Run" << iEvent.id().run() << ".txt";
	  outputDebugFile.open(debugFile.str().c_str(), ios::out);
	  outputDebugFile.close();
	  outputDebugFile.open(debugFile.str().c_str(), ios::app);
	}
    }
  else if (beginTimeOfFit != 0)
    {
      totalHits += HitCounter(iEvent);

      for (vector<Vertex>::const_iterator it3DVx = Vx3DCollection->begin(); it3DVx != Vx3DCollection->end(); it3DVx++) {
	
	if ((it3DVx->isValid() == true) && (it3DVx->isFake() == false) && (it3DVx->ndof() >= minVxDoF))
	  {
	    for (i = 0; i < DIM; i++)
	      {
		for (j = 0; j < DIM; j++)
		  {
		    MyVertex.Covariance[i][j] = it3DVx->covariance(i,j);
		    if (isnan(MyVertex.Covariance[i][j]) == true) break;
		  }
		if (j != DIM) break;
	      }
	    det = fabs(MyVertex.Covariance[0][0])*(fabs(MyVertex.Covariance[1][1])*fabs(MyVertex.Covariance[2][2]) - MyVertex.Covariance[1][2]*MyVertex.Covariance[1][2]) -
	      MyVertex.Covariance[0][1]*(MyVertex.Covariance[0][1]*fabs(MyVertex.Covariance[2][2]) - MyVertex.Covariance[0][2]*MyVertex.Covariance[1][2]) +
	      MyVertex.Covariance[0][2]*(MyVertex.Covariance[0][1]*MyVertex.Covariance[1][2] - MyVertex.Covariance[0][2]*fabs(MyVertex.Covariance[1][1]));
	    if ((i == DIM) && (det > 0.))
	      {
		MyVertex.x = it3DVx->x();
		MyVertex.y = it3DVx->y();
		MyVertex.z = it3DVx->z();
		Vertices.push_back(MyVertex);
	      }
	    else if (internalDebug == true)
	      {
		cout << "Vertex discarded !" << endl;
		for (i = 0; i < DIM; i++)
		  for (j = 0; j < DIM; j++)
		    cout << "(i,j) --> " << i << "," << j << " --> " << MyVertex.Covariance[i][j] << endl;
	      }
	    
	    Vx_X->Fill(it3DVx->x());
	    Vx_Y->Fill(it3DVx->y());
	    Vx_Z->Fill(it3DVx->z());
	    
	    Vx_ZX->Fill(it3DVx->z(), it3DVx->x());
	    Vx_ZY->Fill(it3DVx->z(), it3DVx->y());
	    Vx_XY->Fill(it3DVx->x(), it3DVx->y());
	   
	    Vx_ZX_profile->Fill(it3DVx->z(), it3DVx->x());
	    Vx_ZY_profile->Fill(it3DVx->z(), it3DVx->y());
	  }
      }      
    }
}


unsigned int Vx3DHLTAnalyzer::HitCounter(const Event& iEvent)
{
  edm::Handle<SiPixelRecHitCollection> rechitspixel;
  iEvent.getByLabel("siPixelRecHits",rechitspixel);

  unsigned int counter = 0;
  
  for (SiPixelRecHitCollection::const_iterator j = rechitspixel->begin(); j != rechitspixel->end(); j++)
    for (edmNew::DetSet<SiPixelRecHit>::const_iterator h = j->begin(); h != j->end(); h++) counter += h->cluster()->size();	  
  
  return counter;
}


char* Vx3DHLTAnalyzer::formatTime(const time_t t)
{
  static char ts[] = "yyyy.Mm.dd hh:mm:ss TZN     ";
  strftime(ts, strlen(ts)+1, "%Y.%m.%d %H:%M:%S %Z", gmtime(&t));

#ifdef STRIP_TRAILING_BLANKS_IN_TIMEZONE
  // Strip trailing blanks that would come when the time zone is not as
  // long as the maximum allowed
  unsigned int b = strlen(ts);
  while (ts[--b] == ' ') ts[b] = 0;
#endif
  
  return ts;
}


void Gauss3DFunc(int& /*npar*/, double* /*gin*/, double& fval, double* par, int /*iflag*/)
{
  double K[DIM][DIM]; // Covariance Matrix
  double M[DIM][DIM]; // K^-1
  double coef,det;
  double sumlog = 0.;
  double precision = 1.e-9;

//   par[0] = K(0,0)
//   par[1] = K(1,1)
//   par[2] = K(2,2)
//   par[3] = K(0,1) = K(1,0)
//   par[4] = K(1,2) = K(2,1)
//   par[5] = K(0,2) = K(2,0)
//   par[6] = mean x
//   par[7] = mean y
//   par[8] = mean z

  counterVx = 0;
  for (unsigned int i = 0; i < Vertices.size(); i++)
    {
      if ((sqrt((Vertices[i].x-xPos)*(Vertices[i].x-xPos) + (Vertices[i].y-yPos)*(Vertices[i].y-yPos)) <= maxTransRadius) &&
	  (fabs(Vertices[i].z-zPos) <= maxLongLength))
	{
	  if (considerVxCovariance == true)
	    {
	      K[0][0] = fabs(par[0]) + fabs(Vertices[i].Covariance[0][0]);
	      K[1][1] = fabs(par[1]) + fabs(Vertices[i].Covariance[1][1]);
	      K[2][2] = fabs(par[2]) + fabs(Vertices[i].Covariance[2][2]);
	      K[0][1] = K[1][0] = par[3] + Vertices[i].Covariance[0][1];
	      K[1][2] = K[2][1] = par[4] + Vertices[i].Covariance[1][2];
	      K[0][2] = K[2][0] = par[5] + Vertices[i].Covariance[0][2];
	    }
	  else
	    {
	      K[0][0] = fabs(par[0]);
	      K[1][1] = fabs(par[1]);
	      K[2][2] = fabs(par[2]);
	      K[0][1] = K[1][0] = par[3];
	      K[1][2] = K[2][1] = par[4];
	      K[0][2] = K[2][0] = par[5];
	    }

	  det = K[0][0]*(K[1][1]*K[2][2] - K[1][2]*K[1][2]) -
		K[0][1]*(K[0][1]*K[2][2] - K[0][2]*K[1][2]) +
		K[0][2]*(K[0][1]*K[1][2] - K[0][2]*K[1][1]);

	  M[0][0] = (K[1][1]*K[2][2] - K[1][2]*K[1][2]) / det;
	  M[1][1] = (K[0][0]*K[2][2] - K[0][2]*K[0][2]) / det;
	  M[2][2] = (K[0][0]*K[1][1] - K[0][1]*K[0][1]) / det;
	  M[0][1] = M[1][0] = (K[0][2]*K[1][2] - K[0][1]*K[2][2]) / det;
	  M[1][2] = M[2][1] = (K[0][2]*K[0][1] - K[1][2]*K[0][0]) / det;
	  M[0][2] = M[2][0] = (K[0][1]*K[1][2] - K[0][2]*K[1][1]) / det;

	  coef = 1. / sqrt(powf(2.*pi,DIM)*fabs(det));
	  
	  if ((coef * exp(-1./2. * (M[0][0]*(Vertices[i].x-par[6])*(Vertices[i].x-par[6]) +
				    M[1][1]*(Vertices[i].y-par[7])*(Vertices[i].y-par[7]) +
				    M[2][2]*(Vertices[i].z-par[8])*(Vertices[i].z-par[8]) +
				    2.*M[0][1]*(Vertices[i].x-par[6])*(Vertices[i].y-par[7]) +
				    2.*M[1][2]*(Vertices[i].y-par[7])*(Vertices[i].z-par[8]) +
				    2.*M[0][2]*(Vertices[i].x-par[6])*(Vertices[i].z-par[8])))) >= precision)
// 	    sumlog += double(DIM)*log(2.*pi) + log(fabs(det)) +
	    sumlog += log(fabs(det)) +
	      (M[0][0]*(Vertices[i].x-par[6])*(Vertices[i].x-par[6]) +
	       M[1][1]*(Vertices[i].y-par[7])*(Vertices[i].y-par[7]) +
	       M[2][2]*(Vertices[i].z-par[8])*(Vertices[i].z-par[8]) +
	       2.*M[0][1]*(Vertices[i].x-par[6])*(Vertices[i].y-par[7]) +
	       2.*M[1][2]*(Vertices[i].y-par[7])*(Vertices[i].z-par[8]) +
	       2.*M[0][2]*(Vertices[i].x-par[6])*(Vertices[i].z-par[8]));
	  else sumlog += -2.*log(precision);
 
	  counterVx++;
	}
    }

  fval = sumlog;
}


int Vx3DHLTAnalyzer::MyFit(vector<double>* vals)
{
  // RETURN CODE:
  //  0 == OK
  // -2 == NO OK - not enough "minNentries"
  // Any other number == NO OK
  unsigned int nParams = 9;
 
  if ((vals != NULL) && (vals->size() == nParams*2))
    {
      double nSigmaXY    = 4.;
      double nSigmaZ     = 4.;
      double varFactor   = 2./5.; // Take into account the difference between the RMS and sigma (RMS usually greater than sigma)
      double parDistance = 0.01;
      double det;
      double bestEdm = 1.;
      double deltaMean;
      int bestMovementX = 1;
      int bestMovementY = 1;
      int bestMovementZ = 1;
      int goodData;

      double arglist[2];
      double amin,errdef,edm;
      int nvpar,nparx;
      
      vector<double>::const_iterator it = vals->begin();

      TFitterMinuit* Gauss3D = new TFitterMinuit(nParams);
      if (internalDebug == true) Gauss3D->SetPrintLevel(3);
      else Gauss3D->SetPrintLevel(0);
      // 	  Gauss3D->SetStrategy(0);
      Gauss3D->SetFCN(Gauss3DFunc);
      arglist[0] = 10000; // Max number of function calls
      arglist[1] = 1e-9;  // Tolerance on likelihood

      if (internalDebug == true) cout << "\n@@@ START FITTING @@@" << endl;

      // @@@ Fit at X-deltaMean | X | X+deltaMean @@@
      bestEdm = 1.;
      for (int i = 0; i < 3; i++)
	{
	  deltaMean = (double(i)-1.)*sqrt((*(it+0))*varFactor);
	  if (internalDebug == true) cout << "deltaMean --> " << deltaMean << endl;

	  Gauss3D->Clear();

	  // arg3 - first guess of parameter value
	  // arg4 - step of the parameter
	  Gauss3D->SetParameter(0,"var x ", *(it+0)*varFactor, parDistance, 0., 0.);
	  Gauss3D->SetParameter(1,"var y ", *(it+1)*varFactor, parDistance, 0., 0.);
	  Gauss3D->SetParameter(2,"var z ", *(it+2)*varFactor, parDistance, 0., 0.);
	  Gauss3D->SetParameter(3,"cov xy", *(it+3), parDistance, 0., 0.);
	  Gauss3D->SetParameter(4,"cov yz", *(it+4), parDistance, 0., 0.);
	  Gauss3D->SetParameter(5,"cov xz", *(it+5), parDistance, 0., 0.);
	  Gauss3D->SetParameter(6,"mean x", *(it+6)+deltaMean, parDistance, 0., 0.);
	  Gauss3D->SetParameter(7,"mean y", *(it+7), parDistance, 0., 0.);
	  Gauss3D->SetParameter(8,"mean z", *(it+8), parDistance, 0., 0.);

	  // Set the central positions of the centroid for vertex rejection
	  xPos = Gauss3D->GetParameter(6);
	  yPos = Gauss3D->GetParameter(7);
	  zPos = Gauss3D->GetParameter(8);

	  // Set dimensions of the centroid for vertex rejection
	  maxTransRadius = nSigmaXY * sqrt(fabs(Gauss3D->GetParameter(0)) + fabs(Gauss3D->GetParameter(1)));
	  maxLongLength  = nSigmaZ  * sqrt(fabs(Gauss3D->GetParameter(2)));
	  goodData = Gauss3D->ExecuteCommand("MIGRAD",arglist,2);
	  // Re-set dimensions of the centroid for vertex rejection
	  maxTransRadius = nSigmaXY * sqrt(fabs(Gauss3D->GetParameter(0)) + fabs(Gauss3D->GetParameter(1)));
	  maxLongLength  = nSigmaZ  * sqrt(fabs(Gauss3D->GetParameter(2)));
	  goodData = Gauss3D->ExecuteCommand("MIGRAD",arglist,2);

	  Gauss3D->GetStats(amin, edm, errdef, nvpar, nparx);

	  if (counterVx < minNentries) goodData = -2;
	  else if (isnan(edm) == true) goodData = -1;
	  else for (unsigned int j = 0; j < nParams; j++) if (isnan(Gauss3D->GetParError(j)) == true) { goodData = -1; break; }
	  if (goodData == 0)
	    {
	      det = fabs(Gauss3D->GetParameter(0)) * (fabs(Gauss3D->GetParameter(1))*fabs(Gauss3D->GetParameter(2)) - Gauss3D->GetParameter(4)*Gauss3D->GetParameter(4)) -
		Gauss3D->GetParameter(3) * (Gauss3D->GetParameter(3)*fabs(Gauss3D->GetParameter(2)) - Gauss3D->GetParameter(5)*Gauss3D->GetParameter(4)) +
		Gauss3D->GetParameter(5) * (Gauss3D->GetParameter(3)*Gauss3D->GetParameter(4) - Gauss3D->GetParameter(5)*fabs(Gauss3D->GetParameter(1)));
	      if (det < 0.) goodData = -1;
	    }

	  if ((goodData == 0) && (edm < bestEdm)) { bestEdm = edm; bestMovementX = i; }
	}
      if (internalDebug == true) cout << "Found bestMovementX --> " << bestMovementX << endl;

      // @@@ Fit at Y-deltaMean | Y | Y+deltaMean @@@
      bestEdm = 1.;
      for (int i = 0; i < 3; i++)
	{
	  deltaMean = (double(i)-1.)*sqrt((*(it+1))*varFactor);
	  if (internalDebug == true)
	    {
	      cout << "deltaMean --> " << deltaMean << endl;
	      cout << "deltaMean X --> " << (double(bestMovementX)-1.)*sqrt((*(it+0))*varFactor) << endl;
	    }

	  Gauss3D->Clear();

	  // arg3 - first guess of parameter value
	  // arg4 - step of the parameter
	  Gauss3D->SetParameter(0,"var x ", *(it+0)*varFactor, parDistance, 0., 0.);
	  Gauss3D->SetParameter(1,"var y ", *(it+1)*varFactor, parDistance, 0., 0.);
	  Gauss3D->SetParameter(2,"var z ", *(it+2)*varFactor, parDistance, 0., 0.);
	  Gauss3D->SetParameter(3,"cov xy", *(it+3), parDistance, 0., 0.);
	  Gauss3D->SetParameter(4,"cov yz", *(it+4), parDistance, 0., 0.);
	  Gauss3D->SetParameter(5,"cov xz", *(it+5), parDistance, 0., 0.);
	  Gauss3D->SetParameter(6,"mean x", *(it+6)+(double(bestMovementX)-1.)*sqrt((*(it+0))*varFactor), parDistance, 0., 0.);
	  Gauss3D->SetParameter(7,"mean y", *(it+7)+deltaMean, parDistance, 0., 0.);
	  Gauss3D->SetParameter(8,"mean z", *(it+8), parDistance, 0., 0.);

	  // Set the central positions of the centroid for vertex rejection
	  xPos = Gauss3D->GetParameter(6);
	  yPos = Gauss3D->GetParameter(7);
	  zPos = Gauss3D->GetParameter(8);

	  // Set dimensions of the centroid for vertex rejection
	  maxTransRadius = nSigmaXY * sqrt(fabs(Gauss3D->GetParameter(0)) + fabs(Gauss3D->GetParameter(1)));
	  maxLongLength  = nSigmaZ  * sqrt(fabs(Gauss3D->GetParameter(2)));
	  goodData = Gauss3D->ExecuteCommand("MIGRAD",arglist,2);
	  // Re-set dimensions of the centroid for vertex rejection
	  maxTransRadius = nSigmaXY * sqrt(fabs(Gauss3D->GetParameter(0)) + fabs(Gauss3D->GetParameter(1)));
	  maxLongLength  = nSigmaZ  * sqrt(fabs(Gauss3D->GetParameter(2)));
	  goodData = Gauss3D->ExecuteCommand("MIGRAD",arglist,2);

	  Gauss3D->GetStats(amin, edm, errdef, nvpar, nparx);

	  if (counterVx < minNentries) goodData = -2;
	  else if (isnan(edm) == true) goodData = -1;
	  else for (unsigned int j = 0; j < nParams; j++) if (isnan(Gauss3D->GetParError(j)) == true) { goodData = -1; break; }
	  if (goodData == 0)
	    {
	      det = fabs(Gauss3D->GetParameter(0)) * (fabs(Gauss3D->GetParameter(1))*fabs(Gauss3D->GetParameter(2)) - Gauss3D->GetParameter(4)*Gauss3D->GetParameter(4)) -
		Gauss3D->GetParameter(3) * (Gauss3D->GetParameter(3)*fabs(Gauss3D->GetParameter(2)) - Gauss3D->GetParameter(5)*Gauss3D->GetParameter(4)) +
		Gauss3D->GetParameter(5) * (Gauss3D->GetParameter(3)*Gauss3D->GetParameter(4) - Gauss3D->GetParameter(5)*fabs(Gauss3D->GetParameter(1)));
	      if (det < 0.) goodData = -1;
	    }
	  
	  if ((goodData == 0) && (edm < bestEdm)) { bestEdm = edm; bestMovementY = i; }
	}
      if (internalDebug == true) cout << "Found bestMovementY --> " << bestMovementY << endl;

      // @@@ Fit at Z-deltaMean | Z | Z+deltaMean @@@
      bestEdm = 1.;
      for (int i = 0; i < 3; i++)
	{
	  deltaMean = (double(i)-1.)*sqrt((*(it+2))*varFactor);
	  if (internalDebug == true)
	    {
	      cout << "deltaMean --> " << deltaMean << endl;
	      cout << "deltaMean X --> " << (double(bestMovementX)-1.)*sqrt((*(it+0))*varFactor) << endl;
	      cout << "deltaMean Y --> " << (double(bestMovementY)-1.)*sqrt((*(it+1))*varFactor) << endl;
	    }

	  Gauss3D->Clear();

	  // arg3 - first guess of parameter value
	  // arg4 - step of the parameter
	  Gauss3D->SetParameter(0,"var x ", *(it+0)*varFactor, parDistance, 0., 0.);
	  Gauss3D->SetParameter(1,"var y ", *(it+1)*varFactor, parDistance, 0., 0.);
	  Gauss3D->SetParameter(2,"var z ", *(it+2)*varFactor, parDistance, 0., 0.);
	  Gauss3D->SetParameter(3,"cov xy", *(it+3), parDistance, 0., 0.);
	  Gauss3D->SetParameter(4,"cov yz", *(it+4), parDistance, 0., 0.);
	  Gauss3D->SetParameter(5,"cov xz", *(it+5), parDistance, 0., 0.);
	  Gauss3D->SetParameter(6,"mean x", *(it+6)+(double(bestMovementX)-1.)*sqrt((*(it+0))*varFactor), parDistance, 0., 0.);
	  Gauss3D->SetParameter(7,"mean y", *(it+7)+(double(bestMovementY)-1.)*sqrt((*(it+1))*varFactor), parDistance, 0., 0.);
	  Gauss3D->SetParameter(8,"mean z", *(it+8)+deltaMean, parDistance, 0., 0.);

	  // Set the central positions of the centroid for vertex rejection
	  xPos = Gauss3D->GetParameter(6);
	  yPos = Gauss3D->GetParameter(7);
	  zPos = Gauss3D->GetParameter(8);

	  // Set dimensions of the centroid for vertex rejection
	  maxTransRadius = nSigmaXY * sqrt(fabs(Gauss3D->GetParameter(0)) + fabs(Gauss3D->GetParameter(1)));
	  maxLongLength  = nSigmaZ  * sqrt(fabs(Gauss3D->GetParameter(2)));
	  goodData = Gauss3D->ExecuteCommand("MIGRAD",arglist,2);
	  // Re-set dimensions of the centroid for vertex rejection
	  maxTransRadius = nSigmaXY * sqrt(fabs(Gauss3D->GetParameter(0)) + fabs(Gauss3D->GetParameter(1)));
	  maxLongLength  = nSigmaZ  * sqrt(fabs(Gauss3D->GetParameter(2)));
	  goodData = Gauss3D->ExecuteCommand("MIGRAD",arglist,2);

	  Gauss3D->GetStats(amin, edm, errdef, nvpar, nparx);

	  if (counterVx < minNentries) goodData = -2;
	  else if (isnan(edm) == true) goodData = -1;
	  else for (unsigned int j = 0; j < nParams; j++) if (isnan(Gauss3D->GetParError(j)) == true) { goodData = -1; break; }
	  if (goodData == 0)
	    {
	      det = fabs(Gauss3D->GetParameter(0)) * (fabs(Gauss3D->GetParameter(1))*fabs(Gauss3D->GetParameter(2)) - Gauss3D->GetParameter(4)*Gauss3D->GetParameter(4)) -
		Gauss3D->GetParameter(3) * (Gauss3D->GetParameter(3)*fabs(Gauss3D->GetParameter(2)) - Gauss3D->GetParameter(5)*Gauss3D->GetParameter(4)) +
		Gauss3D->GetParameter(5) * (Gauss3D->GetParameter(3)*Gauss3D->GetParameter(4) - Gauss3D->GetParameter(5)*fabs(Gauss3D->GetParameter(1)));
	      if (det < 0.) goodData = -1;
	    }
	  
	  if ((goodData == 0) && (edm < bestEdm)) { bestEdm = edm; bestMovementZ = i; }
	}
      if (internalDebug == true) cout << "Found bestMovementZ --> " << bestMovementZ << endl;

      Gauss3D->Clear();

      // @@@ FINAL FIT @@@
      // arg3 - first guess of parameter value
      // arg4 - step of the parameter
      Gauss3D->SetParameter(0,"var x ", *(it+0)*varFactor, parDistance, 0., 0.);
      Gauss3D->SetParameter(1,"var y ", *(it+1)*varFactor, parDistance, 0., 0.);
      Gauss3D->SetParameter(2,"var z ", *(it+2)*varFactor, parDistance, 0., 0.);
      Gauss3D->SetParameter(3,"cov xy", *(it+3), parDistance, 0., 0.);
      Gauss3D->SetParameter(4,"cov yz", *(it+4), parDistance, 0., 0.);
      Gauss3D->SetParameter(5,"cov xz", *(it+5), parDistance, 0., 0.);
      Gauss3D->SetParameter(6,"mean x", *(it+6)+(double(bestMovementX)-1.)*sqrt((*(it+0))*varFactor), parDistance, 0., 0.);
      Gauss3D->SetParameter(7,"mean y", *(it+7)+(double(bestMovementY)-1.)*sqrt((*(it+1))*varFactor), parDistance, 0., 0.);
      Gauss3D->SetParameter(8,"mean z", *(it+8)+(double(bestMovementZ)-1.)*sqrt((*(it+2))*varFactor), parDistance, 0., 0.);

      // Set the central positions of the centroid for vertex rejection
      xPos = Gauss3D->GetParameter(6);
      yPos = Gauss3D->GetParameter(7);
      zPos = Gauss3D->GetParameter(8);
      
      // Set dimensions of the centroid for vertex rejection
      maxTransRadius = nSigmaXY * sqrt(fabs(Gauss3D->GetParameter(0)) + fabs(Gauss3D->GetParameter(1)));
      maxLongLength  = nSigmaZ  * sqrt(fabs(Gauss3D->GetParameter(2)));
      goodData = Gauss3D->ExecuteCommand("MIGRAD",arglist,2);
      // Re-set dimensions of the centroid for vertex rejection
      maxTransRadius = nSigmaXY * sqrt(fabs(Gauss3D->GetParameter(0)) + fabs(Gauss3D->GetParameter(1)));
      maxLongLength  = nSigmaZ  * sqrt(fabs(Gauss3D->GetParameter(2)));
      goodData = Gauss3D->ExecuteCommand("MIGRAD",arglist,2);
      
      Gauss3D->GetStats(amin, edm, errdef, nvpar, nparx);
      
      if (counterVx < minNentries) goodData = -2;
      else if (isnan(edm) == true) goodData = -1;
      else for (unsigned int j = 0; j < nParams; j++) if (isnan(Gauss3D->GetParError(j)) == true) { goodData = -1; break; }
      if (goodData == 0)
	{
	  det = fabs(Gauss3D->GetParameter(0)) * (fabs(Gauss3D->GetParameter(1))*fabs(Gauss3D->GetParameter(2)) - Gauss3D->GetParameter(4)*Gauss3D->GetParameter(4)) -
	    Gauss3D->GetParameter(3) * (Gauss3D->GetParameter(3)*fabs(Gauss3D->GetParameter(2)) - Gauss3D->GetParameter(5)*Gauss3D->GetParameter(4)) +
	    Gauss3D->GetParameter(5) * (Gauss3D->GetParameter(3)*Gauss3D->GetParameter(4) - Gauss3D->GetParameter(5)*fabs(Gauss3D->GetParameter(1)));
	  if (det < 0.) goodData = -1;
	}
      if ((goodData != 0) && (goodData != -2))
	{
	  Gauss3D->Clear();
	  
	  if (internalDebug == true) cout << "FIT WITH ENLARGED PARAMETER DISTANCES - STEP 1" << endl;      
	  // @@@ FIT WITH ENLARGED PARAMETER DISTANCES - STEP 1 @@@
	  // arg3 - first guess of parameter value
	  // arg4 - step of the parameter
	  Gauss3D->SetParameter(0,"var x ", *(it+0)*varFactor, parDistance*5., 0, 0);
	  Gauss3D->SetParameter(1,"var y ", *(it+1)*varFactor, parDistance*5., 0, 0);
	  Gauss3D->SetParameter(2,"var z ", *(it+2)*varFactor, parDistance*100., 0, 0);
	  Gauss3D->SetParameter(3,"cov xy", *(it+3), parDistance*5., 0, 0);
	  Gauss3D->SetParameter(4,"cov yz", *(it+4), parDistance*5., 0, 0);
	  Gauss3D->SetParameter(5,"cov xz", *(it+5), parDistance*5., 0, 0);
	  Gauss3D->SetParameter(6,"mean x", *(it+6)+(double(bestMovementX)-1.)*sqrt((*(it+0))*varFactor), parDistance*5., 0, 0);
	  Gauss3D->SetParameter(7,"mean y", *(it+7)+(double(bestMovementY)-1.)*sqrt((*(it+1))*varFactor), parDistance*5., 0, 0);
	  Gauss3D->SetParameter(8,"mean z", *(it+8)+(double(bestMovementZ)-1.)*sqrt((*(it+2))*varFactor), parDistance*50., 0, 0);

	  // Set the central positions of the centroid for vertex rejection
	  xPos = Gauss3D->GetParameter(6);
	  yPos = Gauss3D->GetParameter(7);
	  zPos = Gauss3D->GetParameter(8);

	  // Set dimensions of the centroid for vertex rejection
	  maxTransRadius = nSigmaXY * sqrt(fabs(Gauss3D->GetParameter(0)) + fabs(Gauss3D->GetParameter(1)));
	  maxLongLength  = nSigmaZ  * sqrt(fabs(Gauss3D->GetParameter(2)));
	  goodData = Gauss3D->ExecuteCommand("MIGRAD",arglist,2);
	  // Re-set dimensions of the centroid for vertex rejection
	  maxTransRadius = nSigmaXY * sqrt(fabs(Gauss3D->GetParameter(0)) + fabs(Gauss3D->GetParameter(1)));
	  maxLongLength  = nSigmaZ  * sqrt(fabs(Gauss3D->GetParameter(2)));
	  goodData = Gauss3D->ExecuteCommand("MIGRAD",arglist,2);

	  Gauss3D->GetStats(amin, edm, errdef, nvpar, nparx);
      
	  if (counterVx < minNentries) goodData = -2;
	  else if (isnan(edm) == true) goodData = -1;
	  else for (unsigned int j = 0; j < nParams; j++) if (isnan(Gauss3D->GetParError(j)) == true) { goodData = -1; break; }
	  if (goodData == 0)
	    {
	      det = fabs(Gauss3D->GetParameter(0)) * (fabs(Gauss3D->GetParameter(1))*fabs(Gauss3D->GetParameter(2)) - Gauss3D->GetParameter(4)*Gauss3D->GetParameter(4)) -
		Gauss3D->GetParameter(3) * (Gauss3D->GetParameter(3)*fabs(Gauss3D->GetParameter(2)) - Gauss3D->GetParameter(5)*Gauss3D->GetParameter(4)) +
		Gauss3D->GetParameter(5) * (Gauss3D->GetParameter(3)*Gauss3D->GetParameter(4) - Gauss3D->GetParameter(5)*fabs(Gauss3D->GetParameter(1)));
	      if (det < 0.) goodData = -1;
	    }
	  if ((goodData != 0) && (goodData != -2))
	    {
	      Gauss3D->Clear();
	  
	      if (internalDebug == true) cout << "FIT WITH ENLARGED PARAMETER DISTANCES - STEP 2" << endl;      
	      // @@@ FIT WITH ENLARGED PARAMETER DISTANCES - STEP 2 @@@
	      // arg3 - first guess of parameter value
	      // arg4 - step of the parameter
	      Gauss3D->SetParameter(0,"var x ", *(it+0)*varFactor, parDistance*10., 0, 0);
	      Gauss3D->SetParameter(1,"var y ", *(it+1)*varFactor, parDistance*10., 0, 0);
	      Gauss3D->SetParameter(2,"var z ", *(it+2)*varFactor, parDistance*500., 0, 0);
	      Gauss3D->SetParameter(3,"cov xy", 0.0, parDistance*10., 0, 0);
	      Gauss3D->SetParameter(4,"cov yz", 0.0, parDistance*10., 0, 0);
	      Gauss3D->SetParameter(5,"cov xz", 0.0, parDistance*10., 0, 0);
	      Gauss3D->SetParameter(6,"mean x", *(it+6)+(double(bestMovementX)-1.)*sqrt((*(it+0))*varFactor), parDistance*10., 0, 0);
	      Gauss3D->SetParameter(7,"mean y", *(it+7)+(double(bestMovementY)-1.)*sqrt((*(it+1))*varFactor), parDistance*10., 0, 0);
	      Gauss3D->SetParameter(8,"mean z", *(it+8)+(double(bestMovementZ)-1.)*sqrt((*(it+2))*varFactor), parDistance*100., 0, 0);

	      // Set the central positions of the centroid for vertex rejection
	      xPos = Gauss3D->GetParameter(6);
	      yPos = Gauss3D->GetParameter(7);
	      zPos = Gauss3D->GetParameter(8);
	      
	      // Set dimensions of the centroid for vertex rejection
	      maxTransRadius = nSigmaXY * sqrt(fabs(Gauss3D->GetParameter(0)) + fabs(Gauss3D->GetParameter(1)));
	      maxLongLength  = nSigmaZ  * sqrt(fabs(Gauss3D->GetParameter(2)));
	      goodData = Gauss3D->ExecuteCommand("MIGRAD",arglist,2);
	      // Re-set dimensions of the centroid for vertex rejection
	      maxTransRadius = nSigmaXY * sqrt(fabs(Gauss3D->GetParameter(0)) + fabs(Gauss3D->GetParameter(1)));
	      maxLongLength  = nSigmaZ  * sqrt(fabs(Gauss3D->GetParameter(2)));
	      goodData = Gauss3D->ExecuteCommand("MIGRAD",arglist,2);
	      
	      Gauss3D->GetStats(amin, edm, errdef, nvpar, nparx);
      
	      if (counterVx < minNentries) goodData = -2;
	      else if (isnan(edm) == true) goodData = -1;
	      else for (unsigned int j = 0; j < nParams; j++) if (isnan(Gauss3D->GetParError(j)) == true) { goodData = -1; break; }
	      if (goodData == 0)
		{
		  det = fabs(Gauss3D->GetParameter(0)) * (fabs(Gauss3D->GetParameter(1))*fabs(Gauss3D->GetParameter(2)) - Gauss3D->GetParameter(4)*Gauss3D->GetParameter(4)) -
		    Gauss3D->GetParameter(3) * (Gauss3D->GetParameter(3)*fabs(Gauss3D->GetParameter(2)) - Gauss3D->GetParameter(5)*Gauss3D->GetParameter(4)) +
		    Gauss3D->GetParameter(5) * (Gauss3D->GetParameter(3)*Gauss3D->GetParameter(4) - Gauss3D->GetParameter(5)*fabs(Gauss3D->GetParameter(1)));
		  if (det < 0.) goodData = -1;
		}
	      if ((goodData != 0) && (goodData != -2))
		{
		  Gauss3D->Clear();
	  
		  if (internalDebug == true) cout << "FIT WITH ENLARGED PARAMETER DISTANCES - STEP 3" << endl;      
		  // @@@ FIT WITH ENLARGED PARAMETER DISTANCES - STEP 3 @@@
		  // arg3 - first guess of parameter value
		  // arg4 - step of the parameter
		  Gauss3D->SetParameter(0,"var x ", *(it+0)*varFactor, parDistance*100., 0, 0);
		  Gauss3D->SetParameter(1,"var y ", *(it+1)*varFactor, parDistance*100., 0, 0);
		  Gauss3D->SetParameter(2,"var z ", *(it+2)*varFactor, parDistance*500., 0, 0);
		  Gauss3D->SetParameter(3,"cov xy", 0.0, parDistance*10., 0, 0);
		  Gauss3D->SetParameter(4,"cov yz", 0.0, parDistance*10., 0, 0);
		  Gauss3D->SetParameter(5,"cov xz", 0.0, parDistance*10., 0, 0);
		  Gauss3D->SetParameter(6,"mean x", *(it+6)+(double(bestMovementX)-1.)*sqrt((*(it+0))*varFactor), parDistance*100., 0, 0);
		  Gauss3D->SetParameter(7,"mean y", *(it+7)+(double(bestMovementY)-1.)*sqrt((*(it+1))*varFactor), parDistance*100., 0, 0);
		  Gauss3D->SetParameter(8,"mean z", *(it+8)+(double(bestMovementZ)-1.)*sqrt((*(it+2))*varFactor), parDistance*500., 0, 0);

		  // Set the central positions of the centroid for vertex rejection
		  xPos = Gauss3D->GetParameter(6);
		  yPos = Gauss3D->GetParameter(7);
		  zPos = Gauss3D->GetParameter(8);
		  
		  // Set dimensions of the centroid for vertex rejection
		  maxTransRadius = nSigmaXY * sqrt(fabs(Gauss3D->GetParameter(0)) + fabs(Gauss3D->GetParameter(1)));
		  maxLongLength  = nSigmaZ  * sqrt(fabs(Gauss3D->GetParameter(2)));
		  goodData = Gauss3D->ExecuteCommand("MIGRAD",arglist,2);
		  // Re-set dimensions of the centroid for vertex rejection
		  maxTransRadius = nSigmaXY * sqrt(fabs(Gauss3D->GetParameter(0)) + fabs(Gauss3D->GetParameter(1)));
		  maxLongLength  = nSigmaZ  * sqrt(fabs(Gauss3D->GetParameter(2)));
		  goodData = Gauss3D->ExecuteCommand("MIGRAD",arglist,2);

		  Gauss3D->GetStats(amin, edm, errdef, nvpar, nparx);
		  
		  if (counterVx < minNentries) goodData = -2;
		  else if (isnan(edm) == true) goodData = -1;
		  else for (unsigned int j = 0; j < nParams; j++) if (isnan(Gauss3D->GetParError(j)) == true) { goodData = -1; break; }
		  if (goodData == 0)
		    {
		      det = fabs(Gauss3D->GetParameter(0)) * (fabs(Gauss3D->GetParameter(1))*fabs(Gauss3D->GetParameter(2)) - Gauss3D->GetParameter(4)*Gauss3D->GetParameter(4)) -
			Gauss3D->GetParameter(3) * (Gauss3D->GetParameter(3)*fabs(Gauss3D->GetParameter(2)) - Gauss3D->GetParameter(5)*Gauss3D->GetParameter(4)) +
			Gauss3D->GetParameter(5) * (Gauss3D->GetParameter(3)*Gauss3D->GetParameter(4) - Gauss3D->GetParameter(5)*fabs(Gauss3D->GetParameter(1)));
		      if (det < 0.) goodData = -1;
		    }
		}
	    }
	}
      
      if (goodData == 0)
	for (unsigned int i = 0; i < nParams; i++)
	  {
	    vals->operator[](i) = Gauss3D->GetParameter(i);
	    vals->operator[](i+nParams) = Gauss3D->GetParError(i);
	  }
      
      delete Gauss3D;
      return goodData;
    }
  
  return -1;
}


void Vx3DHLTAnalyzer::reset(string ResetType)
{
  if (ResetType.compare("scratch") == 0)
    {
      runNumber      = 0;
      numberGoodFits = 0;
      numberFits     = 0;
      lastLumiOfFit  = 0;
      
      Vx_X->Reset();
      Vx_Y->Reset();
      Vx_Z->Reset();
      
      Vx_ZX->Reset();
      Vx_ZY->Reset();
      Vx_XY->Reset();
      
      Vx_ZX_profile->Reset();
      Vx_ZY_profile->Reset();

      Vertices.clear();
      
      lumiCounter    = 0;
      totalHits      = 0;
      beginTimeOfFit = 0;
      endTimeOfFit   = 0;
      beginLumiOfFit = 0;
      endLumiOfFit   = 0;
    }
  else if (ResetType.compare("whole") == 0)
    {
      Vx_X->Reset();
      Vx_Y->Reset();
      Vx_Z->Reset();
      
      Vx_ZX->Reset();
      Vx_ZY->Reset();
      Vx_XY->Reset();
      
      Vx_ZX_profile->Reset();
      Vx_ZY_profile->Reset();

      Vertices.clear();
      
      lumiCounter    = 0;
      totalHits      = 0;
      beginTimeOfFit = 0;
      endTimeOfFit   = 0;
      beginLumiOfFit = 0;
      endLumiOfFit   = 0;
    }
  else if (ResetType.compare("partial") == 0)
    {
      Vertices.clear();
      
      lumiCounter    = 0;
      totalHits      = 0;
      beginTimeOfFit = 0;
      endTimeOfFit   = 0;
      beginLumiOfFit = 0;
      endLumiOfFit   = 0;
    }
}


void Vx3DHLTAnalyzer::writeToFile(vector<double>* vals,
				  edm::TimeValue_t BeginTimeOfFit,
				  edm::TimeValue_t EndTimeOfFit,
				  unsigned int BeginLumiOfFit,
				  unsigned int EndLumiOfFit,
				  int dataType)
{
  stringstream BufferString;
  BufferString.precision(5);

  outputFile.open(fileName.c_str(), ios::out);

  if ((outputFile.is_open() == true) && (vals != NULL) && (vals->size() == 8*2))
    {
      vector<double>::const_iterator it = vals->begin();

      outputFile << "Runnumber " << runNumber << endl;
      outputFile << "BeginTimeOfFit " << formatTime(beginTimeOfFit / pow(2,32)) << endl;
      outputFile << "EndTimeOfFit " << formatTime(endTimeOfFit / pow(2,32)) << endl;
      outputFile << "LumiRange " << beginLumiOfFit << " - " << endLumiOfFit << endl;
      outputFile << "Type " << dataType << endl;
      // 3D Vertexing with Pixel Tracks good data = type 3
      // Bad data = type -1

      BufferString << *(it+0);
      outputFile << "X0 " << BufferString.str().c_str() << endl;
      BufferString.str("");

      BufferString << *(it+1);
      outputFile << "Y0 " << BufferString.str().c_str() << endl;
      BufferString.str("");

      BufferString << *(it+2);
      outputFile << "Z0 " << BufferString.str().c_str() << endl;
      BufferString.str("");

      BufferString << *(it+3);
      outputFile << "sigmaZ0 " << BufferString.str().c_str() << endl;
      BufferString.str("");

      BufferString << *(it+4);
      outputFile << "dxdz " << BufferString.str().c_str() << endl;
      BufferString.str("");

      BufferString << *(it+5);
      outputFile << "dydz " << BufferString.str().c_str() << endl;
      BufferString.str("");

      BufferString << *(it+6);
      outputFile << "BeamWidthX " << BufferString.str().c_str() << endl;
      BufferString.str("");

      BufferString << *(it+7);
      outputFile << "BeamWidthY " << BufferString.str().c_str() << endl;
      BufferString.str("");

      outputFile << "Cov(0,j) " << *(it+8) << " 0.0 0.0 0.0 0.0 0.0 0.0" << endl;
      outputFile << "Cov(1,j) 0.0 " << *(it+9) << " 0.0 0.0 0.0 0.0 0.0" << endl;
      outputFile << "Cov(2,j) 0.0 0.0 " << *(it+10) << " 0.0 0.0 0.0 0.0" << endl;
      outputFile << "Cov(3,j) 0.0 0.0 0.0 " << *(it+11) << " 0.0 0.0 0.0" << endl;
      outputFile << "Cov(4,j) 0.0 0.0 0.0 0.0 " << *(it+12) << " 0.0 0.0" << endl;
      outputFile << "Cov(5,j) 0.0 0.0 0.0 0.0 0.0 " << *(it+13) << " 0.0" << endl;
      outputFile << "Cov(6,j) 0.0 0.0 0.0 0.0 0.0 0.0 " << ((*(it+14))+(*(it+15)))/2. << endl;

      outputFile << "EmittanceX 0.0" << endl;
      outputFile << "EmittanceY 0.0" << endl;
      outputFile << "BetaStar 0.0" << endl;
    }
  outputFile.close();

  if ((debugMode == true) && (outputDebugFile.is_open() == true) && (vals != NULL) && (vals->size() == 8*2))
    {
      vector<double>::const_iterator it = vals->begin();
	  
      outputDebugFile << "Runnumber " << runNumber << endl;
      outputDebugFile << "BeginTimeOfFit " << formatTime(beginTimeOfFit / pow(2,32)) << endl;
      outputDebugFile << "EndTimeOfFit " << formatTime(endTimeOfFit / pow(2,32)) << endl;
      outputDebugFile << "LumiRange " << beginLumiOfFit << " - " << endLumiOfFit << endl;
      outputDebugFile << "Type " << dataType << endl;
      // 3D Vertexing with Pixel Tracks good data = type 3
      // Bad data = type -1
	  
      BufferString << *(it+0);
      outputDebugFile << "X0 " << BufferString.str().c_str() << endl;
      BufferString.str("");
	  
      BufferString << *(it+1);
      outputDebugFile << "Y0 " << BufferString.str().c_str() << endl;
      BufferString.str("");
	  
      BufferString << *(it+2);
      outputDebugFile << "Z0 " << BufferString.str().c_str() << endl;
      BufferString.str("");
	  
      BufferString << *(it+3);
      outputDebugFile << "sigmaZ0 " << BufferString.str().c_str() << endl;
      BufferString.str("");
	  
      BufferString << *(it+4);
      outputDebugFile << "dxdz " << BufferString.str().c_str() << endl;
      BufferString.str("");
	  
      BufferString << *(it+5);
      outputDebugFile << "dydz " << BufferString.str().c_str() << endl;
      BufferString.str("");
	  
      BufferString << *(it+6);
      outputDebugFile << "BeamWidthX " << BufferString.str().c_str() << endl;
      BufferString.str("");
	  
      BufferString << *(it+7);
      outputDebugFile << "BeamWidthY " << BufferString.str().c_str() << endl;
      BufferString.str("");
	  
      outputDebugFile << "Cov(0,j) " << *(it+8) << " 0.0 0.0 0.0 0.0 0.0 0.0" << endl;
      outputDebugFile << "Cov(1,j) 0.0 " << *(it+9) << " 0.0 0.0 0.0 0.0 0.0" << endl;
      outputDebugFile << "Cov(2,j) 0.0 0.0 " << *(it+10) << " 0.0 0.0 0.0 0.0" << endl;
      outputDebugFile << "Cov(3,j) 0.0 0.0 0.0 " << *(it+11) << " 0.0 0.0 0.0" << endl;
      outputDebugFile << "Cov(4,j) 0.0 0.0 0.0 0.0 " << *(it+12) << " 0.0 0.0" << endl;
      outputDebugFile << "Cov(5,j) 0.0 0.0 0.0 0.0 0.0 " << *(it+13) << " 0.0" << endl;
      outputDebugFile << "Cov(6,j) 0.0 0.0 0.0 0.0 0.0 0.0 " << ((*(it+14))+(*(it+15)))/2. << endl;
	  
      outputDebugFile << "EmittanceX 0.0" << endl;
      outputDebugFile << "EmittanceY 0.0" << endl;
      outputDebugFile << "BetaStar 0.0" << endl;
    }
}


void Vx3DHLTAnalyzer::beginLuminosityBlock(const LuminosityBlock& lumiBlock, 
					   const EventSetup& iSetup)
{
  if ((lumiCounter == 0) && (lumiBlock.luminosityBlock() > lastLumiOfFit))
    {
      beginTimeOfFit = lumiBlock.beginTime().value();
      beginLumiOfFit = lumiBlock.luminosityBlock();
      lumiCounter++;
    }
    else if ((lumiCounter != 0) && (lumiBlock.luminosityBlock() >= (beginLumiOfFit+lumiCounter))) lumiCounter++;
}


void Vx3DHLTAnalyzer::endLuminosityBlock(const LuminosityBlock& lumiBlock,
					 const EventSetup& iSetup)
{
  unsigned int nParams = 9;
  int goodData;

  if ((lumiCounter%nLumiReset == 0) && (nLumiReset != 0) && (beginTimeOfFit != 0) && (runNumber != 0))
    {            
      endTimeOfFit  = lumiBlock.endTime().value();
      endLumiOfFit  = lumiBlock.luminosityBlock();
      lastLumiOfFit = endLumiOfFit;
      vector<double> vals;
      stringstream histTitle;

      hitCounter->ShiftFillLast(totalHits, sqrt(totalHits), nLumiReset);

      if (dataFromFit == true)
	{
	  double dxdz, dydz;
	  vector<double> fitResults;

	  fitResults.push_back(Vx_X->getTH1()->GetRMS()*Vx_X->getTH1()->GetRMS());
	  fitResults.push_back(Vx_Y->getTH1()->GetRMS()*Vx_Y->getTH1()->GetRMS());
	  fitResults.push_back(Vx_Z->getTH1()->GetRMS()*Vx_Z->getTH1()->GetRMS());
	  fitResults.push_back(0.0);
	  fitResults.push_back(0.0);
	  fitResults.push_back(0.0);
	  fitResults.push_back(Vx_X->getTH1()->GetMean());
	  fitResults.push_back(Vx_Y->getTH1()->GetMean());
	  fitResults.push_back(Vx_Z->getTH1()->GetMean());
	  for (unsigned int i = 0; i < nParams; i++) fitResults.push_back(0.0);
	  
	  goodData = MyFit(&fitResults);	  	      

	  if (internalDebug == true) 
	    {
	      cout << "goodData --> " << goodData << endl;
	      cout << "Used vertices --> " << counterVx << endl;
	      cout << "var x -->  " << fitResults[0] << " +/- " << fitResults[0+nParams] << endl;
	      cout << "var y -->  " << fitResults[1] << " +/- " << fitResults[1+nParams] << endl;
	      cout << "var z -->  " << fitResults[2] << " +/- " << fitResults[2+nParams] << endl;
	      cout << "cov xy --> " << fitResults[3] << " +/- " << fitResults[3+nParams] << endl;
	      cout << "cov yz --> " << fitResults[4] << " +/- " << fitResults[4+nParams] << endl;
	      cout << "cov xz --> " << fitResults[5] << " +/- " << fitResults[5+nParams] << endl;
	      cout << "mean x --> " << fitResults[6] << " +/- " << fitResults[6+nParams] << endl;
	      cout << "mean y --> " << fitResults[7] << " +/- " << fitResults[7+nParams] << endl;
	      cout << "mean z --> " << fitResults[8] << " +/- " << fitResults[8+nParams] << endl;
	    }

	  if (goodData == 0)
	    {		 
	      dxdz = (fitResults[3]*fitResults[4] + fitResults[5]*(fitResults[2] - fitResults[1])) / ((fitResults[2]-fitResults[1])*(fitResults[2]-fitResults[0]) - fitResults[3]*fitResults[3]);
	      dydz = (fitResults[3]*fitResults[5] + fitResults[4]*(fitResults[2] - fitResults[0])) / ((fitResults[2]-fitResults[1])*(fitResults[2]-fitResults[0]) - fitResults[3]*fitResults[3]);
		  
	      vals.push_back(fitResults[6]);
	      vals.push_back(fitResults[7]);
	      vals.push_back(fitResults[8]);
	      vals.push_back(sqrt(fabs(fitResults[2])));
	      vals.push_back(dxdz);
	      vals.push_back(dydz);
	      vals.push_back(sqrt(fabs(fitResults[0])));
	      vals.push_back(sqrt(fabs(fitResults[1])));

	      vals.push_back(powf(fitResults[6+nParams],2.));
	      vals.push_back(powf(fitResults[7+nParams],2.));
	      vals.push_back(powf(fitResults[8+nParams],2.));
	      vals.push_back(powf(fabs(fitResults[2+nParams]) / (2.*sqrt(fabs(fitResults[2]))),2.));
	      vals.push_back(0.0);
	      vals.push_back(0.0);
	      vals.push_back(powf(fabs(fitResults[0+nParams]) / (2.*sqrt(fabs(fitResults[0]))),2.));
	      vals.push_back(powf(fabs(fitResults[1+nParams]) / (2.*sqrt(fabs(fitResults[1]))),2.));
	    }
	  else for (unsigned int i = 0; i < 8*2; i++) vals.push_back(0.0);

	  fitResults.clear();
	}
      else
	{
	  if (Vx_X->getTH1F()->GetEntries() >= minNentries)
	    {
	    goodData = 0;
	    
	    vals.push_back(Vx_X->getTH1F()->GetMean());
	    vals.push_back(Vx_Y->getTH1F()->GetMean());
	    vals.push_back(Vx_Z->getTH1F()->GetMean());
	    vals.push_back(Vx_Z->getTH1F()->GetRMS());
	    vals.push_back(0.0);
	    vals.push_back(0.0);
	    vals.push_back(Vx_X->getTH1F()->GetRMS());
	    vals.push_back(Vx_Y->getTH1F()->GetRMS());
	    
	    vals.push_back(powf(Vx_X->getTH1F()->GetMeanError(),2.));
	    vals.push_back(powf(Vx_Y->getTH1F()->GetMeanError(),2.));
	    vals.push_back(powf(Vx_Z->getTH1F()->GetMeanError(),2.));
	    vals.push_back(powf(Vx_Z->getTH1F()->GetRMSError(),2.));
	    vals.push_back(0.0);
	    vals.push_back(0.0);
	    vals.push_back(powf(Vx_X->getTH1F()->GetRMSError(),2.));
	    vals.push_back(powf(Vx_Y->getTH1F()->GetRMSError(),2.));
	    }
	  else
	    {
	      goodData = -2;
	      for (unsigned int i = 0; i < 8*2; i++) vals.push_back(0.0);
	    }
	}

      // vals[0]  = X0
      // vals[1]  = Y0
      // vals[2]  = Z0
      // vals[3]  = sigmaZ0
      // vals[4]  = dxdz
      // vals[5]  = dydz
      // vals[6]  = BeamWidthX
      // vals[7]  = BeamWidthY

      // vals[8]  = err^2 X0
      // vals[9]  = err^2 Y0
      // vals[10] = err^2 Z0
      // vals[11] = err^2 sigmaZ0
      // vals[12] = err^2 dxdz
      // vals[13] = err^2 dydz
      // vals[14] = err^2 BeamWidthX
      // vals[15] = err^2 BeamWidthY

      // "goodData" CODE:
      //  0 == OK --> Reset
      // -2 == NO OK - not enough "minNentries" --> Wait for more lumisections
      // Any other number == NO OK --> Reset

      numberFits++;
      if (goodData == 0)
	{
	  writeToFile(&vals, beginTimeOfFit, endTimeOfFit, beginLumiOfFit, endLumiOfFit, 3);
	  if ((internalDebug == true) && (outputDebugFile.is_open() == true)) outputDebugFile << "Used vertices: " << counterVx << endl;

	  numberGoodFits++;

	  histTitle << "Fitted Beam Spot [cm] (Lumi start: " << beginLumiOfFit << " - Lumi end: " << endLumiOfFit << ")";

	  reset("whole");
	}
      else
	{
	  writeToFile(&vals, beginTimeOfFit, endTimeOfFit, beginLumiOfFit, endLumiOfFit, -1);
	  if ((internalDebug == true) && (outputDebugFile.is_open() == true)) outputDebugFile << "Used vertices: " << counterVx << endl;

	  if (goodData == -2) { reset("whole"); histTitle << "Fitted Beam Spot [cm] (not enough statistics)"; }
	  else { histTitle << "Fitted Beam Spot [cm] (problems)"; if (lumiCounter == maxLumiIntegration) reset("whole"); }
	}

      reportSummary->Fill(numberFits != 0 ? (double)numberGoodFits/(double)numberFits : 0.0);
      reportSummaryMap->Fill(0.5, 0.5, numberFits != 0 ? (double)numberGoodFits/(double)numberFits : 0.0);

      fitResults->setAxisTitle(histTitle.str().c_str(), 1);
      
      fitResults->setBinContent(1, 9, vals[0]);
      fitResults->setBinContent(1, 8, vals[1]);
      fitResults->setBinContent(1, 7, vals[2]);
      fitResults->setBinContent(1, 6, vals[3]);
      fitResults->setBinContent(1, 5, vals[4]);
      fitResults->setBinContent(1, 4, vals[5]);
      fitResults->setBinContent(1, 3, vals[6]);
      fitResults->setBinContent(1, 2, vals[7]);
      fitResults->setBinContent(1, 1, counterVx);
      
      fitResults->setBinContent(2, 9, sqrt(vals[8]));
      fitResults->setBinContent(2, 8, sqrt(vals[9]));
      fitResults->setBinContent(2, 7, sqrt(vals[10]));
      fitResults->setBinContent(2, 6, sqrt(vals[11]));
      fitResults->setBinContent(2, 5, sqrt(vals[12]));
      fitResults->setBinContent(2, 4, sqrt(vals[13]));
      fitResults->setBinContent(2, 3, sqrt(vals[14]));
      fitResults->setBinContent(2, 2, sqrt(vals[15]));
      fitResults->setBinContent(2, 1, 0.0);

      // Linear fit to the historical plots
      TF1* myLinFit = new TF1("myLinFit", "[0] + [1]*x", mXlumi->getTH1()->GetXaxis()->GetXmin(), mXlumi->getTH1()->GetXaxis()->GetXmax());
      myLinFit->SetLineColor(2);
      myLinFit->SetLineWidth(2);
      myLinFit->SetParName(0,"Intercept");
      myLinFit->SetParName(1,"Slope");

      mXlumi->ShiftFillLast(vals[0], sqrt(vals[8]), nLumiReset);
      myLinFit->SetParameter(0, mXlumi->getTH1()->GetMean(2));
      myLinFit->SetParameter(1, 0.0);
      mXlumi->getTH1()->Fit("myLinFit","QR");

      mYlumi->ShiftFillLast(vals[1], sqrt(vals[9]), nLumiReset);
      myLinFit->SetParameter(0, mYlumi->getTH1()->GetMean(2));
      myLinFit->SetParameter(1, 0.0);
      mYlumi->getTH1()->Fit("myLinFit","QR");

      mZlumi->ShiftFillLast(vals[2], sqrt(vals[10]), nLumiReset);
      myLinFit->SetParameter(0, mZlumi->getTH1()->GetMean(2));
      myLinFit->SetParameter(1, 0.0);
      mZlumi->getTH1()->Fit("myLinFit","QR");

      sXlumi->ShiftFillLast(vals[6], sqrt(vals[14]), nLumiReset);
      myLinFit->SetParameter(0, sXlumi->getTH1()->GetMean(2));
      myLinFit->SetParameter(1, 0.0);
      sXlumi->getTH1()->Fit("myLinFit","QR");

      sYlumi->ShiftFillLast(vals[7], sqrt(vals[15]), nLumiReset);
      myLinFit->SetParameter(0, sYlumi->getTH1()->GetMean(2));
      myLinFit->SetParameter(1, 0.0);
      sYlumi->getTH1()->Fit("myLinFit","QR");

      sZlumi->ShiftFillLast(vals[3], sqrt(vals[11]), nLumiReset);
      myLinFit->SetParameter(0, sZlumi->getTH1()->GetMean(2));
      myLinFit->SetParameter(1, 0.0);
      sZlumi->getTH1()->Fit("myLinFit","QR");

      dxdzlumi->ShiftFillLast(vals[4], (vals[4] != 0.) ? 0.0001 : 0.0, nLumiReset);
      myLinFit->SetParameter(0, dxdzlumi->getTH1()->GetMean(2));
      myLinFit->SetParameter(1, 0.0);
      dxdzlumi->getTH1()->Fit("myLinFit","QR");

      dydzlumi->ShiftFillLast(vals[5], (vals[5] != 0.) ? 0.0001 : 0.0, nLumiReset);
      myLinFit->SetParameter(0, dydzlumi->getTH1()->GetMean(2));
      myLinFit->SetParameter(1, 0.0);
      dydzlumi->getTH1()->Fit("myLinFit","QR");
      
      delete myLinFit;

      vals.clear();
    }
  else if (nLumiReset == 0)
    {
      reportSummaryMap->Fill(0.5, 0.5, 1.0);
      hitCounter->ShiftFillLast(totalHits, sqrt(totalHits), 1);
      reset("partial");
    }
}


void Vx3DHLTAnalyzer::beginJob()
{
  DQMStore* dbe = 0;
  dbe = Service<DQMStore>().operator->();
 
  nBinsHistoricalPlot = 80;

  if ( dbe ) 
    {
      dbe->setCurrentFolder("BeamPixel");

      Vx_X = dbe->book1D("vertex x", "Primary Vertex X Coordinate Distribution", (int)(xRange/xStep), -xRange/2., xRange/2.);
      Vx_Y = dbe->book1D("vertex y", "Primary Vertex Y Coordinate Distribution", (int)(yRange/yStep), -yRange/2., yRange/2.);
      Vx_Z = dbe->book1D("vertex z", "Primary Vertex Z Coordinate Distribution", (int)(zRange/zStep), -zRange/2., zRange/2.);

      Vx_X->setAxisTitle("Primary Vertices X [cm]",1);
      Vx_X->setAxisTitle("Entries [#]",2);
      Vx_Y->setAxisTitle("Primary Vertices Y [cm]",1);
      Vx_Y->setAxisTitle("Entries [#]",2);
      Vx_Z->setAxisTitle("Primary Vertices Z [cm]",1);
      Vx_Z->setAxisTitle("Entries [#]",2);
 
      mXlumi = dbe->book1D("muX vs lumi", "\\mu_{x} vs. Lumisection", nBinsHistoricalPlot, 0.5, (double)nBinsHistoricalPlot+0.5);
      mYlumi = dbe->book1D("muY vs lumi", "\\mu_{y} vs. Lumisection", nBinsHistoricalPlot, 0.5, (double)nBinsHistoricalPlot+0.5);
      mZlumi = dbe->book1D("muZ vs lumi", "\\mu_{z} vs. Lumisection", nBinsHistoricalPlot, 0.5, (double)nBinsHistoricalPlot+0.5);

      mXlumi->setAxisTitle("Lumisection [#]",1);
      mXlumi->setAxisTitle("\\mu_{x} [cm]",2);
      mXlumi->getTH1()->SetOption("E1");
      mYlumi->setAxisTitle("Lumisection [#]",1);
      mYlumi->setAxisTitle("\\mu_{y} [cm]",2);
      mYlumi->getTH1()->SetOption("E1");
      mZlumi->setAxisTitle("Lumisection [#]",1);
      mZlumi->setAxisTitle("\\mu_{z} [cm]",2);
      mZlumi->getTH1()->SetOption("E1");

      sXlumi = dbe->book1D("sigmaX vs lumi", "\\sigma_{x} vs. Lumisection", nBinsHistoricalPlot, 0.5, (double)nBinsHistoricalPlot+0.5);
      sYlumi = dbe->book1D("sigmaY vs lumi", "\\sigma_{y} vs. Lumisection", nBinsHistoricalPlot, 0.5, (double)nBinsHistoricalPlot+0.5);
      sZlumi = dbe->book1D("sigmaZ vs lumi", "\\sigma_{z} vs. Lumisection", nBinsHistoricalPlot, 0.5, (double)nBinsHistoricalPlot+0.5);

      sXlumi->setAxisTitle("Lumisection [#]",1);
      sXlumi->setAxisTitle("\\sigma_{x} [cm]",2);
      sXlumi->getTH1()->SetOption("E1");
      sYlumi->setAxisTitle("Lumisection [#]",1);
      sYlumi->setAxisTitle("\\sigma_{y} [cm]",2);
      sYlumi->getTH1()->SetOption("E1");
      sZlumi->setAxisTitle("Lumisection [#]",1);
      sZlumi->setAxisTitle("\\sigma_{z} [cm]",2);
      sZlumi->getTH1()->SetOption("E1");

      dxdzlumi = dbe->book1D("dxdz vs lumi", "dX/dZ vs. Lumisection", nBinsHistoricalPlot, 0.5, (double)nBinsHistoricalPlot+0.5);
      dydzlumi = dbe->book1D("dydz vs lumi", "dY/dZ vs. Lumisection", nBinsHistoricalPlot, 0.5, (double)nBinsHistoricalPlot+0.5);

      dxdzlumi->setAxisTitle("Lumisection [#]",1);
      dxdzlumi->setAxisTitle("dX/dZ [rad]",2);
      dxdzlumi->getTH1()->SetOption("E1");
      dydzlumi->setAxisTitle("Lumisection [#]",1);
      dydzlumi->setAxisTitle("dY/dZ [rad]",2);
      dydzlumi->getTH1()->SetOption("E1");

      Vx_ZX = dbe->book2D("vertex zx", "Primary Vertex ZX Coordinate Distribution", (int)(zRange/zStep/10.), -zRange/2., zRange/2., (int)(xRange/xStep/10.), -xRange/2., xRange/2.);
      Vx_ZY = dbe->book2D("vertex zy", "Primary Vertex ZY Coordinate Distribution", (int)(zRange/zStep/10.), -zRange/2., zRange/2., (int)(yRange/yStep/10.), -yRange/2., yRange/2.);
      Vx_XY = dbe->book2D("vertex xy", "Primary Vertex XY Coordinate Distribution", (int)(xRange/xStep/10.), -xRange/2., xRange/2., (int)(yRange/yStep/10.), -yRange/2., yRange/2.);

      Vx_ZX->setAxisTitle("Primary Vertices Z [cm]",1);
      Vx_ZX->setAxisTitle("Primary Vertices X [cm]",2);
      Vx_ZX->setAxisTitle("Entries [#]",3);
      Vx_ZY->setAxisTitle("Primary Vertices Z [cm]",1);
      Vx_ZY->setAxisTitle("Primary Vertices Y [cm]",2);
      Vx_ZY->setAxisTitle("Entries [#]",3);
      Vx_XY->setAxisTitle("Primary Vertices X [cm]",1);
      Vx_XY->setAxisTitle("Primary Vertices Y [cm]",2);
      Vx_XY->setAxisTitle("Entries [#]",3);

      Vx_ZX_profile = dbe->bookProfile("zx profile","ZX Profile", (int)(zRange/zStep/20.), -zRange/2., zRange/2., (int)(xRange/xStep/20.), -xRange/2., xRange/2., "");
      Vx_ZX_profile->setAxisTitle("Primary Vertices Z [cm]",1);
      Vx_ZX_profile->setAxisTitle("Primary Vertices X [cm]",2);

      Vx_ZY_profile = dbe->bookProfile("zy profile","ZY Profile", (int)(zRange/zStep/20.), -zRange/2., zRange/2., (int)(yRange/yStep/20.), -yRange/2., yRange/2., "");
      Vx_ZY_profile->setAxisTitle("Primary Vertices Z [cm]",1);
      Vx_ZY_profile->setAxisTitle("Primary Vertices Y [cm]",2);

      hitCounter = dbe->book1D("pixelHits vs lumi", "# Pixel-Hits vs. Lumisection", nBinsHistoricalPlot, 0.5, (double)nBinsHistoricalPlot+0.5);

      hitCounter->setAxisTitle("Lumisection [#]",1);
      hitCounter->setAxisTitle("# Pixel-Hits [#]",2);
      hitCounter->getTH1()->SetOption("E1");

      fitResults = dbe->book2D("fit results","Results of Beam Spot Fit", 2, 0., 2., 9, 0., 9.);
      fitResults->setAxisTitle("Fitted Beam Spot [cm]", 1);
      fitResults->setBinLabel(9, "X0", 2);
      fitResults->setBinLabel(8, "Y0", 2);
      fitResults->setBinLabel(7, "Z0", 2);
      fitResults->setBinLabel(6, "sigmaZ0", 2);
      fitResults->setBinLabel(5, "dX/dZ", 2);
      fitResults->setBinLabel(4, "dY/dZ", 2);
      fitResults->setBinLabel(3, "sigmaX0", 2);
      fitResults->setBinLabel(2, "sigmaY0", 2);
      fitResults->setBinLabel(1, "Vertices", 2);
      fitResults->setBinLabel(1, "Value", 1);
      fitResults->setBinLabel(2, "Stat. Error", 1);
      fitResults->getTH1()->SetOption("text");

      dbe->setCurrentFolder("BeamPixel/EventInfo");
      reportSummary = dbe->bookFloat("reportSummary");
      reportSummary->Fill(0.);
      reportSummaryMap = dbe->book2D("reportSummaryMap","Pixel-Vertices Beam Spot: % Good Fits", 1, 0., 1., 1, 0., 1.);
      reportSummaryMap->Fill(0.5, 0.5, 0.);
      dbe->setCurrentFolder("BeamPixel/EventInfo/reportSummaryContents");

      // Convention for reportSummary:
      // -   0% at the moment of creation of the histogram
      // -  95% if iether not not enough "minNentries" or bad fit
      // - 100% if good fit

      // Convention for reportSummaryMap:
      // - 0%  at the moment of creation of the histogram
      // - n%  numberGoodFits / numberFits
    }

  reset("scratch");
  maxLumiIntegration   = 100;
  minVxDoF             = 4.;
  internalDebug        = false;
  considerVxCovariance = true;

  pi = 3.141592653589793238;
}


void Vx3DHLTAnalyzer::endJob() { reset("scratch"); }


// Define this as a plug-in
DEFINE_FWK_MODULE(Vx3DHLTAnalyzer);
