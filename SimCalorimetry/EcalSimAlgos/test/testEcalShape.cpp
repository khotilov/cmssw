#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EcalSimParameterMap.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EcalShape.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include<iostream>
#include<iomanip>

#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TCanvas.h"

int main() {

  edm::MessageDrop::instance()->debugEnabled = false;

  EcalSimParameterMap parameterMap;

  EBDetId barrel(1,1);
  double thisPhase = parameterMap.simParameters(barrel).timePhase();
  EcalShape theShape(thisPhase);

  std::cout << "Parameters for the ECAL MGPA shape \n" << std::endl;

  std::cout << "Rising time for ECAL shape (timePhase) = " << parameterMap.simParameters(barrel).timePhase() << std::endl;
  std::cout << "Bin of maximum = " << parameterMap.simParameters(barrel).binOfMaximum() << std::endl;

  // standard display of the implemented shape function

  theShape.display();

  double ToM = theShape.computeTimeOfMaximum();
  double T0 = theShape.computeT0();
  double risingTime = theShape.computeRisingTime();

  std::cout << "\n Maximum time from tabulated values = " << std::setprecision(2) << ToM << std::endl;
  std::cout << "\n Tzero from tabulated values        = " << std::setprecision(2) << T0 << std::endl;
  std::cout << "\n Rising time from tabulated values  = " << std::setprecision(2) << risingTime << std::endl;

  // signal used with the nominal parameters and no jitter

  const int csize = 500;
  TCanvas * showShape = new TCanvas("showShape","showShape",2*csize,csize);

  std::cout << "\n computed ECAL pulse shape and its derivative (LHC timePhaseShift = 1) \n" << std::endl;
  double x = risingTime-(parameterMap.simParameters(barrel).binOfMaximum()-1.)*25.;
  double startx = x;

  int nsamp = 250;
  int tconv = 10;
  unsigned int histsiz = nsamp*tconv;

  TH1F* shape2 = new TH1F("shape2","Computed Ecal MGPA shape",nsamp,0.,(float)(nsamp));
  TH1F* deriv2 = new TH1F("deriv2","Computed Ecal MGPA derivative",nsamp,0.,(float)(nsamp));
  double y = 0.;
  double dy = 0.;
  for ( unsigned int i = 0; i < histsiz; ++i ) {
    y = (theShape)(x);
    dy = theShape.derivative(x);
    shape2->Fill((float)(x-startx),(float)y);
    deriv2->Fill((float)(x-startx),(float)dy);
    std::cout << " time (ns) = " << std::fixed << std::setw(6) << std::setprecision(2) << x-startx 
              << " shape = " << std::setw(11) << std::setprecision(8) << y
              << " derivative = " << std::setw(11) << std::setprecision(8) << dy << std::endl;
    x = x+1./(double)tconv;
  }

  showShape->Divide(2,1);
  showShape->cd(1);
  gPad->SetGrid();
  shape2->GetXaxis()->SetNdivisions(10,kFALSE);
  shape2->Draw();
  showShape->cd(2);
  gPad->SetGrid();
  deriv2->GetXaxis()->SetNdivisions(10,kFALSE);
  deriv2->Draw();
  showShape->SaveAs("EcalShapeUsed.jpg");

  delete shape2;
  delete deriv2;
  delete showShape;
  

  return 0;

} 
