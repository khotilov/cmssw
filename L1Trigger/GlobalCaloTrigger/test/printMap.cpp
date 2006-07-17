
#include <iostream>
#include <fstream>
#include <iomanip>



using std::ofstream;
using std::setw;
using std::endl;


#include "DataFormats/L1CaloTrigger/interface/L1CaloRegionDetId.h"

// print the region map and region-associated data (rct/gct card/input numbers)

void makePlot(ofstream& of, int plot);

int main() {


  // GCT card number
  ofstream gctCardPlot("gctCardPlot.txt");
  makePlot(gctCardPlot, 0);

  // GCT region number
  ofstream gctRgnPlot("gctRgnPlot.txt");
  makePlot(gctRgnPlot, 1);

  // RCT crate number
  ofstream rctCratePlot("rctCratePlot.txt");
  makePlot(rctCratePlot, 2);

  // RCT card number
  ofstream rctCardPlot("rctCardPlot.txt");
  makePlot(rctCardPlot, 3);

  // RCT region number
  ofstream rctRgnPlot("rctRgnPlot.txt");
  makePlot(rctRgnPlot, 4);

  // RCT phi
  ofstream rctPhiPlot("rctPhiPlot.txt");
  makePlot(rctPhiPlot, 5);

  // RCT eta
  ofstream rctEtaPlot("rctEtaPlot.txt");
  makePlot(rctEtaPlot, 6);
  
  // forward
  ofstream fwdRgnPlot("fwdRgnPlot.txt");
  makePlot(fwdRgnPlot, 7);



  return 0;

};


void makePlot(ofstream &of, int plot) {

  // print header line
  of << "     ieta->" << endl;
  of << "    : ";
  for (int ieta=0; ieta<22; ieta++) {
    of << setw(2) << ieta << " ";
  }
  of << endl;

  // main loop
  for (int iphi=0; iphi<18; iphi++) {
  
    of << setw(3) << iphi << "  : ";
    
    for (int ieta=0; ieta<22; ieta++) {    

      L1CaloRegionDetId rgn(ieta, iphi);

      switch (plot) {
      case 0 : 
	of << setw(3) << rgn.gctCard() << " "; break;
      case 1 : 
	of << setw(3) << rgn.gctRegion() << " "; break;
      case 2 : 
	of << setw(3) << rgn.rctCrate() << " "; break;
      case 3 : 
	of << setw(3) << rgn.rctCard() << " "; break;
      case 4 : 
	of << setw(3) << rgn.rctRegion() << " "; break;
      case 5 : 
	of << setw(3) << rgn.rctPhi() << " "; break;
      case 6 : 
	of << setw(3) << rgn.rctEta() << " "; break;
      case 7 : 
	of << setw(3) << (rgn.isForward()?1:0) << " "; break;


      default :
	break;
      }

    }

    of << endl;
  }

}
