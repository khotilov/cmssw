#include "L1Trigger/RegionalCaloTrigger/interface/L1RCT.h"
#include "L1Trigger/RegionalCaloTrigger/interface/L1RCTLookupTables.h"
#include "CondFormats/L1TObjects/interface/L1RCTParameters.h"

#include <time.h>
#include <iostream>
using namespace std;
#include <vector>

int main (){
  // For testing use 1:1 LUT
  std::vector<double> eGammaECalScaleFactors(32, 1.0);
  std::vector<double> eGammaHCalScaleFactors(32, 1.0);
  std::vector<double> jetMETECalScaleFactors(32, 1.0);
  std::vector<double> jetMETHCalScaleFactors(32, 1.0);
  L1RCTParameters* rctParameters = 
    new L1RCTParameters(1.0,                       // eGammaLSB
			1.0,                       // jetMETLSB
			3.0,                       // eMinForFGCut
			40.0,                      // eMaxForFGCut
			0.5,                       // hOeCut
			1.0,                       // eMinForHoECut
			50.0,                      // eMaxForHoECut
			2.0,                       // eActivityCut
			3.0,                       // hActivityCut
			3.0,                       // eicIsolationThreshold
			false,                     // ignoreTowerHB
			false,                     // ignoreTowerHEplus
			false,                     // ignoreTowerHEminus
			eGammaECalScaleFactors,
			eGammaHCalScaleFactors,
			jetMETECalScaleFactors,
			jetMETHCalScaleFactors
			);
  L1RCTLookupTables* lut = new L1RCTLookupTables();
  lut->setRCTParameters(rctParameters);  // transcoder and etScale are not used
  L1RCT rct(lut);
  vector<vector<unsigned short> > hf(18,vector<unsigned short>(8));
  vector<vector<vector<unsigned short> > > barrel(18,vector<vector<unsigned short> >(7,vector<unsigned short>(64)));
  time_t s1 = time (NULL);
  for(int i=0;i<10000;i++){
    rct.input(barrel,hf);
    rct.processEvent();
    cout << ".";
    if(!(i % 100)) cout << " " << i << " ";
  }
  time_t s2 = time (NULL);
  double diff =difftime (s2,s1);
  cout << diff << endl;
}
