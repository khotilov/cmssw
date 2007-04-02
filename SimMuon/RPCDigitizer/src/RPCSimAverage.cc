#include "Geometry/RPCGeometry/interface/RPCRoll.h"
#include "Geometry/RPCGeometry/interface/RPCRollSpecs.h"
#include "SimMuon/RPCDigitizer/src/RPCSimAverage.h"
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"

#include <cmath>

#include <CLHEP/Random/RandGaussQ.h>
#include <CLHEP/Random/RandFlat.h>

#include<cstring>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<stdlib.h>

using namespace std;


RPCSimAverage::RPCSimAverage(const edm::ParameterSet& config) : RPCSim(config){
  aveEff = config.getParameter<double>("averageEfficiency");
  aveCls = config.getParameter<double>("averageClusterSize");
  resRPC = config.getParameter<double>("timeResolution");
  timOff = config.getParameter<double>("timingRPCOffset");
  dtimCs = config.getParameter<double>("deltatimeAdjacentStrip");
  resEle = config.getParameter<double>("timeJitter");
  sspeed = config.getParameter<double>("signalPropagationSpeed");
  lbGate = config.getParameter<double>("linkGateWidth");
  rpcdigiprint = config.getParameter<bool>("printOutDigitizer");
  if (rpcdigiprint) {
    std::cout <<"Average Efficiency        = "<<aveEff<<std::endl;
    std::cout <<"Average Cluster Size      = "<<aveCls<<" strips"<<std::endl;
    std::cout <<"RPC Time Resolution       = "<<resRPC<<" ns"<<std::endl;
    std::cout <<"RPC Signal formation time = "<<timOff<<" ns"<<std::endl;
    std::cout <<"RPC adjacent strip delay  = "<<dtimCs<<" ns"<<std::endl;
    std::cout <<"Electronic Jitter         = "<<resEle<<" ns"<<std::endl;
    std::cout <<"Signal propagation time   = "<<sspeed<<" x c"<<std::endl;
    std::cout <<"Link Board Gate Width     = "<<lbGate<<" ns"<<std::endl;
  }

  string ifile="ClSizeTot.dat";

  MyOutput1 = new fstream("clsizeOut1.dat",ios::out);
  MyOutput2 = new fstream("clsizeOut2.dat",ios::out);
  MyOutput3 = new fstream("clsizeOut3.dat",ios::out);

  infile = new ifstream(ifile.c_str(), ios::in);
  if(! *infile) {
    cerr << "error: unable to open input file: "
	 <<  ifile  << endl;
    //    return -1;
  }

  string buffer;
  double sum = 0;
  unsigned int counter = 1;
  int row = 1;
  std::vector<double> sum_clsize;

  while ( *infile >> buffer ) {
    const char *buffer1 = buffer.c_str();
    double dato = atof(buffer1);
    sum += dato;
    sum_clsize.push_back(sum);

    std::cout<<dato<<std::endl;
    
    if(counter == row*20) {

      for(int i = 0; i < sum_clsize.size(); ++i){
	std::cout<<sum_clsize[i]<<std::endl;
      }

      clsMap[row] = sum_clsize;
      row++;
      sum = 0;
      sum_clsize.clear();
    }
    counter++;
  }
  std::cout<<"Somma = "<<(clsMap[1])[0]<<"  "<<(clsMap[2])[0]<<"  "<<(clsMap[3])[0]<<"  "<<(clsMap[4])[0]<<std::endl;

}

int RPCSimAverage::getClSize(float posX)
{

// //     string ifile = getenv("ORCA_DATA_PATH");
// //     int pos = ifile.find_first_of(':');
// //     if (pos != string::npos){
// //       ifile.erase(pos, ifile.size() -1);
// //       ifile=ifile + "/Param_Clust/ClSize_tot_9600.00.dat";
// //     }

  int cnt = 1;
  int min = 1;
  int max = 1;
  double func;
  std::vector<double> sum_clsize;

  double rr_cl = RandFlat::shoot();
  if(0.0 <= posX && posX < 0.2)  {
    func = (clsMap[1])[(clsMap[1]).size()-1]*(rr_cl);
    sum_clsize = clsMap[1];
  }
  if(0.2 <= posX && posX < 0.4) {
    func = (clsMap[2])[(clsMap[2]).size()-1]*(rr_cl);
    sum_clsize = clsMap[2];
  }
  if(0.4 <= posX && posX < 0.6) {
    func = (clsMap[3])[(clsMap[3]).size()-1]*(rr_cl);
    sum_clsize = clsMap[3];
  }
  if(0.6 <= posX && posX < 0.8) {
    func = (clsMap[4])[(clsMap[4]).size()-1]*(rr_cl);
    sum_clsize = clsMap[4];
  }
  if(0.8 <= posX && posX < 1.0)  {
    func = (clsMap[5])[(clsMap[5]).size()-1]*(rr_cl);
    sum_clsize = clsMap[5];
  }
  cout<<"rr_cl = "<<rr_cl<<"  "<<"SumCl = "<<(clsMap[1])[(clsMap[1]).size()-1]<<"  "<<"func = "<<func<<endl;
  
  for(vector<double>::iterator iter = sum_clsize.begin();
      iter != sum_clsize.end(); ++iter){
      cnt++;
      if(func > (*iter)){
        min = cnt;
      }
      else if(func < (*iter)){
        max = cnt;
        break;
      }
  }
  //  *MyOutput<<min<<endl;
  return min;
}


void
RPCSimAverage::simulate(const RPCRoll* roll,
		      const edm::PSimHitContainer& rpcHits )
{
  const Topology& topology=roll->specs()->topology();
  for (edm::PSimHitContainer::const_iterator _hit = rpcHits.begin();
       _hit != rpcHits.end(); ++_hit){

 
    // Here I hould check if the RPC are up side down;
    const LocalPoint& entr=_hit->entryPoint();
    //    const LocalPoint& exit=_hit->exitPoint();

    float posX = roll->strip(_hit->localPosition()) - static_cast<int>(roll->strip(_hit->localPosition()));

    // Effinciecy
    if (RandFlat::shoot() < aveEff) {
      int centralStrip = topology.channel(entr)+1;  
      int fstrip=centralStrip;
      int lstrip=centralStrip;
      // Compute the cluster size
      double w = RandFlat::shoot();
      if (w < 1.e-10) w=1.e-10;
      int clsize = RPCSimAverage::getClSize(posX);//static_cast<int>( -1.*aveCls*log(w)+1.);

      //---------------------VERIFICA------------------------------

      if(clsize == 1) *MyOutput1<<posX<<endl;
      if(clsize == 2) *MyOutput2<<posX<<endl;
      if(clsize == 3) *MyOutput3<<posX<<endl;

      //-----------------------------------------------------------

      std::vector<int> cls;
      cls.push_back(centralStrip);
      if (clsize > 1){
	for (int cl = 0; cl < (clsize-1)/2; cl++)
	  if (centralStrip - cl -1 >= 1  ){
	    fstrip = centralStrip-cl-1;
	    cls.push_back(fstrip);
	  }
	for (int cl = 0; cl < (clsize-1)/2; cl++)
	  if (centralStrip + cl + 1 <= roll->nstrips() ){
	    lstrip = centralStrip+cl+1;
	    cls.push_back(lstrip);
	  }
	if (clsize%2 == 0 ){
	  // insert the last strip according to the 
	  // simhit position in the central strip 
	  double deltaw=roll->centreOfStrip(centralStrip).x()-entr.x();
	  if (deltaw>0.) {
	    if (lstrip < roll->nstrips() ){
	      lstrip++;
	      cls.push_back(lstrip);
	    }
	  }else{
	    if (fstrip > 1 ){
	      fstrip--;
	      cls.push_back(fstrip);
	    }
	  }
	}
      }
      for (std::vector<int>::iterator i=cls.begin(); i!=cls.end();i++){
	// Check the timing of the adjacent strip
	strips.insert(*i);
      }
    }
  }
}


