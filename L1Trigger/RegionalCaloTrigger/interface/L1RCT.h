#ifndef L1RCT_h
#define L1RCT_h

#include <vector>
#include "L1RCTCrate.h"
#include "L1RCTNeighborMap.h"

class L1RCT {

 public:
  
  L1RCT();

  //Should organize the input into the crates and cards then pass on
  //the output.  The data will come into the trigger in the form
  //of two multilayered vectors of vectors.
  //The first is of the actual barrel information.
  //18 crates -> 7 RCs -> 64 unsigned shorts per RC
  //so it should be a vector<vector<vector<unsigned short> > >
  //The second is of the HF regions which is just of type
  //vector<vector<unsigned short> >
  void input(vector<vector<vector<unsigned short> > > barrel,
	     vector<vector<unsigned short> > hf);
  //Should send commands to all crates to send commands to all RCs to
  //process the input data and then send it on to the EICs and then
  //to the JSCs
  void processEvent();

  void fileInput(char* filename);

  void randomInput();

  void print();
  void printCrate(int i){
    crates.at(i).print();
  }
  void printJSC(int i){
    crates.at(i).printJSC();
  }
  void printJSC(){
    for(int i=0;i<18;i++){
      cout << "JSC for Crate " << i << endl;
      crates.at(i).printJSC();
    }
  }
  void printRC(int i, int j){
    crates.at(i).printRC(j);
  }
  void printEIC(int i, int j){
    crates.at(i).printEIC(j);
  }
  void printEICEdges(int i, int j){
    crates.at(i).printEICEdges(j);
  }

  vector<unsigned short> getIsolatedEGObjects(int crate){
    crates.at(crate).getIsolatedEGObjects();
  }
  vector<unsigned short> getNonisolatedEGObjects(int crate){
    crates.at(crate).getNonisolatedEGObjects();
  }
  vector<unsigned short> getJetRegions(int crate){
    crates.at(crate).getJetRegions();
  }
  
  
 private:
  
  //Method called by constructor to set all the neighbors properly.
  //Will make use of the internal neighborMap
  void configureCards();
  void shareNeighbors();

  L1RCTRegion empty;

  //Helper class containing information to set all the neighbors for
  //the receiver cards.  We will use the methods
  //north,south,west,east,se,sw,ne,nw,which each take in the
  //indices of the region in question and then return the indices
  //of the region that is the corresponding neighbor to be used.
  L1RCTNeighborMap neighborMap;

  //Vector of all 18 crates.
  //Will follow numbering convention listed
  //in the CaloTrigger Tower Mapping
  //So 0->8 are eta -5 -> 0
  //While 9-17 are eta 0 -> 5
  //Crate i and crate i+8 are next to each other  
  vector<L1RCTCrate> crates;
};

#endif
