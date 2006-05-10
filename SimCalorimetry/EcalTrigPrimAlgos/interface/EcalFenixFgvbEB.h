#ifndef ECAL_FENIX_FGVB_EB_H
#define ECAL_FENIX_FGVB_EB_H

#include <SimCalorimetry/EcalTrigPrimAlgos/interface/EcalVFgvb.h>
#include <vector>

namespace tpg {

// global type definitions for header defined by Tag entries in ArgoUML
// Result: typedef <typedef_global_header> <tag_value>;


  /** 
   \class EcalFenixFgvbEB
   \brief calculation of Fgvb for Fenix Tcp, format barrel
   *  calculates fgvb for the barrel
   *  
   *  
   *  input: 2X12 bits ( 12 bits Ettot + 12 bits maxof2)
   *  output: 1 bit 
   *  
   *  
   *  makes comparisons between maxof2 and 2 fractions of Ettot and  uses this comparison to decide ---> needs to get some values from outside
   */
class EcalFenixFgvbEB : public EcalVFgvb {


 public:
  EcalFenixFgvbEB() {;}
  int process() {return 0;} //FIXME: find better base methods

  std::vector<int> process( std::vector<int> add_out, std::vector<int> maxof2_out);
  };

} /* End of namespace tpg */

#endif
