#ifndef DDEcalPreshowerAlgo_h
#define DDEcalPreshowerAlgo_h

#include <vector>
#include <string>

#include "DetectorDescription/Algorithm/interface/DDAlgorithm.h"
#include "DetectorDescription/Core/interface/DDMaterial.h"

using namespace std;

class DDEcalPreshowerAlgo : public DDAlgorithm
{
 public:
  
  DDMaterial getMaterial(unsigned int i)   const {return DDMaterial(materials_[i]);}
  DDMaterial getLaddMaterial()  const { return DDMaterial(LaddMaterial_) ; }
  string getLayName(unsigned int i)   const {return layName_[i];}
  string getLadPrefix(unsigned int i)   const {return ladPfx_[i];}

  DDEcalPreshowerAlgo();
  void initialize(const DDNumericArguments & nArgs,
		  const DDVectorArguments & vArgs,
		  const DDMapArguments & mArgs,
		  const DDStringArguments & sArgs,
		  const DDStringVectorArguments & vsArgs);
  void execute();

 private:

  void doLayers();
  void doLadders(); 
  void doSens();
  
  int nmat_;                       // number of preshower layers
  double thickness_;               // overall thickness of the preshower envelope
  vector<string> materials_;       // materials of the presh-layers
  vector<string> layName_;         // names of the presh-layers
  vector<string> ladPfx_ ;         // name prefix for ladders
  string LaddMaterial_;            // ladd material - air
  vector<double> thickLayers_; 
  vector<double> abs1stx;
  vector<double> abs1sty;
  vector<double> abs2ndx;
  vector<double> abs2ndy;
  vector<double> asym_ladd_;
  vector<double> rminVec; 
  vector<double> rmaxVec;
  vector<double> noLaddInCol_;
  vector<double> startOfFirstLadd_;
  vector<string> types_l5_;
  vector<string> types_l4_;
  vector<double> ladd_l5_map_;
  vector<double> ladd_l4_map_;
  vector<string> typeOfLaddRow0;
  vector<string> typeOfLaddRow1;
  vector<string> typeOfLaddRow2;
  vector<string> typeOfLaddRow3;

  double zlead1_, zlead2_, zfoam1_, zfoam2_;
  double waf_intra_col_sep, waf_inter_col_sep, waf_active, wedge_length, wedge_offset, zwedge_ceramic_diff, ywedge_ceramic_diff, wedge_angle, box_thick,dee_separation, In_rad_Abs_Al, In_rad_Abs_Pb;
  double ladder_thick, yladder_1stwedge_diff, ladder_width, ladder_length, micromodule_length;
  double absAlX_X_, absAlX_Y_, absAlX_subtr1_Xshift_, absAlX_subtr1_Yshift_, rMax_Abs_Al_;
  double absAlY_X_, absAlY_Y_, absAlY_subtr1_Xshift_, absAlY_subtr1_Yshift_;
  double LdrBck_Length, LdrFrnt_Length,LdrFrnt_Offset,LdrBck_Offset, ceramic_length, wedge_back_thick;
  
};


#endif // DDEcalPreshowerAlgo_h
