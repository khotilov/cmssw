///////////////////////////////////////////////////////////////////////////////
// File: HcalNumberingFromDDD.h
// Description: Usage of DDD to get to numbering scheme for hadron calorimeter
///////////////////////////////////////////////////////////////////////////////
#ifndef HcalNumberingFromDDD_h
#define HcalNumberingFromDDD_h

#include "Geometry/HcalCommonData/interface/HcalCellType.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DetectorDescription/Core/interface/DDsvalues.h"

#include "CLHEP/Vector/ThreeVector.h"

#include <vector>
#include <string>

class DDCompactView;    
class DDFilteredView;

class HcalNumberingFromDDD {

public:

  HcalNumberingFromDDD(std::string & name, const DDCompactView & cpv);
  ~HcalNumberingFromDDD();
	 
  struct HcalID {
    int subdet, zside, depth, etaR, phi, phis, lay;
    HcalID(int det=0, int zs=0, int d=0, int et=0, int fi=0, int phiskip=0, int ly=-1) :
      subdet(det), zside(zs), depth(d), etaR(et), phi(fi), phis(phiskip), lay(ly) {}
  };

  HcalID         unitID(int det, Hep3Vector pos, int depth, int lay=-1) const;
  HcalID         unitID(double eta, double phi, int depth=1, int lay=-1) const;
  HcalID         unitID(int det, int zside, int depth, int etaR, int phi, 
			int lay=-1) const;
  HcalCellType::HcalCell cell(int det, int zside, int depth, int etaR, 
			      int iphi, bool corr=true) const;
  std::vector<double> getEtaTable() const;
  std::vector<HcalCellType::HcalCellType> HcalCellTypes() const;
  std::vector<HcalCellType::HcalCellType> HcalCellTypes(HcalSubdetector) const;

private:

  double         getEta(int det, int etaR, int zside, int depth=1) const;
  double         getEta(double r, double z) const;
  double         deltaEta(int det, int eta, int depth) const;
  void           initialize(std::string & name, const DDCompactView & cpv);
  void           loadSpecPars(DDFilteredView);
  void           loadGeometry(DDFilteredView);
  std::vector<double> getDDDArray(const std::string &, const DDsvalues_type &,
				  int&) const;
  int            getShift(HcalSubdetector subdet, int depth) const;
  double         getGain (HcalSubdetector subdet, int depth) const;
  unsigned       find (int element, std::vector<int> array) const;

private:

  std::vector<double> phioff;   // Phi offset for barrel, endcap, forward
  std::vector<double> etaTable; // Eta table 
  int                 nEta;     // Number of bins in eta for HB and HE
  std::vector<double> rTable;   // R-table
  int                 nR;       // Number of bins in r
  std::vector<int>    etaMin;   // Minimum eta bin number for HB/HE/HF
  std::vector<int>    etaMax;   // Maximum eta bin number for HB/HE/HF
  std::vector<double> phibin;   // Phi step for all eta bins (HB, HE and HF)
  int                 nPhi;     // Number of bins in dphi
  std::vector<int>    depth1;   // Maximum layer number for depth 1
  std::vector<int>    depth2;   // Maximum layer number for depth 2
  std::vector<int>    depth3;   // Maximum layer number for depth 3
  int                 nDepth;   // Number of bins in depth1/depth2/depth3
  std::vector<double> gainHB;   // Gain factor   for HB
  std::vector<int>    shiftHB;  // Readout shift ..  ..
  std::vector<double> gainHE;   // Gain factor   for HE
  std::vector<int>    shiftHE;  // Readout shift ..  ..
  std::vector<double> gainHF;   // Gain factor   for HF
  std::vector<int>    shiftHF;  // Readout shift ..  ..
  double              zVcal;    // Z-position  of the HF
  double              dzVcal;   // Half length of the HF
  std::vector<int>    nOff;     // Speical eta bin #'s in barrel and endcap
  std::vector<double> rHB;      // Radial positions of HB layers
  std::vector<double> zHE;      // Z-positions of HE layers
  int                 nzHB, nmodHB; // Number of halves and modules in HB
  int                 nzHE, nmodHE; // Number of halves and modules in HE
};

#endif
