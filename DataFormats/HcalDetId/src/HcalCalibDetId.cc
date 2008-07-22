#include "DataFormats/HcalDetId/interface/HcalCalibDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h" 
#include "FWCore/Utilities/interface/Exception.h"

HcalCalibDetId::HcalCalibDetId() : HcalOtherDetId() {
}


HcalCalibDetId::HcalCalibDetId(uint32_t rawid) : HcalOtherDetId(rawid) {
}

HcalCalibDetId::HcalCalibDetId(HcalSubdetector subdet, int ieta, int iphi, int ctype) : HcalOtherDetId(HcalCalibration) {
  id_|=(CalibrationBox<<17); // Calibration Category, bits [17:19] (= "1" for CalibrationBox)
  id_|=(ctype&0xF)           // calibration channel type, bits [0:3]
      |((iphi&0x1F)<<4)      // phi index, bits [4:10]
      |(((ieta+2)&0x7)<<11)     // eta index, bits [11:13]
      |((subdet&0x7)<<14);   // subdetector, bits [14:16]
}

HcalCalibDetId::HcalCalibDetId(int ieta, int iphi) : HcalOtherDetId(HcalCalibration) {
  id_|=(HOCrosstalk<<17); // Calibration Category, bits [17:19] (= "2" for HOX)
  id_|=(iphi&0x3F)               // phi index, bits [0:6]
      |((abs(ieta)&0xF)<<7)     // eta index, bits [7:10]
      |(((ieta > 0)?(1):(0))<<11); // z side, bit [11]
}

HcalCalibDetId::HcalCalibDetId(const DetId& gen) {
  if (!gen.null() && (gen.det()!=Hcal || gen.subdetId()!=HcalOther)) {
    throw cms::Exception("Invalid DetId") << "Cannot initialize HcalCalibDetId from " << std::hex << gen.rawId() << std::dec; 
  }
  id_=gen.rawId();
  if (subdet()!=HcalCalibration) {
    throw cms::Exception("Invalid DetId") << "Cannot initialize HcalCalibDetId from " << std::hex << gen.rawId() << std::dec; 
  }
}

HcalCalibDetId& HcalCalibDetId::operator=(const DetId& gen) {
  if (!gen.null() && (gen.det()!=Hcal || gen.subdetId()!=HcalOther)) {
    throw cms::Exception("Invalid DetId") << "Cannot assign HcalCalibDetId from " << std::hex << gen.rawId() << std::dec; 
  }
  id_=gen.rawId();
  if (subdet()!=HcalCalibration) {
    throw cms::Exception("Invalid DetId") << "Cannot assign HcalCalibDetId from " << std::hex << gen.rawId() << std::dec; 
  }
  return *this;
}

int HcalCalibDetId::cboxChannel() const {
  return (calibFlavor()==CalibrationBox)?(id_&0xF):(0);
}

HcalSubdetector HcalCalibDetId::hcalSubdet() const {  
  return (HcalSubdetector)((calibFlavor()==CalibrationBox)?((id_>>14)&0x7):(0));
}
    
int HcalCalibDetId::ieta() const {
  return (calibFlavor()==CalibrationBox)?(((id_>>11)&0x7)-2):((calibFlavor()==HOCrosstalk)?(((id_>>7)&0xF)*zside()):(0));
}

int HcalCalibDetId::iphi() const {
  return (calibFlavor()==CalibrationBox)?((id_>>4)&0x3F):((calibFlavor()==HOCrosstalk)?(id_&0x3F):(0));
}

int HcalCalibDetId::zside() const {
  return (calibFlavor()==HOCrosstalk)?(((id_>>11)&0x1)?(1):(-1)):(0);
}

std::string HcalCalibDetId::cboxChannelString() const {
  switch (cboxChannel()) {
  case(cbox_MixerHigh): return "Mixer-High";
  case(cbox_MixerLow): return "Mixer-Low";
  case(cbox_LaserMegatile): return "Megatile";
  case(cbox_RadDam_Layer0_RM4): return "RadDam-L0-RM4";
  case(cbox_RadDam_Layer7_RM4): return "RadDam-L7-RM4";
  case(cbox_RadDam_Layer0_RM1): return "RadDam-L0-RM1";
  case(cbox_RadDam_Layer7_RM1): return "RadDam-L7-RM1";
  case(cbox_HOCrosstalkPIN): return "HO-Crosstalk-PIN";
  case(cbox_HF_ScintillatorPIN): return "HF-Scint-PIN";
  default : return "";
  }
}

std::ostream& operator<<(std::ostream& s,const HcalCalibDetId& id) {
  std::string sd;
  switch (id.hcalSubdet()) {
    case(HcalBarrel) : sd="HB"; break;
    case(HcalEndcap) : sd="HE"; break;
    case(HcalOuter) : sd="HO"; break;
    case(HcalForward) : sd="HF"; break;
    default: break;
  }
  switch (id.calibFlavor()) {
  case(HcalCalibDetId::CalibrationBox):
    return s << "(HcalCalibBox " << sd << ' ' << id.ieta() << "," << id.iphi()
	     << ' ' << id.cboxChannelString() << ')';
  case(HcalCalibDetId::HOCrosstalk):
    return s << "(HOCrosstalk "  << id.ieta() << "," << id.iphi() 
	     << ')';
  default: return s;
  };
}
