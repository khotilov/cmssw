#include "EventFilter/CSCRawToDigi/interface/CSCTMBHeader.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDMBHeader.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include <math.h>
#include <string.h> // memcpy

bool CSCTMBHeader::debug = false;
short unsigned int CSCTMBHeader::firmwareVersion=2006;

CSCTMBHeader::CSCTMBHeader() {
  firmwareVersion=2006;
  header2006.nHeaderFrames = 26;
  header2006.e0bline = 0x6E0B;
  header2006.b0cline = 0x6B0C;
  header2006.nTBins = 7; 
  header2006.nCFEBs = 5;
}

void CSCTMBHeader::setEventInformation(const CSCDMBHeader & dmbHeader) {
  header2006.cscID = dmbHeader.dmbID();
  header2006.l1aNumber = dmbHeader.l1a();
  header2006.bxnCount = dmbHeader.bxn();
}

CSCTMBHeader::CSCTMBHeader(const CSCTMBStatusDigi & digi) {
  CSCTMBHeader(digi.header());
}

CSCTMBHeader::CSCTMBHeader(const unsigned short * buf) {
  ///first determine the format
  if (buf[0]==0xDB0C) {
    firmwareVersion=2007;
  }
  else if (buf[0]==0x6B0C) {
    firmwareVersion=2006;
  }
  else {
    edm::LogError("CSCTMBHeader") <<"failed to determine TMB firmware version!!";
  }
    
  ///now fill data
  switch (firmwareVersion) {
  case 2006:
    memcpy(&header2006, buf, header2006.sizeInWords()*2);    
    break;
  case 2007:
    memcpy(&header2007, buf, header2007.sizeInWords()*2);
    break;
  default:
    edm::LogError("CSCTMBHeader")
      <<"coundn't construct: TMB firmware version is bad/not defined!";
    break;
  }
}	
    

std::vector<CSCCLCTDigi> CSCTMBHeader::CLCTDigis() const {
  std::vector<CSCCLCTDigi> result;
  
  switch (firmwareVersion) {
  case 2006: {
    ///fill digis here
    /// for the zeroth clct:
    int shape=0;
    int type=0;
  
    if ( header2006.firmRevCode < 3769 ) { //3769 is may 25 2007 - date of firmware with halfstrip only patterns 
    shape = header2006.clct0_shape;
    type  = header2006.clct0_strip_type;
    }else {//new firmware only halfstrip pattern => stripType==1 and shape is 4 bits 
      shape = ( header2006.clct0_strip_type<<3)+header2006.clct0_shape;
      type = 1;
    }
    CSCCLCTDigi digi0(header2006.clct0_valid, header2006.clct0_quality, shape, type, header2006.clct0_bend, 
		      header2006.clct0_key, (header2006.clct0_cfeb_low)|(header2006.clct0_cfeb_high<<1), 
		      header2006.clct0_bxn, 1);
    digi0.setFullBX(header2006.bxnPreTrigger);
    result.push_back(digi0);
    
    /// for the first clct:
    if ( header2006.firmRevCode < 3769 ) { 
      shape = header2006.clct1_shape;
      type  = header2006.clct1_strip_type;;
    } else {
      shape = (header2006.clct1_strip_type<<3)+header2006.clct1_shape;
      type = 1;
    }
    CSCCLCTDigi digi1(header2006.clct1_valid, header2006.clct1_quality, shape, type, header2006.clct1_bend,
		      header2006.clct1_key, (header2006.clct1_cfeb_low)|(header2006.clct1_cfeb_high<<1),
		      header2006.clct1_bxn, 2);
    digi1.setFullBX(header2006.bxnPreTrigger);
    result.push_back(digi1);
    break;
  }
  case 2007: {
    CSCCLCTDigi digi0(header2007.clct0_valid, header2007.clct0_quality, header2007.clct0_shape, 1, 
		      header2007.clct0_bend, header2007.clct0_key, 
		      (header2007.clct0_cfeb_low)|(header2007.clct0_cfeb_high<<1),
                      header2007.clct0_bxn, 1);
    digi0.setFullBX(header2007.bxnPreTrigger);
    result.push_back(digi0);
    CSCCLCTDigi digi1(header2007.clct1_valid, header2007.clct1_quality, header2007.clct1_shape, 1,
                      header2007.clct1_bend, header2007.clct1_key,
                      (header2007.clct1_cfeb_low)|(header2007.clct1_cfeb_high<<1),
                      header2007.clct1_bxn, 2);
    digi1.setFullBX(header2007.bxnPreTrigger);
    result.push_back(digi1);
    break;
  }
  default:
    edm::LogError("CSCTMBHeader")
      <<"Empty Digis: TMB firmware version is bad/not defined!"; 
    break;
  }
  return result;
}

std::vector<CSCCorrelatedLCTDigi> CSCTMBHeader::CorrelatedLCTDigis() const {
  std::vector<CSCCorrelatedLCTDigi> result;  
  switch (firmwareVersion) {
  case 2006: {
    /// for the zeroth MPC word:
    CSCCorrelatedLCTDigi digi(1, header2006.MPC_Muon0_vpf_, header2006.MPC_Muon0_quality_, 
			      header2006.MPC_Muon0_wire_, header2006.MPC_Muon0_halfstrip_clct_pattern, 
			      header2006.MPC_Muon0_clct_pattern_, header2006.MPC_Muon0_bend_, 
			      header2006.MPC_Muon0_bx_, 0, header2006.MPC_Muon0_bc0_, header2006.MPC_Muon0_SyncErr_, 
			      header2006.MPC_Muon0_cscid_low | (header2006.MPC_Muon0_cscid_bit4<<3) );
    result.push_back(digi);  
    /// for the first MPC word:
    digi = CSCCorrelatedLCTDigi(2, header2006.MPC_Muon1_vpf_, header2006.MPC_Muon1_quality_, 
				header2006.MPC_Muon1_wire_, header2006.MPC_Muon1_halfstrip_clct_pattern,
				header2006.MPC_Muon1_clct_pattern_, header2006.MPC_Muon1_bend_, 
				header2006.MPC_Muon1_bx_, 0, header2006.MPC_Muon1_bc0_, header2006.MPC_Muon1_SyncErr_,
				header2006.MPC_Muon1_cscid_low | (header2006.MPC_Muon1_cscid_bit4<<3) ); 
    result.push_back(digi);
    break;
  }
  case 2007: {
    /// for the zeroth MPC word:
    CSCCorrelatedLCTDigi digi(1, header2007.MPC_Muon0_vpf_, header2007.MPC_Muon0_quality_,
                              header2007.MPC_Muon0_wire_, header2007.MPC_Muon0_halfstrip_clct_pattern,
                              header2007.MPC_Muon0_clct_pattern_, header2007.MPC_Muon0_bend_,
                              header2007.MPC_Muon0_bx_, 0, header2007.MPC_Muon0_bc0_, header2007.MPC_Muon0_SyncErr_,
                              header2007.MPC_Muon0_cscid_low | (header2007.MPC_Muon0_cscid_bit4<<3));
    result.push_back(digi);
    /// for the first MPC word:
    digi = CSCCorrelatedLCTDigi(2, header2007.MPC_Muon1_vpf_, header2007.MPC_Muon1_quality_,
                                header2007.MPC_Muon1_wire_, header2007.MPC_Muon1_halfstrip_clct_pattern,
                                header2007.MPC_Muon1_clct_pattern_, header2007.MPC_Muon1_bend_,
                                header2007.MPC_Muon1_bx_, 0, header2007.MPC_Muon1_bc0_, header2007.MPC_Muon1_SyncErr_,
                                header2007.MPC_Muon1_cscid_low | (header2007.MPC_Muon1_cscid_bit4<<3));
    result.push_back(digi);

    break;
  }
  default:
    edm::LogError("CSCTMBHeader")
      <<"Empty CorrDigis: TMB firmware version is bad/not defined!";
    break;
  }

  return result;
}


std::ostream & operator<<(std::ostream & os, const CSCTMBHeader & hdr) {
  os << "...............TMB Header.................." << std::endl;
  os << std::hex << "BOC LINE " << hdr.header2006.b0cline << " EOB " << hdr.header2006.e0bline << std::endl;
  os << std::dec << "fifoMode = " << hdr.header2006.fifoMode 
     << ", nTBins = " << hdr.header2006.nTBins << std::endl;
  os << "dumpCFEBs = " << hdr.header2006.dumpCFEBs << ", nHeaderFrames = "
     << hdr.header2006.nHeaderFrames << std::endl;
  os << "boardID = " << hdr.header2006.boardID << ", cscID = " << hdr.header2006.cscID << std::endl;
  os << "l1aNumber = " << hdr.header2006.l1aNumber << ", bxnCount = " << hdr.header2006.bxnCount << std::endl;
  os << "preTrigTBins = " << hdr.header2006.preTrigTBins << ", nCFEBs = "<< hdr.header2006.nCFEBs<< std::endl;
  os << "trigSourceVect = " << hdr.header2006.trigSourceVect
     << ", activeCFEBs = " << hdr.header2006.activeCFEBs << std::endl;
  os << "bxnPreTrigger = " << hdr.header2006.bxnPreTrigger << std::endl;
  os << "tmbMatch = " << hdr.header2006.tmbMatch << " alctOnly = " << hdr.header2006.alctOnly
     << " clctOnly = " << hdr.header2006.clctOnly
     << " alctMatchTime = " << hdr.header2006.alctMatchTime << std::endl;
  os << "hs_thresh = " << hdr.header2006.hs_thresh << ", ds_thresh = " << hdr.header2006.ds_thresh
     << std::endl;
  os << "clct0_key = " << hdr.header2006.clct0_key << " clct0_shape = " << hdr.header2006.clct0_shape
     << " clct0_quality = " << hdr.header2006.clct0_quality << std::endl;
  os << "r_buf_nbusy = " << hdr.header2006.r_buf_nbusy << std::endl;

  os << "..................CLCT....................." << std::endl;
  return os;
}
