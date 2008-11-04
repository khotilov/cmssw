#include "CondFormats/L1TObjects/interface/L1MuCSCTFConfiguration.h"
#include <sstream>

edm::ParameterSet L1MuCSCTFConfiguration::parameters(int sp) const {
  edm::ParameterSet pset;
  if(sp>=12) return pset;

  pset.addParameter<int>("CoreLatency",8);
  std::vector<unsigned int> etamin(8), etamax(8), etawin(6);

  int eta_cnt=0;
  std::stringstream conf(registers[sp]);
  while( !conf.eof() ){
    char buff[1024];
    conf.getline(buff,1024);
    std::stringstream line(buff);

    std::string register_;     line>>register_;
    std::string chip_;         line>>chip_;
    std::string muon_;         line>>muon_;
    std::string writeValue_;   line>>writeValue_;
    std::string comments_;     std::getline(line,comments_);

    if( register_=="CSR_REQ" && chip_=="SP" ){
        unsigned int value = strtol(writeValue_.c_str(),'\0',16);
        pset.addParameter<bool>("run_core",       value&0x8000);
        pset.addParameter<bool>("trigger_on_ME1a",value&0x0001);
        pset.addParameter<bool>("trigger_on_ME1b",value&0x0002);
        pset.addParameter<bool>("trigger_on_ME2", value&0x0004);
        pset.addParameter<bool>("trigger_on_ME3", value&0x0008);
        pset.addParameter<bool>("trigger_on_ME4", value&0x0010);
        pset.addParameter<bool>("trigger_on_MB1a",value&0x0100);
        pset.addParameter<bool>("trigger_on_MB1d",value&0x0200);
    }
    if( register_=="DAT_FTR" && chip_=="SP" ){
        unsigned int value = strtol(writeValue_.c_str(),'\0',16);
        pset.addParameter<unsigned int>("singlesTrackPt",value); // 0x1F - rank, 0x60 - Q1,Q0, 0x80 - charge
    }
    if( register_=="CSR_SFC" && chip_=="SP" ){
        unsigned int value = strtol(writeValue_.c_str(),'\0',16);
        pset.addParameter<unsigned int>("singlesTrackOutput",(value&0x3000)>>12);
    }
    if( register_=="CNT_ETA" && chip_=="SP" ){
        unsigned int value = strtol(writeValue_.c_str(),'\0',16);
        eta_cnt = value;
    }
	if( register_=="CSR_SCC" && chip_=="SP" ){
        unsigned int value = strtol(writeValue_.c_str(),'\0',16);
        pset.addParameter<unsigned int>("BXAdepth",      value&0x3     );
        pset.addParameter<bool>        ("AllowALCTonly",(value&0x10)>>4);
        pset.addParameter<bool>        ("AllowCLCTonly",(value&0x20)>>5);
        pset.addParameter<bool>        ("useDT",        (value&0x80)>>8);
        pset.addParameter<unsigned int>("PreTrigger",   (value&0x300)>>8);
    }
    if( register_=="DAT_ETA" && chip_=="SP" ){
        unsigned int value = strtol(writeValue_.c_str(),'\0',16);
        if( eta_cnt< 8                ) etamin[eta_cnt   ] = value;
        if( eta_cnt>=8  && eta_cnt<16 ) etamax[eta_cnt-8 ] = value;
        if( eta_cnt>=16 && eta_cnt<22 ) etawin[eta_cnt-16] = value;
		// 4 line below is just an exaple (need to verify a sequence):
        if( eta_cnt==22 ) pset.addParameter<unsigned int>("mindphip",    value);
        if( eta_cnt==23 ) pset.addParameter<unsigned int>("mindeta_accp",value);
        if( eta_cnt==24 ) pset.addParameter<unsigned int>("maxdeta_accp",value);
        if( eta_cnt==25 ) pset.addParameter<unsigned int>("maxdphi_accp",value);
        eta_cnt++;
    }
    if( register_=="CSR_LQE" && chip_=="F1" && muon_=="M1" )
        pset.addParameter<unsigned int>("QualityEnableME1a",strtol(writeValue_.c_str(),'\0',16));
    if( register_=="CSR_LQE" && chip_=="F1" && muon_=="M2" )
        pset.addParameter<unsigned int>("QualityEnableME1b",strtol(writeValue_.c_str(),'\0',16));
    if( register_=="CSR_LQE" && chip_=="F1" && muon_=="M3" )
        pset.addParameter<unsigned int>("QualityEnableME1c",strtol(writeValue_.c_str(),'\0',16));
    if( register_=="CSR_LQE" && chip_=="F2" && muon_=="M1" )
        pset.addParameter<unsigned int>("QualityEnableME1d",strtol(writeValue_.c_str(),'\0',16));
    if( register_=="CSR_LQE" && chip_=="F2" && muon_=="M2" )
        pset.addParameter<unsigned int>("QualityEnableME1e",strtol(writeValue_.c_str(),'\0',16));
    if( register_=="CSR_LQE" && chip_=="F2" && muon_=="M3" )
        pset.addParameter<unsigned int>("QualityEnableME1f",strtol(writeValue_.c_str(),'\0',16));
    if( register_=="CSR_LQE" && chip_=="F3" && muon_=="M1" )
        pset.addParameter<unsigned int>("QualityEnableME2a",strtol(writeValue_.c_str(),'\0',16));
    if( register_=="CSR_LQE" && chip_=="F3" && muon_=="M2" )
        pset.addParameter<unsigned int>("QualityEnableME2b",strtol(writeValue_.c_str(),'\0',16));
    if( register_=="CSR_LQE" && chip_=="F3" && muon_=="M3" )
        pset.addParameter<unsigned int>("QualityEnableME2c",strtol(writeValue_.c_str(),'\0',16));
    if( register_=="CSR_LQE" && chip_=="F4" && muon_=="M1" )
        pset.addParameter<unsigned int>("QualityEnableME3a",strtol(writeValue_.c_str(),'\0',16));
    if( register_=="CSR_LQE" && chip_=="F4" && muon_=="M2" )
        pset.addParameter<unsigned int>("QualityEnableME3b",strtol(writeValue_.c_str(),'\0',16));
    if( register_=="CSR_LQE" && chip_=="F4" && muon_=="M3" )
        pset.addParameter<unsigned int>("QualityEnableME3c",strtol(writeValue_.c_str(),'\0',16));
    if( register_=="CSR_LQE" && chip_=="F5" && muon_=="M1" )
        pset.addParameter<unsigned int>("QualityEnableME4a",strtol(writeValue_.c_str(),'\0',16));
    if( register_=="CSR_LQE" && chip_=="F5" && muon_=="M2" )
        pset.addParameter<unsigned int>("QualityEnableME4b",strtol(writeValue_.c_str(),'\0',16));
    if( register_=="CSR_LQE" && chip_=="F5" && muon_=="M3" )
        pset.addParameter<unsigned int>("QualityEnableME4c",strtol(writeValue_.c_str(),'\0',16));

    if( register_=="CSR_KFL" )//&& chip_=="SP" && muon_=="MA" )
        pset.addParameter<unsigned int>("kill_fiber",strtol(writeValue_.c_str(),'\0',16));
  }

  if( eta_cnt     ) pset.addParameter< std::vector<unsigned int> >("EtaMin",etamin);
  if( eta_cnt>=8  ) pset.addParameter< std::vector<unsigned int> >("EtaMax",etamax);
  if( eta_cnt>=16 ) pset.addParameter< std::vector<unsigned int> >("EtaWindows",etawin);

  return pset;
}
