#include <L1Trigger/CSCTrackFinder/interface/CSCTFSectorProcessor.h>
#include <L1Trigger/CSCTrackFinder/interface/CSCTrackFinderDataTypes.h>
#include <DataFormats/MuonDetId/interface/CSCTriggerNumbering.h>

#include <FWCore/MessageLogger/interface/MessageLogger.h>

#include "CondFormats/L1TObjects/interface/L1MuCSCTFConfiguration.h"
#include "CondFormats/DataRecord/interface/L1MuCSCTFConfigurationRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <stdlib.h>
#include <sstream>

const std::string CSCTFSectorProcessor::FPGAs[5] = {"F1","F2","F3","F4","F5"};

CSCTFSectorProcessor::CSCTFSectorProcessor(const unsigned& endcap,
					   const unsigned& sector,
					   const edm::ParameterSet& pset)
{
  m_endcap = endcap;
  m_sector = sector;
  m_bxa_on = pset.getParameter<bool>("UseBXA");
  m_extend_length = pset.getParameter<unsigned>("BXAExtendLength");
  m_latency = pset.getParameter<unsigned>("CoreLatency");
  m_minBX = pset.getParameter<int>("MinBX");
  m_maxBX = pset.getParameter<int>("MaxBX");

  std::vector<unsigned> etawins = pset.getParameter<std::vector<unsigned> >("EtaWindows");
  std::vector<unsigned>::const_iterator iter = etawins.begin();

  int i=0; 
  for(; iter!=etawins.end() && i<6; iter++,i++) m_etawin[i] = *iter;

  std::vector<unsigned> etamins = pset.getParameter<std::vector<unsigned> >("EtaMin");
  for(iter=etamins.begin(),i=0; iter!=etamins.end() && i<8; iter++,i++) m_etamin[i] = *iter;

  std::vector<unsigned> etamaxs = pset.getParameter<std::vector<unsigned> >("EtaMax");
  for(iter=etamaxs.begin(),i=0; iter!=etamaxs.end() && i<8; iter++,i++) m_etamax[i] = *iter;

  m_mindphip     = pset.getParameter<unsigned>("mindphip");
  m_mindeta_accp = pset.getParameter<unsigned>("mindeta_accp");
  m_maxdeta_accp = pset.getParameter<unsigned>("maxdeta_accp");
  m_maxdphi_accp = pset.getParameter<unsigned>("maxdphi_accp");

  try {
    edm::ParameterSet srLUTset = pset.getParameter<edm::ParameterSet>("SRLUT");
    for(i = 1; i <= 4; ++i)
      {
        if(i == 1)
          for(int j = 0; j < 2; j++)
          {
            srLUTs_[FPGAs[j]] = new CSCSectorReceiverLUT(endcap, sector, j+1, i, srLUTset);
	      }
          else
            srLUTs_[FPGAs[i]] = new CSCSectorReceiverLUT(endcap, sector, 0, i, srLUTset);
      }
      LogDebug("CSCTFSectorProcessor") << "Using stand-alone SR LUT for endcap="<<m_endcap<<", sector="<<m_sector;
  } catch(...) {
    LogDebug("CSCTFSectorProcessor") << "Using SR LUT from EventSetup for endcap="<<m_endcap<<", sector="<<m_sector;
  }

  core_ = new CSCTFSPCoreLogic();

  try {
    edm::ParameterSet ptLUTset = pset.getParameter<edm::ParameterSet>("PTLUT");
	ptLUT_ = new CSCTFPtLUT(ptLUTset);
    LogDebug("CSCTFSectorProcessor") << "Using stand-alone PT LUT for endcap="<<m_endcap<<", sector="<<m_sector;
  } catch(...){
    ptLUT_=0;
    LogDebug("CSCTFSectorProcessor") << "Using PT LUT from EventSetup for endcap="<<m_endcap<<", sector="<<m_sector;
  }
}

void CSCTFSectorProcessor::initialize(const edm::EventSetup& c){
  for(int i = 1; i <= 4; ++i)
    {
      if(i == 1)
	for(int j = 0; j < 2; j++)
	  {
		  if(!srLUTs_[FPGAs[j]]){
			  LogDebug("CSCTFSectorProcessor") << "Initializing SR LUT for endcap="<<m_endcap<<", station=1, sector="<<m_sector<<", sub_sector="<<j<<" from EventSetup";
			  srLUTs_[FPGAs[j]] = new CSCSectorReceiverLUT(m_endcap, m_sector, j+1, i, c);
		  }
	  }
      else
		  if(!srLUTs_[FPGAs[i]]){
			  LogDebug("CSCTFSectorProcessor") << "Initializing SR LUT for endcap="<<m_endcap<<", station="<<i<<", sector="<<m_sector<<" from EventSetup";
			  srLUTs_[FPGAs[i]] = new CSCSectorReceiverLUT(m_endcap, m_sector, 0, i, c);
		  }
    }

  if(!ptLUT_){
	  LogDebug("CSCTFSectorProcessor") << "Initializing PT LUT from EventSetup";
	  ptLUT_ = new CSCTFPtLUT(c);
  }

  edm::ESHandle<L1MuCSCTFConfiguration> config;
  c.get<L1MuCSCTFConfigurationRcd>().get(config);
  std::stringstream conf(config.product()->parameters());
  int eta_cnt=0;
  while( !conf.eof() ){
    char buff[1024];
    conf.getline(buff,1024);
    std::stringstream line(buff);

    std::string register_;     line>>register_;
    std::string chip_ ;        line>>chip_;
    std::string muon_;         line>>muon_;
    std::string writeValue_;   line>>writeValue_;
    std::string comments_;     std::getline(line,comments_);

    if( register_=="CNT_ETA" && chip_=="SP" ){
        unsigned int value = strtol(writeValue_.c_str(),'\0',16);
        eta_cnt = value;
    }
    if( register_=="DAT_ETA" && chip_=="SP" ){
        unsigned int value = strtol(writeValue_.c_str(),'\0',16);
        if( eta_cnt< 8                ) m_etamin[eta_cnt   ] = value;
        if( eta_cnt>=8  && eta_cnt<16 ) m_etamax[eta_cnt-8 ] = value;
		if( eta_cnt>=16 && eta_cnt<22 ) m_etawin[eta_cnt-16] = value;
        eta_cnt++;
    }
    if( register_=="CSR_SCC" && chip_=="SP" ){
        unsigned int value = strtol(writeValue_.c_str(),'\0',16);
        m_bxa_on = value&0x1;
    }
    //uint32 mindphip     = 70;
    //uint32 mindeta_accp = 4;
    //uint32 maxdeta_accp = 16;
    //uint32 maxdphi_accp = 128;
  }

}

CSCTFSectorProcessor::~CSCTFSectorProcessor()
{
  for(int i = 0; i < 5; ++i)
    {
      if(srLUTs_[FPGAs[i]]) delete srLUTs_[FPGAs[i]]; // delete the pointer
      srLUTs_[FPGAs[i]] = NULL; // point it at a safe place
    }

  delete core_;
  core_ = NULL;

  if(ptLUT_) delete ptLUT_;
  ptLUT_ = NULL;
}

bool CSCTFSectorProcessor::run(const CSCTriggerContainer<csctf::TrackStub>& stubs)
{

  if( !ptLUT_ )
    throw cms::Exception("Initialize CSC TF LUTs first (missed call to CSCTFTrackProducer::beginJob?)")<<"CSCTFSectorProcessor::run";
  for(int i = 0; i < 5; ++i)
    if(!srLUTs_[FPGAs[i]])
		throw cms::Exception("Initialize CSC TF LUTs first (missed call to CSCTFTrackProducer::beginJob?)")<<"CSCTFSectorProcessor::run";


  l1_tracks.clear();
  dt_stubs.clear();

  /** STEP ONE
   *  We take stubs from the MPC and assign their eta and phi
   *  coordinates using the SR Lookup tables.
   *  This is independent of what BX we are on so we can
   *  process one large vector of stubs.
   *  After this we append the stubs gained from the DT system.
   */

  std::vector<csctf::TrackStub> stub_vec = stubs.get();
  std::vector<csctf::TrackStub>::iterator itr = stub_vec.begin();
  std::vector<csctf::TrackStub>::const_iterator end = stub_vec.end();

  for(; itr != end; itr++)
    {
      if(itr->station() != 5)
	{
	  CSCDetId id(itr->getDetId().rawId());
	  unsigned fpga = (id.station() == 1) ? CSCTriggerNumbering::triggerSubSectorFromLabels(id) - 1 : id.station();

          lclphidat lclPhi;
          try {
            lclPhi = srLUTs_[FPGAs[fpga]]->localPhi(itr->getStrip(), itr->getPattern(), itr->getQuality(), itr->getBend());
          } catch(...) { bzero(&lclPhi,sizeof(lclPhi)); }
          gblphidat gblPhi;
          try {
            gblPhi = srLUTs_[FPGAs[fpga]]->globalPhiME(lclPhi.phi_local, itr->getKeyWG(), itr->cscid());
          } catch(...) { bzero(&gblPhi,sizeof(gblPhi)); }
          gbletadat gblEta;
          try {
            gblEta = srLUTs_[FPGAs[fpga]]->globalEtaME(lclPhi.phi_bend_local, lclPhi.phi_local, itr->getKeyWG(), itr->cscid());
          } catch(...) { bzero(&gblEta,sizeof(gblEta)); }

	  itr->setEtaPacked(gblEta.global_eta);
	  itr->setPhiPacked(gblPhi.global_phi);

	  if(itr->station() == 1) dt_stubs.push_back(*itr); // send stubs to DT

	  LogDebug("CSCTFSectorProcessor:run()") << "LCT found, processed by FPGA: " << FPGAs[fpga] << std::endl
						 << " LCT now has (eta, phi) of: (" << itr->etaValue() << "," << itr->phiValue() <<")\n";
	}
    }

  CSCTriggerContainer<csctf::TrackStub> processedStubs(stub_vec);

  /** STEP TWO
   *  We take the stubs filled by the SR LUTs and load them
   *  for processing into the SP core logic.
   *  After loading we run and then retrieve any tracks generated.
   */

  std::vector<csc::L1Track> tftks;

  core_->loadData(processedStubs, m_endcap, m_sector, m_minBX, m_maxBX);

  if( core_->run(m_endcap, m_sector, m_latency, m_etawin[0],
		 m_etamin[0], m_etamin[1], m_etamin[2], m_etamin[3],
		 m_etamin[4], m_etamin[5], m_etamin[6], m_etamin[7],
		 m_etamax[0], m_etamax[1], m_etamax[2], m_etamax[3],
		 m_etamax[4], m_etamax[5], m_etamax[6], m_etamax[7],
		 m_etawin[1], m_etawin[2], m_etawin[3],
		 m_etawin[4], m_etawin[5],
		 m_mindphip,  m_mindeta_accp,  m_maxdeta_accp, m_maxdphi_accp,
		 m_bxa_on,    m_extend_length, m_minBX, m_maxBX) )
    {
      l1_tracks = core_->tracks();
    }

  tftks = l1_tracks.get();

  /** STEP THREE
   *  Now that we have the found tracks from the core,
   *  we must assign their Pt.
   */

  std::vector<csc::L1Track>::iterator titr = tftks.begin();

  for(; titr != tftks.end(); titr++)
    {
      ptadd thePtAddress(titr->ptLUTAddress());
      ptdat thePtData = ptLUT_->Pt(thePtAddress);

      if(thePtAddress.track_fr)
	{
	  titr->setRank(thePtData.front_rank);
	  titr->setChargeValidPacked(thePtData.charge_valid_front);
	}
      else
	{
	  titr->setRank(thePtData.rear_rank);
	  titr->setChargeValidPacked(thePtData.charge_valid_rear);
	}
    }

  l1_tracks = tftks;

  return (l1_tracks.get().size() > 0);
}

