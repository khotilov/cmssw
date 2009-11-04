#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/SiStripCommon/interface/SiStripFedKey.h"
#include "DataFormats/SiStripCommon/interface/ConstantsForHardwareSystems.h"

#include "DQM/SiStripMonitorHardware/interface/FEDHistograms.hh"



FEDHistograms::FEDHistograms()
{
  dqm_ = 0;
  histogramConfig_.clear();
}

FEDHistograms::~FEDHistograms()
{
}
  
void FEDHistograms::initialise(const edm::ParameterSet& iConfig,
			       std::ostringstream* pDebugStream
			       )
{
  getConfigForHistogram("DataPresent",iConfig,pDebugStream);
  getConfigForHistogram("AnyFEDErrors",iConfig,pDebugStream);
  getConfigForHistogram("AnyDAQProblems",iConfig,pDebugStream);
  getConfigForHistogram("AnyFEProblems",iConfig,pDebugStream);
  getConfigForHistogram("CorruptBuffers",iConfig,pDebugStream);
  getConfigForHistogram("BadChannelStatusBits",iConfig,pDebugStream);
  getConfigForHistogram("BadActiveChannelStatusBits",iConfig,pDebugStream);
  
  getConfigForHistogram("FEOverflows",iConfig,pDebugStream);
  getConfigForHistogram("FEMissing",iConfig,pDebugStream);
  getConfigForHistogram("BadMajorityAddresses",iConfig,pDebugStream);
  
  getConfigForHistogram("DataMissing",iConfig,pDebugStream);
  getConfigForHistogram("BadIDs",iConfig,pDebugStream);
  getConfigForHistogram("BadDAQPacket",iConfig,pDebugStream);
  getConfigForHistogram("InvalidBuffers",iConfig,pDebugStream);
  getConfigForHistogram("BadDAQCRCs",iConfig,pDebugStream);
  getConfigForHistogram("BadFEDCRCs",iConfig,pDebugStream);
  
  getConfigForHistogram("FEOverflowsDetailed",iConfig,pDebugStream);
  getConfigForHistogram("FEMissingDetailed",iConfig,pDebugStream);
  getConfigForHistogram("BadMajorityAddressesDetailed",iConfig,pDebugStream);
  getConfigForHistogram("BadAPVStatusBitsDetailed",iConfig,pDebugStream);
  getConfigForHistogram("APVErrorBitsDetailed",iConfig,pDebugStream);
  getConfigForHistogram("APVAddressErrorBitsDetailed",iConfig,pDebugStream);
  getConfigForHistogram("UnlockedBitsDetailed",iConfig,pDebugStream);
  getConfigForHistogram("OOSBitsDetailed",iConfig,pDebugStream);
  
  getConfigForHistogram("nFEDErrors",iConfig,pDebugStream);
  getConfigForHistogram("nFEDDAQProblems",iConfig,pDebugStream);
  getConfigForHistogram("nFEDsWithFEProblems",iConfig,pDebugStream);
  getConfigForHistogram("nFEDCorruptBuffers",iConfig,pDebugStream);
  getConfigForHistogram("nBadChannelStatusBits",iConfig,pDebugStream);
  getConfigForHistogram("nBadActiveChannelStatusBits",iConfig,pDebugStream);
  getConfigForHistogram("nFEDsWithFEOverflows",iConfig,pDebugStream);
  getConfigForHistogram("nFEDsWithMissingFEs",iConfig,pDebugStream);
  getConfigForHistogram("nFEDsWithFEBadMajorityAddresses",iConfig,pDebugStream);

  getConfigForHistogram("nFEDErrorsvsTime",iConfig,pDebugStream);
  getConfigForHistogram("nFEDCorruptBuffersvsTime",iConfig,pDebugStream);
  getConfigForHistogram("nFEDsWithFEProblemsvsTime",iConfig,pDebugStream);

  getConfigForHistogram("nUnconnectedChannels",iConfig,pDebugStream);

  getConfigForHistogram("nTotalBadChannels",iConfig,pDebugStream);
  getConfigForHistogram("nTotalBadActiveChannels",iConfig,pDebugStream);

  getConfigForHistogram("nTotalBadChannelsvsTime",iConfig,pDebugStream);
  getConfigForHistogram("nTotalBadActiveChannelsvsTime",iConfig,pDebugStream);

  getConfigForHistogram("nAPVStatusBit",iConfig,pDebugStream);
  getConfigForHistogram("nAPVError",iConfig,pDebugStream);
  getConfigForHistogram("nAPVAddressError",iConfig,pDebugStream);
  getConfigForHistogram("nUnlocked",iConfig,pDebugStream);
  getConfigForHistogram("nOutOfSync",iConfig,pDebugStream);

  getConfigForHistogram("nAPVStatusBitvsTime",iConfig,pDebugStream);
  getConfigForHistogram("nAPVErrorvsTime",iConfig,pDebugStream);
  getConfigForHistogram("nAPVAddressErrorvsTime",iConfig,pDebugStream);
  getConfigForHistogram("nUnlockedvsTime",iConfig,pDebugStream);
  getConfigForHistogram("nOutOfSyncvsTime",iConfig,pDebugStream);

  tkMapConfigName_ = "TkHistoMap";
  getConfigForHistogram(tkMapConfigName_,iConfig,pDebugStream);


}

void FEDHistograms::fillCountersHistograms(const FEDErrors::FEDCounters & fedLevelCounters, 
					   const FEDErrors::ChannelCounters & chLevelCounters, 
					   const double aTime )
{
  fillHistogram(nFEDErrors_,fedLevelCounters.nFEDErrors);
  fillHistogram(nFEDDAQProblems_,fedLevelCounters.nDAQProblems);
  fillHistogram(nFEDsWithFEProblems_,fedLevelCounters.nFEDsWithFEProblems);
  fillHistogram(nFEDCorruptBuffers_,fedLevelCounters.nCorruptBuffers);
  fillHistogram(nFEDsWithFEOverflows_,fedLevelCounters.nFEDsWithFEOverflows);
  fillHistogram(nFEDsWithFEBadMajorityAddresses_,fedLevelCounters.nFEDsWithFEBadMajorityAddresses);
  fillHistogram(nFEDsWithMissingFEs_,fedLevelCounters.nFEDsWithMissingFEs);
  fillHistogram(nBadChannelStatusBits_,fedLevelCounters.nBadChannels);
  fillHistogram(nBadActiveChannelStatusBits_,fedLevelCounters.nBadActiveChannels);

  fillHistogram(nFEDErrorsvsTime_,aTime,fedLevelCounters.nFEDErrors);
  fillHistogram(nFEDCorruptBuffersvsTime_,aTime,fedLevelCounters.nCorruptBuffers);
  fillHistogram(nFEDsWithFEProblemsvsTime_,aTime,fedLevelCounters.nFEDsWithFEProblems);

  fillHistogram(nUnconnectedChannels_,chLevelCounters.nNotConnected);

  fillHistogram(nTotalBadChannels_,fedLevelCounters.nTotalBadChannels);
  fillHistogram(nTotalBadActiveChannels_,fedLevelCounters.nTotalBadActiveChannels);

  fillHistogram(nTotalBadChannelsvsTime_,aTime,fedLevelCounters.nTotalBadChannels);
  fillHistogram(nTotalBadActiveChannelsvsTime_,aTime,fedLevelCounters.nTotalBadActiveChannels);
  
  fillHistogram(nAPVStatusBit_,chLevelCounters.nAPVStatusBit);
  fillHistogram(nAPVError_,chLevelCounters.nAPVError);
  fillHistogram(nAPVAddressError_,chLevelCounters.nAPVAddressError);
  fillHistogram(nUnlocked_,chLevelCounters.nUnlocked);
  fillHistogram(nOutOfSync_,chLevelCounters.nOutOfSync);

  fillHistogram(nAPVStatusBitvsTime_,aTime,chLevelCounters.nAPVStatusBit);
  fillHistogram(nAPVErrorvsTime_,aTime,chLevelCounters.nAPVError);
  fillHistogram(nAPVAddressErrorvsTime_,aTime,chLevelCounters.nAPVAddressError);
  fillHistogram(nUnlockedvsTime_,aTime,chLevelCounters.nUnlocked);
  fillHistogram(nOutOfSyncvsTime_,aTime,chLevelCounters.nOutOfSync);

}

void FEDHistograms::fillFEDHistograms(FEDErrors & aFedErr, 
				      bool lFullDebug)
{
  const FEDErrors::FEDLevelErrors & lFedLevelErrors = aFedErr.getFEDLevelErrors();
  const unsigned int lFedId = aFedErr.fedID();

  if (lFedLevelErrors.DataPresent) fillHistogram(dataPresent_,lFedId);

  if (lFedLevelErrors.HasCabledChannels && lFedLevelErrors.DataMissing) fillHistogram(dataMissing_,lFedId);
  
  if (lFedLevelErrors.InvalidBuffers) fillHistogram(invalidBuffers_,lFedId);
  else if (lFedLevelErrors.BadFEDCRCs) fillHistogram(badFEDCRCs_,lFedId);
  else if (lFedLevelErrors.BadDAQCRCs) fillHistogram(badDAQCRCs_,lFedId);
  else if (lFedLevelErrors.BadIDs) fillHistogram(badIDs_,lFedId);
  else if (lFedLevelErrors.BadDAQPacket) fillHistogram(badDAQPacket_,lFedId);
  else if (lFedLevelErrors.CorruptBuffer) fillHistogram(corruptBuffers_,lFedId);

  if (aFedErr.anyFEDErrors()) fillHistogram(anyFEDErrors_,lFedId);
  if (lFedLevelErrors.HasCabledChannels && aFedErr.anyDAQProblems()) fillHistogram(anyDAQProblems_,lFedId);
  if (aFedErr.anyFEProblems()) fillHistogram(anyFEProblems_,lFedId);

  if (lFedLevelErrors.FEsOverflow) fillHistogram(feOverflows_,lFedId);
  if (lFedLevelErrors.FEsMissing) fillHistogram(feMissing_,lFedId);
  if (lFedLevelErrors.FEsBadMajorityAddress) fillHistogram(badMajorityAddresses_,lFedId);

  if (lFedLevelErrors.BadChannelStatusBit) fillHistogram(badChannelStatusBits_,lFedId);
  if (lFedLevelErrors.BadActiveChannelStatusBit) fillHistogram(badActiveChannelStatusBits_,lFedId);

  std::vector<FEDErrors::FELevelErrors> & lFeVec = aFedErr.getFELevelErrors();
  for (unsigned int iFe(0); iFe<lFeVec.size(); iFe++){
    fillFEHistograms(lFedId,lFeVec.at(iFe));
  }

  std::vector<FEDErrors::ChannelLevelErrors> & lChVec = aFedErr.getChannelLevelErrors();
  for (unsigned int iCh(0); iCh < lChVec.size(); iCh++){
    fillChannelsHistograms(lFedId,lChVec.at(iCh),lFullDebug);
  }

  std::vector<FEDErrors::APVLevelErrors> & lAPVVec = aFedErr.getAPVLevelErrors();
  for (unsigned int iApv(0); iApv < lAPVVec.size(); iApv++){
    fillAPVsHistograms(lFedId,lAPVVec.at(iApv),lFullDebug);
  }


}

//fill a histogram if the pointer is not NULL (ie if it has been booked)
void FEDHistograms::fillFEHistograms(const unsigned int aFedId, 
				     const FEDErrors::FELevelErrors & aFeLevelErrors)
{
  const unsigned short lFeId = aFeLevelErrors.FeID;
  bookFEDHistograms(aFedId);
  if (aFeLevelErrors.Overflow) fillHistogram(feOverflowDetailed_[aFedId],lFeId);
  else if (aFeLevelErrors.Missing) fillHistogram(feMissingDetailed_[aFedId],lFeId);
  else if (aFeLevelErrors.BadMajorityAddress) fillHistogram(badMajorityAddressDetailed_[aFedId],lFeId);
}

//fill a histogram if the pointer is not NULL (ie if it has been booked)
void FEDHistograms::fillChannelsHistograms(const unsigned int aFedId, 
					   const FEDErrors::ChannelLevelErrors & aChErr, 
					   bool fullDebug)
{
  unsigned int lChId = aChErr.ChannelID;
  bookFEDHistograms(aFedId,fullDebug);
  if (aChErr.Unlocked) {
    fillHistogram(unlockedDetailed_[aFedId],lChId);
  }
  if (aChErr.OutOfSync) {
    fillHistogram(outOfSyncDetailed_[aFedId],lChId);    
  }
}


void FEDHistograms::fillAPVsHistograms(const unsigned int aFedId, 
				       const FEDErrors::APVLevelErrors & aAPVErr, 
				       bool fullDebug)
{
  unsigned int lChId = aAPVErr.APVID;
  bookFEDHistograms(aFedId,fullDebug);
  if (aAPVErr.APVStatusBit) fillHistogram(badStatusBitsDetailed_[aFedId],lChId);
  if (aAPVErr.APVError) fillHistogram(apvErrorDetailed_[aFedId],lChId);
  if (aAPVErr.APVAddressError) fillHistogram(apvAddressErrorDetailed_[aFedId],lChId);
}

void FEDHistograms::bookTopLevelHistograms(DQMStore* dqm)
{
  //get FED IDs
  const unsigned int siStripFedIdMin = FEDNumbering::MINSiStripFEDID;
  const unsigned int siStripFedIdMax = FEDNumbering::MAXSiStripFEDID;

  //get the pointer to the dqm object
  dqm_ = dqm;

  //book FED level histograms
  histosBooked_.resize(siStripFedIdMax+1,false);
  debugHistosBooked_.resize(siStripFedIdMax+1,false);


  //book histos
  dataPresent_ = bookHistogram("DataPresent","DataPresent",
                               "Number of events where the data from a FED is seen",
                               siStripFedIdMax-siStripFedIdMin+1,
                               siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");

  dataMissing_ = bookHistogram("DataMissing","DataMissing",
                               "Number of events where the data from a FED with cabled channels is missing",
                               siStripFedIdMax-siStripFedIdMin+1,
                               siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");

  anyFEDErrors_ = bookHistogram("AnyFEDErrors","AnyFEDErrors",
				"Number of buffers with any FED error (excluding bad channel status bits, FE problems except overflows) per FED",
				siStripFedIdMax-siStripFedIdMin+1,
				siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");

  corruptBuffers_ = bookHistogram("CorruptBuffers","CorruptBuffers",
                                  "Number of corrupt FED buffers per FED",
                                  siStripFedIdMax-siStripFedIdMin+1,
                                  siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");

  invalidBuffers_ = bookHistogram("InvalidBuffers","InvalidBuffers",
                                  "Number of invalid FED buffers per FED",
                                  siStripFedIdMax-siStripFedIdMin+1,
                                  siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");

  anyDAQProblems_ = bookHistogram("AnyDAQProblems","AnyDAQProblems",
                                  "Number of buffers with any problems flagged in DAQ header (including CRC)",
                                  siStripFedIdMax-siStripFedIdMin+1,
                                  siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");

  badIDs_ = bookHistogram("BadIDs","BadIDs",
                          "Number of buffers with non-SiStrip source IDs in DAQ header",
                          siStripFedIdMax-siStripFedIdMin+1,
                          siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");

  badChannelStatusBits_ = bookHistogram("BadChannelStatusBits","BadChannelStatusBits",
                                        "Number of buffers with one or more enabled channel with bad status bits",
                                        siStripFedIdMax-siStripFedIdMin+1,
                                        siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");

  badActiveChannelStatusBits_ = bookHistogram("BadActiveChannelStatusBits","BadActiveChannelStatusBits",
                                              "Number of buffers with one or more active channel with bad status bits",
                                              siStripFedIdMax-siStripFedIdMin+1,
                                              siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");

  anyFEProblems_ = bookHistogram("AnyFEProblems","AnyFEProblems",
				 "Number of buffers with any FE unit problems",
				 siStripFedIdMax-siStripFedIdMin+1,
				 siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");
  
  badDAQCRCs_ = bookHistogram("BadDAQCRCs","BadDAQCRCs",
                              "Number of buffers with bad CRCs from the DAQ",
                              siStripFedIdMax-siStripFedIdMin+1,
                              siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");

  badFEDCRCs_ = bookHistogram("BadFEDCRCs","BadFEDCRCs",
                              "Number of buffers with bad CRCs from the FED",
                              siStripFedIdMax-siStripFedIdMin+1,
                              siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");

  badDAQPacket_ = bookHistogram("BadDAQPacket","BadDAQPacket",
                               "Number of buffers with (non-CRC) problems flagged in DAQ header/trailer",
                               siStripFedIdMax-siStripFedIdMin+1,
                               siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");

  feOverflows_ = bookHistogram("FEOverflows","FEOverflows",
                               "Number of buffers with one or more FE overflow",
                               siStripFedIdMax-siStripFedIdMin+1,
                               siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");

  badMajorityAddresses_ = bookHistogram("BadMajorityAddresses","BadMajorityAddresses",
                                        "Number of buffers with one or more FE with a bad majority APV address",
                                        siStripFedIdMax-siStripFedIdMin+1,
                                        siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");

  feMissing_ = bookHistogram("FEMissing","FEMissing",
                             "Number of buffers with one or more FE unit payload missing",
                             siStripFedIdMax-siStripFedIdMin+1,
                             siStripFedIdMin-0.5,siStripFedIdMax+0.5,"FED-ID");
  
  nFEDErrors_ = bookHistogram("nFEDErrors",
			      "nFEDErrors",
                              "Number of FEDs with errors (exclusing channel status bits) per event",
			      "# FEDErrors");

  nFEDDAQProblems_ = bookHistogram("nFEDDAQProblems",
				   "nFEDDAQProblems",
                                   "Number of FEDs with DAQ problems per event",
				   "# FEDDAQProblems");

  nFEDsWithFEProblems_ = bookHistogram("nFEDsWithFEProblems",
				       "nFEDsWithFEProblems",
                                       "Number of FEDs with FE problems per event",
				       "# FEDs with FE problems");

  nFEDCorruptBuffers_ = bookHistogram("nFEDCorruptBuffers",
				      "nFEDCorruptBuffers",
                                      "Number of FEDs with corrupt buffers per event",
				      "# FEDs with corrupt buffer");

  nBadChannelStatusBits_ = bookHistogram("nBadChannelStatusBits",
					 "nBadChannelStatusBits",
					 "Number of channels with bad status bits per event",
					 "# bad enabled channels");

  nBadActiveChannelStatusBits_ = bookHistogram("nBadActiveChannelStatusBits",
					       "nBadActiveChannelStatusBits",
                                               "Number of active channels with bad status bits per event",
					       "# bad active channels");

  nFEDsWithFEOverflows_ = bookHistogram("nFEDsWithFEOverflows",
					"nFEDsWithFEOverflows",
                                        "Number FEDs with FE units which overflowed per event",
					"# FEDs with FE overflows");

  nFEDsWithFEBadMajorityAddresses_ = bookHistogram("nFEDsWithFEBadMajorityAddresses",
						   "nFEDsWithFEBadMajorityAddresses",
                                                   "Number of FEDs with FE units with a bad majority address per event",
						   "# FEDs with bad address");

  nFEDsWithMissingFEs_ = bookHistogram("nFEDsWithMissingFEs",
				       "nFEDsWithMissingFEs",
                                       "Number of FEDs with missing FE unit payloads per event",
				       "# FEDs with missing FEs");

  nUnconnectedChannels_ = bookHistogram("nUnconnectedChannels",
					"nUnconnectedChannels",
					"Number of channels not connected per event",
					"# unconnected channels");

  nAPVStatusBit_ = bookHistogram("nAPVStatusBit",
			       "nAPVStatusBit",
			       "Number of APVs with APVStatusBit error per event",
			       "# APVs with APVStatusBit error");

  nAPVError_ = bookHistogram("nAPVError",
			   "nAPVError",
			   "Number of APVs with APVError per event",
			   "#APVs with APVError");

  nAPVAddressError_ = bookHistogram("nAPVAddressError",
				  "nAPVAddressError",
				  "Number of APVs with APVAddressError per event",
				  "#APVs with APVAddressError");

  nUnlocked_ = bookHistogram("nUnlocked",
			     "nUnlocked",
			     "Number of channels Unlocked per event",
			     "# channels unlocked");

  nOutOfSync_ = bookHistogram("nOutOfSync",
			      "nOutOfSync",
			      "Number of channels OutOfSync per event",
			      "# channels out-of-sync");

  nTotalBadChannels_ = bookHistogram("nTotalBadChannels",
				     "nTotalBadChannels",
				     "Number of channels with any error",
				     "Total # bad enabled channels");

  nTotalBadActiveChannels_ = bookHistogram("nTotalBadActiveChannels",
					   "nTotalBadActiveChannels",
					   "Number of active channels with any error",
					   "Total # bad active channels");


  nTotalBadChannelsvsTime_ = bookProfile("nTotalBadChannelsvsTime",
					 "nTotalBadChannelsvsTime",
					 "Number of channels with any error vs time",
					 0,
					 42241, //total number of channels
					 "Time",
					 "Total # bad enabled channels"
					 );


  nTotalBadActiveChannelsvsTime_  = bookProfile("nTotalBadActiveChannelsvsTime",
						"nTotalBadActiveChannelsvsTime",
						"Number of active channels with any error vs time",
						0,
						42241, //total number of channels
						"Time",
						"Total # bad active channels"
						);


  nFEDErrorsvsTime_ = bookProfile("nFEDErrorsvsTime",
				  "nFEDErrorsvsTime",
				  "Number of FEDs with any error vs time",
				  0,
				  42241, //total number of channels
				  "Time",
				  "# FEDErrors"
				  );

  nFEDCorruptBuffersvsTime_ = bookProfile("nFEDCorruptBuffersvsTime",
					  "nFEDCorruptBuffersvsTime",
					  "Number of FEDs with corrupt buffer vs time",
					  0,
					  42241, //total number of channels
					  "Time",
					  "# FEDCorruptBuffer"
					  );

  nFEDsWithFEProblemsvsTime_ = bookProfile("nFEDsWithFEProblemsvsTime",
					   "nFEDsWithFEProblemsvsTime",
					   "Number of FEDs with any FE error vs time",
					   0,
					   42241, //total number of channels
					   "Time",
					   "# FEDsWithFEProblems"
					   );

  nAPVStatusBitvsTime_ = bookProfile("nAPVStatusBitvsTime",
				     "nAPVStatusBitvsTime",
				     "Number of APVs with APVStatusBit error vs time",
				     0,
				     42241, //total number of channels
				     "Time",
				     "# APVs with APVStatusBit error"
				     );

  nAPVErrorvsTime_ = bookProfile("nAPVErrorvsTime",
				 "nAPVErrorvsTime",
				 "Number of APVs with APVError vs time",
				 0,
				 42241, //total number of channels
				 "Time",
				 "# APVs with APVError"
				 );

  nAPVAddressErrorvsTime_ = bookProfile("nAPVAddressErrorvsTime",
					"nAPVAddressErrorvsTime",
					"Number of APVs with APVAddressError vs time",
					0,
					42241, //total number of channels
					"Time",
					"# APVs with APVAddressError"
					);

  nUnlockedvsTime_ = bookProfile("nUnlockedvsTime",
				 "nUnlockedvsTime",
				 "Number of channels Unlocked vs time",
				 0,
				 42241, //total number of channels
				 "Time",
				 "# channels unlocked "
				 );

  nOutOfSyncvsTime_ = bookProfile("nOutOfSyncvsTime",
				  "nOutOfSyncvsTime",
				  "Number of channels OutOfSync vs time",
				  0,
				  42241, //total number of channels
				  "Time",
				  "# channels out-of-sync"
				  );


  //book map after, as it creates a new folder...
  if (histogramConfig_[tkMapConfigName_].enabled){
    //const std::string dqmPath = dqm_->pwd();
    tkmapFED_ = new TkHistoMap("SiStrip/TkHisto","TkHMap_FractionOfBadChannels",0.,1);
  }
  else tkmapFED_ = 0;

}

void FEDHistograms::bookFEDHistograms(unsigned int fedId,
				      bool fullDebugMode
				      )
{
  if (!histosBooked_[fedId]) {
    //will do that only once
    SiStripFedKey fedKey(fedId,0,0,0);
    std::stringstream fedIdStream;
    fedIdStream << fedId;
    dqm_->setCurrentFolder(fedKey.path());
    feOverflowDetailed_[fedId] = bookHistogram("FEOverflowsDetailed",
                                               "FEOverflowsForFED"+fedIdStream.str(),
                                               "FE overflows per FE unit for FED ID "+fedIdStream.str(),
                                               sistrip::FEUNITS_PER_FED,0,sistrip::FEUNITS_PER_FED,
                                               "FE-Index");
    badMajorityAddressDetailed_[fedId] = bookHistogram("BadMajorityAddressesDetailed",
                                                       "BadMajorityAddressesForFED"+fedIdStream.str(),
                                                       "Bad majority APV addresses per FE unit for FED ID "+fedIdStream.str(),
                                                       sistrip::FEUNITS_PER_FED,0,sistrip::FEUNITS_PER_FED,
                                                       "FE-Index");
    feMissingDetailed_[fedId] = bookHistogram("FEMissingDetailed",
                                              "FEMissingForFED"+fedIdStream.str(),
                                              "Buffers with FE Unit payload missing per FE unit for FED ID "+fedIdStream.str(),
                                              sistrip::FEUNITS_PER_FED,0,sistrip::FEUNITS_PER_FED,
                                              "FE-Index");
    badStatusBitsDetailed_[fedId] = bookHistogram("BadAPVStatusBitsDetailed",
                                                  "BadAPVStatusBitsForFED"+fedIdStream.str(),
                                                  "Bad apv status bits for FED ID "+fedIdStream.str(),
                                                  sistrip::APVS_PER_FED,0,sistrip::APVS_PER_FED,
                                                  "APV-Index");
     histosBooked_[fedId] = true;
  }
  if (fullDebugMode && !debugHistosBooked_[fedId]) {
    //will do that only once
    SiStripFedKey fedKey(fedId,0,0,0);
    std::stringstream fedIdStream;
    fedIdStream << fedId;
    dqm_->setCurrentFolder(fedKey.path());

    apvErrorDetailed_[fedId] = bookHistogram("APVErrorBitsDetailed",
                                             "APVErrorBitsForFED"+fedIdStream.str(),
                                             "APV errors for FED ID "+fedIdStream.str(),
                                             sistrip::APVS_PER_FED,0,sistrip::APVS_PER_FED,
                                             "APV-Index");
    apvAddressErrorDetailed_[fedId] = bookHistogram("APVAddressErrorBitsDetailed",
                                                    "APVAddressErrorBitsForFED"+fedIdStream.str(),
                                                    "Wrong APV address errors for FED ID "+fedIdStream.str(),
                                                    sistrip::APVS_PER_FED,0,sistrip::APVS_PER_FED,
                                                    "APV-Index");
    unlockedDetailed_[fedId] = bookHistogram("UnlockedBitsDetailed",
                                             "UnlockedBitsForFED"+fedIdStream.str(),
                                             "Unlocked channels for FED ID "+fedIdStream.str(),
                                             sistrip::FEDCH_PER_FED,0,sistrip::FEDCH_PER_FED,
                                             "Channel-Index");
    outOfSyncDetailed_[fedId] = bookHistogram("OOSBitsDetailed",
                                              "OOSBitsForFED"+fedIdStream.str(),
                                              "Out of sync channels for FED ID "+fedIdStream.str(),
                                              sistrip::FEDCH_PER_FED,0,sistrip::FEDCH_PER_FED,
                                              "Channel-Index");
    debugHistosBooked_[fedId] = true;
  }
}

void FEDHistograms::bookAllFEDHistograms()
{
  //get FED IDs
  const unsigned int siStripFedIdMin = FEDNumbering::MINSiStripFEDID;
  const unsigned int siStripFedIdMax = FEDNumbering::MAXSiStripFEDID;
  //book them
  for (unsigned int iFed = siStripFedIdMin; iFed <= siStripFedIdMax; iFed++) {
    bookFEDHistograms(iFed,true);
  }
}

std::string FEDHistograms::tkHistoMapName(unsigned int aIndex){
  return tkMapConfigName_;
}

TkHistoMap * FEDHistograms::tkHistoMapPointer(unsigned int aIndex){
  return tkmapFED_;
}
