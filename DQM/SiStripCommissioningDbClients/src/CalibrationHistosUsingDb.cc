// Last commit: $Id: CalibrationHistosUsingDb.cc,v 1.3 2008/03/06 13:30:52 delaer Exp $

#include "DQM/SiStripCommissioningDbClients/interface/CalibrationHistosUsingDb.h"
#include "CondFormats/SiStripObjects/interface/CalibrationAnalysis.h"
#include "DataFormats/SiStripCommon/interface/SiStripConstants.h"
#include "DataFormats/SiStripCommon/interface/SiStripFecKey.h"
#include "DQM/SiStripCommon/interface/ExtractTObject.h"
#include <iostream>

using namespace sistrip;

// -----------------------------------------------------------------------------
/** */
CalibrationHistosUsingDb::CalibrationHistosUsingDb( DQMOldReceiver* mui,
					            const DbParams& params,
						    const sistrip::RunType& task )
  : CommissioningHistosUsingDb( params ),
    CalibrationHistograms( mui, task )
{
  LogTrace(mlDqmClient_) 
    << "[CalibrationHistosUsingDb::" << __func__ << "]"
    << " Constructing object...";
}

// -----------------------------------------------------------------------------
/** */
CalibrationHistosUsingDb::CalibrationHistosUsingDb( DQMOldReceiver* mui,
					      SiStripConfigDb* const db,
					      const sistrip::RunType& task )
  : CommissioningHistograms( mui, task ),
    CommissioningHistosUsingDb( db, mui, task),
    CalibrationHistograms( mui, task )
{
  LogTrace(mlDqmClient_) 
    << "[CalibrationHistosUsingDb::" << __func__ << "]"
    << " Constructing object...";
  // Load and dump the current ISHA/VFS values. This is used by the standalone analysis script
  const SiStripConfigDb::DeviceDescriptions & apvDescriptions = db->getDeviceDescriptions(APV25);
  for(SiStripConfigDb::DeviceDescriptions::const_iterator apv = apvDescriptions.begin();apv!=apvDescriptions.end();++apv) {
    apvDescription* desc = dynamic_cast<apvDescription*>( *apv );
    if ( !desc ) { continue; }
    // Retrieve device addresses from device description
    const SiStripConfigDb::DeviceAddress& addr = db->deviceAddress(*desc);
    std::stringstream bin;
    bin        << std::setw(1) << std::setfill('0') << addr.fecCrate_;
    bin << "." << std::setw(2) << std::setfill('0') << addr.fecSlot_;
    bin << "." << std::setw(1) << std::setfill('0') << addr.fecRing_;
    bin << "." << std::setw(3) << std::setfill('0') << addr.ccuAddr_;
    bin << "." << std::setw(2) << std::setfill('0') << addr.ccuChan_;
    bin << "." << desc->getAddress();
    LogTrace(mlDqmClient_) << "Present values for ISHA/VFS of APV " 
      << bin.str() << " : " 
      << static_cast<uint16_t>(desc->getIsha()) << " " << static_cast<uint16_t>(desc->getVfs());
  }
  // Load the histograms with the results
  std::string pwd = bei()->pwd();
  std::string ishaPath = pwd.substr(0,pwd.find(sistrip::root_ + "/")+sistrip::root_.size()+1);
  ishaPath += "ControlView/isha";
  ishaHistogram_ = ExtractTObject<TH1F>().extract( bei()->get(ishaPath) );
  std::string vfsPath = pwd.substr(0,pwd.find(sistrip::root_ + "/")+sistrip::root_.size()+1);
  vfsPath += "ControlView/vfs";
  vfsHistogram_ = ExtractTObject<TH1F>().extract( bei()->get(vfsPath) );

}

// -----------------------------------------------------------------------------
/** */
CalibrationHistosUsingDb::CalibrationHistosUsingDb( DQMStore* bei,
					      SiStripConfigDb* const db,
					      const sistrip::RunType& task ) 
  : CommissioningHistosUsingDb( db, task ),
    CalibrationHistograms( bei, task )
{
  LogTrace(mlDqmClient_) 
    << "[CalibrationHistosUsingDb::" << __func__ << "]"
    << " Constructing object...";
  // Load and dump the current ISHA/VFS values. This is used by the standalone analysis script
  const SiStripConfigDb::DeviceDescriptions & apvDescriptions = db->getDeviceDescriptions(APV25);
  for(SiStripConfigDb::DeviceDescriptions::const_iterator apv = apvDescriptions.begin();apv!=apvDescriptions.end();++apv) {
    apvDescription* desc = dynamic_cast<apvDescription*>( *apv );
    if ( !desc ) { continue; }
    // Retrieve device addresses from device description
    const SiStripConfigDb::DeviceAddress& addr = db->deviceAddress(*desc);
    std::stringstream bin;
    bin        << std::setw(1) << std::setfill('0') << addr.fecCrate_;
    bin << "." << std::setw(2) << std::setfill('0') << addr.fecSlot_;
    bin << "." << std::setw(1) << std::setfill('0') << addr.fecRing_;
    bin << "." << std::setw(3) << std::setfill('0') << addr.ccuAddr_;
    bin << "." << std::setw(2) << std::setfill('0') << addr.ccuChan_;
    bin << "." << desc->getAddress();
    LogTrace(mlDqmClient_) << "Present values for ISHA/VFS of APV " 
      << bin.str() << " : " 
      << static_cast<uint16_t>(desc->getIsha()) << " " << static_cast<uint16_t>(desc->getVfs());
  }
  // Load the histograms with the results
  std::string pwd = bei->pwd();
  std::string ishaPath = pwd.substr(0,pwd.find(sistrip::root_ + "/")+sistrip::root_.size()+1);
  ishaPath += "ControlView/isha";
  ishaHistogram_ = ExtractTObject<TH1F>().extract( bei->get(ishaPath) );
  std::string vfsPath = pwd.substr(0,pwd.find(sistrip::root_ + "/")+sistrip::root_.size()+1);
  vfsPath += "ControlView/vfs";
  vfsHistogram_ = ExtractTObject<TH1F>().extract( bei->get(vfsPath) );
  
}

// -----------------------------------------------------------------------------
/** */
CalibrationHistosUsingDb::~CalibrationHistosUsingDb() {
  LogTrace(mlDqmClient_) 
    << "[CalibrationHistosUsingDb::" << __func__ << "]"
    << " Destructing object...";
}

// -----------------------------------------------------------------------------
/** */
void CalibrationHistosUsingDb::uploadConfigurations() {
  
  LogTrace(mlDqmClient_)
    << "[CalibrationHistosUsingDb::" << __func__ << "]";

  if(!ishaHistogram_ && !vfsHistogram_) return;

  if ( !db() ) {
    edm::LogWarning(mlDqmClient_) 
      << "[CalibrationHistosUsingDb::" << __func__ << "]"
      << " NULL pointer to SiStripConfigDb interface!"
      << " Aborting upload...";
    return;
  }
 
  // Update all APV device descriptions with new ISHA and VFS settings
  const SiStripConfigDb::DeviceDescriptions& devices = db()->getDeviceDescriptions();
  update( const_cast<SiStripConfigDb::DeviceDescriptions&>(devices) );
  if ( doUploadConf() ) {
    edm::LogVerbatim(mlDqmClient_)
      << "[CalibrationHistosUsingDb::" << __func__ << "]"
      << " Uploading ISHA/VFS settings to DB...";
    db()->uploadDeviceDescriptions(true);
    edm::LogVerbatim(mlDqmClient_)
      << "[CalibrationHistosUsingDb::" << __func__ << "]"
      << " Uploaded ISHA/VFS settings to DB!";
  } else {
    edm::LogWarning(mlDqmClient_)
      << "[CalibrationHistosUsingDb::" << __func__ << "]"
      << " TEST only! No ISHA/VFS settings will be uploaded to DB...";
  }

  LogTrace(mlDqmClient_)
    << "[CalibrationHistosUsingDb::" << __func__ << "]"
    << " Upload of ISHA/VFS settings to DB finished!";

}

// -----------------------------------------------------------------------------
/** */
void CalibrationHistosUsingDb::update( SiStripConfigDb::DeviceDescriptions& devices ) {

  if(!ishaHistogram_ || !vfsHistogram_) return;

  // Iterate through devices and update device descriptions
  SiStripConfigDb::DeviceDescriptions::iterator idevice;
  for ( idevice = devices.begin(); idevice != devices.end(); idevice++ ) {

    // Check device type
    if ( (*idevice)->getDeviceType() != APV25 ) { continue; }

    // Cast to retrieve appropriate description object
    apvDescription* desc = dynamic_cast<apvDescription*>( *idevice );
    if ( !desc ) { continue; }

    // Retrieve the device address from device description
    const SiStripConfigDb::DeviceAddress& addr = db()->deviceAddress(*desc);

    // Construct the string for that address
    std::stringstream bin;
    bin        << std::setw(1) << std::setfill('0') << addr.fecCrate_;
    bin << "." << std::setw(2) << std::setfill('0') << addr.fecSlot_;
    bin << "." << std::setw(1) << std::setfill('0') << addr.fecRing_;
    bin << "." << std::setw(3) << std::setfill('0') << addr.ccuAddr_;
    bin << "." << std::setw(2) << std::setfill('0') << addr.ccuChan_;
    bin << "." << desc->getAddress();

    // Iterate over the histo bins and find the right one
    for(int i = 1;i <= ishaHistogram_->GetNbinsX(); ++i) {
      std::string label = ishaHistogram_->GetXaxis()->GetBinLabel(i);
      if(label == bin.str()) {
        desc->setIsha( (int)round(ishaHistogram_->GetBinContent(i)) );
      }
    }
    for(int i = 1;i <= vfsHistogram_->GetNbinsX(); ++i) {
      std::string label = vfsHistogram_->GetXaxis()->GetBinLabel(i);
      if(label == bin.str()) {
        desc->setVfs( (int)round(vfsHistogram_->GetBinContent(i)) );
      }
    }
    
  }

}

// -----------------------------------------------------------------------------
/** */
void CalibrationHistosUsingDb::create( SiStripConfigDb::AnalysisDescriptions& desc,
				     Analysis analysis) {

#ifdef USING_NEW_DATABASE_MODEL

  CalibrationAnalysis* anal = dynamic_cast<CalibrationAnalysis*>( analysis->second );
  if ( !anal ) { return; }

  SiStripFecKey fec_key( anal->fecKey() );
  SiStripFedKey fed_key( anal->fedKey() );

  for ( uint16_t iapv = 0; iapv < 2; ++iapv ) {

    // Create description
    //TODO: should also store calchan, isha and vfs.
    // the values can be accessed from the CalibrationHistogram members.
    int calchan = calchan_;
    int isha = isha_;
    int vfs = vfs_;
    CalibrationAnalysisDescription *tmp;
    tmp = new CalibrationAnalysisDescription(anal->amplitudeMean()[iapv],
                                             anal->tailMean()[iapv],
  					     anal->riseTimeMean()[iapv],
  					     anal->timeConstantMean()[iapv],
  					     anal->smearingMean()[iapv],
  					     anal->chi2Mean()[iapv],
  					     anal->deconvMode(),
  					     fec_key.fecCrate(),
  					     fec_key.fecSlot(),
  					     fec_key.fecRing(),
  					     fec_key.ccuAddr(),
  					     fec_key.ccuChan(),
  					     SiStripFecKey::i2cAddr( fec_key.lldChan(), !iapv ),
  					     db()->dbParams().partition_,
  					     db()->dbParams().runNumber_,
  					     anal->isValid(),
  					     "",
  					     fed_key.fedId(),
  					     fed_key.feUnit(),
  					     fed_key.feChan(),
  					     fed_key.fedApv() );
  
    // Add comments
    typedef std::vector<std::string> Strings;
    Strings errors = anal->getErrorCodes();
    Strings::const_iterator istr = errors.begin();
    Strings::const_iterator jstr = errors.end();
    for ( ; istr != jstr; ++istr ) { tmp->addComments( *istr ); }
  
    // Store description
    desc.push_back( tmp );
  }
  
#endif

}

