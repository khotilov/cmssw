// Last commit: $Id: SiStripConfigDb.cc,v 1.41 2007/12/05 17:05:21 bainbrid Exp $

#include "OnlineDB/SiStripConfigDb/interface/SiStripConfigDb.h"
#include "DataFormats/SiStripCommon/interface/SiStripConstants.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace sistrip;

// -----------------------------------------------------------------------------
// 
uint32_t SiStripConfigDb::cntr_ = 0;

// -----------------------------------------------------------------------------
// 
SiStripConfigDb::SiStripConfigDb( const edm::ParameterSet& pset,
				  const edm::ActivityRegistry& activity ) :
  factory_(0), 
  dbParams_(),
  // Local cache
  devices_(), 
  feds_(), 
  connections_(), 
  dcuDetIdMap_(), 
  fedIds_(),
  // Misc
  usingStrips_(true),
  openConnection_(false)
{
  cntr_++;
  LogTrace(mlConfigDb_)
    << "[SiStripConfigDb::" << __func__ << "]"
    << " Constructing object using Service..."
    << " (Class instance: " << cntr_ << ")";

  // Set DB connection parameters
  dbParams_.reset();
  dbParams_.setParams( pset );
  
  stringstream ss;
  ss << "[SiStripConfigDb::" << __func__ << "]"
     << " Database connection parameters:" << endl
     << dbParams_ << endl;
  edm::LogVerbatim(mlConfigDb_) << ss.str();
  
  // Open connection
  openDbConnection();
  
}

// -----------------------------------------------------------------------------
// 
SiStripConfigDb::SiStripConfigDb( string confdb, 
				  string partition,
				  uint32_t major,
				  uint32_t minor ) :
  factory_(0), 
  dbParams_(),
  // Local cache
  devices_(), 
  feds_(), 
  connections_(), 
  dcuDetIdMap_(), 
  fedIds_(),
  // Misc
  usingStrips_(true),
  openConnection_(false)
{
  cntr_++;
  LogTrace(mlConfigDb_)
    << "[SiStripConfigDb::" << __func__ << "]"
    << " Constructing object..."
    << " (Class instance: " << cntr_ << ")";

  dbParams_.reset();
  dbParams_.usingDb_ = true; 
  dbParams_.confdb( confdb );
  dbParams_.partition_ = partition; 
  dbParams_.cabMajor_ = major;
  dbParams_.cabMinor_ = minor;
  dbParams_.fedMajor_ = major; //@@
  dbParams_.fedMinor_ = minor; //@@
  dbParams_.fecMajor_ = major; //@@
  dbParams_.fecMinor_ = minor; //@@
  dbParams_.dcuMajor_ = major; //@@
  dbParams_.dcuMinor_ = minor; //@@
  dbParams_.calMajor_ = major; //@@
  dbParams_.calMinor_ = minor; //@@

  stringstream ss;
  ss << "[SiStripConfigDb::" << __func__ << "]"
     << " Database connection parameters:" << endl
     << dbParams_ << endl;
  edm::LogVerbatim(mlConfigDb_) << ss.str();

}

// -----------------------------------------------------------------------------
// 
SiStripConfigDb::SiStripConfigDb( string user, 
				  string passwd, 
				  string path,
				  string partition,
				  uint32_t major,
				  uint32_t minor ) :
  factory_(0), 
  dbParams_(),
  // Local cache
  devices_(), 
  feds_(), 
  connections_(), 
  dcuDetIdMap_(), 
  fedIds_(), 
  // Misc
  usingStrips_(true),
  openConnection_(false)
{
  cntr_++;
  LogTrace(mlConfigDb_)
    << "[SiStripConfigDb::" << __func__ << "]"
    << " Constructing object..."
    << " (Class instance: " << cntr_ << ")";
  
  dbParams_.reset();
  dbParams_.usingDb_ = true; 
  dbParams_.confdb( user, passwd, path );
  dbParams_.partition_ = partition; 
  dbParams_.cabMajor_ = major;
  dbParams_.cabMinor_ = minor;
  dbParams_.fedMajor_ = major; //@@
  dbParams_.fedMinor_ = minor; //@@
  dbParams_.fecMajor_ = major; //@@
  dbParams_.fecMinor_ = minor; //@@
  dbParams_.calMajor_ = major; //@@
  dbParams_.calMinor_ = minor; //@@

  stringstream ss;
  ss << "[SiStripConfigDb::" << __func__ << "]"
     << " Database connection parameters:" << endl
     << dbParams_ << endl;
  edm::LogVerbatim(mlConfigDb_) << ss.str();

}

// -----------------------------------------------------------------------------
// 
SiStripConfigDb::SiStripConfigDb( string input_module_xml,
				  string input_dcuinfo_xml,
				  vector<string> input_fec_xml,
				  vector<string> input_fed_xml,
				  string output_module_xml,
				  string output_dcuinfo_xml,
				  string output_fec_xml,
				  string output_fed_xml ) : 
  factory_(0), 
  dbParams_(),
  // Local cache
  devices_(), 
  feds_(), 
  connections_(), 
  dcuDetIdMap_(),
  fedIds_(), 
  // Misc
  usingStrips_(true),
  openConnection_(false)
{
  cntr_++;
  LogTrace(mlConfigDb_)
    << "[SiStripConfigDb::" << __func__ << "]"
    << " Constructing object..."
    << " (Class instance: " << cntr_ << ")";

  dbParams_.reset();
  dbParams_.usingDb_ = false; 
  dbParams_.inputModuleXml_ = input_module_xml; 
  dbParams_.inputDcuInfoXml_ = input_dcuinfo_xml; 
  dbParams_.inputFecXml_ = input_fec_xml; 
  dbParams_.inputFedXml_ = input_fed_xml; 
  dbParams_.outputModuleXml_ = output_module_xml; 
  dbParams_.outputDcuInfoXml_ = output_dcuinfo_xml; 
  dbParams_.outputFecXml_ = output_fec_xml; 
  dbParams_.outputFedXml_ = output_fed_xml; 

  stringstream ss;
  ss << "[SiStripConfigDb::" << __func__ << "]"
     << " Database connection parameters:" << endl
     << dbParams_ << endl;
  edm::LogVerbatim(mlConfigDb_) << ss.str();

}

// -----------------------------------------------------------------------------
//
SiStripConfigDb::~SiStripConfigDb() {
  if ( openConnection_ ) { closeDbConnection(); }
  LogTrace(mlConfigDb_)
    << "[SiStripConfigDb::" << __func__ << "]"
    << " Destructing object...";
  if ( cntr_ ) { cntr_--; }
}

// -----------------------------------------------------------------------------
// 
SiStripConfigDb::DbParams::DbParams() :
  usingDb_(false),
  user_(null_),
  passwd_(null_),
  path_(null_),
  partition_(null_), 
  runNumber_(0),
  cabMajor_(0),
  cabMinor_(0),
  fedMajor_(0),
  fedMinor_(0),
  fecMajor_(0),
  fecMinor_(0),
  calMajor_(0),
  calMinor_(0),
  dcuMajor_(0),
  dcuMinor_(0),
  force_(true),
  inputModuleXml_(""),
  inputDcuInfoXml_(""),
  inputFecXml_(),
  inputFedXml_(),
  inputDcuConvXml_(""),
  outputModuleXml_("/tmp/module.xml"),
  outputDcuInfoXml_("/tmp/dcuinfo.xml"),
  outputFecXml_("/tmp/fec.xml"),
  outputFedXml_("/tmp/fed.xml")
{
  reset();
}

// -----------------------------------------------------------------------------
// 
SiStripConfigDb::DbParams::~DbParams() {
  inputFecXml_.clear();
  inputFedXml_.clear();
}

// -----------------------------------------------------------------------------
// 
void SiStripConfigDb::DbParams::reset() {
  usingDb_ = true;
  confdb_ = null_;
  confdb( confdb_ );
  partition_ = null_;
  runNumber_ = 0;
  cabMajor_ = 0;
  cabMinor_ = 0;
  fedMajor_ = 0;
  fedMinor_ = 0;
  fecMajor_ = 0;
  fecMinor_ = 0;
  calMajor_ = 0;
  calMinor_ = 0;
  dcuMajor_ = 0;
  dcuMinor_ = 0;
  force_ = true;
  inputModuleXml_ = "";
  inputDcuInfoXml_ = "";
  inputFecXml_ = vector<string>(1,"");
  inputFedXml_ = vector<string>(1,"");
  inputDcuConvXml_ = "";
  outputModuleXml_ = "";
  outputDcuInfoXml_ = "";
  outputFecXml_ = "";
  outputFedXml_ = "";
}

// -----------------------------------------------------------------------------
// 
void SiStripConfigDb::DbParams::setParams( const edm::ParameterSet& pset ) {
  reset();
  usingDb_ = pset.getUntrackedParameter<bool>("UsingDb",true); 
  confdb( pset.getUntrackedParameter<string>("ConfDb","") );
  partition_ = pset.getUntrackedParameter<string>("Partition","");
  runNumber_ = pset.getUntrackedParameter<unsigned int>("RunNumber",0);
  cabMajor_ = pset.getUntrackedParameter<unsigned int>("MajorVersion",0);
  cabMinor_ = pset.getUntrackedParameter<unsigned int>("MinorVersion",0);
  fedMajor_ = pset.getUntrackedParameter<unsigned int>("FedMajorVersion",0);
  fedMinor_ = pset.getUntrackedParameter<unsigned int>("FedMinorVersion",0);
  fecMajor_ = pset.getUntrackedParameter<unsigned int>("FecMajorVersion",0);
  fecMinor_ = pset.getUntrackedParameter<unsigned int>("FecMinorVersion",0);
  dcuMajor_ = pset.getUntrackedParameter<unsigned int>("DcuDetIdMajorVersion",0);
  dcuMinor_ = pset.getUntrackedParameter<unsigned int>("DcuDetIdMinorVersion",0);
  calMajor_ = pset.getUntrackedParameter<unsigned int>("CalibMajorVersion",0);
  calMinor_ = pset.getUntrackedParameter<unsigned int>("CalibMinorVersion",0);
  force_ = pset.getUntrackedParameter<bool>("ForceDcuDetIdVersions",true);
  inputModuleXml_ = pset.getUntrackedParameter<string>("InputModuleXml","");
  inputDcuInfoXml_ = pset.getUntrackedParameter<string>("InputDcuInfoXml",""); 
  inputFecXml_ = pset.getUntrackedParameter< vector<string> >( "InputFecXml", vector<string>(1,"") ); 
  inputFedXml_ = pset.getUntrackedParameter< vector<string> >( "InputFedXml", vector<string>(1,"") ); 
  inputDcuConvXml_ = pset.getUntrackedParameter<string>( "InputDcuConvXml","" );
  outputModuleXml_ = pset.getUntrackedParameter<string>("OutputModuleXml","/tmp/module.xml");
  outputDcuInfoXml_ = pset.getUntrackedParameter<string>("OutputDcuInfoXml","/tmp/dcuinfo.xml");
  outputFecXml_ = pset.getUntrackedParameter<string>( "OutputFecXml", "/tmp/fec.xml" );
  outputFedXml_ = pset.getUntrackedParameter<string>( "OutputFedXml", "/tmp/fed.xml" );
}

// -----------------------------------------------------------------------------
// 
void SiStripConfigDb::DbParams::confdb( const string& confdb ) {
  confdb_ = confdb;
  uint32_t ipass = confdb.find("/");
  uint32_t ipath = confdb.find("@");
  if ( ipass != string::npos && 
       ipath != string::npos ) {
    user_   = confdb.substr(0,ipass); 
    passwd_ = confdb.substr(ipass+1,ipath-ipass-1); 
    path_   = confdb.substr(ipath+1,confdb.size());
  } else {
    user_   = null_;
    passwd_ = null_;
    path_   = null_;
  }
}

// -----------------------------------------------------------------------------
// 
void SiStripConfigDb::DbParams::confdb( const string& user,
					const string& passwd,
					const string& path ) {
  if ( user != "" && passwd != "" && path != "" &&
       user != null_ && passwd != null_ && path != null_ ) {
    user_   = user;
    passwd_ = passwd;
    path_   = path;
  } else {
    user_   = null_;
    passwd_ = null_;
    path_   = null_;
  }
  confdb_ = user_ + "/" + passwd_ + "@" + path_;
}

// -----------------------------------------------------------------------------
// 
void SiStripConfigDb::DbParams::print( stringstream& ss ) const {
  ss << " Using database            : " << std::boolalpha << usingDb_ << std::noboolalpha << endl
     << " ConfDb                    : " << confdb_ << endl
     << " User/Passwd@Path          : " << user_ << "/" << passwd_ << "@" << path_ << endl
     << " Partition                 : " << partition_ << endl
     << " Run number                : " << runNumber_ << endl
     << " Cabling major/minor vers  : " << cabMajor_ << "." << cabMinor_ << endl
     << " FED major/minor vers      : " << fedMajor_ << "." << fedMinor_ << endl
     << " FEC major/minor vers      : " << fecMajor_ << "." << fecMinor_ << endl
     << " Calibration maj/min vers  : " << calMajor_ << "." << calMinor_ << endl
     << " DCU-DetId maj/min vers    : " << dcuMajor_ << "." << dcuMinor_;
  if ( force_ ) { ss << " (version not overridden by run number)"; }
  ss << endl;
  // Input
  ss << " Input \"module.xml\" file   : " << inputModuleXml_ << endl
     << " Input \"dcuinfo.xml\" file  : " << inputDcuInfoXml_ << endl
     << " Input \"fec.xml\" file(s)   : ";
  vector<string>::const_iterator ifec = inputFecXml_.begin();
  for ( ; ifec != inputFecXml_.end(); ifec++ ) { ss << *ifec << ", "; }
  ss << endl;
  ss << " Input \"fed.xml\" file(s)   : ";
  vector<string>::const_iterator ifed = inputFedXml_.begin();
  for ( ; ifed != inputFedXml_.end(); ifed++ ) { ss << *ifed << ", "; }
  ss << endl;
  // Output 
  ss << " Output \"module.xml\" file  : " << outputModuleXml_ << endl
     << " Output \"dcuinfo.xml\" file : " << outputDcuInfoXml_ << endl
     << " Output \"fec.xml\" file(s)  : " << outputFecXml_ << endl
     << " Output \"fed.xml\" file(s)  : " << outputFedXml_ << endl;
}

// -----------------------------------------------------------------------------
// 
ostream& operator<< ( ostream& os, const SiStripConfigDb::DbParams& params ) {
  stringstream ss;
  params.print(ss);
  os << ss.str();
  return os;
}

// -----------------------------------------------------------------------------
// 
void SiStripConfigDb::openDbConnection() {

  edm::LogVerbatim(mlConfigDb_) 
    << "[SiStripConfigDb::" << __func__ << "]"
    << " Opening connection to database...";

  // Check
  if ( openConnection_ ) {
    edm::LogWarning(mlConfigDb_) 
      << "[SiStripConfigDb::" << __func__ << "]"
      << " Connection already open!";
    return;
  }
  openConnection_ = true;
  
  // Establish database connection
  if ( dbParams_.usingDb_ ) { 
    usingDatabase(); 
  } else { 
    usingXmlFiles(); 
  }

  devices_.clear();
  feds_.clear();
  connections_.clear();
  dcuDetIdMap_.clear();
  fedIds_.clear();
  
}

// -----------------------------------------------------------------------------
//
void SiStripConfigDb::closeDbConnection() {

  edm::LogVerbatim(mlConfigDb_) 
    << "[SiStripConfigDb::" << __func__ << "]"
    << " Closing connection to database...";

  // Check
  if ( !openConnection_ ) {
    edm::LogWarning(mlConfigDb_) 
      << "[SiStripConfigDb::" << __func__ << "]"
      << " No connection open!";
    return;
  }
  openConnection_ = true;

  try { 
    if ( factory_ ) { delete factory_; }
  } catch (...) { handleException( __func__, "Attempting to close database connection..." ); }
  factory_ = 0; 

}

// -----------------------------------------------------------------------------
//
DeviceFactory* const SiStripConfigDb::deviceFactory( string method_name ) const { 
  if ( factory_ ) { return factory_; }
  else { 
    if ( method_name != "" ) { 
      stringstream ss;
      ss << "[SiStripConfigDb::" << method_name << "]"
	 << " NULL pointer to DeviceFactory!";
      edm::LogWarning(mlConfigDb_) << ss.str();
    }
    return 0;
  }
}

// -----------------------------------------------------------------------------
//
void SiStripConfigDb::usingDatabase() {
  
  // Retrieve connection params from CONFDB env. var. and override .cfg values 
  string user = "";
  string passwd = "";
  string path = "";
  DbAccess::getDbConfiguration( user, passwd, path );
  if ( user != "" && passwd != "" && path != "" ) {
    stringstream ss;
    ss << "[SiStripConfigDb::" << __func__ << "]"
       << " Setting \"user/passwd@path\" to \""
       << user << "/" << passwd << "@" << path
       << "\" using 'CONFDB' environmental variable";
    if ( dbParams_.user_ != null_ || 
	 dbParams_.passwd_ != null_ || 
	 dbParams_.path_ != null_ ) { 
      ss << " (Overwriting existing value of \""
	 << dbParams_.user_ << "/" 
	 << dbParams_.passwd_ << "@" 
	 << dbParams_.path_ 
	 << "\" read from .cfg file)";
    }
    LogTrace(mlConfigDb_) << ss.str() << endl;
    dbParams_.confdb( user, passwd, path );
  } else if ( dbParams_.user_ != null_ && 
	      dbParams_.passwd_ != null_ && 
	      dbParams_.path_ != null_ ) { 
    stringstream ss;
    ss << "[SiStripConfigDb::" << __func__ << "]"
       << " Setting \"user/passwd@path\" to \""
       << dbParams_.user_ << "/" 
       << dbParams_.passwd_ << "@" 
       << dbParams_.path_ 
       << "\" using 'ConfDb' configurable read from .cfg file";
    LogTrace(mlConfigDb_) << ss.str();
  } else {
    edm::LogWarning(mlConfigDb_)
      << "[SiStripConfigDb::" << __func__ << "]"
      << " Unable to retrieve 'user/passwd@path' parameters"
      << " from 'CONFDB' environmental variable or .cfg file"
      << " (present value is \"" 
      << user << "/" 
      << passwd << "@" 
      << path 
      << "\"). Aborting connection to database...";
    return;
  }
  
  // Retrieve partition name from ENV_CMS_TK_PARTITION env. var. and override .cfg value
  string partition = "ENV_CMS_TK_PARTITION";
  if ( getenv(partition.c_str()) != NULL ) { 
    stringstream ss;
    ss << "[SiStripConfigDb::" << __func__ << "]"
       << " Setting \"PARTITION\" to \""
       << getenv( partition.c_str() )
       << "\" using 'ENV_CMS_TK_PARTITION' environmental variable";
    if ( dbParams_.partition_ != "" ) { 
      ss << " (Overwriting existing value of \""
	 << dbParams_.partition_
	 << "\" read from .cfg file)";
    }
    LogTrace(mlConfigDb_) << ss.str() << endl;
    dbParams_.partition_ = getenv( partition.c_str() ); 
  } else if ( dbParams_.partition_ != "" ) {
    std::stringstream ss;
    ss << "[SiStripConfigDb::" << __func__ << "]"
       << " Setting \"partition\" to \""
       << dbParams_.partition_
       << "\" using 'Partition' configurable read from .cfg file";
    LogTrace(mlConfigDb_) << ss.str();
  } else { 
    edm::LogWarning(mlConfigDb_)
      << "[SiStripConfigDb::" << __func__ << "]"
      << " Unable to retrieve 'partition' parameter"
      << " from 'CONFDB' environmental variable or .cfg file!"
      << " Aborting connection to database...";
    return;
  } 

  // Check TNS_ADMIN environmental variable
  std::string pattern = "TNS_ADMIN";
  std::string tns_admin = "/afs/cern.ch/project/oracle/admin";
  if ( getenv( pattern.c_str() ) != NULL ) { 
    tns_admin = getenv( pattern.c_str() ); 
    LogTrace(mlConfigDb_)
      << "[SiStripConfigDb::" << __func__ << "]"
      << " TNS_ADMIN is set to: \"" 
      << tns_admin << "\"";
  } else {
    edm::LogWarning(mlConfigDb_)
      << "[SiStripConfigDb::" << __func__ << "]"
      << " TNS_ADMIN is not set!"
      << " Trying to use /afs and setting to: \"" 
      << tns_admin << "\"";
  }
  
  // Remove trailing slash and set TNS_ADMIN
  if ( tns_admin.empty() ) { tns_admin = "."; }
  std::string slash = tns_admin.substr( tns_admin.size()-1, 1 ); 
  if ( slash == sistrip::dir_ ) { tns_admin = tns_admin.substr( 0, tns_admin.size()-1 ); }
  setenv( pattern.c_str(), tns_admin.c_str(), 1 ); 
  
  // Check if database is found in tnsnames.ora file
  std::string filename( tns_admin + "/tnsnames.ora" ); 
  std::ifstream tnsnames_ora( filename.c_str() );
  bool ok = false;
  if ( tnsnames_ora.is_open() ) {
    std::string line;
    while ( !tnsnames_ora.eof() ) {
      getline( tnsnames_ora, line );
      if ( !dbParams_.path_.empty() && 
	   line.find( dbParams_.path_ ) != std::string::npos ) { ok = true; }
    }
  } else {
    edm::LogWarning(mlConfigDb_)
      << "[SiStripConfigDb::" << __func__ << "]"
      << " Cannot open file \""
      << filename << "\"";
  }

  if ( ok ) {
    LogTrace(mlConfigDb_)
      << "[SiStripConfigDb::" << __func__ << "]"
      << " Found database account \"" 
      << dbParams_.path_ << "\" in file \""
      << filename << "\"!";
  } else {
    edm::LogWarning(mlConfigDb_)
      << "[SiStripConfigDb::" << __func__ << "]"
      << " Cannot find database account \"" 
      << dbParams_.path_ << "\" in file \""
      << filename << "\""
      << " Aborting connection to database...";
    return; 
  }
  
  // Create device factory object
  try { 
    edm::LogVerbatim(mlConfigDb_)
      << "[SiStripConfigDb::" << __func__ << "]"
      << " Creating DeviceFactory object...";
    factory_ = new DeviceFactory( dbParams_.user_, 
				  dbParams_.passwd_, 
				  dbParams_.path_ ); 
    edm::LogVerbatim(mlConfigDb_)
      << "[SiStripConfigDb::" << __func__ << "]"
      << " Created DeviceFactory object!";
  } catch (...) { 
    stringstream ss; 
    ss << "Failed to connect to database using parameters '" 
       << dbParams_.user_ << "/" 
       << dbParams_.passwd_ << "@" 
       << dbParams_.path_ 
       << "' and partition '" 
       << dbParams_.partition_ << "'";
    handleException( __func__, ss.str() );
  }
  
  // Check for valid pointer to DeviceFactory
  if ( deviceFactory(__func__) ) { 
    stringstream ss;
    ss << "[SiStripConfigDb::" << __func__ << "]"
       << " DeviceFactory created at address 0x" 
       << hex << setw(8) << setfill('0') << factory_ << dec
       << ", using database account with parameters '" 
       << dbParams_.user_ << "/" 
       << dbParams_.passwd_ << "@" 
       << dbParams_.path_
       << "' and partition '" 
       << dbParams_.partition_ << "'";
    LogTrace(mlConfigDb_) << ss.str();
  } else {
    edm::LogError(mlConfigDb_)
      << "[SiStripConfigDb::" << __func__ << "]"
      << " NULL pointer to DeviceFactory!"
      << " Unable to connect to database using connection parameters '" 
      << dbParams_.user_ << "/" 
      << dbParams_.passwd_ << "@" 
      << dbParams_.path_
      << "' and partition '" 
      << dbParams_.partition_ << "'";
    return; 
  }
  
  try { 
    deviceFactory(__func__)->setUsingDb( dbParams_.usingDb_ ); 
  } catch (...) { 
    handleException( __func__, "Attempted to 'setUsingDb'" );
  }
  
  // Retrieve versioning from DB for a non-zero run numbers
  if ( dbParams_.runNumber_ ) {
    TkRun* run = deviceFactory(__func__)->getRun( dbParams_.partition_, 
						  dbParams_.runNumber_ );
    if ( run ) {
      if ( run->getRunNumber() ) {
#ifdef USING_NEW_DATABASE_MODEL
	dbParams_.cabMajor_ = run->getConnectionVersionMajorId();
	dbParams_.cabMinor_ = run->getConnectionVersionMinorId();
#endif
	dbParams_.fecMajor_ = run->getFecVersionMajorId();
	dbParams_.fecMinor_ = run->getFecVersionMinorId();
	dbParams_.fedMajor_ = run->getFedVersionMajorId();
	dbParams_.fedMinor_ = run->getFedVersionMinorId();
#ifdef USING_NEW_DATABASE_MODEL
	if ( !dbParams_.force_ ) { //@@ check if forcing versions specified in .cfi
	  dbParams_.dcuMajor_ = run->getDcuInfoVersionMajorId();
	  dbParams_.dcuMinor_ = run->getDcuInfoVersionMinorId();
	}
#endif
	std::stringstream ss;
	ss << "[SiStripConfigDb::" << __func__ << "]"
	   << " Versions overridden using values retrieved for specific run number: "
	   << std::endl << dbParams_;
	edm::LogVerbatim(mlConfigDb_) << ss.str();
      } else {
	edm::LogWarning(mlConfigDb_)
	  << "[SiStripConfigDb::" << __func__ << "]"
	  << " NULL run number returned!";
      }
    } else {
      edm::LogWarning(mlConfigDb_)
	<< "[SiStripConfigDb::" << __func__ << "]"
	<< " NULL pointer to TkRun object!"
	<< " Unable to retrieve versions for run number "
	<< dbParams_.runNumber_;
    }
  }
  
  // DCU-DetId 
  try { 
#ifdef USING_NEW_DATABASE_MODEL
    deviceFactory(__func__)->addDetIdPartition( dbParams_.partition_,
						dbParams_.dcuMajor_, 
						dbParams_.dcuMinor_ );
#else
    deviceFactory(__func__)->addDetIdPartition( dbParams_.partition_ );
#endif
  } catch (...) { 
    stringstream ss;
    ss << "Attempted to 'addDetIdPartition; for partition: " << dbParams_.partition_;
    handleException( __func__, ss.str() );
  }

#ifndef USING_NEW_DATABASE_MODEL
  // FED-FEC connections
  try { 
    deviceFactory(__func__)->createInputDBAccess();
  } catch (...) { 
    handleException( __func__, "Attempted to 'createInputDBAccess' for FED-FEC connections!" );
  }
  
  try {
    deviceFactory(__func__)->setInputDBVersion( dbParams_.partition_,
						dbParams_.cabMajor_,
						dbParams_.cabMinor_ );
  } catch (...) { 
    stringstream ss;
    ss << "Attempted to 'setInputDBVersion; for partition: " << dbParams_.partition_;
    handleException( __func__, ss.str() ); 
  }
#endif
  
}

// -----------------------------------------------------------------------------
//
void SiStripConfigDb::usingXmlFiles() {
  LogTrace(mlConfigDb_)
    << "[SiStripConfigDb::" << __func__ << "]"
    << " Using XML description files...";

  // Create device factory object
  try { 
    factory_ = new DeviceFactory(); 
  } catch (...) { 
    handleException( __func__, "Attempting to create DeviceFactory for use with xml files" );
  }
  
 // Check for valid pointer to DeviceFactory
  if ( deviceFactory(__func__) ) { 
    stringstream ss;
    ss << "[SiStripConfigDb::" << __func__ << "]"
       << " DeviceFactory created at address 0x" 
       << hex << setw(8) << setfill('0') << factory_ << dec
       << ", using XML description files";
    LogTrace(mlConfigDb_) << ss.str();
    cout << ss.str() << endl;
  } else {    
    edm::LogError(mlConfigDb_)
      << "[SiStripConfigDb::" << __func__ << "]"
      << " NULL pointer to DeviceFactory!"
      << " Unable to connect to database!";
    return; 
  }
  
  try { 
    deviceFactory(__func__)->setUsingDb( dbParams_.usingDb_ );
  } catch (...) { 
    handleException( __func__, "Attempted to 'setUsingDb'" );
  }

  try { 
#ifndef USING_NEW_DATABASE_MODEL
    deviceFactory(__func__)->createInputFileAccess();
#endif
  } catch (...) { 
    handleException( __func__, "Attempted to 'createInputFileAccess'" ); 
  }
  
  // Input module.xml file
  if ( dbParams_.inputModuleXml_ == "" ) {
    edm::LogWarning(mlConfigDb_)
      << "[SiStripConfigDb::" << __func__ << "]"
      << " NULL path to input 'module.xml' file!";
  } else {
    if ( checkFileExists( dbParams_.inputModuleXml_ ) ) { 
      try { 
#ifdef USING_NEW_DATABASE_MODEL
	deviceFactory(__func__)->setConnectionInputFileName( dbParams_.inputModuleXml_ ); 
#else
	deviceFactory(__func__)->setFedFecConnectionInputFileName( dbParams_.inputModuleXml_ ); 
#endif
      } catch (...) { 
	handleException( __func__ ); 
      }
      LogTrace(mlConfigDb_)
	<< "[SiStripConfigDb::" << __func__ << "]"
	<< " Added input 'module.xml' file: " << dbParams_.inputModuleXml_;
    } else {
      edm::LogWarning(mlConfigDb_)
	<< "[SiStripConfigDb::" << __func__ << "]"
	<< " No 'module.xml' file found at " << dbParams_.inputModuleXml_;
      dbParams_.inputModuleXml_ = ""; 
    }
  }
  
  // Input dcuinfo.xml file
  if ( dbParams_.inputDcuInfoXml_ == "" ) {
    edm::LogWarning(mlConfigDb_)
      << "[SiStripConfigDb::" << __func__ << "]"
      << " NULL path to input 'dcuinfo.xml' file!";
  } else { 
    if ( checkFileExists( dbParams_.inputDcuInfoXml_ ) ) { 
      try { 
	deviceFactory(__func__)->setTkDcuInfoInputFileName( dbParams_.inputDcuInfoXml_ ); 
      } catch (...) { 
	handleException( __func__ ); 
      }
      LogTrace(mlConfigDb_)
	<< "[SiStripConfigDb::" << __func__ << "]"
	<< " Added 'dcuinfo.xml' file: " << dbParams_.inputDcuInfoXml_;
    } else {
      edm::LogWarning(mlConfigDb_)
	<< "[SiStripConfigDb::" << __func__ << "]"
	<< " No 'dcuinfo.xml' file found at " << dbParams_.inputDcuInfoXml_;
      dbParams_.inputDcuInfoXml_ = ""; 
    } 
  }

  // Input FEC xml files
  if ( dbParams_.inputFecXml_.empty() ) {
    edm::LogWarning(mlConfigDb_) 
      << "[SiStripConfigDb::" << __func__ << "]"
      << " NULL paths to input 'fec.xml' files!";
  } else {
    vector<string>::iterator iter = dbParams_.inputFecXml_.begin();
    for ( ; iter != dbParams_.inputFecXml_.end(); iter++ ) {
      if ( *iter == "" ) {
	edm::LogWarning(mlConfigDb_)
	  << "[SiStripConfigDb::" << __func__ << "]"
	  << " NULL path to input 'fec.xml' file!";
      } else {
	if ( checkFileExists( *iter ) ) { 
	  try { 
	    if ( dbParams_.inputFecXml_.size() == 1 ) {
	      deviceFactory(__func__)->setFecInputFileName( *iter ); 
	    } else {
	      deviceFactory(__func__)->addFecFileName( *iter ); 
	    }
	  } catch (...) { handleException( __func__ ); }
	  LogTrace(mlConfigDb_) 
	    << "[SiStripConfigDb::" << __func__ << "]"
	    << " Added 'fec.xml' file: " << *iter;
	} else {
	  edm::LogWarning(mlConfigDb_) 
	    << "[SiStripConfigDb::" << __func__ << "]"
	    << " No 'fec.xml' file found at " << *iter;
	  *iter = ""; 
	} 
      }
    }
  }
    
  // Input FED xml files
  if ( dbParams_.inputFedXml_.empty() ) {
    edm::LogWarning(mlConfigDb_) 
      << "[SiStripConfigDb::" << __func__ << "]"
      << " NULL paths to input 'fed.xml' files!";
  } else {
    vector<string>::iterator iter = dbParams_.inputFedXml_.begin();
    for ( ; iter != dbParams_.inputFedXml_.end(); iter++ ) {
      if ( *iter == "" ) {
	edm::LogWarning(mlConfigDb_) 
	  << "[SiStripConfigDb::" << __func__ << "]"
	  << " NULL path to input 'fed.xml' file!";
      } else {
	if ( checkFileExists( *iter ) ) { 
	  try { 
	    if ( dbParams_.inputFedXml_.size() == 1 ) {
	      deviceFactory(__func__)->setFedInputFileName( *iter ); 
	    } else {
	      deviceFactory(__func__)->addFedFileName( *iter ); 
	    }
	  } catch (...) { 
	    handleException( __func__ ); 
	  }
	  LogTrace(mlConfigDb_) 
	    << "[SiStripConfigDb::" << __func__ << "]"
	    << " Added 'fed.xml' file: " << *iter;
	} else {
	  edm::LogWarning(mlConfigDb_) 
	    << "[SiStripConfigDb::" << __func__ << "]"
	    << " No 'fed.xml' file found at " << *iter;
	  *iter = ""; 
	} 
      }
    }
  }

  // Output module.xml file
  if ( dbParams_.outputModuleXml_ == "" ) { 
    edm::LogWarning(mlConfigDb_) 
      << "[SiStripConfigDb::" << __func__ << "]"
      << " NULL path to output 'module.xml' file!"
      << " Setting to '/tmp/module.xml'...";
    dbParams_.outputModuleXml_ = "/tmp/module.xml"; 
  } else {
    try { 
#ifdef USING_NEW_DATABASE_MODEL
      ConnectionFactory* factory = deviceFactory(__func__);
#else
      FedFecConnectionDeviceFactory* factory = deviceFactory(__func__);
#endif
      factory->setOutputFileName( dbParams_.outputModuleXml_ ); 
    } catch (...) { 
      handleException( __func__, "Problems setting output 'module.xml' file!" ); 
    }
  }

  // Output dcuinfo.xml file
  if ( dbParams_.outputDcuInfoXml_ == "" ) { 
    edm::LogWarning(mlConfigDb_) 
      << "[SiStripConfigDb::" << __func__ << "]"
      << " NULL path to output 'dcuinfo.xml' file!"
      << " Setting to '/tmp/dcuinfo.xml'...";
    dbParams_.outputModuleXml_ = "/tmp/dcuinfo.xml"; 
  } else {
    try { 
      TkDcuInfoFactory* factory = deviceFactory(__func__);
      factory->setOutputFileName( dbParams_.outputDcuInfoXml_ ); 
    } catch (...) { 
      handleException( __func__, "Problems setting output 'dcuinfo.xml' file!" ); 
    }
  }

  // Output fec.xml file
  if ( dbParams_.outputFecXml_ == "" ) {
    edm::LogWarning(mlConfigDb_) 
      << "[SiStripConfigDb::" << __func__ << "]"
      << " NULL path to output 'fec.xml' file!"
      << " Setting to '/tmp/fec.xml'...";
    dbParams_.outputFecXml_ = "/tmp/fec.xml";
  } else {
    try { 
      FecDeviceFactory* factory = deviceFactory(__func__);
      factory->setOutputFileName( dbParams_.outputFecXml_ ); 
    } catch (...) { 
      handleException( __func__, "Problems setting output 'fec.xml' file!" ); 
    }
  }

  // Output fed.xml file
  if ( dbParams_.outputFedXml_ == "" ) {
    edm::LogWarning(mlConfigDb_) 
      << "[SiStripConfigDb::" << __func__ << "]"
      << " NULL path to output 'fed.xml' file!"
      << " Setting to '/tmp/fed.xml'...";
    dbParams_.outputFedXml_ = "/tmp/fed.xml";
  } else {
    try { 
      Fed9U::Fed9UDeviceFactory* factory = deviceFactory(__func__);
      factory->setOutputFileName( dbParams_.outputFedXml_ ); 
    } catch (...) { 
      handleException( __func__, "Problems setting output 'fed.xml' file!" ); 
    }
  }

}

// -----------------------------------------------------------------------------
// 
void SiStripConfigDb::createPartition( const string& partition_name,
				       const SiStripFecCabling& fec_cabling,
				       const DcuDetIdMap& dcu_detid_map ) {
  
  // Set partition name and version
  dbParams_.partition_ = partition_name;
  dbParams_.fecMajor_ = 0;
  dbParams_.fecMinor_ = 0;

  LogTrace(mlConfigDb_)
    << "[SiStripConfigDb::" << __func__ << "]"
    << " Creating partition " << dbParams_.partition_;

  // Create new partition based on device descriptions
  createDeviceDescriptions( fec_cabling );
  if ( !devices_.empty() ) {
    try {
      stringstream ss; 
      ss << "/tmp/fec_" << dbParams_.partition_ << ".xml";
      FecDeviceFactory* factory = deviceFactory(__func__);
      factory->setOutputFileName( ss.str() );
#ifdef USING_NEW_DATABASE_MODEL
      deviceFactory(__func__)->createPartition( devices_,
						&dbParams_.fecMajor_, 
						&dbParams_.fecMinor_, 
						dbParams_.partition_ );
#else
      deviceFactory(__func__)->createPartition( devices_,
						&dbParams_.fecMajor_, 
						&dbParams_.fecMinor_, 
						dbParams_.partition_,
						dbParams_.partition_ ); 
#endif
    } catch (...) { 
      stringstream ss; 
      ss << "Failed to create new partition with name "
	 << dbParams_.partition_ << " and version " 
	 << dbParams_.fecMajor_ << "." << dbParams_.fecMinor_;
      handleException( __func__, ss.str() );
    } 
  }
  
  // Create and upload FED descriptions
  createFedDescriptions( fec_cabling );
  if ( !feds_.empty() ) {
    try {
      stringstream ss; 
      ss << "/tmp/fed_" << dbParams_.partition_ << ".xml";
      Fed9U::Fed9UDeviceFactory* factory = deviceFactory(__func__);
      factory->setOutputFileName( ss.str() );
      uploadFedDescriptions( true );
    } catch(...) {
      stringstream ss; 
      ss << "Failed to create and upload FED descriptions"
	 << " to partition with name "
	 << dbParams_.partition_ << " and version " 
	 << dbParams_.fedMajor_ << "." << dbParams_.fedMinor_;
      handleException( __func__, ss.str() );
    }
  }

  // Create and upload FED connections
  createFedConnections( fec_cabling );
  if ( !connections_.empty() ) {
    try {
      uploadFedConnections( true );
    } catch(...) {
      stringstream ss; 
      ss << "Failed to add FedChannelConnectionDescription!";
      handleException( __func__, ss.str() );
    }
    try {
      stringstream ss; 
      ss << "/tmp/module_" << dbParams_.partition_ << ".xml";
#ifdef USING_NEW_DATABASE_MODEL
      ConnectionFactory* factory = deviceFactory(__func__);
#else
      FedFecConnectionDeviceFactory* factory = deviceFactory(__func__);
#endif
      factory->setOutputFileName( ss.str() );
#ifndef USING_NEW_DATABASE_MODEL
      deviceFactory(__func__)->write(); //@@ corresponding method in new model???
#endif
    } catch(...) {
      stringstream ss; 
      ss << "Failed to create and upload ConnectionDescriptions"
	 << " to partition with name "
	 << dbParams_.partition_ << " and version " 
	 << dbParams_.cabMajor_ << "." << dbParams_.cabMinor_;
      handleException( __func__, ss.str() );
    }
  }

  // Create and upload DCU-DetId map
  dcuDetIdMap_.clear();
  dcuDetIdMap_ = dcu_detid_map;
  if ( !dcuDetIdMap_.empty() ) {
    uploadDcuDetIdMap(); 
  }
  
  LogTrace(mlConfigDb_)
    << "[SiStripConfigDb::" << __func__ << "]"
    << " Finished creating partition " << dbParams_.partition_;
  
}

// -----------------------------------------------------------------------------
// 
void SiStripConfigDb::handleException( const string& method_name,
				       const string& extra_info ) { //throw (cms::Exception) {

  stringstream ss;
  try {
    //throw; // rethrow caught exception to be dealt with below
  } 

  catch ( const cms::Exception& e ) { 
    ss << " Caught cms::Exception in method "
       << method_name << " with message: " << endl 
       << e.what();
    if ( extra_info != "" ) { ss << "Additional info: " << extra_info << endl; }
    //throw e; // rethrow cms::Exception
  }
  
  catch ( const oracle::occi::SQLException& e ) { 
    ss << " Caught oracle::occi::SQLException in method "
       << method_name << " with message: " << endl 
       << e.getMessage();
    if ( extra_info != "" ) { ss << "Additional info: " << extra_info << endl; }
    //throw cms::Exception(mlConfigDb_) << ss.str() << endl;
  }

  catch ( const FecExceptionHandler& e ) {
    ss << " Caught FecExceptionHandler exception in method "
       << method_name << " with message: " << endl 
       << const_cast<FecExceptionHandler&>(e).getMessage();
    if ( extra_info != "" ) { ss << "Additional info: " << extra_info << endl; }
    //throw cms::Exception(mlConfigDb_) << ss.str() << endl;
  }

  catch ( const ICUtils::ICException& e ) {
    ss << " Caught ICUtils::ICException in method "
       << method_name << " with message: " << endl 
       << e.what();
    if ( extra_info != "" ) { ss << "Additional info: " << extra_info << endl; }
    //throw cms::Exception(mlConfigDb_) << ss.str() << endl;
  }

  catch ( const exception& e ) {
    ss << " Caught std::exception in method "
       << method_name << " with message: " << endl 
       << e.what();
    if ( extra_info != "" ) { ss << "Additional info: " << extra_info << endl; }
    //throw cms::Exception(mlConfigDb_) << ss.str() << endl;
  }

  catch (...) {
    ss << " Caught unknown exception in method "
       << method_name << " (No message) " << endl;
    if ( extra_info != "" ) { ss << "Additional info: " << extra_info << endl; }
    //throw cms::Exception(mlConfigDb_) << ss.str() << endl;
  }
  
  // Message
  edm::LogError(mlConfigDb_) << ss.str();
  
}

// -----------------------------------------------------------------------------
//
bool SiStripConfigDb::checkFileExists( const std::string& path ) {
  fstream fs; 
  fs.open( path.c_str(), ios::in );
  if( !fs.is_open() ) { return false; }
  fs.close();
  return true;
}






  
