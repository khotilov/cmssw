/*
 * \file DQMStoreStats.cc
 * \author Andreas Meyer
 * Last Update:
 * $Date: 2009/12/15 08:59:50 $
 * $Revision: 1.9 $
 * $Author: dellaric $
 *
 * Description: Print out statistics of histograms in DQMStore
*/

#include "DQMServices/Components/src/DQMStoreStats.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

using namespace std;
using namespace edm;

//==================================================================//
//================= Constructor and Destructor =====================//
//==================================================================//
DQMStoreStats::DQMStoreStats( const edm::ParameterSet& ps )
  : subsystem_ (""),
    subfolder_ (""),
    nbinsglobal_ (0),
    nbinssubsys_ (0),
    nmeglobal_ (0),
    nmesubsys_ (0),
    maxbinsglobal_ (0),
    maxbinssubsys_ (0),
    maxbinsmeglobal_ (""),
    maxbinsmesubsys_ (""),
    statsdepth_ (1),
    pathnamematch_ ("*"),
    verbose_ (0)
{
  parameters_ = ps;
  pathnamematch_ = ps.getUntrackedParameter<std::string>( "pathNameMatch", pathnamematch_ );
  statsdepth_ = ps.getUntrackedParameter<int>( "statsDepth", statsdepth_ );
  verbose_ = ps.getUntrackedParameter<int>( "verbose", verbose_ );
  dumpMemHistory_ = ps.getUntrackedParameter<bool>( "dumpMemoryHistory", false );
  runonendrun_    = ps.getUntrackedParameter<bool>( "runOnEndRun", true );
  runonendjob_    = ps.getUntrackedParameter<bool>( "runOnEndJob", false );
  runonendlumi_   = ps.getUntrackedParameter<bool>( "runOnEndLumi", false );
  runineventloop_ = ps.getUntrackedParameter<bool>( "runInEventLoop", false );

  startingTime_ = time( 0 );

}

DQMStoreStats::~DQMStoreStats(){
}


///
/// do the stats here and produce output;
///
/// mode is coded in DQMStoreStats::statMode enum
/// (select subsets of ME, e.g. those with getLumiFlag() == true)
///
int DQMStoreStats::calcstats( int mode = DQMStoreStats::considerAllME ) {

  ////---- initialise Event and LS counters
  nbinsglobal_ = 0; 
  nbinssubsys_ = 0; 
  maxbinsglobal_ = 0; 
  maxbinssubsys_ = 0; 
  std::string path = "";
  std::string subsystemname = "";
  std::string subfoldername = "";
  size_t subsysStringEnd = 0, subfolderStringEnd  = 0;


  std::vector<MonitorElement*> melist;
  melist = dbe_->getMatchingContents( pathnamematch_ );

  DQMStoreStatsTopLevel dqmStoreStatsTopLevel;

  // loop all ME
  typedef std::vector <MonitorElement*>::iterator meIt;
  for(meIt it = melist.begin(); it != melist.end(); ++it) {

    // consider only ME with getLumiFlag() == true ?
    if( mode == DQMStoreStats::considerOnlyLumiProductME && 
	!( (*it)->getLumiFlag() ) ) continue;
    
    // figure out subsystem/subfolder names
    std::string path = (*it)->getPathname();

    // protection against ghost ME with empty paths
    if( 0 == path.size() ) continue;

    subsysStringEnd = path.find( '/', 0 );
    if( std::string::npos == subsysStringEnd ) subsysStringEnd = path.size(); // no subfolder

    // new subsystem?
    if( path.substr( 0, subsysStringEnd ) != subsystemname ) {
      DQMStoreStatsSubsystem aSubsystem;
      subsystemname = path.substr( 0, subsysStringEnd );
      aSubsystem.subsystemName_ = subsystemname;
      dqmStoreStatsTopLevel.push_back( aSubsystem );
    }

    // get subfolder name (if there is one..)
    if( path.size() == subsysStringEnd ) {
      // no subfolders in subsystem, make dummy
      DQMStoreStatsSubfolder aSubfolder;
      aSubfolder.subfolderName_ = subsystemname; // <-- for tagging this case
      dqmStoreStatsTopLevel.back().push_back( aSubfolder );
    }

    else {

      // there is a subfolder, get its name
      subfolderStringEnd = path.find( '/', subsysStringEnd + 1 );
      if( std::string::npos == subfolderStringEnd ) subfolderStringEnd = path.size();

      // new subfolder?
      if( path.substr( subsysStringEnd + 1, subfolderStringEnd - subsysStringEnd - 1 ) != subfoldername ) {
	subfoldername = path.substr( subsysStringEnd + 1, subfolderStringEnd - subsysStringEnd - 1 );
	DQMStoreStatsSubfolder aSubfolder;
	aSubfolder.subfolderName_ = subfoldername;
	dqmStoreStatsTopLevel.back().push_back( aSubfolder );
      }

    }

    // shortcut
    DQMStoreStatsSubfolder& currentSubfolder = dqmStoreStatsTopLevel.back().back();

    switch( (*it)->kind() ) {
      
      // one-dim ME
      case MonitorElement::DQM_KIND_TH1F: currentSubfolder.AddBinsF( (*it)->getNbinsX() ); break;
      case MonitorElement::DQM_KIND_TH1S: currentSubfolder.AddBinsS( (*it)->getNbinsX() ); break;
      case MonitorElement::DQM_KIND_TH1D: currentSubfolder.AddBinsD( (*it)->getNbinsX() ); break;
      case MonitorElement::DQM_KIND_TPROFILE: currentSubfolder.AddBinsD( (*it)->getNbinsX() ); break;

      // two-dim ME
      case MonitorElement::DQM_KIND_TH2F: currentSubfolder.AddBinsF( (*it)->getNbinsX() * (*it)->getNbinsY() ); break;
      case MonitorElement::DQM_KIND_TH2S: currentSubfolder.AddBinsS( (*it)->getNbinsX() * (*it)->getNbinsY() ); break;
      case MonitorElement::DQM_KIND_TH2D: currentSubfolder.AddBinsD( (*it)->getNbinsX() * (*it)->getNbinsY() ); break;
      case MonitorElement::DQM_KIND_TPROFILE2D: currentSubfolder.AddBinsD( (*it)->getNbinsX() * (*it)->getNbinsY() ); break;
 
      // three-dim ME
      case MonitorElement::DQM_KIND_TH3F: 
	currentSubfolder.AddBinsF( (*it)->getNbinsX() * (*it)->getNbinsY() * (*it)->getNbinsZ() ); break;

      default: {}
	// here we have a DQM_KIND_INVALID, DQM_KIND_INT, DQM_KIND_REAL or DQM_KIND_STRING
	// which we don't care much about. Alternatively:

	//   std::cerr << "[DQMStoreStats::calcstats] ** WARNING: monitor element of kind: " 
	// 	       << (*it)->kind() << ", name: \"" << (*it)->getName() << "\"\n"
	// 	       << "  in path: \"" << path << "\" not considered." << std::endl;
    }
      
  } 



  // OUTPUT

  std::cout << endl;
  std::cout << "===========================================================================================" << std::endl;
  std::cout << "[DQMStoreStats::calcstats] -- Dumping stats results ";
  if( mode == DQMStoreStats::considerAllME ) std::cout << "FOR ALL ME" << std::endl;
  else if( mode == DQMStoreStats::considerOnlyLumiProductME ) std::cout << "FOR LUMI PRODUCTS ONLY" << std::endl;
  std::cout << "===========================================================================================" << std::endl;
  std::cout << endl;

  std::cout << "------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "Configuration:" << std::endl;
  std::cout << "------------------------------------------------------------------------------------------" << std::endl;
  std::cout << " > running ";
  if (runonendrun_) std::cout << "on run end." << std::endl;
  if (runonendlumi_) std::cout << "on lumi end." << std::endl;
  if (runonendjob_) std::cout << "on job end." << std::endl;
  if (runineventloop_) std::cout << "in event loop." << std::endl;
  std::cout << " > pathNameMatch = \"" << pathnamematch_ << "\"" << std::endl;
  std::cout << std::endl;

  // dump folder structure
  std::cout << "------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "Top level folder tree:" << std::endl;
  std::cout << "------------------------------------------------------------------------------------------" << std::endl;
  for( DQMStoreStatsTopLevel::const_iterator it0 = dqmStoreStatsTopLevel.begin(); it0 < dqmStoreStatsTopLevel.end(); ++it0 ) {
    std::cout << it0->subsystemName_ << " (subsystem)" << std::endl;
    
    for( DQMStoreStatsSubsystem::const_iterator it1 = it0->begin(); it1 < it0->end(); ++it1 ) {
      std::cout << "  |--> " << it1->subfolderName_ << " (subfolder)" << std::endl;
    }
    
  }

  // dump mem/bin table

  unsigned int overallNHistograms = 0, overallNBins = 0, overallNBytes = 0;

  std::cout << std::endl;
  std::cout << "------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "Detailed ressource usage information ";
  if( mode == DQMStoreStats::considerAllME ) std::cout << "FOR ALL ME" << std::endl;
  else if( mode == DQMStoreStats::considerOnlyLumiProductME ) std::cout << "FOR LUMI PRODUCTS ONLY" << std::endl;
  std::cout << "------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "subsystem/folder                  histograms     bins      bins per      MB        kB per" << std::endl;
  std::cout << "                                   (total)     (total)    histogram    (total)   histogram  " << std::endl;
  std::cout << "------------------------------------------------------------------------------------------" << std::endl;
  for( DQMStoreStatsTopLevel::const_iterator it0 = dqmStoreStatsTopLevel.begin(); it0 < dqmStoreStatsTopLevel.end(); ++it0 ) {
    std::cout << it0->subsystemName_ << std::endl;
    
    unsigned int nHistograms = 0, nBins = 0, nBytes = 0;

    for( DQMStoreStatsSubsystem::const_iterator it1 = it0->begin(); it1 < it0->end(); ++it1 ) {

      // fixed-size working copy
      std::string thisSubfolderName( it1->subfolderName_ );
      if( thisSubfolderName.size() > 30 ) {
	thisSubfolderName.resize( 30 );
	thisSubfolderName.replace( thisSubfolderName.size() - 3, 3, 3, '.' );
      }

      std::cout << " -> " << std::setw( 30 ) << std::left << thisSubfolderName;
      std::cout << std::setw( 7 ) << std::right << it1->totalHistos_;
      std::cout << std::setw( 12 ) << std::right << it1->totalBins_;

      // bins/histogram, need to catch nan if histos=0
      if( it1->totalHistos_ ) {
	std::cout << std::setw( 12 ) << std::right << std::setprecision( 3 ) << it1->totalBins_ / float( it1->totalHistos_ );
      } 
      else std::cout << std::setw( 12 ) << std::right << "-";

      std::cout << std::setw( 12 ) << std::right << std::setprecision( 3 ) << it1->totalMemory_ / 1024. / 1000.;

      // mem/histogram, need to catch nan if histos=0
      if( it1->totalHistos_ ) {
	std::cout << std::setw( 12 ) << std::right << std::setprecision( 3 ) << it1->totalMemory_ / 1024. / it1->totalHistos_;
      }
      else std::cout << std::setw( 12 ) << std::right << "-";

      std::cout << std::endl;

      // collect totals
      nHistograms += it1->totalHistos_; 
      nBins       += it1->totalBins_;   
      nBytes      += it1->totalMemory_; 

    }




    overallNHistograms += nHistograms;
    overallNBins       += nBins;
    overallNBytes      += nBytes;

    // display totals
    std::cout << "    " << std::setw( 30 ) << std::left << "SUBSYSTEM TOTAL";
    std::cout << std::setw( 7 ) << std::right << nHistograms;
    std::cout << std::setw( 12 ) << std::right << nBins;
    std::cout << std::setw( 12 ) << std::right << std::setprecision( 3 ) << nBins / float( nHistograms );
    std::cout << std::setw( 12 ) << std::right << std::setprecision( 3 ) << nBytes / 1024. / 1000.;
    std::cout << std::setw( 12 ) << std::right << std::setprecision( 3 ) << nBytes / 1024. / nHistograms;
    std::cout << std::endl;
      
    std::cout << ".........................................................................................." << std::endl;

  }


  // dump total
  std::cout << std::endl;
  std::cout << "------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "Grand total ";
  if( mode == DQMStoreStats::considerAllME ) std::cout << "FOR ALL ME:" << std::endl;
  else if( mode == DQMStoreStats::considerOnlyLumiProductME ) std::cout << "FOR LUMI PRODUCTS ONLY:" << std::endl;
  std::cout << "------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "Number of subsystems: " << dqmStoreStatsTopLevel.size() << std::endl;
  std::cout << "Total number of histograms: " << overallNHistograms << " with: " << overallNBins << " bins alltogether" << std::endl;
  std::cout << "Total memory occupied by histograms (excl. overhead): " << overallNBytes / 1024. / 1000. << " MB" << std::endl;



  std::cout << endl;
  std::cout << "===========================================================================================" << std::endl;
  std::cout << "[DQMStoreStats::calcstats] -- End of output ";
  if( mode == DQMStoreStats::considerAllME ) std::cout << "FOR ALL ME." << std::endl;
  else if( mode == DQMStoreStats::considerOnlyLumiProductME ) std::cout << "FOR LUMI PRODUCTS ONLY." << std::endl;
  std::cout << "===========================================================================================" << std::endl;
  std::cout << endl;

  return 0;

}



///
///
///
void DQMStoreStats::dumpMemoryProfile( void ) {

  std::cout << std::endl;
  std::cout << "------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "Memory profile:" << std::endl;
  std::cout << "------------------------------------------------------------------------------------------" << std::endl;

  // determine virtual memory maximum
  std::pair<time_t, unsigned int> maxItem( 0, 0 );
  for( std::vector<std::pair<time_t, unsigned int> >::const_iterator it = memoryHistoryVector_.begin();
       it < memoryHistoryVector_.end(); ++it ) {
    if( it->second > maxItem.second ) {
      maxItem = *it;
    }
  }

  std::stringstream rootOutputFileName;
  rootOutputFileName << "dqmStoreStats_memProfile_" << getpid() << ".root";

  // dump memory history to root file
  if( dumpMemHistory_ && isOpenProcFileSuccessful_ ) {

    TFile outputFile( rootOutputFileName.str().c_str(), "RECREATE" );

    int aTime;
    float aMb;

    TTree memHistoryTree( "dqmstorestats_memhistory", "memory history" );
    memHistoryTree.Branch( "seconds", &aTime, "seconds/I" );
    memHistoryTree.Branch( "megabytes", &aMb, "megabytes/F" );
    for( std::vector<std::pair<time_t, unsigned int> >::const_iterator it = memoryHistoryVector_.begin();
	 it < memoryHistoryVector_.end(); ++it ) {
      aTime = it->first - startingTime_;
      aMb = it->second / 1000.;
      memHistoryTree.Fill();
    }

    outputFile.Write();
    outputFile.Close();

  }

  std::cout << "Approx. maximum total virtual memory size of job: ";
  if( isOpenProcFileSuccessful_ && memoryHistoryVector_.size() ) {
    std::cout << maxItem.second / 1000.
	      << " MB (reached " << maxItem.first - startingTime_ << " sec. after constructor called)," << std::endl;
    std::cout << " memory history written to: " << rootOutputFileName.str() << " (" << memoryHistoryVector_.size() << " samples)" << std::endl;
  } else {
    std::cout << "(could not be determined)" << std::endl;
  }

  std::cout << std::endl << std::endl;
  
}



///
///
///
void DQMStoreStats::print(){
  // subsystem info printout
  std::cout << " ---------- " << subsystem_ << " ---------- " << std::endl;
  std::cout <<  "  " << subfolder_ << ": " ;
  std::cout <<  nmesubsys_ << " histograms with " 
            <<  nbinssubsys_  << " bins. " ; 
  if (nmesubsys_ > 0) std::cout <<  nbinssubsys_/nmesubsys_ << " bins/histogram " ;
  std::cout << std::endl;
  std::cout <<  "  Largest histogram: " << maxbinsmesubsys_ << " with " <<
		                         maxbinssubsys_ << " bins." <<  std::endl;
}




///
/// read virtual memory size from /proc/<pid>/status file
///
std::pair<unsigned int, unsigned int> DQMStoreStats::readMemoryEntry( void ) const {
  
  // see if initial test reading was successful
  if( isOpenProcFileSuccessful_ ) {

    std::ifstream procFile( procFileName_.str().c_str(), ios::in );

    std::string readBuffer( "" );
    unsigned int memSize = 0;

    // scan procfile
    while( !procFile.eof() ) {
      procFile >> readBuffer;
      if( std::string( "VmSize:" ) == readBuffer ) {
	procFile >> memSize;
	break;
      }
    }

    procFile.close();
    return std::pair<time_t, unsigned int>( time( 0 ), memSize );
  }

  return std::pair<time_t, unsigned int>( 0, 0 );

}



//==================================================================//
//========================= beginJob ===============================//
//==================================================================//
void DQMStoreStats::beginJob() {

  ////---- get DQM store interface
  dbe_ = Service<DQMStore>().operator->();

  // access the proc/ folder for memory information
  procFileName_ << "/proc/" << getpid() << "/status";

  // open for a test
  std::ifstream procFile( procFileName_.str().c_str(), ios::in );

  if( procFile.good() ) {
    isOpenProcFileSuccessful_ = true;
  }
  else {
    std::cerr << " [DQMStoreStats::beginJob] ** WARNING: could not open file: " << procFileName_.str() << std::endl;
    std::cerr << "  Total memory profile will not be available." << std::endl;
    isOpenProcFileSuccessful_ = false;
  }

  procFile.close();

}

//==================================================================//
//========================= beginRun ===============================//
//==================================================================//
void DQMStoreStats::beginRun(const edm::Run& r, const EventSetup& context) {
}


//==================================================================//
//==================== beginLuminosityBlock ========================//
//==================================================================//
void DQMStoreStats::beginLuminosityBlock(const LuminosityBlock& lumiSeg,
					    const EventSetup& context) {
}


//==================================================================//
//==================== analyse (takes each event) ==================//
//==================================================================//
void DQMStoreStats::analyze(const Event& iEvent, const EventSetup& iSetup) {

  //now read virtual memory size from proc folder
  memoryHistoryVector_.push_back( readMemoryEntry() );

  if (runineventloop_) {
    calcstats( DQMStoreStats::considerAllME );
    calcstats( DQMStoreStats::considerOnlyLumiProductME );
    dumpMemoryProfile();
  }

}


//==================================================================//
//========================= endLuminosityBlock =====================//
//==================================================================//
void DQMStoreStats::endLuminosityBlock(const LuminosityBlock& lumiSeg,
					  const EventSetup& context) {
  if (runonendlumi_) { 
    calcstats( DQMStoreStats::considerAllME );
    calcstats( DQMStoreStats::considerOnlyLumiProductME );
    dumpMemoryProfile();
  }

}

//==================================================================//
//============================= endRun =============================//
//==================================================================//
void DQMStoreStats::endRun(const Run& r, const EventSetup& context) {

  if (runonendrun_) {
    calcstats( DQMStoreStats::considerAllME );
    calcstats( DQMStoreStats::considerOnlyLumiProductME );
    dumpMemoryProfile();
  }

}

//==================================================================//
//============================= endJob =============================//
//==================================================================//
void DQMStoreStats::endJob() {

  if (runonendjob_) {
    calcstats( DQMStoreStats::considerAllME );
    calcstats( DQMStoreStats::considerOnlyLumiProductME );
    dumpMemoryProfile();
  }

}
