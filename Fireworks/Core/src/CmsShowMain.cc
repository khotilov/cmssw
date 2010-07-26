// -*- C++ -*-
//
// Package:     Core
// Class  :     CmsShowMain
//
// Implementation:
//     <Notes on implementation>
//
// Original Author:
//         Created:  Mon Dec  3 08:38:38 PST 2007
// $Id: CmsShowMain.cc,v 1.172 2010/07/23 08:35:03 eulisse Exp $
//

// system include files
#include <sstream>
#include <boost/bind.hpp>
#include <boost/program_options.hpp>
#include <string.h>

#include "TSystem.h"
#include "TGLWidget.h"
#include "TTimer.h"
#include "TROOT.h"
#include "TGFileDialog.h"
#include "TGMsgBox.h"
#include "TMonitor.h"
#include "TServerSocket.h"
#include "TEveLine.h"
#include "TEveManager.h"
#include "TFile.h"

#include "Fireworks/Core/src/CmsShowMain.h"

#include "Fireworks/Core/interface/FWEveViewManager.h"

#include "Fireworks/Core/interface/FWTableViewManager.h"
#include "Fireworks/Core/interface/FWL1TriggerTableViewManager.h"
#include "Fireworks/Core/interface/FWTriggerTableViewManager.h"
#include "Fireworks/Core/interface/FWEventItemsManager.h"
#include "Fireworks/Core/interface/FWViewManagerManager.h"
#include "Fireworks/Core/interface/FWGUIManager.h"
#include "Fireworks/Core/interface/FWLiteJobMetadataManager.h"
#include "Fireworks/Core/interface/FWModelChangeManager.h"
#include "Fireworks/Core/interface/FWColorManager.h"
#include "Fireworks/Core/src/FWColorSelect.h"
#include "Fireworks/Core/interface/FWSelectionManager.h"
#include "Fireworks/Core/interface/FWConfigurationManager.h"
#include "Fireworks/Core/interface/FWMagField.h"
#include "Fireworks/Core/interface/Context.h"

#include "Fireworks/Core/interface/CmsShowNavigator.h"
#include "Fireworks/Core/interface/CSGAction.h"
#include "Fireworks/Core/interface/CSGContinuousAction.h"
#include "Fireworks/Core/interface/FWLiteJobMetadataUpdateRequest.h"

#include "Fireworks/Core/interface/ActionsList.h"

#include "Fireworks/Core/src/CmsShowTaskExecutor.h"
#include "Fireworks/Core/interface/CmsShowMainFrame.h"
#include "Fireworks/Core/interface/CmsShowSearchFiles.h"

#include "Fireworks/Core/interface/fwLog.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

//
// constants, enums and typedefs
//

static const char* const kInputFilesOpt        = "input-files";
static const char* const kInputFilesCommandOpt = "input-files,i";
static const char* const kConfigFileOpt        = "config-file";
static const char* const kConfigFileCommandOpt = "config-file,c";
static const char* const kGeomFileOpt          = "geom-file";
static const char* const kGeomFileCommandOpt   = "geom-file,g";
static const char* const kNoConfigFileOpt      = "noconfig";
static const char* const kNoConfigFileCommandOpt = "noconfig,n";
static const char* const kPlayOpt              = "play";
static const char* const kPlayCommandOpt       = "play,p";
static const char* const kLoopOpt              = "loop";
static const char* const kLoopCommandOpt       = "loop";
static const char* const kDebugOpt             = "debug";
static const char* const kDebugCommandOpt      = "debug,d";
static const char* const kLogLevelCommandOpt   = "log";
static const char* const kLogLevelOpt          = "log";
static const char* const kEveOpt               = "eve";
static const char* const kEveCommandOpt        = "eve";
static const char* const kAdvancedRenderOpt        = "shine";
static const char* const kAdvancedRenderCommandOpt = "shine,s";
static const char* const kHelpOpt        = "help";
static const char* const kHelpCommandOpt = "help,h";
static const char* const kSoftCommandOpt = "soft";
static const char* const kPortCommandOpt = "port";
static const char* const kPlainRootCommandOpt = "prompt";
static const char* const kRootInteractiveCommandOpt = "root-interactive,r";
static const char* const kChainCommandOpt = "chain";
static const char* const kLiveCommandOpt  = "live";
static const char* const kFieldCommandOpt = "field";
static const char* const kFreePaletteCommandOpt = "free-palette";
static const char* const kAutoSaveAllViews = "auto-save-all-views";


//
// constructors and destructor
//
CmsShowMain::CmsShowMain(int argc, char *argv[]) 
   : CmsShowMainBase(),
     m_navigator(new CmsShowNavigator(*this)),
     m_metadataManager(new FWLiteJobMetadataManager()),
     m_context(new fireworks::Context(changeManager(),
                                      selectionManager(),
                                      eiManager(),
                                      colorManager(),
                                      m_metadataManager.get()))
{
   eiManager()->setContext(m_context.get());

   try {
      std::string descString(argv[0]);
      descString += " [options] <data file>\nAllowed options";

      namespace po = boost::program_options;
      po::options_description desc(descString);
      desc.add_options()
         (kInputFilesCommandOpt, po::value< std::vector<std::string> >(),   "Input root files")
         (kConfigFileCommandOpt, po::value<std::string>(),   "Include configuration file")
         (kGeomFileCommandOpt,   po::value<std::string>(),   "Include geometry file")
         (kNoConfigFileCommandOpt,                           "Don't load any configuration file")
         (kPlayCommandOpt, po::value<float>(),               "Start in play mode with given interval between events in seconds")
         (kPortCommandOpt, po::value<unsigned int>(),        "Listen to port for new data files to open")
         (kEveCommandOpt,                                    "Show Eve browser to help debug problems")
         (kLoopCommandOpt,                                   "Loop events in play mode")
         (kPlainRootCommandOpt,                              "Plain ROOT without event display")
         (kRootInteractiveCommandOpt,                        "Enable root interactive prompt")
         (kDebugCommandOpt,                                  "Start the display from a debugger and producer a crash report")
         (kLogLevelCommandOpt, po::value<unsigned int>(),    "Set log level starting from 0 to 4 : kDebug, kInfo, kWarning, kError")
         (kAdvancedRenderCommandOpt,                         "Use advance options to improve rendering quality       (anti-alias etc)")
         (kSoftCommandOpt,                                   "Try to force software rendering to avoid problems with bad hardware drivers")
         (kChainCommandOpt, po::value<unsigned int>(),       "Chain up to a given number of recently open files. Default is 1 - no chain")
         (kLiveCommandOpt,                                   "Enforce playback mode if a user is not using display")
         (kFieldCommandOpt, po::value<double>(),             "Set magnetic field value explicitly. Default is auto-field estimation")
         (kFreePaletteCommandOpt,                            "Allow free color selection (requires special configuration!)")
         (kAutoSaveAllViews, po::value<std::string>(),       "Auto-save all views with given prefix (run_event_lumi_view.png is appended)")
         (kHelpCommandOpt,                                   "Display help message");
      po::positional_options_description p;
      p.add(kInputFilesOpt, -1);

      int newArgc = argc;
      char **newArgv = argv;
      po::variables_map vm;
      //po::store(po::parse_command_line(newArgc, newArgv, desc), vm);
      //po::notify(vm);
      po::store(po::command_line_parser(newArgc, newArgv).
                options(desc).positional(p).run(), vm);
      po::notify(vm);
      if(vm.count(kHelpOpt)) {
         std::cout << desc <<std::endl;
         exit(0);
      }
      
      if(vm.count(kLogLevelOpt)) {
         fwlog::LogLevel level = (fwlog::LogLevel)(vm[kLogLevelOpt].as<unsigned int>());
         fwlog::setPresentLogLevel(level);
      }

      if(vm.count(kPlainRootCommandOpt)) {
         fwLog(fwlog::kInfo) << "Plain ROOT prompt requested" << std::endl;
         return;
      }

      const char* cmspath = gSystem->Getenv("CMSSW_BASE");
      if(0 == cmspath) {
         throw std::runtime_error("CMSSW_BASE environment variable not set");
      }

      // input file
      if (vm.count(kInputFilesOpt)) {
         m_inputFiles = vm[kInputFilesOpt].as< std::vector<std::string> >();
      }

      if (!m_inputFiles.size())
         fwLog(fwlog::kInfo) << "No data file given." << std::endl;
      else if (m_inputFiles.size() == 1)
         fwLog(fwlog::kInfo) << "Input " << m_inputFiles.front() << std::endl;
      else
         fwLog(fwlog::kInfo) << m_inputFiles.size() << " input files; first: " << m_inputFiles.front() << ", last: " << m_inputFiles.back() << std::endl;

      // configuration file
      if (vm.count(kConfigFileOpt)) {
         setConfigFilename(vm[kConfigFileOpt].as<std::string>());
         if (access(configFilename(), R_OK) == -1)
         {
            fwLog(fwlog::kError) << "Specified configuration file does not exist. Quitting.\n";
            exit(1);
         }
      } else {
         if (vm.count(kNoConfigFileOpt)) {
            fwLog(fwlog::kInfo) << "No configuration is loaded, show everything.\n";
            setConfigFilename("");
         } else
            setConfigFilename("default.fwc");
      }
      fwLog(fwlog::kInfo) << "Config "  <<  configFilename() << std::endl;

      // geometry
      if (vm.count(kGeomFileOpt)) {
         m_geomFileName = vm[kGeomFileOpt].as<std::string>();
      } else {
         fwLog(fwlog::kInfo) << "No geom file name.  Choosing default.\n";
         m_geomFileName.append("cmsGeom10.root");
      }
      fwLog(fwlog::kInfo) << "Geom " <<  m_geomFileName.c_str() << std::endl;

      // Free-palette palette
      if (vm.count(kFreePaletteCommandOpt)) {
         FWColorPopup::EnableFreePalette();
         fwLog(fwlog::kInfo) << "Palette restriction removed on user request!\n";
      }
      bool eveMode = vm.count(kEveOpt);

      //Delay creating guiManager and enabling autoloading until here so that if we have a 'help' request we don't
      // open any graphics or build dictionaries
      AutoLibraryLoader::enable();

      TEveManager::Create(kFALSE, "FIV");

      setup(m_navigator.get(), m_context.get(), m_metadataManager.get());

      if ( vm.count(kAdvancedRenderOpt) ) 
      {
         TEveLine::SetDefaultSmooth(kTRUE);
      }

      //figure out where to find macros
      //tell ROOT where to find our macros
      CmsShowTaskExecutor::TaskFunctor f;
      // first check if port is not occupied
      if (vm.count(kPortCommandOpt)) { 	 
         f=boost::bind(&CmsShowMain::setupSocket, this, vm[kPortCommandOpt].as<unsigned int>()); 	 
         startupTasks()->addTask(f); 	 
      }
    
      f=boost::bind(&CmsShowMain::loadGeometry,this);
      startupTasks()->addTask(f);
      f=boost::bind(&CmsShowMain::setupViewManagers,this);
      startupTasks()->addTask(f);
      f=boost::bind(&CmsShowMainBase::setupConfiguration,this);
      startupTasks()->addTask(f);
      f=boost::bind(&CmsShowMain::setupDataHandling,this);
      startupTasks()->addTask(f);

      if (vm.count(kLoopOpt))
         setPlayLoop();

      gSystem->IgnoreSignal(kSigSegmentationViolation, true);
      if (eveMode) {
         f = boost::bind(&CmsShowMainBase::setupDebugSupport,this);
         startupTasks()->addTask(f);
      }
      if(vm.count(kChainCommandOpt)) {
         f = boost::bind(&CmsShowNavigator::setMaxNumberOfFilesToChain, m_navigator.get(), vm[kChainCommandOpt].as<unsigned int>());
         startupTasks()->addTask(f);
      }
      if (vm.count(kPlayOpt)) {
         f = boost::bind(&CmsShowMainBase::setupAutoLoad, this, vm[kPlayOpt].as<float>());
         startupTasks()->addTask(f);
      }

      if(vm.count(kLiveCommandOpt))
      {
         f = boost::bind(&CmsShowMainBase::setLiveMode, this);
         startupTasks()->addTask(f);
      }
      
      if(vm.count(kFieldCommandOpt)) 
      {
         m_context->getField()->setSource(FWMagField::kUser);
         m_context->getField()->setUserField(vm[kFieldCommandOpt].as<double>());
      }
      if(vm.count(kAutoSaveAllViews)) {
         m_autoSaveAllViewsFormat  = vm[kAutoSaveAllViews].as<std::string>();
         m_autoSaveAllViewsFormat += "%d_%d_%d_%s.png";
      }
      startupTasks()->startDoingTasks();
   } catch(std::exception& iException) {
      std::cerr <<"CmsShowMain caught exception "<<iException.what()<<std::endl;
      throw;
   }
}

//
// Destruction
//

CmsShowMain::~CmsShowMain()
{}

class DieTimer : public TTimer
{
protected:
   CmsShowMain* fApp;
public:
   DieTimer(CmsShowMain* app) : TTimer(), fApp(app)
   {
      Start(0, kTRUE);
   }

   virtual Bool_t Notify()
   {
      TurnOff();
      fApp->doExit();
      delete this;
      return kFALSE;
   }
};

void CmsShowMain::quit()
{
   new DieTimer(this);
}

void CmsShowMain::doExit()
{
   // pre terminate eve
   m_context->deleteEveElements();
   guiManager()->evePreTerminate();

   // sleep at least 150 ms
   // windows in ROOT GUI are destroyed in 150 ms timeout after
   gSystem->Sleep(151);
   gSystem->ProcessEvents();

   gSystem->ExitLoop();
}

//
// assignment operators
//
// const CmsShowMain& CmsShowMain::operator=(const CmsShowMain& rhs)
// {
//   //An exception safe implementation is
//   CmsShowMain temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//

const fwlite::Event* 
CmsShowMain::getCurrentEvent() const
{
   if (m_navigator.get())
     return static_cast<const fwlite::Event*>(m_navigator->getCurrentEvent());
   return 0;
}

void
CmsShowMain::fileChangedSlot(const TFile *file)
{
   m_openFile = file;
   if (file)
      guiManager()->titleChanged(m_openFile->GetName());
   m_metadataManager->update(new FWLiteJobMetadataUpdateRequest(getCurrentEvent(), m_openFile));
}

void
CmsShowMain::eventChangedSlot()
{
   m_metadataManager->update(new FWLiteJobMetadataUpdateRequest(getCurrentEvent(), m_openFile));
}

void CmsShowMain::resetInitialization() {
   //printf("Need to reset\n");
}

void CmsShowMain::openData()
{
   const char* kRootType[] = {"ROOT files","*.root", 0, 0};
   TGFileInfo fi;
   fi.fFileTypes = kRootType;
   /* this is how things used to be done:
      fi.fIniDir = ".";
      this is bad because the destructor calls delete[] on fIniDir.
    */
   fi.fIniDir = new char[10];
   strcpy(fi.fIniDir, ".");
   new TGFileDialog(gClient->GetDefaultRoot(), guiManager()->getMainFrame(), kFDOpen, &fi);
   guiManager()->updateStatus("loading file ...");
   if (fi.fFilename) {
      m_navigator->openFile(fi.fFilename);
      m_loadedAnyInputFile = true;
      m_navigator->firstEvent();
      checkPosition();
      draw();
   }
   guiManager()->clearStatus();
}

void CmsShowMain::appendData()
{
   const char* kRootType[] = {"ROOT files","*.root", 0, 0};
   TGFileInfo fi;
   fi.fFileTypes = kRootType;
   /* this is how things used to be done:
      fi.fIniDir = ".";
      this is bad because the destructor calls delete[] on fIniDir.
    */
   fi.fIniDir = new char[10];
   strcpy(fi.fIniDir, ".");
   new TGFileDialog(gClient->GetDefaultRoot(), guiManager()->getMainFrame(), kFDOpen, &fi);
   guiManager()->updateStatus("loading file ...");
   if (fi.fFilename) {
      m_navigator->appendFile(fi.fFilename, false, false);
      m_loadedAnyInputFile = true;
      checkPosition();
      draw();
   }
   guiManager()->clearStatus();
}

void
CmsShowMain::openDataViaURL()
{
   if (m_searchFiles.get() == 0) {
      m_searchFiles = std::auto_ptr<CmsShowSearchFiles>(new CmsShowSearchFiles("",
                                                                               "Open Remote Data Files",
                                                                               guiManager()->getMainFrame(),
                                                                               500, 400));
      m_searchFiles->CenterOnParent(kTRUE,TGTransientFrame::kBottomRight);
   }
   std::string chosenFile = m_searchFiles->chooseFileFromURL();
   if(!chosenFile.empty()) {
      guiManager()->updateStatus("loading file ...");
      if(m_navigator->openFile(chosenFile.c_str())) {
         m_navigator->firstEvent();
         checkPosition();
         draw();
         guiManager()->clearStatus();
      } else {
         guiManager()->updateStatus("failed to load data file");
      }
   }
}

//
// const member functions
//

//STARTUP TASKS

void
CmsShowMain::loadGeometry()
{   // prepare geometry service
   // ATTN: this should be made configurable
   try {
      guiManager()->updateStatus("Loading geometry...");
      m_detIdToGeo.loadGeometry( m_geomFileName.c_str() );
      m_detIdToGeo.loadMap( m_geomFileName.c_str() );
      m_context->setGeom(&m_detIdToGeo);
   }
   catch (const std::runtime_error& iException)
   {
      fwLog(fwlog::kError) << "CmsShowMain::loadGeometry() caught exception: \n"
                           << iException.what() << std::endl;
      exit(0);
   }
}

void
CmsShowMain::setupViewManagers()
{
   guiManager()->updateStatus("Setting up view manager...");

   boost::shared_ptr<FWViewManagerBase> eveViewManager(new FWEveViewManager(guiManager()));
   eveViewManager->setContext(m_context.get());
   viewManager()->add(eveViewManager);

   boost::shared_ptr<FWTableViewManager> tableViewManager(new FWTableViewManager(guiManager()));
   configurationManager()->add(std::string("Tables"), tableViewManager.get());
   viewManager()->add(tableViewManager);
   eiManager()->goingToClearItems_.connect(boost::bind(&FWTableViewManager::removeAllItems, tableViewManager.get()));

   boost::shared_ptr<FWTriggerTableViewManager> triggerTableViewManager(new FWTriggerTableViewManager(guiManager()));
   configurationManager()->add(std::string("TriggerTables"), triggerTableViewManager.get());
   viewManager()->add( triggerTableViewManager );

   boost::shared_ptr<FWL1TriggerTableViewManager> l1TriggerTableViewManager(new FWL1TriggerTableViewManager(guiManager()));
   configurationManager()->add(std::string("L1TriggerTables"), l1TriggerTableViewManager.get());
   viewManager()->add( l1TriggerTableViewManager );
   
   // Unfortunately, due to the plugin mechanism, we need to delay
   // until here the creation of the FWJobMetadataManager, because
   // otherwise the supportedTypesAndRepresentations map is empty.
   // FIXME: should we have a signal for whenever the above mentioned map
   //        changes? Can that actually happer (maybe if we add support
   //        for loading plugins on the fly??).
   m_metadataManager->initReps(viewManager()->supportedTypesAndRepresentations());
}

//_______________________________________________________________________________
void 
CmsShowMain::autoLoadNewEvent()
{
   stopAutoLoadTimer();
   
   // case when start with no input file
   if (!m_loadedAnyInputFile)
   {
      if (m_monitor.get()) 
         startAutoLoadTimer();
      return;
   }

   bool reachedEnd = (forward() && m_navigator->isLastEvent()) || (!forward() && m_navigator->isFirstEvent());

   if (loop() && reachedEnd)
   {
      forward() ? m_navigator->firstEvent() : m_navigator->lastEvent();
      draw();
   }
   else if (!reachedEnd)
   {
      forward() ? m_navigator->nextEvent() : m_navigator->previousEvent();
      draw();
   }

   // stop loop in case no loop or monitor mode
   if (reachedEnd && (loop() || m_monitor.get()) == kFALSE)
   {
      if (forward() && m_navigator->isLastEvent())
      {
         guiManager()->enableActions();
         checkPosition();
      }

      if ((!forward()) && m_navigator->isFirstEvent())
      {
         guiManager()->enableActions();
         checkPosition();
      }
   }
   else
      startAutoLoadTimer();
}

//______________________________________________________________________________

void 
CmsShowMain::checkPosition()
{
   if ((m_monitor.get() || loop() ) && isPlaying())
      return;
   
   guiManager()->getMainFrame()->enableNavigatorControls();

   if (m_navigator->isFirstEvent())
      guiManager()->disablePrevious();

   if (m_navigator->isLastEvent())
   {
      guiManager()->disableNext();
      // force enable play events action in --port mode
      if (m_monitor.get() && !guiManager()->playEventsAction()->isEnabled())
         guiManager()->playEventsAction()->enable();
   }
}

//==============================================================================
void
CmsShowMain::setupDataHandling()
{
   guiManager()->updateStatus("Setting up data handling...");

   // Event / file change actions which require different response for different
   // implementations. 
   m_navigator->newEvent_.connect(boost::bind(&CmsShowMain::eventChangedSlot, this));
   m_navigator->fileChanged_.connect(boost::bind(&CmsShowMain::fileChangedSlot, this, _1));

   // navigator filtering  ->
   m_navigator->editFiltersExternally_.connect(boost::bind(&FWGUIManager::updateEventFilterEnable, guiManager(), _1));
   m_navigator->filterStateChanged_.connect(boost::bind(&CmsShowMain::navigatorChangedFilterState, this, _1));
   m_navigator->postFiltering_.connect(boost::bind(&CmsShowMain::postFiltering, this));

   // navigator fitlering <-
   guiManager()->showEventFilterGUI_.connect(boost::bind(&CmsShowNavigator::showEventFilterGUI, m_navigator.get(),_1));
   guiManager()->filterButtonClicked_.connect(boost::bind(&CmsShowMain::filterButtonClicked,this));

   // Data handling. File related and therefore not in the base class.
   if (guiManager()->getAction(cmsshow::sOpenData)    != 0) 
      guiManager()->getAction(cmsshow::sOpenData)->activated.connect(sigc::mem_fun(*this, &CmsShowMain::openData));
   if (guiManager()->getAction(cmsshow::sAppendData)  != 0) 
      guiManager()->getAction(cmsshow::sAppendData)->activated.connect(sigc::mem_fun(*this, &CmsShowMain::appendData));
   if (guiManager()->getAction(cmsshow::sSearchFiles) != 0)
      guiManager()->getAction(cmsshow::sSearchFiles)->activated.connect(sigc::mem_fun(*this, &CmsShowMain::openDataViaURL));

   setupActions();
   // init data from  CmsShowNavigator configuration, can do this with signals since there were not connected yet
   guiManager()->setFilterButtonIcon(m_navigator->getFilterState());

   for (unsigned int ii = 0; ii < m_inputFiles.size(); ++ii)
   {
      const std::string& fname = m_inputFiles[ii];
      if (fname.empty())
         continue;
      guiManager()->updateStatus("loading data file ...");
      if (!m_navigator->appendFile(fname, false, false))
      {
         guiManager()->updateStatus("failed to load data file");
         openData();
      }
      else
         m_loadedAnyInputFile = true;
   }

   if (m_loadedAnyInputFile)
   {
      m_navigator->firstEvent();
      checkPosition();
      draw();
   }
   else if (m_monitor.get() == 0)
      openData();
}

void
CmsShowMain::setupSocket(unsigned int iSocket)
{
   m_monitor = std::auto_ptr<TMonitor>(new TMonitor);
   TServerSocket* server = new TServerSocket(iSocket,kTRUE);
   if (server->GetErrorCode())
   {
      fwLog(fwlog::kError) << "CmsShowMain::setupSocket, can't create socket on port "<< iSocket << "." << std::endl;
      exit(0);
   }
   m_monitor->Connect("Ready(TSocket*)","CmsShowMain",this,"notified(TSocket*)");
   m_monitor->Add(server);
}

void
CmsShowMain::notified(TSocket* iSocket)
{
   TServerSocket* server = dynamic_cast<TServerSocket*> (iSocket);
   if (server)
   {
      TSocket* connection = server->Accept();
      if (connection)
      {
         m_monitor->Add(connection);
         std::stringstream s;
         s << "received connection from "<<iSocket->GetInetAddress().GetHostName();
         guiManager()->updateStatus(s.str().c_str());
      }
   }
   else
   {
      char buffer[4096];
      memset(buffer,0,sizeof(buffer));
      if (iSocket->RecvRaw(buffer, sizeof(buffer)) <= 0)
      {
         m_monitor->Remove(iSocket);
         //std::stringstream s;
         //s << "closing connection to "<<iSocket->GetInetAddress().GetHostName();
         //m_guiManager->updateStatus(s.str().c_str());
         delete iSocket;
         return;
      }
      std::string fileName(buffer);
      std::string::size_type lastNonSpace = fileName.find_last_not_of(" \n\t");
      if (lastNonSpace != std::string::npos)
      {
         fileName.erase(lastNonSpace+1);
      }

      std::stringstream s;
      s <<"New file notified '"<<fileName<<"'";
      guiManager()->updateStatus(s.str().c_str());

      bool appended = m_navigator->appendFile(fileName, true, live());

      if (appended)
      {
         if (live() && isPlaying() && forward())
            m_navigator->activateNewFileOnNextEvent();
         else if (!isPlaying())
            checkPosition();

         // bootstrap case: --port  and no input file
         if (!m_loadedAnyInputFile)
         {
            m_loadedAnyInputFile = true;
            m_navigator->firstEvent();
            if (!isPlaying())
               draw();
         }

         std::stringstream sr;
         sr <<"New file registered '"<<fileName<<"'";
         guiManager()->updateStatus(sr.str().c_str());
      }
      else
      {
         std::stringstream sr;
         sr <<"New file NOT registered '"<<fileName<<"'";
         guiManager()->updateStatus(sr.str().c_str());
      }
   }
}

void
CmsShowMain::stopPlaying()
{
   stopAutoLoadTimer();
   if (live())
      m_navigator->resetNewFileOnNextEvent();
   setIsPlaying(false);
   guiManager()->enableActions();
   checkPosition();
}

void
CmsShowMain::navigatorChangedFilterState(int state)
{
   guiManager()->setFilterButtonIcon(state);
   if (m_navigator->filesNeedUpdate() == false)
   {
      guiManager()->setFilterButtonText(m_navigator->filterStatusMessage());
      checkPosition();
   }
}

void
CmsShowMain::filterButtonClicked()
{
   if (m_navigator->getFilterState() == CmsShowNavigator::kWithdrawn )
      guiManager()->showEventFilterGUI();
   else
      m_navigator->toggleFilterEnable();
}

void
CmsShowMain::preFiltering()
{
   // called only if filter has changed
   guiManager()->updateStatus("Filtering events");
}

void
CmsShowMain::postFiltering()
{
   // called only filter is changed
   guiManager()->clearStatus();
   draw();
   checkPosition();
   guiManager()->setFilterButtonText(m_navigator->filterStatusMessage());
}

