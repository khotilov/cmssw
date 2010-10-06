#include "TauAnalysis/Core/plugins/GenericAnalyzer.h"

#include "TauAnalysis/DQMTools/interface/dqmAuxFunctions.h"
#include "TauAnalysis/DQMTools/interface/generalAuxFunctions.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "TauAnalysis/Core/interface/SysUncertaintyService.h"
#include "TauAnalysis/Core/interface/sysUncertaintyAuxFunctions.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include <algorithm>
#include <limits.h>
#include <iostream>

unsigned GenericAnalyzer::analysisSequenceEntry_filter::filterId_ = 1;

GenericAnalyzer::analysisSequenceEntry::analysisSequenceEntry(const std::string& name)
  : name_(name)
{}

GenericAnalyzer::analysisSequenceEntry::~analysisSequenceEntry()
{
// nothing to be done yet...
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

GenericAnalyzer::analysisSequenceEntry_filter::analysisSequenceEntry_filter(const std::string& name, const std::string& title,
									    const edm::ParameterSet& cfgFilter, 
									    bool estimateSysUncertainties, int& cfgError)
  : analysisSequenceEntry(name),
    estimateSysUncertainties_(estimateSysUncertainties)
{
  //std::cout << "<analysisSequenceEntry_filter::analysisSequenceEntry_filter>:" << std::endl;
  //std::cout << " name = " << name_ << std::endl;
  
  std::string filterType = cfgFilter.getParameter<std::string>("pluginType");

  EventSelectorBase* filterPlugin_cumulative = 0;
  EventSelectorBase* filterPlugin_individual = 0;

  if ( cfgFilter.exists("src_cumulative") &&
       cfgFilter.exists("src_individual") ) {
    edm::ParameterSet cfgFilter_cumulative = cfgFilter;
    cfgFilter_cumulative.addParameter<edm::InputTag>("src", cfgFilter.getParameter<edm::InputTag>("src_cumulative"));
    filterPlugin_cumulative = EventSelectorPluginFactory::get()->create(filterType, cfgFilter_cumulative);

    edm::ParameterSet cfgFilter_individual = cfgFilter;
    cfgFilter_individual.addParameter<edm::InputTag>("src", cfgFilter.getParameter<edm::InputTag>("src_individual"));
    filterPlugin_individual = EventSelectorPluginFactory::get()->create(filterType, cfgFilter_individual);
  } else {
    filterPlugin_cumulative = EventSelectorPluginFactory::get()->create(filterType, cfgFilter);

    filterPlugin_individual = EventSelectorPluginFactory::get()->create(filterType, cfgFilter);
  }

  const std::string& nameCentralValue = SysUncertaintyService::getNameCentralValue();

  filterPlugins_cumulative_.insert(std::pair<std::string, EventSelectorBase*>(nameCentralValue, filterPlugin_cumulative));
  filterPlugins_individual_.insert(std::pair<std::string, EventSelectorBase*>(nameCentralValue, filterPlugin_individual));

//--- handle systematic uncertainties:
//    set InputTag to name of collection with systematic shifts applied
//
//    NOTE: for estimating systematic uncertainties
//          need to take into account only "cumulative" event selection
//          ("individual" event selection has effect on filter statistics table only,
//           not on histogram filling or computation of binning results)
//
  if ( estimateSysUncertainties && cfgFilter.exists("systematics") ) {
    typedef std::vector<std::string> vstring;
    vstring systematics = cfgFilter.getParameter<vstring>("systematics");
    for ( vstring::const_iterator sysName = systematics.begin();
	  sysName != systematics.end(); ++sysName ) {
      //std::cout << " sysName = " << (*sysName) << std::endl;

      edm::InputTag src;
      if ( cfgFilter.exists("src_cumulative") ) {
	src = cfgFilter.getParameter<edm::InputTag>("src_cumulative");
      } else { 
	src = cfgFilter.getParameter<edm::InputTag>("src");
      }

      //std::cout << " src = " << src.label() << std::endl;

//--- compose name of InputTag
//
// NOTE: string concatenation rule (first character of systematics name is capitalized) for composition of InputTag
//       needs to match implementation in (python) function composeModuleName
//       defined in TauAnalysis/CandidateTools/python/tools/composeModuleName.py
//
      std::string src_label = src.label();
      std::string src_instance = src.instance();

      if ( sysName->length() >= 1 ) src_label += toupper(sysName->at(0));
      if ( sysName->length() >= 2 ) src_label += std::string(*sysName, 1);

      //std::cout << "src_label = " << src_label << std::endl; 

      edm::ParameterSet cfgFilter_systematic(cfgFilter);
      cfgFilter_systematic.addParameter<edm::InputTag>("src", edm::InputTag(src_label, src_instance));
      EventSelectorBase* filterPlugin_systematic = EventSelectorPluginFactory::get()->create(filterType, cfgFilter_systematic);

      filterPlugins_cumulative_.insert(std::pair<std::string, EventSelectorBase*>(*sysName, filterPlugin_systematic));
    }
  }

  ++filterId_;

  //print();
}

GenericAnalyzer::analysisSequenceEntry_filter::~analysisSequenceEntry_filter()
{
  for ( std::map<std::string, EventSelectorBase*>::const_iterator it = filterPlugins_cumulative_.begin();
	it != filterPlugins_cumulative_.end(); ++it ) {
    delete it->second;
  }
  for ( std::map<std::string, EventSelectorBase*>::const_iterator it = filterPlugins_individual_.begin();
	it != filterPlugins_individual_.end(); ++it ) {
    delete it->second;
  }
}

void GenericAnalyzer::analysisSequenceEntry_filter::print() const
{
  std::cout << "<GenericAnalyzer::analysisSequenceEntry_filter::print>:" << std::endl; 
  std::cout << " name = " << name_ << std::endl;
  std::cout << " filterPlugins_cumulative:" << std::endl;
  for ( std::map<std::string, EventSelectorBase*>::const_iterator filterPlugin_cumulative = filterPlugins_cumulative_.begin();
	filterPlugin_cumulative != filterPlugins_cumulative_.end(); ++filterPlugin_cumulative ) {
    std::cout << "  " << filterPlugin_cumulative->first << ": " << filterPlugin_cumulative->second << std::endl;
  }
  std::cout << " filterPlugins_individual:" << std::endl;
  for ( std::map<std::string, EventSelectorBase*>::const_iterator filterPlugin_individual = filterPlugins_individual_.begin();
	filterPlugin_individual != filterPlugins_individual_.end(); ++filterPlugin_individual ) {
    std::cout << "  " << filterPlugin_individual->first << ": " << filterPlugin_individual->second << std::endl;
  }
  std::cout << " filterId = " << filterId_ << std::endl;
}

bool GenericAnalyzer::analysisSequenceEntry_filter::filter(const edm::Event& evt, const edm::EventSetup& es, 
							   const SysUncertaintyService* sysUncertaintyService,
							   const std::map<std::string, EventSelectorBase*>& filterPlugins)
{
//--- check if specific event selector module is defined for current systematic uncertainty;
//    if so evaluate that event selector module, 
//    else evaluate the event selector module for no systematic shifts applied

  //std::cout << "<analysisSequenceEntry_filter::filter>:" << std::endl;
  //std::cout << " name = " << name_ << std::endl;

  std::map<std::string, EventSelectorBase*>::const_iterator it = ( estimateSysUncertainties_ ) ?
    filterPlugins.find(sysUncertaintyService->getCurrentSystematic()) : filterPlugins.end();
  if ( it == filterPlugins.end() ) it = filterPlugins.find(SysUncertaintyService::getNameCentralValue());
  if ( it == filterPlugins.end() ) {
    edm::LogError ("filter") << " No event selector plugin defined for central value !!";
    return false;
  };
  
  bool filterDecision = true;
  try { 
    edm::Event* evt_nonConst = const_cast<edm::Event*>(&evt);
    filterDecision = (*it->second)(*evt_nonConst, es);
    //std::cout << " filterDecision = " << filterDecision << std::endl;
  } catch ( cms::Exception& e ) {
    edm::LogError("filter") 
      << " Filter plugin name = " << name_ << " caused exception --> rethrowing !!";
    throw e;
  }
  
  return filterDecision;
}

bool GenericAnalyzer::analysisSequenceEntry_filter::filter_cumulative(const edm::Event& evt, const edm::EventSetup& es, 
								      const SysUncertaintyService* sysUncertaintyService)
{
  return filter(evt, es, sysUncertaintyService, filterPlugins_cumulative_);
}

bool GenericAnalyzer::analysisSequenceEntry_filter::filter_individual(const edm::Event& evt, const edm::EventSetup& es, 
								      const SysUncertaintyService* sysUncertaintyService)
{
  return filter(evt, es, sysUncertaintyService, filterPlugins_individual_);
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

GenericAnalyzer::analysisSequenceEntry_analyzer::analysisSequenceEntry_analyzer(const std::string& name, 
										const std::list<edm::ParameterSet>& cfgAnalyzers)
  : analysisSequenceEntry(name)
{
  for ( std::list<edm::ParameterSet>::const_iterator cfgAnalyzer = cfgAnalyzers.begin(); 
	cfgAnalyzer != cfgAnalyzers.end(); ++cfgAnalyzer ) {
    analyzerPluginEntry analyzerEntry;
    std::string analyzerType = cfgAnalyzer->getParameter<std::string>("pluginType");
    analyzerEntry.plugin_ = AnalyzerPluginFactory::get()->create(analyzerType, *cfgAnalyzer);
    analyzerEntry.supportsSystematics_ = cfgAnalyzer->getParameter<bool>("supportsSystematics");
    analyzerPlugins_.push_back(analyzerEntry);
  }

  //print();
}

GenericAnalyzer::analysisSequenceEntry_analyzer::~analysisSequenceEntry_analyzer()
{
  for ( std::list<analyzerPluginEntry>::const_iterator analyzer = analyzerPlugins_.begin();
	analyzer != analyzerPlugins_.end(); ++analyzer ) {
    delete analyzer->plugin_;
  }
}

void GenericAnalyzer::analysisSequenceEntry_analyzer::print() const
{
  std::cout << "<GenericAnalyzer::analysisSequenceEntry_analyzer::print>:" << std::endl;
}

void GenericAnalyzer::analysisSequenceEntry_analyzer::beginJob()
{
  for ( std::list<analyzerPluginEntry>::const_iterator analyzer = analyzerPlugins_.begin();
	analyzer != analyzerPlugins_.end(); ++analyzer ) {
    analyzer->plugin_->beginJob();
  }
}

void GenericAnalyzer::analysisSequenceEntry_analyzer::analyze(const edm::Event& evt, const edm::EventSetup& es, 
							      double evtWeight, bool isSystematicApplied)
{
  for ( std::list<analyzerPluginEntry>::const_iterator analyzer = analyzerPlugins_.begin();
	analyzer != analyzerPlugins_.end(); ++analyzer ) {
    
    if ( analyzer->supportsSystematics_ == false && isSystematicApplied ) continue;

    try { 
      analyzer->plugin_->analyze(evt, es, evtWeight);
    } catch ( cms::Exception& e ) {
      edm::LogError("analyze") 
	<< " Analyzer plugin name = " << analyzer->plugin_->name() << " caused exception --> rethrowing !!";
      throw e;
    }
  }
}

void GenericAnalyzer::analysisSequenceEntry_analyzer::endJob()
{
  for ( std::list<analyzerPluginEntry>::const_iterator analyzer = analyzerPlugins_.begin();
	analyzer != analyzerPlugins_.end(); ++analyzer ) {
    analyzer->plugin_->endJob();
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void GenericAnalyzer::addFilter(const std::string& filterName, const vstring& saveRunEventNumbers)
{
  std::map<std::string, edm::ParameterSet>::const_iterator it = cfgFilters_.find(filterName);
  if ( it != cfgFilters_.end() ) {
    edm::ParameterSet cfgFilter = it->second;
    
    std::string filterTitle = ( cfgFilter.exists("title") ) ? cfgFilter.getParameter<std::string>("title") : filterName;
    
    analysisSequenceEntry_filter* entry 
      = new analysisSequenceEntry_filter(filterName, filterTitle, cfgFilter, estimateSysUncertainties_, cfgError_);
    analysisSequence_.push_back(entry);
  } else {
    edm::LogError("GenericAnalyzer::addFilter") 
      << " Failed to access configuration parameter for filter = " << filterName << " --> skipping !!";
    cfgError_ = 1;
  }
}

std::string capitalizeFirstCharacter(const std::string& input)
{
  std::string output = input;
  if ( output != "" ) std::transform(output.begin(), ++output.begin(), output.begin(), ::toupper);
  return output;
}

std::string dqmSubDirectoryName_filter(const std::string& afterFilterName, const std::string& beforeFilterName)
{
  std::string subDirName;
  if ( afterFilterName != "" ) subDirName.append("after").append(capitalizeFirstCharacter(afterFilterName));
  if ( afterFilterName != "" && beforeFilterName != "" ) subDirName.append("_");
  if ( beforeFilterName != "" ) subDirName.append("before").append(capitalizeFirstCharacter(beforeFilterName));
  return subDirName;
}

void GenericAnalyzer::addAnalyzers(const vstring& analyzerNames, 
				   const std::string& afterFilterName, const std::string& beforeFilterName,
				   const vstring& replaceCommands)
{
  std::list<edm::ParameterSet> cfgAnalyzers;

  for ( vstring::const_iterator analyzerName = analyzerNames.begin();
	analyzerName != analyzerNames.end(); ++analyzerName ) {

    std::map<std::string, edm::ParameterSet>::const_iterator it = cfgAnalyzers_.find(*analyzerName);
    if ( it != cfgAnalyzers_.end() ) {
      edm::ParameterSet cfgAnalyzer = it->second;
      
      for ( vstring::const_iterator replaceCommand = replaceCommands.begin();
	    replaceCommand != replaceCommands.end(); ++replaceCommand ) {
	int errorFlag;
	std::string replaceCommand_noWhiteSpace = replace_string(*replaceCommand, " ", "", 0, UINT_MAX, errorFlag);

	std::string keyword = std::string(*analyzerName).append(".");
	if ( replaceCommand_noWhiteSpace.find(keyword) != std::string::npos ) {
	  std::string replaceCommand_noKeyword = replace_string(replaceCommand_noWhiteSpace, keyword, "", 0, UINT_MAX, errorFlag);

	  size_t posAssignmentOperator = replaceCommand_noKeyword.find("=");
	  if ( posAssignmentOperator == std::string::npos ||
	       !(posAssignmentOperator >= 1 && posAssignmentOperator < (replaceCommand_noKeyword.length() - 1)) ) {
	    edm::LogError("GenericAnalyzer::addAnalyzers") 
	      << " Failed to parse replaceCommand = " << (*replaceCommand) << " --> skipping !!";
	    cfgError_ = 1;
	    continue;
	  }
	
	  std::string leftHandSide = std::string(replaceCommand_noKeyword, 0, posAssignmentOperator);
	  std::string rightHandSide = std::string(replaceCommand_noKeyword, posAssignmentOperator + 1);

//--- check that parameter that is to be replaced
//    exists in configuration parameter set 
//    and is defined to be a supported type
//    (so far, only std::string and edm::InputTag are supported types
//     for configuration parameter replacements)
	  if ( cfgAnalyzer.existsAs<std::string>(leftHandSide) ) {
	    cfgAnalyzer.addParameter<std::string>(leftHandSide, rightHandSide);
	  } else if ( cfgAnalyzer.existsAs<edm::InputTag>(leftHandSide) ) {
	    edm::InputTag rightHandSide_inputTag = rightHandSide;
	    cfgAnalyzer.addParameter<edm::InputTag>(leftHandSide, rightHandSide_inputTag);
	  } else if ( cfgAnalyzer.existsAs<bool>(leftHandSide) ) {
	    bool rightHandSide_bool;
	    if ( rightHandSide == "True" || rightHandSide == "true" ) {
	      rightHandSide_bool = true;
	    } else if ( rightHandSide == "False" || rightHandSide == "false" ) {
	      rightHandSide_bool = false;
	    } else {
	      edm::LogError("GenericAnalyzer::addAnalyzers") 
		<< " Invalid value = " << rightHandSide << " for replacing"
		<< " Configuration parameter = " << leftHandSide << " of boolean type;"
		<< " valid values for Replacement = { 'True', 'true', 'False', 'false' }" << " --> skipping !!";
	      cfgError_ = 1;
	      continue;
	    }
	    cfgAnalyzer.addParameter<bool>(leftHandSide, rightHandSide_bool);
	  } else {
	    edm::LogError("GenericAnalyzer::addAnalyzers") 
	      << " Configuration parameter to be replaced = " << leftHandSide 
	      << " does either not exist or is not of a supported type" << " --> skipping !!";
	    cfgError_ = 1;
	    continue;
	  }

//--- store modified configuration ParameterSet
//    (NOTE: by storing the modified configuration ParameterSet,
//           replace statements become cumulative)
	  cfgAnalyzers_[it->first] = cfgAnalyzer;
	}
      }

      cfgAnalyzer.addParameter<std::string>("afterFilterName", afterFilterName);
      cfgAnalyzer.addParameter<std::string>("beforeFilterName", beforeFilterName);

      if ( cfgAnalyzer.exists("dqmDirectory_store") ) {
	std::string dqmDirectory_store = cfgAnalyzer.getParameter<std::string>("dqmDirectory_store");
	std::string dirName = dqmDirectoryName(name_).append(dqmSubDirectoryName_filter(afterFilterName, beforeFilterName));
	std::string dirName_full = dqmDirectoryName(dirName).append(dqmDirectory_store);
	cfgAnalyzer.addParameter<std::string>("dqmDirectory_store", dirName_full);
      }

      cfgAnalyzers.push_back(cfgAnalyzer);
    } else {
      edm::LogError("GenericAnalyzer::addAnalyzers") 
	<< " Failed to access configuration parameter for analyzer = " << (*analyzerName) << " --> skipping !!";
      cfgError_ = 1;
      return;
    }
  }

  std::string entryName = dqmDirectoryName(name_).append(dqmSubDirectoryName_filter(afterFilterName, beforeFilterName));
  analysisSequenceEntry_analyzer* entry = new analysisSequenceEntry_analyzer(entryName, cfgAnalyzers);
  analysisSequence_.push_back(entry);
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

GenericAnalyzer::GenericAnalyzer(const edm::ParameterSet& cfg)
{
  std::cout << "<GenericAnalyzer::GenericAnalyzer>:" << std::endl;

  cfgError_ = 0;

  name_ = cfg.getParameter<std::string>("name");

  eventWeightSrc_ = ( cfg.exists("eventWeightSource") ) ? cfg.getParameter<vInputTag>("eventWeightSource") : vInputTag();

  estimateSysUncertainties_ = cfg.exists("estimateSysUncertainties") ?
    cfg.getParameter<bool>("estimateSysUncertainties") : false;

  typedef std::vector<edm::ParameterSet> vParameterSet;

//--- store configuration parameters for filters
  //std::cout << "--> storing configuration parameters for filters..." << std::endl;
  if ( cfg.exists("filters") ) {
    vParameterSet cfgFilters = cfg.getParameter<vParameterSet>("filters");
    for ( vParameterSet::const_iterator cfgFilter = cfgFilters.begin(); 
	  cfgFilter != cfgFilters.end(); ++cfgFilter ) {
      std::string cfgFilterName = cfgFilter->getParameter<std::string>("pluginName");
      cfgFilters_.insert(std::pair<std::string, edm::ParameterSet>(cfgFilterName, *cfgFilter)); 
    }
  }

//--- store configuration parameters for analyzers
//    (with and without support for estimating systematic uncertainties)
  //std::cout << "--> storing configuration parameters for analyzers..." << std::endl;
  if ( cfg.exists("analyzers") ) {
    vParameterSet cfgAnalyzers = cfg.getParameter<vParameterSet>("analyzers");
    for ( vParameterSet::iterator cfgAnalyzer = cfgAnalyzers.begin(); 
	  cfgAnalyzer != cfgAnalyzers.end(); ++cfgAnalyzer ) {
      cfgAnalyzer->addParameter<bool>("supportsSystematics", false);
      std::string cfgAnalyzerName = cfgAnalyzer->getParameter<std::string>("pluginName");
      cfgAnalyzers_.insert(std::pair<std::string, edm::ParameterSet>(cfgAnalyzerName, *cfgAnalyzer));
    }
  }

  if ( estimateSysUncertainties_ && cfg.exists("analyzers_systematic") ) {
    vParameterSet cfgAnalyzers_systematic = cfg.getParameter<vParameterSet>("analyzers_systematic");
    for ( vParameterSet::iterator cfgAnalyzer_systematic = cfgAnalyzers_systematic.begin(); 
	  cfgAnalyzer_systematic != cfgAnalyzers_systematic.end(); ++cfgAnalyzer_systematic ) {
      cfgAnalyzer_systematic->addParameter<bool>("supportsSystematics", true);
      std::string cfgAnalyzerName = cfgAnalyzer_systematic->getParameter<std::string>("pluginName");
      cfgAnalyzers_.insert(std::pair<std::string, edm::ParameterSet>(cfgAnalyzerName, *cfgAnalyzer_systematic));
    }
  }

//--- configure names of systematic uncertainties
//    to be taken into account in analysis
  if ( estimateSysUncertainties_ ) {
    vstring cfgSystematics = cfg.getParameter<vstring>("systematics");
    for ( vstring::const_iterator sysName = cfgSystematics.begin();
	  sysName != cfgSystematics.end(); ++sysName ) {
      vstring sysNames_expanded = expandSysName(*sysName);      
      systematics_.insert(systematics_.end(), sysNames_expanded.begin(), sysNames_expanded.end());
    }
  }

  systematics_.insert(systematics_.begin(), SysUncertaintyService::getNameCentralValue());
  std::cout << " systematics = " << format_vstring(systematics_) << std::endl;

//--- configure analysisSequence
  //std::cout << "--> configuring analysisSequence..." << std::endl;
  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgAnalysisSequenceEntries = cfg.getParameter<vParameterSet>("analysisSequence");
  std::string lastFilterName = "";
  for ( vParameterSet::const_iterator cfgAnalysisSequenceEntry = cfgAnalysisSequenceEntries.begin();
	cfgAnalysisSequenceEntry != cfgAnalysisSequenceEntries.end(); ++cfgAnalysisSequenceEntry ) {
//--- check that analysisSequenceEntry is either a filter or an analyzer
    if ( !(cfgAnalysisSequenceEntry->exists("filter")    ||
	   cfgAnalysisSequenceEntry->exists("analyzers")) ) {
      edm::LogError("GenericAnalyzer") << " Either filter or analyzers must be specified for sequenceEntries !!";
      cfgError_ = 1;
      continue;
    }
 
    if ( cfgAnalysisSequenceEntry->exists("filter")      &&
	 cfgAnalysisSequenceEntry->exists("analyzers") ) {
      edm::LogError("GenericAnalyzer") << " Must not specify filter and analyzers for same sequenceEntry !!";
      cfgError_ = 1;
      continue;
    }

    if ( cfgAnalysisSequenceEntry->exists("filter") ) {
      std::string filterName = cfgAnalysisSequenceEntry->getParameter<std::string>("filter");

      vstring saveRunEventNumbers = cfgAnalysisSequenceEntry->exists("saveRunEventNumbers") ? 
	cfgAnalysisSequenceEntry->getParameter<vstring>("saveRunEventNumbers") : vstring();
      
      addFilter(filterName, saveRunEventNumbers);
      
      lastFilterName = filterName;
    }

    if ( cfgAnalysisSequenceEntry->exists("analyzers") ) {
      vstring analyzerNames = cfgAnalysisSequenceEntry->getParameter<vstring>("analyzers");

      if ( analyzerNames.size() == 0 ) {
	edm::LogWarning("GenericAnalyzer") << " List of analyzers is empty !!";
      }

      std::string nextFilterName = "";
      for ( vParameterSet::const_iterator cfgAnalysisSequenceEntry_nextFilter = cfgAnalysisSequenceEntry;
	    cfgAnalysisSequenceEntry_nextFilter != cfgAnalysisSequenceEntries.end(); ++cfgAnalysisSequenceEntry_nextFilter ) {
	if (  cfgAnalysisSequenceEntry_nextFilter->exists("filter") ) {
	  nextFilterName = cfgAnalysisSequenceEntry_nextFilter->getParameter<std::string>("filter");
	  break;
	}
      }

      vstring replaceCommands = ( cfgAnalysisSequenceEntry->exists("replace") ) ? 
	cfgAnalysisSequenceEntry->getParameter<vstring>("replace") : vstring();
      //std::cout << " replaceCommands = " << format_vstring(replaceCommands) << std::endl;
      
      addAnalyzers(analyzerNames, lastFilterName, nextFilterName, replaceCommands);
    }
  }

  if ( cfgError_ ) return;

//--- configure filter statistics service
//    and create filter statistics table
  filterStatisticsService_ = new FilterStatisticsService();
  filterStatisticsService_dqmDirectory_ = dqmDirectoryName(name_).append("FilterStatistics");

  edm::ParameterSet cfgFilterStatisticsTable;
  cfgFilterStatisticsTable.addParameter<std::string>("name", std::string(name_).append("-").append("FilterStatisticsService"));
  vParameterSet filterStatisticsTable_config;
  for ( vParameterSet::const_iterator cfgAnalysisSequenceEntry = cfgAnalysisSequenceEntries.begin();
	cfgAnalysisSequenceEntry != cfgAnalysisSequenceEntries.end(); ++cfgAnalysisSequenceEntry ) {
    if ( cfgAnalysisSequenceEntry->exists("filter") ) {
      std::string filterName = cfgAnalysisSequenceEntry->getParameter<std::string>("filter");
      std::string filterTitle = ( cfgAnalysisSequenceEntry->exists("title") ) ? 
	cfgAnalysisSequenceEntry->getParameter<std::string>("title") : filterName;
          
      edm::ParameterSet config;
      config.addParameter<std::string>("filterName", filterName);
      config.addParameter<std::string>("filterTitle", filterTitle);
      
      filterStatisticsTable_config.push_back(config);
    }
  }
  cfgFilterStatisticsTable.addParameter<vParameterSet>("config", filterStatisticsTable_config);
  filterStatisticsTable_ = filterStatisticsService_->createFilterStatisticsTable(cfgFilterStatisticsTable);
 
//--- configure run & event number service
  edm::ParameterSet cfgRunEventNumberService;
  cfgRunEventNumberService.addParameter<std::string>("name", std::string(name_).append("-").append("RunEventNumberService"));
  cfgRunEventNumberService.addParameter<std::string>("dqmDirectory_store", dqmDirectoryName(name_).append("FilterStatistics"));
  vParameterSet runEventNumberService_config;
  for ( vParameterSet::const_iterator cfgAnalysisSequenceEntry = cfgAnalysisSequenceEntries.begin();
	cfgAnalysisSequenceEntry != cfgAnalysisSequenceEntries.end(); ++cfgAnalysisSequenceEntry ) {
    if ( cfgAnalysisSequenceEntry->exists("filter") ) {
      std::string filterName = cfgAnalysisSequenceEntry->getParameter<std::string>("filter");
      vstring saveRunEventNumbers = cfgAnalysisSequenceEntry->exists("saveRunEventNumbers") ? 
	cfgAnalysisSequenceEntry->getParameter<vstring>("saveRunEventNumbers") : vstring();
      
      edm::ParameterSet config;
      config.addParameter<std::string>("filterName", filterName);
      config.addParameter<vstring>("saveRunEventNumbers", saveRunEventNumbers);
      
      runEventNumberService_config.push_back(config);
    }
  }
  cfgRunEventNumberService.addParameter<vParameterSet>("config", runEventNumberService_config);
  runEventNumberService_ = new RunEventNumberService(cfgRunEventNumberService);

//--- configure eventDumps
  if ( cfg.exists("eventDumps") ) {
    //std::cout << "--> configuring eventDumps..." << std::endl;
    vParameterSet cfgEventDumps = cfg.getParameter<vParameterSet>("eventDumps");
    for ( vParameterSet::const_iterator cfgEventDump = cfgEventDumps.begin(); 
	  cfgEventDump != cfgEventDumps.end(); ++cfgEventDump ) {
      std::string eventDumpName = cfgEventDump->getParameter<std::string>("pluginName");
      std::string eventDumpType = cfgEventDump->getParameter<std::string>("pluginType");
 
      EventDumpBase* entry = EventDumpPluginFactory::get()->create(eventDumpType, *cfgEventDump);
      eventDumps_.push_back(entry);
    }
  }
}

GenericAnalyzer::~GenericAnalyzer()
{
  for ( std::list<analysisSequenceEntry*>::const_iterator entry = analysisSequence_.begin();
	entry != analysisSequence_.end(); ++entry ) {
    delete (*entry);
  }

  delete filterStatisticsService_;
  delete filterStatisticsTable_;
  delete runEventNumberService_;
}

void GenericAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{  
  //std::cout << "<GenericAnalyzer::analyze>:" << std::endl; 

//--- check that configuration parameters contain no errors
  if ( cfgError_ ) {
    edm::LogError("GenericAnalyzer::analyze") << " Error in Configuration ParameterSet --> skipping !!";
    return;
  }

  double eventWeight = 1.;
  for ( vInputTag::const_iterator eventWeightSrc_i = eventWeightSrc_.begin();
	eventWeightSrc_i != eventWeightSrc_.end(); ++eventWeightSrc_i ) {
    edm::Handle<double> eventWeight_i;
    evt.getByLabel(*eventWeightSrc_i, eventWeight_i);

    eventWeight *= (*eventWeight_i);
  }

  //std::cout << " eventWeight = " << eventWeight << std::endl;

  SysUncertaintyService* sysUncertaintyService = 0;
  if ( edm::Service<SysUncertaintyService>().isAvailable() ) {
    sysUncertaintyService = &(*edm::Service<SysUncertaintyService>());
    sysUncertaintyService->update(SysUncertaintyService::getNameCentralValue(), evt, es);
  } else if ( estimateSysUncertainties_ ) {
    edm::LogError ("GenericAnalyzer::analyze") << " Failed to access SysUncertaintyService --> skipping !!";
    return;
  }

//--- estimate systematic uncertainties:
//    iterate over names of systematics uncertainties,
//    set name as currently "active" systematic in SysUncertaintyService;
//    filter and analyzer plugins will "pick-up" name from SysUncertaintyService
//    and take care of correct handling of systematic internally
//
//    NOTE: there are two distinct "types" of systematic uncertainties;
//          the first is handled by reweighting the event 
//          (used e.g. for estimation of PDF, ISR/FSR uncertainties),
//          the second type is handled by shifting values of observables and reapplying the event selection
//          (used e.g. for estimation of uncertainties due to imprecise knowledge of tau-jet energy scale)
//
  for ( vstring::const_iterator sysName = systematics_.begin();
	sysName != systematics_.end(); ++sysName ) {
    //std::cout << " sysName = " << (*sysName) << std::endl;

    double eventWeight_systematic = eventWeight;
    if ( estimateSysUncertainties_ ) {
      sysUncertaintyService->update(*sysName, evt, es);
      eventWeight_systematic *= sysUncertaintyService->getWeight();
    }

    //std::cout << " eventWeight_systematic = " << eventWeight_systematic << std::endl;

    bool isSystematicApplied = ( (*sysName) == SysUncertaintyService::getNameCentralValue() ) ? false : true;
    //std::cout << " isSystematicApplied = " << isSystematicApplied << std::endl;

//--- call analyze method of each analyzerPlugin
//    (fill histograms, compute binning results,...)
    typedef std::vector<std::pair<std::string, bool> > filterResults_type;
    filterResults_type filterResults_cumulative;
    bool previousFiltersPassed = true;
    filterResults_type filterResults_individual;
    for ( std::list<analysisSequenceEntry*>::iterator entry = analysisSequence_.begin();
	  entry != analysisSequence_.end(); ++entry ) {
      
      //std::cout << " analysisSequenceEntry: name = " << (*entry)->name_ << std::endl;

      if ( (*entry)->type() == analysisSequenceEntry::kFilter ) {
	bool filterPassed_cumulative = (*entry)->filter_cumulative(evt, es, sysUncertaintyService);
	//std::cout << " filter: passed_cumulative = " << filterPassed_cumulative << std::endl;
	
	filterResults_cumulative.push_back(std::pair<std::string, bool>((*entry)->name_, filterPassed_cumulative));
	
	if ( !filterPassed_cumulative ) previousFiltersPassed = false;
	
	bool filterPassed_individual = (*entry)->filter_individual(evt, es, sysUncertaintyService);
	//std::cout << " filter: passed_individual = " << filterPassed_individual << std::endl;
	
	filterResults_individual.push_back(std::pair<std::string, bool>((*entry)->name_, filterPassed_individual));
      }

      if ( (*entry)->type() == analysisSequenceEntry::kAnalyzer ) {
	//std::cout << " analyzer: previousFiltersPassed = " << previousFiltersPassed << std::endl;
	if ( previousFiltersPassed ) (*entry)->analyze(evt, es, eventWeight_systematic, isSystematicApplied);
      }
    }

    if ( !isSystematicApplied ) {

//--- update filter statistics table
      filterStatisticsTable_->update(filterResults_cumulative, filterResults_individual, eventWeight);

//--- save run & event numbers
      runEventNumberService_->update(evt.id().run(), evt.id().event(), evt.luminosityBlock(),
				     filterResults_cumulative, filterResults_individual, eventWeight);

//--- if requested, dump event information 
      for ( std::list<EventDumpBase*>::const_iterator it = eventDumps_.begin();
	    it != eventDumps_.end(); ++it ) {
	EventDumpBase* eventDump = (*it);
	eventDump->analyze(evt, es, filterResults_cumulative, filterResults_individual, eventWeight);
      }
    }
  }
}

void GenericAnalyzer::beginJob()
{
  //std::cout << "<GenericAnalyzer::beginJob>:" << std::endl;

  benchmark_.Start(std::string("GenericAnalyzer").append("-").append(name_).data());

//--- call beginJob method of each EDFilter/EDAnalyzer
  for ( std::list<analysisSequenceEntry*>::iterator entry = analysisSequence_.begin();
	entry != analysisSequence_.end(); ++entry ) {
    (*entry)->beginJob();
  }
}

void GenericAnalyzer::endJob()
{
  //std::cout << "<GenericAnalyzer::endJob>:" << std::endl;

  if ( cfgError_ ) return;

//--- call endJob method of each EDFilter/EDAnalyzer
  for ( std::list<analysisSequenceEntry*>::iterator entry = analysisSequence_.begin();
	entry != analysisSequence_.end(); ++entry ) {
    (*entry)->endJob();
  }

  benchmark_.Show(std::string("GenericAnalyzer").append("-").append(name_).data());

//--- print filter statistics (cut-flow) table
//    and save results in DQM file
  filterStatisticsTable_->print(std::cout);
  filterStatisticsService_->saveFilterStatisticsTable(filterStatisticsService_dqmDirectory_, filterStatisticsTable_);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GenericAnalyzer);

