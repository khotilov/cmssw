#ifndef DQMSERVICES_COMPONEntS_DQMFILESAVER_H
# define DQMSERVICES_COMPONEntS_DQMFILESAVER_H

# include "FWCore/Framework/interface/EDAnalyzer.h"
# include "DQMServices/Core/interface/MonitorElement.h"
# include <sys/time.h>
# include <string>

class DQMStore;
class DQMFileSaver : public edm::EDAnalyzer
{
public:
  DQMFileSaver(const edm::ParameterSet &ps);

protected:
  virtual void beginJob(const edm::EventSetup &);
  virtual void beginRun(const edm::Run &, const edm::EventSetup &);
  virtual void beginLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup &);
  virtual void analyze(const edm::Event &e, const edm::EventSetup &);
  virtual void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup &);
  virtual void endRun(const edm::Run &, const edm::EventSetup &);
  virtual void endJob(void);

private:
  void saveForOffline(const std::string &workflow, int run, int lumi=0);
  void saveForOnline(const std::string &suffix, const std::string &rewrite);

  enum Convention
  {
    Online,
    Offline
  };

  Convention	convention_;
  std::string	workflow_;
  std::string	producer_;
  std::string	dirName_;
  int        	version_;
  bool		runIsComplete_;

  int		saveByLumiSection_;
  int		saveByEvent_;
  int		saveByMinute_;
  int		saveByTime_;
  int		saveByRun_;
  bool		saveAtJobEnd_;
  int		saveReference_;
  int		saveReferenceQMin_;
  int		forceRunNumber_;

  std::string	fileBaseName_;
  DQMStore	*dbe_;
  int		irun_;
  int		ilumi_;
  int		ilumiprev_;
  int		ievent_;
  int		nrun_;
  int		nlumi_;
  int		nevent_;
  timeval	start_;
  timeval	saved_;

  int			 numKeepSavedFiles_;
  std::list<std::string> pastSavedFiles_;
  
  MonitorElement * versCMSSW_ ;
  MonitorElement * versDataset_ ;
  MonitorElement * versTaglist_ ;
  MonitorElement * versGlobaltag_ ;
  MonitorElement * hostName_;          ///Hostname of the local machine
  MonitorElement * processName_;       ///DQM "name" of the job (eg, Hcal or DT)
  MonitorElement * workingDir_;        ///Current working directory of the job
  MonitorElement * processId_;         ///The PID associated with this job
  MonitorElement * isComplete_;
  MonitorElement * fileVersion_;
  

  
};

#endif // DQMSERVICES_COMPONEntS_DQMFILESAVER_H
