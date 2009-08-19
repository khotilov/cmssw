#ifndef DQM_SiPixelHistoricInfoClient_SiPixelHistoryDQMService_H
#define DQM_SiPixelHistoricInfoClient_SiPixelHistoryDQMService_H

#include "DQMServices/Diagnostic/interface/DQMHistoryServiceBase.h"
#include "DQM/SiPixelHistoricInfoClient/interface/SiPixelSummary.h"


/**
  @author Dean Andrew Hidas <dhidas@cern.ch>
 */

class SiPixelHistoryDQMService : public DQMHistoryServiceBase {
  public:

    explicit SiPixelHistoryDQMService(const edm::ParameterSet&,const edm::ActivityRegistry&);
    ~SiPixelHistoryDQMService();


  private:
    //Methods to be specified by each subdet
    uint32_t returnDetComponent(const MonitorElement* ME);
    bool setDBLabelsForUser  (std::string& keyName, std::vector<std::string>& userDBContent);
    bool setDBValuesForUser(std::vector<MonitorElement*>::const_iterator iterMes, HDQMSummary::InputVector& values  );
    //std::string userTranslator (int);

    edm::ParameterSet iConfig_;
};

#endif //DQM_SiPixelHistoricInfoClient_SiPixelHistoryDQMService_H

