#ifndef ODLTSCONFIG_H
#define ODLTSCONFIG_H

#include <map>
#include <stdexcept>

#include "OnlineDB/EcalCondDB/interface/IODConfig.h"

class ODLTSConfig : public IODConfig {
 public:
  friend class EcalCondDBInterface;
  ODLTSConfig();
  ~ODLTSConfig();

  // User data methods
  inline std::string getTable() { return "ECAL_LTS_CONFIGURATION"; }

  inline void setId(int id) { m_ID = id; }
  inline int getId() const { return m_ID; }

  inline void setTriggerType(std::string x) { m_trg_type = x; }
  inline std::string  getTriggerType() const { return m_trg_type; }

  inline void setNumberOfEvents(int x) { m_num = x; }
  inline int getNumberOfEvents() const { return m_num; }

  inline void setRate(int x) { m_rate = x; }
  inline int getRate() const { return m_rate; }

  inline void setTrigLocL1Delay(int x) { m_delay = x; }
  inline int getTrigLocL1Delay() const { return m_delay; }
  void setParameters(std::map<string,string> my_keys_map);

 private:
  int fetchNextId() throw(std::runtime_error);
  void prepareWrite()  throw(std::runtime_error);
  void writeDB()       throw(std::runtime_error);
  void clear();
  void fetchData(ODLTSConfig * result)     throw(std::runtime_error);
  int fetchID()  throw(std::runtime_error);

  // User data
  int m_ID;
  std::string m_trg_type;
  int m_num;
  int m_rate;
  int m_delay;
  
};

#endif
