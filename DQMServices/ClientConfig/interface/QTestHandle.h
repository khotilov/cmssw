#ifndef QTestHandle_H
#define QTestHandle_H

/** \class QTestHandle
 * *
 *  Handles quality tests (configuring, attaching to ME's, 
 *
 *  $Date: 2007/09/06 13:21:57 $
 *  $Revision: 1.4 $
 *  \author Ilaria Segoni
  */
  
#include<string>
#include<vector>
#include<map>

class DaqMonitorBEInterface;
class QTestConfigurationParser;
class QTestConfigure;
class QTestStatusChecker;

class QTestHandle{
  public:
	///Creator
	QTestHandle();
	///Destructor
	~QTestHandle();
	///Parses Config File and configures the quality tests
	bool configureTests(std::string configFile, DaqMonitorBEInterface * bei);
	///Attaches the quality tests to the MonitorElement
	void attachTests(DaqMonitorBEInterface * bei);
	///Checks global status of Quality Tests
	std::pair<std::string,std::string> checkGlobalQTStatus(DaqMonitorBEInterface * bei) const;
	///Checks alarms for single MonitorElements
	std::map< std::string, std::vector<std::string> > checkDetailedQTStatus(DaqMonitorBEInterface * bei) const;
  
  private:

	QTestConfigurationParser * qtParser;
	QTestConfigure * qtConfigurer;
	QTestStatusChecker * qtChecker;
	bool testsConfigured;


};


#endif
