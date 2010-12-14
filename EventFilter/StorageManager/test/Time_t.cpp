#include <iostream>
#include <iomanip>
#include <ctime>

#include "boost/date_time/posix_time/posix_time.hpp"

#include "Utilities/Testing/interface/CppUnit_testdriver.icpp"
#include "cppunit/extensions/HelperMacros.h"

#include "EventFilter/StorageManager/interface/Utils.h"

using namespace stor;

class testTime :  public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(testTime);
  CPPUNIT_TEST(testDifference);
  CPPUNIT_TEST(testTimeStamp);

  CPPUNIT_TEST_SUITE_END();

public:

  void testDifference();
  void testTimeStamp();

};


void testTime::testDifference()
{
  utils::time_point_t utilsTime = utils::getCurrentTime();
  struct tm utils_tm = boost::posix_time::to_tm(utilsTime);
  time_t raw = time(0);
  struct tm* raw_tm = gmtime(&raw);
  time_t rawTime = mktime(raw_tm);
  double timeDiff = difftime(rawTime,mktime(&utils_tm));

  std::ostringstream msg;
  msg << std::setiosflags(std::ios::fixed)
    << "Difference: rawtime " << rawTime
    << "\t utilstime: " << utilsTime
    << std::resetiosflags(std::ios::fixed)
    << "\t difference: " << timeDiff;
  CPPUNIT_ASSERT_MESSAGE(msg.str(), timeDiff<1);
}

void testTime::testTimeStamp()
{
  char dateString[80]; 
  time_t rawTime = time(0);
  struct tm * timeinfo = localtime( &rawTime);
  strftime(dateString,80,"%d/%m/%Y:%H/%M/%S",timeinfo);
  std::string utilsDateString = utils::timeStamp(utils::getCurrentTime());

  std::ostringstream msg;
  msg << "TimeStamp: system: " << dateString
    << "\t utils: " << utilsDateString;
  CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str(), 0, utilsDateString.compare(dateString));
}


// This macro writes the 'main' for this test.
CPPUNIT_TEST_SUITE_REGISTRATION(testTime);


/// emacs configuration
/// Local Variables: -
/// mode: c++ -
/// c-basic-offset: 2 -
/// indent-tabs-mode: nil -
/// End: -
