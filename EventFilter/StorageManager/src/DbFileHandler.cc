// $Id: DbFileHandler.cc,v 1.5 2010/04/12 12:05:01 mommsen Exp $
/// @file: DbFileHandler.cc

#include <EventFilter/StorageManager/interface/DbFileHandler.h>
#include <EventFilter/StorageManager/interface/Exception.h>

#include <iomanip>

using namespace stor;
using namespace std;


DbFileHandler::DbFileHandler() :
_runNumber(0)
{}


void DbFileHandler::writeOld(const utils::time_point_t& timestamp, const string& str)
{
  std::ofstream outputFile;
  openFile(outputFile, timestamp);
  outputFile << str.c_str();
  outputFile.close();
}


void DbFileHandler::write(const string& str)
{
  const utils::time_point_t timestamp = utils::getCurrentTime();

  std::ofstream outputFile;
  openFile(outputFile, timestamp);
  addReportHeader(outputFile, timestamp);
  outputFile << str.c_str();
  outputFile << endl;
  outputFile.close();
}


void DbFileHandler::configure(const unsigned int runNumber, const DiskWritingParams& params)
{ 
  _dwParams = params;
  _runNumber = runNumber;

  write("BoR");
}



void DbFileHandler::openFile
(
  std::ofstream& outputFile,
  const utils::time_point_t& timestamp
) const
{
  utils::checkDirectory(_dwParams._dbFilePath);

  ostringstream dbfilename;
  dbfilename << _dwParams._dbFilePath << "/"
             << utils::dateStamp(timestamp)
             << "-" << _dwParams._hostName
             << "-" << _dwParams._smInstanceString
             << ".log";

  outputFile.open( dbfilename.str().c_str(), ios_base::ate | ios_base::out | ios_base::app );
  if (! outputFile.is_open() )
  {
    std::ostringstream msg;
    msg << "Failed to open db log file " << dbfilename.str();
    XCEPT_RAISE(stor::exception::DiskWriting, msg.str());
  }
}


void DbFileHandler::addReportHeader
(
  std::ostream& msg,
  const utils::time_point_t& timestamp
) const
{
  msg << "Timestamp:" << static_cast<int>(timestamp)
    << "\trun:" << _runNumber << "\t";
}



/// emacs configuration
/// Local Variables: -
/// mode: c++ -
/// c-basic-offset: 2 -
/// indent-tabs-mode: nil -
/// End: -
