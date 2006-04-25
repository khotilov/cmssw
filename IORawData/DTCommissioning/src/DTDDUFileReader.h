#ifndef DaqSource_DTDDUFileReader_h
#define DaqSource_DTDDUFileReader_h

/** \class DTDDUFileReader
 *  Read DT ROS8 raw data files
 *
 *  $Date: 2006/04/10 12:15:28 $
 *  $Revision: 1.1 $
 *  \author M. Zanetti - INFN Padova
 */

#include <IORawData/DaqSource/interface/DaqBaseReader.h>
#include <DataFormats/Common/interface/EventID.h>

#include <ostream>
#include <fstream>
#include <boost/cstdint.hpp>

class DTDDUFileReader : public DaqBaseReader {
 public:
  /// Constructor
  DTDDUFileReader(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~DTDDUFileReader();

  /// Generate and fill FED raw data for a full event
  virtual bool fillRawData(edm::EventID& eID,
			   edm::Timestamp& tstamp, 
			   FEDRawDataCollection& data);


  /// check for a 64 bits word to be a DDU header
  bool isHeader(uint64_t word, bool dataTag);

  /// check for a 64 bits word to be a DDU trailer
  bool isTrailer(uint64_t word, bool dataTag);

  /// pre-unpack the data if read via DMA
  //  std::pair<uint64_t,bool> dmaUnpack();
  uint64_t dmaUnpack(bool & isData);


  /// swapping the lsBits with the msBits
  void swap(uint64_t & word);
 

  virtual bool checkEndOfFile();

 private:


  std::ifstream inputFile;

  edm::RunNumber_t runNumber;
  edm::EventNumber_t eventNumber;

  bool readFromDMA;

  static const int dduWordLenght = 8;

};
#endif



