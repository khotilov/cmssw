#ifndef EventFilter_DTDigiToRawModule_h
#define EventFilter_DTDigiToRawModule_h


#include <FWCore/Framework/interface/EDProducer.h>



class DTDigiToRaw;

class DTDigiToRawModule : public edm::EDProducer {
public:
  /// Constructor
  DTDigiToRawModule(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~DTDigiToRawModule();

  // Operations
  virtual void produce( edm::Event&, const edm::EventSetup& );

private:
  DTDigiToRaw * packer;
  
  int dduID;
  bool debug;
  bool digibyType;
  std::string digicoll;

};
#endif

