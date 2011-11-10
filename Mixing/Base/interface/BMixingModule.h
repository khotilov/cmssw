#ifndef BMixingModule_h
#define BMixingModule_h

/** \class BMixingModule
 *
 * BMixingModule is the EDProducer subclass 
 * which fills the CrossingFrame object
 * It is the baseclass for all modules mnixing events 
 *
 * \author Ursula Berthon, LLR Palaiseau, Bill Tanenbaum
 *
 * \version   1st Version June 2005
 * \version   2nd Version Sep 2005

 *
 ************************************************************/

#include "boost/shared_ptr.hpp"

#include "FWCore/Framework/interface/EventPrincipal.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "Mixing/Base/interface/PileUp.h"


namespace edm {
  class BMixingModule : public edm::EDProducer {
    public:
      /** standard constructor*/
      explicit BMixingModule(const edm::ParameterSet& ps);

      /**Default destructor*/
      virtual ~BMixingModule();

      /**Cumulates the pileup events onto this event*/
      virtual void produce(edm::Event& e1, const edm::EventSetup& c);

      // Should 'averageNumber' return 0 or 1 if there is no mixing? It is the average number of
      // *crossings*, including the hard scatter, or the average number of overlapping events?
      // We have guessed 'overlapping events'.
      double averageNumber() const {return inputSources_[0] ? inputSources_[0]->averageNumber() : 0.0; }
      // Should 'poisson' return 0 or 1 if there is no mixing? See also averageNumber above.
      bool poisson() const {return inputSources_[0] ? inputSources_[0]->poisson() : 0.0 ;}

      virtual void createnewEDProduct() {std::cout << "BMixingModule::createnewEDProduct must be overwritten!" << std::endl;}
      virtual void checkSignal(const edm::Event &e) {std::cout << "BMixingModule::checkSignal must be overwritten!" << std::endl;}
      virtual void addSignals(const edm::Event &e,const edm::EventSetup& c) {;}
      virtual void addPileups(const int bcr, EventPrincipal *ep, unsigned int eventId,unsigned int worker, const edm::EventSetup& c) {;}
      virtual void setBcrOffset () {std::cout << "BMixingModule::setBcrOffset must be overwritten!" << std::endl;} //FIXME: LogWarning
      virtual void setSourceOffset (const unsigned int s) {std::cout << "BMixingModule::setSourceOffset must be overwritten!" << std::endl;}
      virtual void put(edm::Event &e,const edm::EventSetup& c) {;}
      virtual void doPileUp(edm::Event &e, const edm::EventSetup& c) {std::cout << "BMixingModule::doPileUp must be overwritten!" << std::endl;}
      virtual void setEventStartInfo(const unsigned int s) {;} //to be set in CF
      virtual void getEventStartInfo(edm::Event & e,const unsigned int source) {;} //to be set locally

  protected:
      void dropUnwantedBranches(std::vector<std::string> const& wantedBranches);
      virtual void endJob();
      int bunchSpace_;
      static int vertexoffset;
      bool checktof_;
      int const minBunch_;
      int const maxBunch_;
      bool const mixProdStep1_;	       	
      bool const mixProdStep2_;
      	
      bool playback_;
      const static unsigned int maxNbSources_;
      std::vector<std::string> sourceNames_;
      bool doit_[4];//FIXME
      std::vector< float > TrueNumInteractions_;

      unsigned int eventId_;

      // input, cosmics, beamhalo_plus, beamhalo_minus
      std::vector<boost::shared_ptr<PileUp> > inputSources_;

  };

}//edm

#endif
