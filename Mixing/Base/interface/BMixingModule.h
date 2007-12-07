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

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "Mixing/Base/interface/PileUp.h"


namespace edm {
  class BMixingModule : public edm::EDProducer {
    public:
      typedef PileUp::EventPrincipalVector EventPrincipalVector;

      /** standard constructor*/
      explicit BMixingModule(const edm::ParameterSet& ps);

      /**Default destructor*/
      virtual ~BMixingModule();

      /**Cumulates the pileup events onto this event*/
      virtual void produce(edm::Event& e1, const edm::EventSetup& c);

      // Should 'averageNumber' return 0 or 1 if there is no mixing? It is the average number of
      // *crossings*, including the hard scatter, or the average number of overlapping events?
      // We have guessed 'overlapping events'.
      double averageNumber() const {return input_ ? input_->averageNumber() : 0.0;}
      // Should 'poisson' return 0 or 1 if there is no mixing? See also averageNumber above.
      //bool poisson() const {return input_.poisson();}
      bool poisson() const {return input_ ? input_->poisson() : 0.0 ;}

      virtual void createnewEDProduct() {std::cout << "BMixingModule::createnewEDProduct must be overwritten!" << std::endl;}
      void merge(const int bcr, const EventPrincipalVector& vec);
      virtual void addSignals(const edm::Event &e) {;}
      virtual void addPileups(const int bcr, edm::Event*, unsigned int eventId) {;}
      virtual void setBcrOffset () {std::cout << "BMixingModule::setBcrOffset must be overwritten!" << std::endl;} //FIXME: LogWarning
      virtual void setSourceOffset (const unsigned int s) {std::cout << "BMixingModule::setSourceOffset must be overwritten!" << std::endl;}
      virtual void put(edm::Event &e) {;}
      virtual void setEventStartInfo(edm::EventID &id,int fileNr,const unsigned int source) = 0; //to be set in CF
      virtual void getEventStartInfo(edm::Event & e,const unsigned int source) = 0; //to be set locally

  protected:
      int bunchSpace_;
      static int vertexoffset;
      bool checktof_;
      int const minBunch_;
      int const maxBunch_;

      // for playback
      edm::EventID id_;
      int fileNr_;
      bool playback_;

    private:

      boost::shared_ptr<PileUp> input_;
      boost::shared_ptr<PileUp> cosmics_;
      boost::shared_ptr<PileUp> beamHalo_p_;
      boost::shared_ptr<PileUp> beamHalo_m_;
      ModuleDescription md_;

      unsigned int eventId_;

      const static unsigned int maxNbSources;
  };

}//edm

#endif
