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
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
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

      //int minBunch() const {return input_.minBunch();}
      int minBunch() const {return input_ ? input_->minBunch() : 0 ;}
      //int maxBunch() const {return input_.maxBunch();}
      int maxBunch() const {return input_ ? input_->maxBunch() : 0 ;}
      //double averageNumber() const {return input_.averageNumber();}
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
      virtual void put(edm::Event &e) {;}

    protected:
      int bunchSpace_;
      static int trackoffset;
      static int vertexoffset;
      bool checktof_;

    private:
      boost::shared_ptr<PileUp> input_;
      ModuleDescription md_;
      unsigned int eventId_;
    };
}//edm

#endif
