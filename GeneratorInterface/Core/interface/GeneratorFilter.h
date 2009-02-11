// -*- C++ -*-
//
//

// class template GeneratorFilter<HAD> provides an EDFilter which uses
// the hadronizer type HAD to generate partons, hadronize them, and
// decay the resulting particles, in the CMS framework.

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"


//#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenInfoProduct.h"

namespace edm
{
  template <class HAD> class GeneratorFilter : public EDFilter
  {
  public:
    typedef HAD Hadronizer;

    // The given ParameterSet will be passed to the contained
    // Hadronizer object.
    explicit GeneratorFilter(ParameterSet const& ps);

    virtual ~GeneratorFilter();

    virtual bool filter(Event& e, EventSetup const& es);
    virtual void beginJob(EventSetup const&);
    virtual void endJob();
    virtual bool beginRun(Run &, EventSetup const&);
    virtual bool endRun(Run &, EventSetup const&);
    virtual bool beginLuminosityBlock(LuminosityBlock &, EventSetup const&);
    virtual bool endLuminosityBlock(LuminosityBlock &, EventSetup const&);
    virtual void respondToOpenInputFile(FileBlock const& fb);
    virtual void respondToCloseInputFile(FileBlock const& fb);
    virtual void respondToOpenOutputFiles(FileBlock const& fb);
    virtual void respondToCloseOutputFiles(FileBlock const& fb);

  private:
    Hadronizer hadronizer_;
  };

  //------------------------------------------------------------------------
  //
  // Implementation

  template <class HAD>
  GeneratorFilter<HAD>::GeneratorFilter(ParameterSet const& ps) :
    EDFilter(),
    hadronizer_(ps)
  {
    // TODO:
    // Put the list of types produced by the filters here.
    // The current design calls for:
    //   * LHEGeneratorInfo
    //   * LHEEvent
    //   * HepMCProduct
    // But I can not find the LHEGeneratorInfo class; it might need to
    // be invented.

    // Commented out because of compilation failures; inability to
    // find headers.
    //
    produces<edm::HepMCProduct>();
    produces<edm::GenInfoProduct, edm::InRun>();
  }

  template <class HAD>
  GeneratorFilter<HAD>::~GeneratorFilter()
  { }

  template <class HAD>
  bool
  GeneratorFilter<HAD>::filter(Event& ev, EventSetup const& /* es */)
  {
    std::auto_ptr<HepMCProduct> bare_product(new HepMCProduct());

    // hadronizer_.generatePartons();
    // hadronizer_.hadronize();
    if ( !hadronizer_.generatePartonsAndHadronize() ) return false;

    // When the external decay driver is added to the system, it
    // should be called here.
    if ( !hadronizer_.decay() ) return false;

    HepMC::GenEvent* event = hadronizer_.getGenEvent();
    if ( !event ) return false;

    bare_product->addHepMCData(event);
    ev.put(bare_product);

    return true;
  }

  template <class HAD>
  void
  GeneratorFilter<HAD>::beginJob(EventSetup const&)
  { 

    if (! hadronizer_.declareStableParticles())
      throw edm::Exception(errors::Configuration)
	<< "Failed to declare stable particles in hadronizer "
	<< hadronizer_.classname()
	<< "\n";
  }
  
  template <class HAD>
  void
  GeneratorFilter<HAD>::endJob()
  { }

  template <class HAD>
  bool
  GeneratorFilter<HAD>::beginRun(Run &, EventSetup const&)
  {
    // Create the LHEGeneratorInfo product describing the run
    // conditions here, and insert it into the Run object.

    if (! hadronizer_.initializeForInternalPartons())
      throw edm::Exception(errors::Configuration) 
	<< "Failed to initialize hadronizer "
	<< hadronizer_.classname()
	<< " for internal parton generation\n";

    return true;
  }

  template <class HAD>
  bool
  GeneratorFilter<HAD>::endRun(Run& r, EventSetup const&)
  {
    // If relevant, record the integrated luminosity for this run
    // here.  To do so, we would need a standard function to invoke on
    // the contained hadronizer that would report the integrated
    // luminosity.

    hadronizer_.statistics();
    
    std::auto_ptr<edm::GenInfoProduct> giproduct(new edm::GenInfoProduct(hadronizer_.getGenInfoProduct()));
    r.put(giproduct);

    return true;
  }

  template <class HAD>
  bool
  GeneratorFilter<HAD>::beginLuminosityBlock(LuminosityBlock &, EventSetup const&)
  {
    return true;
  }

  template <class HAD>
  bool
  GeneratorFilter<HAD>::endLuminosityBlock(LuminosityBlock &, EventSetup const&)
  {
    // If relevant, record the integration luminosity of this
    // luminosity block here.  To do so, we would need a standard
    // function to invoke on the contained hadronizer that would
    // report the integrated luminosity.
    return true;
  }

  template <class HAD>
  void
  GeneratorFilter<HAD>::respondToOpenInputFile(FileBlock const& fb)
  { }

  template <class HAD>
  void
  GeneratorFilter<HAD>::respondToCloseInputFile(FileBlock const& fb)
  { }

  template <class HAD>
  void
  GeneratorFilter<HAD>::respondToOpenOutputFiles(FileBlock const& fb)
  { }

  template <class HAD>
  void
  GeneratorFilter<HAD>::respondToCloseOutputFiles(FileBlock const& fb)
  { }

}
