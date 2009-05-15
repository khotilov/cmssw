#include <memory>

#include <HepMC/GenEvent.h>
#include <HepMC/IO_BaseClass.h>

#include <ThePEG/Repository/Repository.h>
#include <ThePEG/EventRecord/Event.h>
#include <ThePEG/Config/ThePEG.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "GeneratorInterface/Core/interface/BaseHadronizer.h"
#include "GeneratorInterface/Core/interface/GeneratorFilter.h"
#include "GeneratorInterface/Core/interface/HadronizerFilter.h"

#include "GeneratorInterface/ThePEGInterface/interface/ThePEGInterface.h"

class ThePEGHadronizer : public ThePEGInterface, public gen::BaseHadronizer {
    public:
	ThePEGHadronizer(const edm::ParameterSet &params);
	virtual ~ThePEGHadronizer();

	bool initializeForInternalPartons();
//	bool initializeForExternalPartons();
	bool declareStableParticles(const std::vector<int> &pdgIds);

	void statistics();

	bool generatePartonsAndHadronize();
//	bool hadronize();
	bool decay();
	bool residualDecay();
	void finalizeEvent();

	const char *classname() const { return "ThePEGHadronizer"; }

    private:
	unsigned int			eventsToPrint;

	ThePEG::EventPtr		thepegEvent;
};

ThePEGHadronizer::ThePEGHadronizer(const edm::ParameterSet &pset) :
	ThePEGInterface(pset),
	BaseHadronizer(pset),
	eventsToPrint(pset.getUntrackedParameter<unsigned int>("eventsToPrint", 0))
{  
	initRepository(pset);

}

ThePEGHadronizer::~ThePEGHadronizer()
{
}

bool ThePEGHadronizer::initializeForInternalPartons()
{
	initGenerator();
	return true;
}

bool ThePEGHadronizer::declareStableParticles(const std::vector<int> &pdgIds)
{
	return false;
}

void ThePEGHadronizer::statistics()
{
	runInfo().setInternalXSec(GenRunInfoProduct::XSec(
				eg_->integratedXSec() / ThePEG::picobarn));
}

bool ThePEGHadronizer::generatePartonsAndHadronize()
{
	edm::LogInfo("Generator|ThePEGHadronizer") << "Start production";

	flushRandomNumberGenerator();

	thepegEvent = eg_->shoot();
	if (!thepegEvent) {
		edm::LogWarning("Generator|ThePEGHadronizer") << "thepegEvent not initialized";
		return false;
	}

	event() = convert(thepegEvent);
	if (!event().get()) {
		edm::LogWarning("Generator|ThePEGHadronizer") << "genEvent not initialized";
		return false;
	}

	return true;
}

void ThePEGHadronizer::finalizeEvent()
{
	HepMC::PdfInfo pdf;
	clearAuxiliary(event().get(), &pdf);
	fillAuxiliary(event().get(), &pdf, thepegEvent);
	if (usePthatEventScale)
		setPthatEventScale(event().get(), thepegEvent);
	event()->set_pdf_info(pdf);

	if (eventsToPrint) {
		eventsToPrint--;
		event()->print();
	}

	if (iobc_.get())
		iobc_->write_event(event().get());

	edm::LogInfo("Generator|ThePEGHadronizer") << "Event produced";
}

bool ThePEGHadronizer::decay()
{
	return true;
}

bool ThePEGHadronizer::residualDecay()
{
	return true;
}

typedef edm::GeneratorFilter<ThePEGHadronizer> ThePEGGeneratorFilter;
DEFINE_FWK_MODULE(ThePEGGeneratorFilter);

#if 0
typedef edm::HadronizerFilter<ThePEGHadronizer> ThePEGHadronizerFilter;
DEFINE_FWK_MODULE(ThePEGHadronizerFilter);
#endif
