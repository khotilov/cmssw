// -*- C++ -*-
//
// Package:     MVAComputer
// Class  :     VarProcessor
// 

// Implementation:
//     Base class for variable processors. Basically only passes calls
//     through to virtual methods in the actual implementation daughter class.
//
// Author:      Christophe Saout
// Created:     Sat Apr 24 15:18 CEST 2007
// $Id: VarProcessor.cc,v 1.5 2009/05/11 16:01:17 saout Exp $
//

#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/PluginManager/interface/PluginManager.h"
#include "FWCore/PluginManager/interface/PluginFactory.h"

#include "PhysicsTools/MVAComputer/interface/VarProcessor.h"
#include "PhysicsTools/MVAComputer/interface/Calibration.h"
#include "PhysicsTools/MVAComputer/interface/BitSet.h"

// #define DEBUG_DERIV

#ifdef DEBUG_DERIV
#	include <Reflex/Tools.h>
#endif

EDM_REGISTER_PLUGINFACTORY(PhysicsTools::VarProcessor::PluginFactory,
                           "PhysicsToolsMVAComputer");

namespace PhysicsTools {

VarProcessor::VarProcessor(const char *name,
                           const Calibration::VarProcessor *calib,
                           const MVAComputer *computer) :
	computer(computer),
	inputVars(Calibration::convert(calib->inputVars)),
	nInputVars(inputVars.bits())
{
}

VarProcessor::~VarProcessor()
{
	inputVars = BitSet(0);
	nInputVars = 0;
}

void VarProcessor::configure(ConfigCtx &config)
{
	ConfigCtx::size_type pos = config.size();
	if (pos != inputVars.size())
		return;

	ConfIterator iter(inputVars.iter(), config);
	configure(iter, nInputVars);

	VarProcessor *loop = config.loop ? config.loop : this;
	ConfigCtx::Context *ctx =
		loop->configureLoop(config.ctx, config.begin(),
		                    config.begin() + pos, config.end());

	if (ctx != config.ctx) {
		delete config.ctx;
		config.ctx = ctx;
	}

	if (config.loop && !ctx)
		config.loop = 0;
	else if (!config.loop && ctx)
		config.loop = this;
}

VarProcessor::ConfigCtx::ConfigCtx(std::vector<Variable::Flags> flags) :
	loop(0), ctx(0)
{
	for(std::vector<Variable::Flags>::const_iterator iter = flags.begin();
	    iter != flags.end(); ++iter)
		configs.push_back(Config(*iter, 1));
}

VarProcessor::ConfigCtx::Context *
VarProcessor::configureLoop(ConfigCtx::Context *ctx, ConfigCtx::iterator begin,
                            ConfigCtx::iterator cur, ConfigCtx::iterator end)
{
	return 0;
}

template<>
VarProcessor *ProcessRegistry<VarProcessor, Calibration::VarProcessor,
                              const MVAComputer>::Factory::create(
	        const char *name, const Calibration::VarProcessor *calib,
		const MVAComputer *parent)
{
	VarProcessor *result = ProcessRegistry::create(name, calib, parent);
	if (!result) {
		// try to load the shared library and retry
		try {
			delete VarProcessor::PluginFactory::get()->create(
					std::string("VarProcessor/") + name);
			result = ProcessRegistry::create(name, calib, parent);
		} catch(const cms::Exception &e) {
			// caller will have to deal with the null pointer
			// in principle this will just give a slightly more
			// descriptive error message (and will rethrow anyhow)
		}
	}
	return result;
}

void VarProcessor::deriv(double *input, int *conf, double *output,
                         int *outConf, int *loop, unsigned int offset,
                         unsigned int in, unsigned int out_,
                         std::vector<double> &deriv) const
{
	ValueIterator iter(inputVars.iter(), input, conf,
	                   output, outConf, loop, offset);

	eval(iter, nInputVars);

	std::vector<double> matrix = this->deriv(iter, nInputVars);

	unsigned int out = outConf[out_] - outConf[0];
	unsigned int size = 0;
	while(iter)
		size += (iter++).size();

#ifdef DEBUG_DERIV
	if (!matrix.empty()) {
		std::cout << "---------------- "
		          << ROOT::Reflex::Tools::Demangle(typeid(*this))
		          << std::endl;
	for(unsigned int i = 0; i < out; i++) {
		for(unsigned int j = 0; j < size; j++)
			std::cout << matrix[i*size+j] << "\t";
		std::cout << std::endl;
	}
	std::cout << "----------------" << std::endl;
}
#endif

	unsigned int sz = deriv.size();
	deriv.resize(sz + out * in);

#ifdef DEBUG_DERIV
	std::cout << "======= in = " << in << ", size = " << size
	          << ", out = " << out << std::endl;
#endif
	if (matrix.empty())
		return;
	else if (matrix.size() != out * size)
		throw cms::Exception("VarProcessor")
			<< "Derivative matrix implausible size in "
			<< typeid(*this).name() << "."
			<< std::endl;

	double *begin = &deriv.front() + sz;
	double *end = begin + out * in;

	double *m0 = &matrix.front();
	BitSet::Iterator cur = inputVars.iter();
	for(unsigned int i = 0; i < nInputVars; i++, ++cur) {
#ifdef DEBUG_DERIV
		std::cout << " inputvar " << i << std::endl;
#endif
		int *curConf = conf + cur();
		unsigned int pos = *curConf;
#ifdef DEBUG_DERIV
		std::cout << " -> cur = " << cur() << ", pos = "
		          << pos << std::endl;
#endif
		if (loop && curConf >= loop) {
			pos += offset;
			loop = 0;
		}

		unsigned int n = loop ? (curConf[1] - curConf[0]) : 1;
		for(unsigned int j = 0; j < n; m0++, j++, pos++) {
#ifdef DEBUG_DERIV
			std::cout << "  multip " << j << std::endl;
#endif
			double *p = begin;
			if (pos >= in) {
#ifdef DEBUG_DERIV
				std::cout << "   deriv " << (pos - in)
				          << std::endl;
#endif
				const double *q = &deriv.front() +
				                  (pos - in) * in;
				for(const double *m = m0; p < end;
				    m += size, q -= size)
					for(unsigned int k = 0; k < in; k++)
						*p++ += *m * *q++;
			} else {
#ifdef DEBUG_DERIV
				std::cout << "   in " << pos << std::endl;
#endif
				for(const double *m = m0; p < end;
				    m += size, p += in)
					p[pos] += *m;
			}
		}
	}

#ifdef DEBUG_DERIV
	std::cout << "================" << std::endl;
	for(const double *p = &deriv.front();
	    p != &deriv.front() + deriv.size();) {
		for(unsigned int j = 0; j < in; j++)
			std::cout << *p++ << "\t";
		std::cout << std::endl;
	}
	std::cout << "================" << std::endl;
#endif
}

} // namespace PhysicsTools
