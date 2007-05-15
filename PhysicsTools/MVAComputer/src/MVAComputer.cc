#include <iostream>
#include <algorithm>
#include <cstddef>
#include <vector>
#include <set>

#include "FWCore/Utilities/interface/Exception.h"

#include "PhysicsTools/MVAComputer/interface/MVAComputer.h"
#include "PhysicsTools/MVAComputer/interface/VarProcessor.h"
#include "PhysicsTools/MVAComputer/interface/Calibration.h"
#include "PhysicsTools/MVAComputer/interface/Variable.h"
#include "PhysicsTools/MVAComputer/interface/AtomicId.h"

// #define DEBUG_EVAL

#ifdef DEBUG_EVAL
#	include <Reflex/Tools.h>
#endif

namespace PhysicsTools {

MVAComputer::MVAComputer(
			const Calibration::MVAComputer *calib) :
	nVars(0), output(0)
{
	setup(calib);
}

void
MVAComputer::setup(const Calibration::MVAComputer *calib)
{
	nVars = calib->inputSet.size();
	output = calib->output;

	VarProcessor::Config_t config(nVars,
				VarProcessor::Config(Variable::FLAG_ALL, 1));
	std::vector<Calibration::VarProcessor*> processors =
							calib->getProcessors();
	for(std::vector<Calibration::VarProcessor*>::const_iterator iter =
							processors.begin();
	    iter != processors.end(); ++iter) {
		std::string name = (*iter)->getInstanceName();
		VarProcessor *processor =
			VarProcessor::create(name.c_str(), *iter, this);
		if (!processor)
			throw cms::Exception("UnknownProcessor")
				<< name << " could not be instantiated."
				<< std::endl;

		VarProcessor::Config_t::iterator::difference_type pos =
						config.end() - config.begin();
		processor->configure(config);
		unsigned int nOutput = (config.end() - config.begin()) - pos;
		if (!nOutput)
			throw cms::Exception("InvalidProcessor")
				<< name << " rejected input variable "
				"configuration" << std::endl;

		varProcessors.push_back(Processor(processor, nOutput));
	}

	for(VarProcessor::Config_t::iterator iter = config.begin() + nVars;
	    iter != config.end(); iter++) {
		if (!(iter->mask & Variable::FLAG_MULTIPLE))
			continue;

		VarProcessor::Config *origin = &config[iter->origin];
		iter->mask = (Variable::Flags)(iter->mask & origin->mask);
		if (iter->origin >= nVars)
			iter->origin = origin->origin;

		if (iter->mask & Variable::FLAG_MULTIPLE)
			config[iter->origin].origin++;
	}

	nVars = config.size();

	if (output >= nVars)
			// FIXME || config[output].mask != Variable::FLAG_NONE)
		throw cms::Exception("InvalidOutput")
			<< "Output variable at index " << output
			<< " invalid." << std::endl;

	std::set<InputVar> variables;
	unsigned int i = 0;
	for(std::vector<Calibration::Variable>::const_iterator iter =
			calib->inputSet.begin(); iter != calib->inputSet.end();
	    ++iter, i++) {
		InputVar var;
		var.var = Variable(iter->name, config[i].mask);
		var.index = i;
		var.multiplicity = config[i].origin;
		variables.insert(var);
	}

	inputVariables.resize(i);
	std::copy(variables.begin(), variables.end(),
	          inputVariables.begin());
}

MVAComputer::~MVAComputer()
{
}

unsigned int MVAComputer::getVariableId(AtomicId name) const
{
	std::vector<InputVar>::const_iterator pos =
		std::lower_bound(inputVariables.begin(), inputVariables.end(),
		                 name);

	if (pos == inputVariables.end() || pos->var.getName() != name)
		throw cms::Exception("InvalidVariable")
			<< "Input variable " << (const char*)name
			<< " not found."  << std::endl;

	return pos->index;
}

void MVAComputer::eval(double *values, int *conf,
                                 unsigned int n) const
{
	double *output = values + n;
	int *outConf = conf + inputVariables.size();

#ifdef DEBUG_EVAL
	std::cout << "Input" << std::endl;
	double *v = values;
	for(int *o = conf; o < outConf; o++) {
		std::cout << "\tVar " << (o - conf) << std::endl;
		for(int i = o[0]; i < o[1]; i++)
			std::cout << "\t\t" << *v++ << std::endl;
	}
#endif

	for(std::vector<Processor>::const_iterator iter = varProcessors.begin();
	    iter != varProcessors.end(); iter++) {
#ifdef DEBUG_EVAL
		std::cout << ROOT::Reflex::Tools::Demangle(
				typeid(*iter->processor)) << std::endl;
#endif

		iter->processor->eval(values, conf, output, outConf);

#ifdef DEBUG_EVAL
		for(unsigned int i = 0; i < iter->nOutput; i++, outConf++) {
			std::cout << "\tVar " << (outConf - conf) << std::endl;
			for(int j = outConf[0]; j < outConf[1]; j++)
				std::cout << "\t\t" << *output++ << std::endl;
		}
#else
		int orig = *outConf;
		outConf += iter->nOutput;
		output += *outConf - orig;
#endif
	}
}

} // namespace PhysicsTools
