// -*- C++ -*-
//
// Package:     MVAComputer
// Class  :     ProcTMVA
// 

// Implementation:
//     TMVA wrapper, needs n non-optional, non-multiple input variables
//     and outputs one result variable. All TMVA algorithms can be used,
//     calibration data is passed via stream and extracted from a zipped
//     buffer.
//
// Author:      Christophe Saout
// Created:     Sat Apr 24 15:18 CEST 2007
// $Id: ProcTMVA.cc,v 1.1 2009/05/11 16:54:21 saout Exp $
//

#include <sstream>
#include <string>
#include <vector>
#include <memory>

// ROOT version magic to support TMVA interface changes in newer ROOT
#include <RVersion.h>

#include <TMVA/DataSet.h>
#include <TMVA/Types.h>
#include <TMVA/MethodBase.h>
#include <TMVA/Methods.h>

#include "PhysicsTools/MVAComputer/interface/memstream.h"
#include "PhysicsTools/MVAComputer/interface/zstream.h"

#include "PhysicsTools/MVAComputer/interface/VarProcessor.h"
#include "PhysicsTools/MVAComputer/interface/Calibration.h"

using namespace PhysicsTools;

namespace { // anonymous

class ProcTMVA : public VarProcessor {
    public:
	typedef VarProcessor::Registry::Registry<ProcTMVA,
					Calibration::ProcExternal> Registry;

	ProcTMVA(const char *name,
	         const Calibration::ProcExternal *calib,
	         const MVAComputer *computer);
	virtual ~ProcTMVA() {}

	virtual void configure(ConfIterator iter, unsigned int n);
	virtual void eval(ValueIterator iter, unsigned int n) const;

    private:
	mutable TMVA::DataSet		data;
	std::auto_ptr<TMVA::MethodBase>	method;
	unsigned int			nVars;
};

static ProcTMVA::Registry registry("ProcTMVA");

#define SWITCH_METHOD(name)					\
	case (TMVA::Types::k##name):				\
		return new TMVA::Method##name(*data, "");

static TMVA::MethodBase *methodInst(TMVA::DataSet *data, TMVA::Types::EMVA type)
{
	switch(type) {
		SWITCH_METHOD(Cuts)
		SWITCH_METHOD(SeedDistance)
		SWITCH_METHOD(Likelihood)
		SWITCH_METHOD(PDERS)
		SWITCH_METHOD(HMatrix)
		SWITCH_METHOD(Fisher)
		SWITCH_METHOD(CFMlpANN)
		SWITCH_METHOD(TMlpANN)
		SWITCH_METHOD(BDT)
		SWITCH_METHOD(RuleFit)
		SWITCH_METHOD(SVM)
		SWITCH_METHOD(MLP)
		SWITCH_METHOD(BayesClassifier)
		SWITCH_METHOD(FDA)
		SWITCH_METHOD(Committee)
	    default:
		return 0;
	}
}

#undef SWITCH_METHOD

ProcTMVA::ProcTMVA(const char *name,
                   const Calibration::ProcExternal *calib,
                   const MVAComputer *computer) :
	VarProcessor(name, calib, computer)
{
	ext::imemstream is(
		reinterpret_cast<const char*>(&calib->store.front()),
	        calib->store.size());
	ext::izstream izs(&is);

	std::string methodName;
	std::getline(izs, methodName);
	std::string tmp;
	std::getline(izs, tmp);
	std::istringstream iss(tmp);
	iss >> nVars;
	for(unsigned int i = 0; i < nVars; i++) {
		std::getline(izs, tmp);
		data.AddVariable(tmp.c_str());
	}

	TMVA::Types::EMVA methodType =
			TMVA::Types::Instance().GetMethodType(methodName);

	method = std::auto_ptr<TMVA::MethodBase>(
					methodInst(&data, methodType));

	method->ReadStateFromStream(izs);
}

void ProcTMVA::configure(ConfIterator iter, unsigned int n)
{
	if (n != nVars)
		return;

	for(unsigned int i = 0; i < n; i++)
		iter++(Variable::FLAG_NONE);

	iter << Variable::FLAG_NONE;
}

void ProcTMVA::eval(ValueIterator iter, unsigned int n) const
{
	for(unsigned int i = 0; i < n; i++)
		data.GetEvent().SetVal(i, *iter++);

	method->GetVarTransform().GetEventRaw().CopyVarValues(data.GetEvent());
	method->GetVarTransform().ApplyTransformation(TMVA::Types::kSignal);
	iter(method->GetMvaValue());
}

} // anonymous namespace

MVA_COMPUTER_DEFINE_PLUGIN(ProcTMVA);
