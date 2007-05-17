// -*- C++ -*-
//
// Package:     Discriminator
// Class  :     ProcLinear
// 

// Implementation:
//     Variable processor to compute a simple linear discriminant using
//     coefficients for each input variable.
//
// Author:      Christophe Saout
// Created:     Sat Apr 24 15:18 CEST 2007
// $Id: ProcLinear.cc,v 1.1 2007/05/07 18:30:55 saout Exp $
//

#include <vector>

#include "PhysicsTools/MVAComputer/interface/VarProcessor.h"
#include "PhysicsTools/MVAComputer/interface/Calibration.h"

using namespace PhysicsTools;

namespace { // anonymous

class ProcLinear : public VarProcessor {
    public:
	typedef VarProcessor::Registry::Registry<ProcLinear,
					Calibration::ProcLinear> Registry;

	ProcLinear(const char *name,
	           const Calibration::ProcLinear *calib,
	           const MVAComputer *computer);
	virtual ~ProcLinear() {}

	virtual void configure(ConfIterator iter, unsigned int n);
	virtual void eval(ValueIterator iter, unsigned int n) const;

    private:
	std::vector<double>	coeffs;
	double			offset;
};

static ProcLinear::Registry registry("ProcLinear");

ProcLinear::ProcLinear(const char *name,
                       const Calibration::ProcLinear *calib,
                       const MVAComputer *computer) :
	VarProcessor(name, calib, computer),
	coeffs(calib->coeffs),
	offset(calib->offset)
{
}

void ProcLinear::configure(ConfIterator iter, unsigned int n)
{
	while(iter)
		iter++(Variable::FLAG_NONE);

	iter << Variable::FLAG_NONE;
}

void ProcLinear::eval(ValueIterator iter, unsigned int n) const
{
	double sum = offset;

	for(std::vector<double>::const_iterator coeff = coeffs.begin();
	    coeff != coeffs.end(); coeff++, ++iter)
		sum += *coeff * *iter;

	iter(sum);
}

} // anonymous namespace
