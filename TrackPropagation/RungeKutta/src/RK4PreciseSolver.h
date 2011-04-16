#ifndef RK4PreciseSolver_H
#define RK4PreciseSolver_H

#include "RKSolver.h"
#include "Utilities/UI/interface/SimpleConfigurable.h"

template <typename T, int N>
class RK4PreciseSolver : public RKSolver<T,N> {
public:

    typedef RKSolver<T,N>                       Base;
    typedef typename Base::Scalar               Scalar;
    typedef typename Base::Vector               Vector;

    virtual Vector operator()( Scalar startPar, const Vector& startState,
			       Scalar step, const RKDerivative<T,N>& deriv,
			       const RKDistance<T,N>& dist,
			       Scalar eps);

    std::pair< Vector, T> 
    stepWithAccuracy( Scalar startPar, const Vector& startState,
		      const RKDerivative<T,N>& deriv,
		      const RKDistance<T,N>& dist, Scalar step);

protected:

    bool verbose() const {
	static bool verb = SimpleConfigurable<bool>(false,"RKSolver:verbose").value();
	return verb;
    }

};

#include "TrackPropagation/RungeKutta/src/RK4PreciseSolver.icc"

#endif
