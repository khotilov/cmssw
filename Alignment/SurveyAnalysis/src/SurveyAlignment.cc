#include <iostream>

#include "Alignment/CommonAlignment/interface/Alignable.h"
#include "Alignment/CommonAlignment/interface/AlignmentParameters.h"
#include "Alignment/SurveyAnalysis/interface/SurveyOutput.h"

#include "Alignment/SurveyAnalysis/interface/SurveyAlignment.h"

using namespace align;

SurveyAlignment::SurveyAlignment(const std::vector<Alignable*>& sensors):
  theSensors(sensors)
{
}

void SurveyAlignment::shiftSensors()
{
  unsigned int nSensor = theSensors.size();

  for (unsigned int i = 0; i < nSensor; ++i)
  {
    Alignable* ali = theSensors[i];

    const AlignableSurface& surf = ali->surface();
    const AlgebraicVector&  pars = ali->alignmentParameters()->parameters();

    EulerAngles angles(3);

    angles(1) = pars[3]; angles(2) = pars[4]; angles(3) = pars[5];

    RotationType rot = surf.toGlobal( toMatrix(angles) );

    rectify(rot); // correct for rounding errors

    ali->move( surf.toGlobal( align::LocalVector(pars[0], pars[1], pars[2]) ) );
    ali->rotateInGlobalFrame(rot);
  }
}

void SurveyAlignment::iterate(unsigned int nIteration,
			      const std::string& fileName)
{
  static const double tolerance = 1e-4; // convergence criteria

  SurveyOutput out(theSensors, fileName);

  for (unsigned int i = 0; i < nIteration; ++i)
  {
    std::cout << "***** Iteration " << i << " *****\n";
    findAlignPars();
    out.write(i);
    shiftSensors();

  // Check convergence

    double parChi2 = 0.;

    unsigned int nSensor = theSensors.size();

    for (unsigned int j = 0; j < nSensor; ++j)
    {
      AlignmentParameters* alignPar = theSensors[j]->alignmentParameters();

      const AlgebraicVector&    par = alignPar->parameters();
      const AlgebraicSymMatrix& cov = alignPar->covariance();

      int dummy;

      parChi2 += cov.inverse(dummy).similarity(par);
    }

    parChi2 /= static_cast<double>(nSensor);
    std::cout << "chi2 = " << std::sqrt(parChi2) << std::endl;
    if (parChi2 < tolerance) break; // converges, so exit loop
  }
}
