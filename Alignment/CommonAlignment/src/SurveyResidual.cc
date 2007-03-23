#include "Alignment/CommonAlignment/interface/Alignable.h"
#include "Alignment/CommonAlignment/interface/SurveyDet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "Alignment/CommonAlignment/interface/SurveyResidual.h"

using namespace align;

SurveyResidual::SurveyResidual(const Alignable& ali,
			       ObjectId type,
			       bool bias):
  theSurface( ali.surface() )
{
// Find mother matching given type

  const Alignable* dau = &ali; // start finding from sensor

  while (dau)
  {
    const Alignable* mom = dau->mother();

    if (!mom)
    {
      throw cms::Exception("ConfigError")
	<< "Alignable (id = " << ali.geomDetId().rawId()
	<< ") does not belong to a composite of type " << type
	<< ". Abort!"
	<< std::endl;
    }

    if (mom->alignableObjectId() == type) break; // found

    dau = mom; // move up a level
  }

  findSisters(dau, bias);

  if (theSisters.size() == 0)
  {
    throw cms::Exception("ConfigError")
      << "You are finding an unbiased residual of an alignable "
      << " (id = " << ali.geomDetId().rawId()
      << ") which has no sister. Abort!"
      << std::endl;
  }

  calculate(ali);
}

AlgebraicVector SurveyResidual::sensorResidual() const
{
  align::LocalVector deltaR = theSurface.toLocal(theCurrentVs[0] - theNominalVs[0]);

// Match the centers of current and nominal surfaces to find the angular
// displacement about the center.

  GlobalVectors nominalVs = theNominalVs;
  GlobalVectors currentVs = theCurrentVs;

  for (unsigned int j = 0; j < nominalVs.size(); ++j)
  {
    nominalVs[j] -= theNominalVs[0]; // move to nominal pos
    currentVs[j] -= theCurrentVs[0]; // move to current pos
  }

  RotationType rot = diffRot(nominalVs, currentVs); // frame rotation

  EulerAngles deltaW = toAngles( theSurface.toLocal(rot) );

  AlgebraicVector deltaRW(6); // (deltaR, deltaW)

  deltaRW(1) = deltaR.x();
  deltaRW(2) = deltaR.y();
  deltaRW(3) = deltaR.z();
  deltaRW(4) = deltaW(1);
  deltaRW(5) = deltaW(2);
  deltaRW(6) = deltaW(3);

  return deltaRW;
}

LocalVectors SurveyResidual::pointsResidual() const
{
  LocalVectors residuals;

  unsigned int nPoint = theNominalVs.size();

  residuals.reserve(nPoint);

  for (unsigned int j = 0; j < nPoint; ++j)
  {
    residuals.push_back( theSurface.toLocal(theCurrentVs[j] - theNominalVs[j]) );
  }

  return residuals;
}

void SurveyResidual::findSisters(const Alignable* ali,
				 bool bias)
{
  theSisters.clear();
  theSisters.reserve(1000);

  if (bias) ali->deepComponents(theSisters);

  const Alignable* mother = ali->mother();

  if (mother)
  {
    const std::vector<Alignable*>& comp = mother->components();

    unsigned int nComp = comp.size();

    for (unsigned int i = 0; i < nComp; ++i)
    {
      const Alignable* dau = comp[i];

      if (dau != ali) dau->deepComponents(theSisters);
    }
  }
}

void SurveyResidual::calculate(const Alignable& ali)
{
  unsigned int nSister = theSisters.size();

// First get sisters' positions

  std::vector<const PositionType*> nominalSisPos; // nominal sisters' pos
  std::vector<const PositionType*> currentSisPos; // current sisters' pos

  nominalSisPos.reserve(nSister);
  currentSisPos.reserve(nSister);

  for (unsigned int i = 0; i < nSister; ++i)
  {
    const Alignable* sis    = theSisters[i];
    const SurveyDet* survey = sis->survey();

    if (!survey)
    {
      throw cms::Exception("ConfigError")
	<< "No survey info is found for Alignable "
	<< " (id = " << sis->geomDetId().rawId() << "). Abort!"
	<< std::endl;
    }

    nominalSisPos.push_back( &survey->position() );
    currentSisPos.push_back( &sis->globalPosition() );
  }

// Then find mother's position using sisters' positions

  PositionType nominalMomPos = motherPosition(nominalSisPos);
  PositionType currentMomPos = motherPosition(currentSisPos);

// Now find rotation from nominal mother to current mother

  GlobalVectors nominalSisVs; // nominal sisters' pos from mother's pos
  GlobalVectors currentSisVs; // current sisters' pos from mother's pos

  for (unsigned int i = 0; i < nSister; ++i)
  {
    const Alignable* sis = theSisters[i];

    const GlobalPoints& nominalSisPoints = sis->survey()->globalPoints();
    const GlobalPoints& currentSisPoints = sis->surface().toGlobal( sis->survey()->localPoints() );

    for (unsigned int j = 0; j < nominalSisPoints.size(); ++j)
    {
      nominalSisVs.push_back(nominalSisPoints[j] - *nominalSisPos[i]);
      currentSisVs.push_back(currentSisPoints[j] - *currentSisPos[i]);
//       nominalSisVs.push_back(nominalSisPoints[j] - nominalMomPos);
//       currentSisVs.push_back(currentSisPoints[j] - currentMomPos);
    }
  }

  RotationType toCurrent = diffRot(currentSisVs, nominalSisVs);

// Finally shift and rotate nominal sensor to current sensor

  const SurveyDet* survey = ali.survey();

  if (!survey)
  {
    throw cms::Exception("ConfigError")
      << "No survey info is found for Alignable "
      << " (id = " << ali.geomDetId().rawId() << "). Abort!"
      << std::endl;
  }

  const GlobalPoints& nominalPoints = survey->globalPoints();
  const GlobalPoints& currentPoints = theSurface.toGlobal( survey->localPoints() );

  for (unsigned int j = 0; j < nominalPoints.size(); ++j)
  {
    align::GlobalVector nv = nominalPoints[j] - nominalMomPos;

    theNominalVs.push_back( align::GlobalVector( toCurrent * nv.basicVector() ) );
    theCurrentVs.push_back(currentPoints[j] - currentMomPos);
  }
}
