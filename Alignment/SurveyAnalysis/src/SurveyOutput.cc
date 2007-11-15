#include <sstream>

#include "TNtupleD.h"

#include "Alignment/CommonAlignment/interface/Alignable.h"
// #include "Alignment/CommonAlignment/interface/AlignmentParameters.h"
#include "Alignment/CommonAlignment/interface/SurveyDet.h"

#include "Alignment/SurveyAnalysis/interface/SurveyOutput.h"

using namespace align;

SurveyOutput::SurveyOutput(const std::vector<Alignable*>& alignables,
			   const std::string& fileName):
  theAlignables(alignables),
  theFile(fileName.c_str(), "RECREATE")
{
}

void SurveyOutput::write(unsigned int iter)
{
  std::ostringstream o;

  o << 't' << iter;

  TNtupleD* nt = new TNtupleD(o.str().c_str(), "", "id:x:y:z:a:b:g");

  unsigned int N = theAlignables.size();

  for (unsigned int i = 0; i < N; ++i)
  {
    const Alignable* ali = theAlignables[i];

    const align::GlobalVector& shifts = ali->displacement();

    EulerAngles angles = toAngles( ali->rotation() );

    nt->Fill( ali->id(), shifts.x(), shifts.y(), shifts.z(),
	      angles(1), angles(2), angles(3) );
//     const AlgebraicVector& pars = ali->alignmentParameters()->parameters();

//     nt->Fill(pars[0], pars[1], pars[2], pars[3], pars[4], pars[5]);
  }

  theFile.Write();

  delete nt;
}
