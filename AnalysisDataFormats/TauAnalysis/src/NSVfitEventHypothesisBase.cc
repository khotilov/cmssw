#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisBase.h"

NSVfitEventHypothesisBase::NSVfitEventHypothesisBase(const NSVfitEventHypothesisBase& bluePrint)
  : name_(bluePrint.name_),
    barcode_(bluePrint.barcode_),
    met_(bluePrint.met_),
    eventVertexPosition_(bluePrint.eventVertexPosition_),
    eventVertexPositionErr_(bluePrint.eventVertexPositionErr_),
    eventVertexPositionShift_(bluePrint.eventVertexPositionShift_),
    eventVertexIsValid_(bluePrint.eventVertexIsValid_),
    resonances_(bluePrint.resonances_)
{}

NSVfitEventHypothesisBase& NSVfitEventHypothesisBase::operator=(const NSVfitEventHypothesisBase& bluePrint)
{
  name_ = bluePrint.name_;
  barcode_ = bluePrint.barcode_;
  met_ = bluePrint.met_;
  eventVertexPosition_ = bluePrint.eventVertexPosition_;
  eventVertexPositionErr_ = bluePrint.eventVertexPositionErr_;
  eventVertexPositionShift_ = bluePrint.eventVertexPositionShift_;
  eventVertexIsValid_ = bluePrint.eventVertexIsValid_;
  resonances_ = bluePrint.resonances_;
  return (*this);
}

//
//-------------------------------------------------------------------------------
//

bool operator<(const NSVfitEventHypothesisBase& hypothesis1, const NSVfitEventHypothesisBase& hypothesis2)
{
  return (hypothesis1.barcode() < hypothesis2.barcode());
}

