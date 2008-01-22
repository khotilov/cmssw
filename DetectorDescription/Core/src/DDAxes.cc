

#include "DetectorDescription/Core/interface/DDAxes.h"

namespace DDI { } using namespace DDI;


AxesNames::AxesNames() {
  axesmap_["x"] = x;
  axesmap_["y"] = y;
  axesmap_["z"] = z;
  axesmap_["rho"] = rho;
  axesmap_["radial3D"] = radial3D;
  axesmap_["phi"] = phi;
  axesmap_["undefined"] = undefined;
}

AxesNames::~AxesNames() { }

const std::string AxesNames::name(const DDAxes& s) 
{
  static std::map<std::string, DDAxes>::const_iterator it;

  for (it = axesmap_.begin(); it != axesmap_.end(); ++it)
    {
      if (it->second == s)
	break;
    }
  return it->first;
}

const DDAxes AxesNames::index(const std::string & s)
{
  return axesmap_[s];
}

const std::string DDAxesNames::name(const DDAxes& s) 
{
  return instance().name(s);
}

const DDAxes DDAxesNames::index(const std::string & s) 
{
  return instance().index(s);
}
