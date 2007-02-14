#include "EventFilter/Utilities/interface/ModuleWebRegistry.h"
#include "EventFilter/Utilities/interface/ModuleWeb.h"

#include <iostream>

namespace evf{
ModuleWebRegistry::ModuleWebRegistry(const edm::ParameterSet &ps){
}



bool ModuleWebRegistry::checkWeb(const std::string &name){return clm_.find(name) != clm_.end();}

void ModuleWebRegistry::invoke(xgi::Input *in, xgi::Output *out, const std::string &name)
{
  idct i = clm_.find(name);
  if(i != clm_.end())
    {
      try{
	(*i).second->defaultWebPage(in,out);
      }
      catch(...)
	{
	  std::cout << "exception caught when calling moduleweb page for " << name << std::endl;
	}
    }
}

void ModuleWebRegistry::publish(xdata::InfoSpace *is)
{
    idct i = clm_.begin();
    while (i != clm_.end())
      {
	(*i).second->publish(is);
	i++;
      }

}
void ModuleWebRegistry::clear(){clm_.clear();}


} //end namespace evf
