
#include "CondFormats//SiStripObjects/interface/SiStripLorentzAngle.h"

#include "CondCore/Utilities/interface/PayLoadInspector.h"
#include "CondCore/Utilities/interface/InspectorPythonWrapper.h"

namespace cond {

  template<>
  class ValueExtractor<SiStripLorentzAngle>: public  BaseValueExtractor<SiStripLorentzAngle> {
  public:

    typedef SiStripLorentzAngle Class;
    typedef ExtractWhat<Class> What;
    static What what() { return What();}

    ValueExtractor(){}
    ValueExtractor(What const & what, std::vector<int> const& which)
    {
      // here one can make stuff really complicated...
    }
    void compute(Class const & it){
    }
  private:
  
  };


  template<>
  std::string
  PayLoadInspector<SiStripLorentzAngle>::dump() const {
    std::stringstream ss;
    return ss.str();
    
  }
  
  template<>
  std::string PayLoadInspector<SiStripLorentzAngle>::summary() const {
    std::stringstream ss;
    return ss.str();
  }
  

  template<>
  std::string PayLoadInspector<SiStripLorentzAngle>::plot(std::string const & filename,
						   std::string const &, 
						   std::vector<int> const&, 
						   std::vector<float> const& ) const {
    std::string fname = filename + ".png";
    std::ofstream f(fname.c_str());
    return fname;
  }


}

PYTHON_WRAPPER(SiStripLorentzAngle,SiStripLorentzAngle);
