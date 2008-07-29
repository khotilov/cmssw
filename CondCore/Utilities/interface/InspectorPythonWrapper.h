#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "CondCore/DBCommon/interface/ClassID.h"


using namespace boost::python;

namespace condPython {
  template<typename T>
  void defineWhat() {
    typedef cond::ExtractWhat<T> What;
    class_<What>("What",init<>());
  }

}

namespace {
  
  template<typename Wrapper>
  void define() {
    typedef typename Wrapper::Extractor Extractor;

    condPython::defineWhat<typename Extractor::Class>();

    class_<Extractor>("Extractor", init<>())
      .def(init<std::string, std::vector<int> >())
      .def("what",Extractor::what)
      .def("values",&Extractor::values, return_value_policy<copy_const_reference>())
      ;

    class_<Wrapper>("Object",init<>()) 
      .def(init<cond::IOVElement>())
      .def("dump",&Wrapper::dump)
      .def("plot",&Wrapper::plot)
      .def("summary",&Wrapper::summary) 
      .def("extract",&Wrapper::extract)
      ; 
  }
}


#define PYTHON_WRAPPER(_class,_name) \
namespace { typedef cond::PayLoadInspector< _class > PythonWrapper;} \
 BOOST_PYTHON_MODULE(plugin ## _name ## PyInterface) { define<PythonWrapper>(); } \
namespace { const char * pluginName_="plugin"  #_name "PyInterface"; }\
PYTHON_ID(PythonWrapper::Class, pluginName_)


