#ifndef WORKERMAKER_H
#define WORKERMAKER_H

#include "FWCore/Framework/src/Worker.h"
#include "FWCore/Framework/interface/ModuleDescription.h"
#include "FWCore/Framework/src/WorkerParams.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// The following includes are temporary until a better
// solution can be found.  Placing these includes here
// leads to more physical coupling than is probably necessary.
// Another solution is to build a typeid lookup table in the 
// implementation file (one every for each XXXWorker) and
// then include all the relevent worker headers in the 
// implementation file only.
#include "FWCore/Framework/src/ProducerWorker.h"
#include "FWCore/Framework/src/FilterWorker.h"
#include "FWCore/Framework/src/AnalyzerWorker.h"
#include "FWCore/Framework/src/OutputWorker.h"

#include <memory>
#include <string>

namespace edm {
  
  class Maker
  {
  public:
    virtual ~Maker();
    virtual std::auto_ptr<Worker> makeWorker(const WorkerParams&) const = 0;
  };

  template <class T>
  class WorkerMaker : public Maker
  {
  public:
    typedef T worker_type;
    explicit WorkerMaker();
    std::auto_ptr<Worker> makeWorker(const WorkerParams&) const;
  };

  template <class T>
  WorkerMaker<T>::WorkerMaker()
  {
  }

  template <class T>
  std::auto_ptr<Worker> WorkerMaker<T>::makeWorker(const WorkerParams& p) const
  {
    typedef T UserType;
    typedef typename UserType::ModuleType ModuleType;
    typedef typename WorkerType<ModuleType>::worker_type  worker_type;

    const ParameterSet& conf = *p.pset_;
    ModuleDescription md;
    md.pid = PS_ID("oink"); // conf.id();
    md.module_name = conf.template getParameter<std::string>("module_type");
    md.module_label = conf.template getParameter<std::string>("module_label");
    md.version_number = p.version_number_;
    md.process_name = p.process_name_;;
    md.pass = p.pass_; 

    std::auto_ptr<ModuleType> module(new UserType(conf));
    std::auto_ptr<Worker> worker(new worker_type(module, md, p));
    return worker;
  }

}

#endif
