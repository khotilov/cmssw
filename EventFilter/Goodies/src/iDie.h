#ifndef EVENTFILTER_GOODIES_IDIE_H
#define EVENTFILTER_GOODIES_IDIE_H

#include "EventFilter/Utilities/interface/Exception.h"

#include "xdata/String.h"
#include "xdata/UnsignedInteger32.h"
#include "xdata/ActionListener.h"

#include "xoap/MessageReference.h"
#include "xoap/MessageFactory.h"
#include "xoap/Method.h"

#include "xgi/Utils.h"
#include "xgi/Input.h"
#include "xgi/Output.h"
#include "xgi/Method.h"

#include "xdaq/Application.h"

#include "toolbox/net/URN.h"
#include "toolbox/fsm/exception/Exception.h"


#include <vector>

#include <sys/time.h>


namespace evf {


  namespace internal{
   struct fu{
      time_t tstamp;
      unsigned int ccount;
      std::vector<pid_t> cpids;
      std::vector<std::string> signals;
      std::vector<std::string> stacktraces;
    };
  }
  typedef std::map<std::string,internal::fu> fmap;
  typedef fmap::iterator ifmap;
  
  class iDie : public xdaq::Application,
    public xdata::ActionListener
  {
  public:
    //
    // xdaq instantiator macro
    //
    XDAQ_INSTANTIATOR();
  
    
    //
    // construction/destruction
    //
    iDie(xdaq::ApplicationStub *s);
    virtual ~iDie();
    //UI
    void defaultWeb(xgi::Input *in,xgi::Output *out)
      throw (xgi::exception::Exception);
    void summaryTable(xgi::Input *in,xgi::Output *out)
      throw (xgi::exception::Exception);
    void detailsTable(xgi::Input *in,xgi::Output *out)
      throw (xgi::exception::Exception);
    void dumpTable(xgi::Input *in,xgi::Output *out)
      throw (xgi::exception::Exception);
    void updater(xgi::Input *in,xgi::Output *out)
      throw (xgi::exception::Exception);
    void iChoke(xgi::Input *in,xgi::Output *out)
      throw (xgi::exception::Exception);
    //AI
    void postEntry(xgi::Input*in,xgi::Output*out)
      throw (xgi::exception::Exception);
    void postEntryiChoke(xgi::Input*in,xgi::Output*out)
      throw (xgi::exception::Exception);
    
    // *fake* fsm soap command callback
    xoap::MessageReference fsmCallback(xoap::MessageReference msg)
      throw (xoap::exception::Exception);

    // xdata:ActionListener interface
    void actionPerformed(xdata::Event& e);


  private:

    struct sorted_indices{
      sorted_indices(const std::vector<int> &arr) : arr_(arr)
      {
	ind_.resize(arr_.size(),0);
	unsigned int i = 0;
	while(i<ind_.size()) ind_[i] = i++;
	std::sort(ind_.rbegin(),ind_.rend(),*this);
      }
      int operator[](size_t ind) const {return arr_[ind_[ind]];}
      
      bool operator()(const size_t a, const size_t b) const
      {
	return arr_[a]<arr_[b];
      }
      int ii(size_t ind){return ind_[ind];}
      std::vector<int> ind_;
      const std::vector<int> &arr_;
    };
    //
    // private member functions
    //
    
    void reset();
    void parseLegenda(std::string);
    void parseHisto(const char *, unsigned int);
    //
    // member data
    //

    // message logger
    Logger                          log_;
		
    // monitored parameters
    xdata::String                   url_;
    xdata::String                   class_;
    xdata::UnsignedInteger32        instance_;
    xdata::String                   hostname_;
    xdata::UnsignedInteger32        runNumber_;
    xdata::String                   configString_;
    fmap                            fus_;
    
    unsigned int                    totalCores_;
    unsigned int                    nstates_;   
    std::vector<int>                cpuentries_;
    std::vector<std::vector<int> >  cpustat_;
    std::vector<std::string>        mapmod_;
    unsigned int                    last_ls_;

  }; // class iDie


} // namespace evf


#endif
