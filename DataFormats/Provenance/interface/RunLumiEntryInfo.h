#ifndef DataFormats_Provenance_RunLumiEntryInfo_h
#define DataFormats_Provenance_RunLumiEntryInfo_h

/*----------------------------------------------------------------------
  
RunLumiEntryInfo: The event dependent portion of the description of a product
and how it came into existence, plus the product identifier and the status.

----------------------------------------------------------------------*/
#include <iosfwd>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "DataFormats/Provenance/interface/BranchID.h"
#include "DataFormats/Provenance/interface/ModuleDescriptionID.h"
#include "DataFormats/Provenance/interface/ProductStatus.h"
#include "DataFormats/Provenance/interface/ProvenanceFwd.h"

/*
  RunLumiEntryInfo
*/
#include "DataFormats/Provenance/interface/ProductID.h"

namespace edm {
  class RunLumiEntryInfo {
  public:
    RunLumiEntryInfo();
    explicit RunLumiEntryInfo(BranchID const& bid);
    RunLumiEntryInfo(BranchID const& bid,
		    ProductStatus status);
    RunLumiEntryInfo(BranchID const& bid,
		    ProductStatus status,
		    ModuleDescriptionID const& mid,
		    ProductID const& pid = ProductID(),
		    std::vector<ProductID> const& parents = std::vector<ProductID>());

    RunLumiEntryInfo(BranchID const& bid,
		    ProductStatus status,
		    ProductID const& pid,
		    EntryDescriptionID const& edid);

    ~RunLumiEntryInfo() {}

    void write(std::ostream& os) const;

    ProductID const& productID() const {assert(0 && "Run and lumi products do not have productID's"); return *new ProductID;}
    BranchID const& branchID() const {return branchID_;}
    ProductStatus const& productStatus() const {return productStatus_;}
    ModuleDescriptionID const& moduleDescriptionID() const {return moduleDescriptionID_;}
    void setStatus(ProductStatus status) {productStatus_ = status;}
    void setModuleDescriptionID(ModuleDescriptionID const& mdid) {moduleDescriptionID_ = mdid;}
    void setPresent();
    void setNotPresent();

  private:
    BranchID branchID_;
    ProductStatus productStatus_;
    ModuleDescriptionID moduleDescriptionID_;
  };

  inline
  bool
  operator < (RunLumiEntryInfo const& a, RunLumiEntryInfo const& b) {
    return a.branchID() < b.branchID();
  }
  
  inline
  std::ostream&
  operator<<(std::ostream& os, RunLumiEntryInfo const& p) {
    p.write(os);
    return os;
  }

  // Only the 'salient attributes' are testing in equality comparison.
  bool operator==(RunLumiEntryInfo const& a, RunLumiEntryInfo const& b);
  inline bool operator!=(RunLumiEntryInfo const& a, RunLumiEntryInfo const& b) { return !(a==b); }

  typedef RunLumiEntryInfo LumiEntryInfo;
  typedef RunLumiEntryInfo RunEntryInfo;
  typedef std::vector<RunLumiEntryInfo> RunLumiEntryInfoVector;
  typedef std::vector<RunLumiEntryInfo> LumiEntryInfoVector;
  typedef std::vector<RunLumiEntryInfo> RunEntryInfoVector;
}
#endif
#ifndef DataFormats_Provenance_EventEntryInfo_h
#define DataFormats_Provenance_EventEntryInfo_h

/*----------------------------------------------------------------------
  
EventEntryInfo: The event dependent portion of the description of a product
and how it came into existence, plus the product identifier and the status.

----------------------------------------------------------------------*/
#include <iosfwd>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "DataFormats/Provenance/interface/BranchID.h"
#include "DataFormats/Provenance/interface/EntryDescriptionID.h"
#include "DataFormats/Provenance/interface/ProductID.h"
#include "DataFormats/Provenance/interface/ProductStatus.h"
#include "DataFormats/Provenance/interface/ProvenanceFwd.h"

/*
  EventEntryInfo
*/

namespace edm {
  class EventEntryInfo {
  public:
    EventEntryInfo();
    explicit EventEntryInfo(BranchID const& bid);
    EventEntryInfo(BranchID const& bid,
		    ProductStatus status,
		    ProductID const& pid);
    EventEntryInfo(BranchID const& bid,
		    ProductStatus status,
		    ProductID const& pid,
		    boost::shared_ptr<EntryDescription> edPtr);
    EventEntryInfo(BranchID const& bid,
		    ProductStatus status,
		    ProductID const& pid,
		    EntryDescriptionID const& edid);

    EventEntryInfo(BranchID const& bid,
		   ProductStatus status,
		   ModuleDescriptionID const& mdid,
		   ProductID const& pid = ProductID(),
		   std::vector<ProductID> const& parents = std::vector<ProductID>());

    ~EventEntryInfo() {}

    void write(std::ostream& os) const;

    BranchID const& branchID() const {return branchID_;}
    ProductID const& productID() const {return productID_;}
    ProductStatus const& productStatus() const {return productStatus_;}
    EntryDescriptionID const& entryDescriptionID() const {return entryDescriptionID_;}
    EntryDescription const& entryDescription() const;
    void setStatus(ProductStatus status) {productStatus_ = status;}
    void setPresent();
    void setNotPresent();

  private:
    BranchID branchID_;
    ProductID productID_;
    ProductStatus productStatus_;
    EntryDescriptionID entryDescriptionID_;
    mutable boost::shared_ptr<EntryDescription> entryDescriptionPtr_;
  };

  inline
  bool
  operator < (EventEntryInfo const& a, EventEntryInfo const& b) {
    return a.branchID() < b.branchID();
  }
  
  inline
  std::ostream&
  operator<<(std::ostream& os, EventEntryInfo const& p) {
    p.write(os);
    return os;
  }

  // Only the 'salient attributes' are testing in equality comparison.
  bool operator==(EventEntryInfo const& a, EventEntryInfo const& b);
  inline bool operator!=(EventEntryInfo const& a, EventEntryInfo const& b) { return !(a==b); }

  typedef std::vector<EventEntryInfo> EventEntryInfoVector;
}
#endif
