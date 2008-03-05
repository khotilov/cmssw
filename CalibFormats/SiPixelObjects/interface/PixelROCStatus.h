#ifndef PixelROCStatus_h
#define PixelROCStatus_h
//
// This class keeps the possible non-standard
// status a ROC can have.
//
//
//

#include <stdint.h>
#include <set>
#include <string>

namespace pos{

  class PixelROCStatus {


  private:

    uint32_t bits_;

  public:

    //Insert new status before nStatus
    enum ROCstatus {off=0, noHits, nStatus};

    PixelROCStatus();
    PixelROCStatus(const std::set<ROCstatus>& stat);
    virtual ~PixelROCStatus();

    std::string statusName(ROCstatus stat) const;
    std::string statusName() const;

    void set(ROCstatus stat);
    void clear(ROCstatus stat);
    void set(ROCstatus stat, bool mode);
    void set(const std::string& statName);
    bool get(ROCstatus stat) const ;

    // Added by Dario (March 4th 2008)
    void reset(void) ;
 
  };
}
#endif
