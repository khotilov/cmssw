// Last commit: $Id: $

#ifndef DataFormats_SiStripCommon_SiStripKey_H
#define DataFormats_SiStripCommon_SiStripKey_H

#include "DataFormats/SiStripCommon/interface/Constants.h"
#include "DataFormats/SiStripCommon/interface/ConstantsForGranularity.h"
#include <boost/cstdint.hpp>
#include <ostream>
#include <string>

class SiStripKey;

/** Debug info for SiStripKey class. */
std::ostream& operator<< ( std::ostream&, const SiStripKey& );

/**
   @class SiStripKey
   @author R.Bainbridge

   @brief Base utility class that identifies a position within a
   logical structure of the strip tracker.
*/
class SiStripKey {
  
 public:
  
  // ---------- Constructors ----------
  
  /** Constructor using 32-bit "key". */
  SiStripKey( const uint32_t& key );
  
  /** Constructor using directory path. */
  SiStripKey( const std::string& directory_path );
  
  /** Default constructor. */
  SiStripKey();

  /** Virtual destructor. */
  virtual ~SiStripKey() {;}
  
  // ---------- Public interface to member data ----------
  
  /** Returns 32-bit key. */
  inline const uint32_t& key() const;

  /** Returns directory path. */
  inline const std::string& path() const;

  /** Returns granularity to which key is unambiguous. */
  inline const sistrip::Granularity& granularity() const;

  /** Returns channel for key granularity. */
  inline const uint16_t& channel() const;

  // ---------- Pure virtual utility methods ---------- 
  
  /** Identifies key objects with identical member data. */
  virtual bool isEqual( const SiStripKey& ) const = 0;
  
  /** "Consistent" means identical and/or null (ie, "all") data. */
  virtual bool isConsistent( const SiStripKey& ) const = 0;

  /** Identifies all member data as being "valid" or null ("all"). */
  virtual bool isValid() const = 0;
  
  /** All member data to level of "Granularity" are valid. If
      sistrip::Granularity is "undefined", returns false. */
  virtual bool isValid( const sistrip::Granularity& ) const = 0;
  
  /** Identifies all member data as being invalid. */
  virtual bool isInvalid() const = 0;

  /** All member data to level of "Granularity" are invalid. If
      sistrip::Granularity is "undefined", returns true.  */
  virtual bool isInvalid( const sistrip::Granularity& ) const = 0;

 protected: 

  // ---------- Protected methods ----------

  virtual void initFromValue() = 0;
  virtual void initFromKey() = 0;
  virtual void initFromPath() = 0;
  virtual void initGranularity() = 0;
  
  inline void key( const uint32_t& );
  inline void path( const std::string& );
  inline void granularity( const sistrip::Granularity& );
  inline void channel( const uint16_t& );

 private: 
  
  // ---------- Private member data ----------
  
  /** 32-bit key. */
  uint32_t key_; 

  /** Directory path. */
  std::string path_;

  /** Granularity to which FED key is unambiguous. */
  sistrip::Granularity granularity_;
  
  /** Channel of key granularity. */
  uint16_t channel_;
  
};

// ---------- Inline methods ----------

const uint32_t& SiStripKey::key() const { return key_; }
const std::string& SiStripKey::path() const { return path_; }
const sistrip::Granularity& SiStripKey::granularity() const { return granularity_; }
const uint16_t& SiStripKey::channel() const { return channel_; }

void SiStripKey::key( const uint32_t& key ) { key_ = key; }
void SiStripKey::path( const std::string& path ) { path_ = path; }
void SiStripKey::granularity( const sistrip::Granularity& gran ) { granularity_ = gran; }
void SiStripKey::channel( const uint16_t& chan ) { channel_ = chan; }

#endif // DataFormats_SiStripCommon_SiStripKey_H


