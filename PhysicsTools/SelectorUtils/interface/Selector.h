#ifndef Selector_h
#define Selector_h

/**
  \class    Selector Selector.h "CommonTools/Utils/interface/Selector.h"
  \brief    Implements a string-indexed bit_vector

  class template that implements an interface to Selector utilities. This
  allows the user to access individual cuts based on a string index.
  The user can then turn individual cuts on and off at will. 

  \author Salvatore Rappoccio
  \version  $Id: Selector.h,v 1.4 2010/03/10 17:26:55 srappocc Exp $
*/


#include "PhysicsTools/SelectorUtils/interface/strbitset.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <fstream>
#include <functional>

/// Functor that operates on <T>
template<class T>
class Selector : public std::binary_function<T, std::strbitset, bool>  {
  
 public:
  typedef T                                            data_type;
  typedef std::binary_function<T,std::strbitset,bool>  base_type;
  typedef std::pair<std::string, size_t>               cut_flow_item;
  typedef std::vector<cut_flow_item>                   cut_flow_map;
  typedef std::map<std::string, int>                   int_map;
  typedef std::map<std::string, double>                double_map;

  /// Constructor clears the bits
  Selector( ) { 
    bits_.clear(); 
    intCuts_.clear(); 
    doubleCuts_.clear();
    cutFlow_.clear();
    retInternal_ = getBitTemplate();
  }
  virtual ~Selector() {}

  /// This is the registration of an individual cut string
  virtual void push_back( std::string s) {
    bits_.push_back(s);
    // don't need to check to see if the key is already there,
    // bits_ does that.
    cutFlow_.push_back( cut_flow_item(s,0) );
  }


  /// This is the registration of an individual cut string, with an int cut value
  virtual void push_back( std::string s, int cut) {
    bits_.push_back(s);
    intCuts_[s] = cut;
    // don't need to check to see if the key is already there,
    // bits_ does that.
    cutFlow_.push_back( cut_flow_item(s,0) );
  }

  /// This is the registration of an individual cut string, with a double cut value
  virtual void push_back( std::string s, double cut) {
    bits_.push_back(s);
    doubleCuts_[s] = cut;
    // don't need to check to see if the key is already there,
    // bits_ does that.
    cutFlow_.push_back( cut_flow_item(s,0) );
  }

  /// This provides the interface for base classes to select objects
  virtual bool operator()( T const & t, std::strbitset & ret ) = 0;

  /// This provides an alternative signature without the second ret
  virtual bool operator()( T const & t )
  {
    retInternal_.set(false);
    operator()(t, retInternal_);
    setIgnored(retInternal_);
    return (bool)retInternal_;
  }
  
  /// Set a given selection cut, on or off
  void set(std::string s, bool val = true) {
    bits_[s] = val;
  }

  /// Set a given selection cut, on or off, and reset int cut value
  void set(std::string s, int cut, bool val = true) {
    bits_[s] = val;
    intCuts_[s] = cut;
  }

  /// Set a given selection cut, on or off, and reset int cut value
  void set(std::string s, double cut, bool val = true) {
    bits_[s] = val;
    doubleCuts_[s] = cut;
  }

  /// Turn off a given selection cut. 
  void clear(std::string s) {
    bits_[s] = false;
  }

  /// Access the selector cut at index "s".
  /// "true" means to consider the cut.
  /// "false" means to ignore the cut.
  bool operator[] ( std::string s ) const {
    return bits_[s];
  }

  /// consider the cut at index "s"
  bool considerCut( std::string s ) const {
    return bits_[s] == true;
  }

  /// ignore the cut at index "s"
  bool ignoreCut( std::string s ) const {
    return bits_[s] == false;
  }

  /// set the bits to ignore from a vector
  void setIgnoredCuts( std::vector<std::string> const & bitsToIgnore ) {
    for ( std::vector<std::string>::const_iterator ignoreBegin = bitsToIgnore.begin(),
	    ignoreEnd = bitsToIgnore.end(), ibit = ignoreBegin;
	  ibit != ignoreEnd; ++ibit ) {
      set(*ibit, false );
    }
  }

  /// Passing cuts
  void passCut( std::strbitset & ret, std::string const & s ) {
    ret[s] = true;
    cut_flow_map::iterator found = cutFlow_.end();
    for ( cut_flow_map::iterator cutsBegin = cutFlow_.begin(),
	    cutsEnd = cutFlow_.end(), icut = cutsBegin;
	  icut != cutsEnd && found == cutsEnd; ++icut ) {
      if ( icut->first == s ) {
	found = icut;
      }
    }
    ++(found->second);
  }

  /// Access the int cut values at index "s"
  int cut( std::string s, int val ) const {
    return intCuts_.find( s )->second;
  };
  /// Access the double cut values at index "s"
  double cut( std::string s, double val ) const {
    return doubleCuts_.find( s )->second;
  };

  /// Get an empty bitset with the proper names
  std::strbitset getBitTemplate() const { 
    std::strbitset ret = bits_; 
    ret.set(false);
    for ( cut_flow_map::const_iterator cutsBegin = cutFlow_.begin(),
	    cutsEnd = cutFlow_.end(), icut = cutsBegin;
	  icut != cutsEnd; ++icut ) {
      if ( ignoreCut(icut->first) ) ret[icut->first] = true;
    }     
    return ret;
  }

  /// set ignored bits
  void setIgnored( std::strbitset & ret ) {
    for ( cut_flow_map::const_iterator cutsBegin = cutFlow_.begin(),
	    cutsEnd = cutFlow_.end(), icut = cutsBegin;
	  icut != cutsEnd; ++icut ) {
      if ( ignoreCut(icut->first) ) ret[icut->first] = true;
    } 
  }

  /// Print the cut flow
  void print(std::ostream & out) const { 
    for ( cut_flow_map::const_iterator cutsBegin = cutFlow_.begin(),
	    cutsEnd = cutFlow_.end(), icut = cutsBegin;
	  icut != cutsEnd; ++icut ) {
      char buff[1000];
      if ( considerCut( icut->first ) ) {
	sprintf(buff, "%6lu : %20s %10lu", 
		static_cast<unsigned long>(icut - cutsBegin),
		icut->first.c_str(),
		static_cast<unsigned long>(icut->second) );
      } else {
	sprintf(buff, "%6lu : %20s %10s", 
		static_cast<unsigned long>(icut - cutsBegin),
		icut->first.c_str(),
		"off" );
      }
      out << buff << std::endl;
    } 
  }

 protected:
  std::strbitset bits_;        //!< the bitset indexed by strings
  std::strbitset retInternal_; //!< internal ret if users don't care about return bits
  int_map        intCuts_;     //!< the int-value cut map
  double_map     doubleCuts_;  //!< the double-value cut map
  cut_flow_map   cutFlow_;     //!< map of cut flows in "human" order
};

#endif
