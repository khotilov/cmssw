#ifndef __PFBlockElement__
#define __PFBlockElement__

#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"


#include <iostream>


namespace reco {
  class PFBlockElementCluster;
  class PFBlockElementTrack;
  
  
  /// \brief Base element of a PFBlock (track, cluster...)
  /// 
  /// this class is essentially wraps a PFRecTrackRef of a 
  /// PFClusterRef, depending on the type of the element
  class PFBlockElement {
  public:
    
    /// number of element types
    static const unsigned nTypes_;

    /// possible types for the element
    enum Type { 
      NONE=0,
      TRACK, 
      PS1, 
      PS2, 
      ECAL, 
      HCAL, 
      MUON
    };
    

    /// default constructor 
    PFBlockElement() :  
      type_( NONE ), 
      locked_(false), 
      index_( static_cast<unsigned>(-1) ) {
    }

    /// standard constructor
    PFBlockElement(Type type) :  
      type_(type), 
      locked_(false),
      index_( static_cast<unsigned>(-1) ) {
    }


    /// destructor
    virtual ~PFBlockElement() {}
  
    /// print the object inside the element
    virtual void Dump(std::ostream& out=std::cout, 
		      const char* tab=" " ) const;

    /// necessary to have the edm::OwnVector<PFBlockElement> working
    virtual PFBlockElement* clone() const = 0;
      
    /// lock element
    void lock() {locked_ = true;}

    /// unlock element
    void unLock() {locked_ = false;}

    /// \return type
    Type type() const { return type_; }

    /// locked ? 
    bool    locked() const {return locked_;}
    
    /// set index 
    void     setIndex(unsigned index) { index_ = index; }

    /// \return index
    unsigned index() const {return index_;} 

    virtual PFRecTrackRef trackRef()  const {return PFRecTrackRef();}
    virtual PFClusterRef  clusterRef() const {return PFClusterRef();}


    friend std::ostream& operator<<( std::ostream& out, 
				     const PFBlockElement& element );
    
  protected:  
  
    /// type
    Type     type_;
  
    /// locked ? can probably be transient. should be replaced by a 
    // "remaining energy"
    bool       locked_;
    
    /// index in block vector 
    unsigned   index_;

  };
}
#endif
