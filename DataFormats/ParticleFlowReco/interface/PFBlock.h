#ifndef RecoParticleFlow_PFAlgo_PFBlock_h
#define RecoParticleFlow_PFAlgo_PFBlock_h 

#include <map>
#include <iostream>

/* #include "boost/graph/adjacency_matrix.hpp" */


// #include "DataFormats/ParticleFlowReco/interface/PFBlockLink.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
// #include "DataFormats/ParticleFlowReco/interface/PFBlockParticle.h"
#include "DataFormats/Common/interface/OwnVector.h"



namespace reco {

  
  /// \brief Block of elements
  /*!
    \author Colin Bernet
    \date January 2006

    A PFBlock is: 
    - a set of topologically connected elements.
    - a set of links between these elements
  */
  
  class PFBlock {


  public:

    typedef edm::OwnVector< reco::PFBlockElement >::const_iterator IE;
    /*     typedef std::vector< reco::PFBlockLink >::const_iterator IL; */
    
    // typedef std::vector< std::vector<double> > LinkData;
    typedef std::vector< std::vector<std::pair<float,float> > > LinkData;
    
    enum LinkTest {
      LINKTEST_CHI2=0,
      LINKTEST_RECHIT,
      LINKTEST_TANGENT,
      LINKTEST_NLINKTEST,
      LINKTEST_ALL
    };

    PFBlock() {}
    // PFBlock(const PFBlock& other);

    /// add an element to the current PFBlock
    /// the block will keep a copy.
    void addElement( reco::PFBlockElement* element );
    
    void bookLinkData();

    /// makes the correspondance between a 2d element matrix and 
    /// the 1D vector which is the most compact way to store the matrix
    bool matrix2vector(unsigned i, unsigned j, unsigned& index) const;

    /// set a link between elements of indices i1 and i2, of "distance" chi2
    /// the link is set in the linkData vector provided as an argument.
    /// As indicated by the 'const' statement, 'this' is not modified.
    void setLink(unsigned i1, unsigned i2, 
		 double chi2, double dist, 
                 LinkData& linkData, 
		 LinkTest  test=LINKTEST_CHI2 ) const;

    /// lock an element ( unlink it from the others )
    /// Colin: this function is misleading
    /// void lock(unsigned i, LinkData& linkData ) const;


    /// fills a map with the elements associated to element i.
    /// elements are sorted by increasing chi2.
    /// if specified, only the elements of type "type" will be considered
    /// if specified, only the link calculated from a certain "test" will 
    /// be considered: chi2 test, rechit test, tangent test etc..
    void associatedElements( unsigned i,
                             const LinkData& linkData, 
                             std::multimap<double, unsigned>& sortedAssociates,
                             reco::PFBlockElement::Type type = PFBlockElement::NONE,
			     LinkTest test=LINKTEST_CHI2 ) const; 
      

    /// \return chi2 of link
    double chi2( unsigned ie1, unsigned ie2, 
                 const LinkData& linkData, 
		 LinkTest  test=LINKTEST_CHI2 ) const;

    /// \return distance of link
    double dist( unsigned ie1, unsigned ie2, 
                 const LinkData& linkData, 
		 LinkTest  test=LINKTEST_CHI2 ) const;

    /// \return elements
    const edm::OwnVector< reco::PFBlockElement >& elements() const 
      {return elements_;}

    /// \return link data
    const LinkData& linkData() const 
      {return linkData_;}

    /// \return link data
    LinkData& linkData()  
      {return linkData_;}

    friend std::ostream& operator<<( std::ostream& out, const PFBlock& co );

  private:
    
    /// \return size of linkData_, calculated from the number of elements
    unsigned linkDataSize() const;
    
    /// link data (permanent)
    LinkData                                        linkData_;
     
    /// all elements 
    edm::OwnVector< reco::PFBlockElement >          elements_;
        
  };
}

#endif


  
