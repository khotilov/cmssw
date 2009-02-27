#ifndef CastorReco_CastorTower_h
#define CastorReco_CastorTower_h
/** \class reco::CastorTower CastorTower.h DataFormats/CastorReco/CastorTower.h
 *  
 * Class for Castor towers
 *
 * \author Hans Van Haevermaet, University of Antwerp
 *
 * \version $Id: CastorTower.h,v 1.2 2008/11/24 22:19:10 hvanhaev Exp $
 *
 */
#include <vector>
#include <memory>
#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/CastorReco/interface/CastorCell.h"

namespace reco {

  class CastorTower {
  public:

    // default constructor. Set energy and position to zero
 
    CastorTower() : energy_(0.), position_(ROOT::Math::XYZPoint(0.,0.,0.)), emEnergy_(0.), hadEnergy_(0.), fem_(0.), depth_(0.), fhot_(0.) { }

    // constructor from values
    CastorTower(const double energy, const ROOT::Math::XYZPoint& position, const double emEnergy, const double hadEnergy, const double fem,
		const double depth, const double fhot, const CastorCellRefVector& usedCells);

    /// destructor
    virtual ~CastorTower();

    /// tower energy
    double energy() const { return energy_; }

    /// tower centroid position
    ROOT::Math::XYZPoint position() const { return position_; }
    
    /// tower em energy
    double emEnergy() const { return emEnergy_; }
    
    /// tower had energy
    double hadEnergy() const { return hadEnergy_; }
    
    /// tower em/tot ratio
    double fem() const { return fem_; }
    
    /// tower depth in z
    double depth() const { return depth_; }
    
    /// tower  hotcell/tot ratio
    double fhot() const { return fhot_; }

    /// vector of used Cells
    CastorCellRefVector getUsedCells() const { return usedCells_; }

    /// fist iterator over CastorCell constituents
    CastorCell_iterator cellsBegin() const { return usedCells_.begin(); }

    /// last iterator over CastorCell constituents
    CastorCell_iterator cellsEnd() const { return usedCells_.end(); }

    /// number of CastorCell constituents
    size_t cellsSize() const { return usedCells_.size(); }

    /// add reference to constituent CastorCell
    void add( const CastorCellRef & cell ) { usedCells_.push_back( cell ); }

    /// comparison >= operator
    bool operator >=(const CastorTower& rhs) const { return (energy_>=rhs.energy_); }

    /// comparison > operator
    bool operator > (const CastorTower& rhs) const { return (energy_> rhs.energy_); }

    /// comparison <= operator
    bool operator <=(const CastorTower& rhs) const { return (energy_<=rhs.energy_); }

    /// comparison <= operator
    bool operator < (const CastorTower& rhs) const { return (energy_< rhs.energy_); }

    /// pseudorapidity of tower centroid
    double eta() const { return position_.eta(); }

    /// azimuthal angle of tower centroid
    double phi() const { return position_.phi(); }

    /// x of tower centroid
    double x() const { return position_.x(); }

    /// y of tower centroid
    double y() const { return position_.y(); }

    /// rho of tower centroid
    double rho() const { return position_.rho(); }

  private:

    /// tower energy
    double energy_;
    
    /// tower centroid position
    ROOT::Math::XYZPoint position_;
    
    /// tower em energy
    double emEnergy_;
    
    /// tower had energy
    double hadEnergy_;
    
    /// tower em/tot Ratio
    double fem_;
    
    /// tower depth
    double depth_;

    /// tower  hotcell/tot ratio
    double fhot_;

    /// references to CastorCell constituents
    CastorCellRefVector usedCells_;
  };
  
  /// collection of CastorTower objects
  typedef std::vector<CastorTower> CastorTowerCollection;

  // persistent reference to CastorTower objects
  typedef edm::Ref<CastorTowerCollection> CastorTowerRef;

  /// vector of references to CastorTower objects all in the same collection
  typedef edm::RefVector<CastorTowerCollection> CastorTowerRefVector;

  /// iterator over a vector of references to CastorTower objects all in the same collection
  typedef CastorTowerRefVector::iterator CastorTower_iterator;

}

#endif
