#ifndef Geometry_CSCGeometry_CSCChamber_H
#define Geometry_CSCGeometry_CSCChamber_H

/** \class CSCChamber
 *
 * Describes the geometry of the second-level detector unit 
 * modelled by a C++ object in the endcap muon CSC system.
 * A CSCChamber is composed of 6 CSCLayer's and is,
 * of course, a Cathode Strip Chamber Chamber!
 *
 * \author Tim Cox
 */

#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <Geometry/CommonDetUnit/interface/GeomDetType.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include <Geometry/CSCGeometry/interface/CSCChamberSpecs.h>
//#include <Geometry/CSCGeometry/interface/CSCLayer.h>

class CSCLayer;

class CSCChamber : public GeomDet {

public:

  CSCChamber( BoundPlane* bp, CSCDetId id, CSCChamberSpecs* specs ) :
  GeomDet( bp ), theId( id ), theChamberSpecs( specs ), 
    theComponents(6,(const CSCLayer*)0) {}

  ~CSCChamber();

  const GeomDetType& type() const { return *(specs()); }

  DetId geographicalId() const { return theId; } //@@ Slices base

  /// Get the (concrete) DetId.
  CSCDetId id() const { return theId; }

  // Which subdetector
  virtual SubDetector subDetector() const {return GeomDetEnumerators::CSC;}

  const CSCChamberSpecs* specs() const { return theChamberSpecs; }

  /// Return the layers in this chamber
  virtual std::vector< const GeomDet* > components() const;

  /// Return the layer with a given id in this chamber
  virtual const GeomDet* component(DetId id) const;


  // Extension of the interface

  /// Add a layer
  void addComponent( int n, const CSCLayer* gd );

  /// Return all layers
  const std::vector< const CSCLayer* >& layers() const { return theComponents; }

  /// Return the layer corresponding to the given id 
  const CSCLayer* layer(CSCDetId id) const;
  
  /// Return the given layer.
  /// Layers are numbered 1-6.
  const CSCLayer* layer(int ilay) const;

private:

  CSCDetId theId;
  CSCChamberSpecs* theChamberSpecs;
  std::vector< const CSCLayer* > theComponents; // the 6 CSCLayers comprising a CSCChamber; are owned by this class
};

#endif // Geometry_CSCGeometry_CSCChamber_H
