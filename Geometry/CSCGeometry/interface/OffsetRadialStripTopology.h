#ifndef ABS_OFFSET_RADIAL_STRIP_TOPOLOGY_H
#define ABS_OFFSET_RADIAL_STRIP_TOPOLOGY_H

/** \class OffsetRadialStripTopology
 *  ABC defining  RadialStripTopology with shifted offset so that it
 *  is not centred on local y (of parent chamber.)
 *  The offset is specified as a fraction of the strip angular width.
 *
 *  \author Tim Cox
 * 
 */

#include "Geometry/CommonTopologies/interface/RadialStripTopology.h"
#include <iosfwd>

class OffsetRadialStripTopology : public RadialStripTopology
{
public:

  OffsetRadialStripTopology( int numberOfStrips, float stripPhiPitch,
     float detectorHeight, float radialDistance, float stripOffset);

  virtual ~OffsetRadialStripTopology(){};

  /** Fraction of a strip offset of layer relative to
   *  symmetry axis (local y). (This is an ANGULAR value)
   */
  virtual float stripOffset( void ) const { return theStripOffset; }

   /** LocalPoint for a given MeasurementPoint <BR>
   * What's a MeasurementPoint?  <BR>
   * A MeasurementPoint is a 2-dim object.<BR>
   * The first dimension measures the
   * angular position wrt central line of symmetry of detector,
   * in units of strip (angular) widths (range 0 to total angle subtended
   * by a detector).<BR>
   * The second dimension measures
   * the fractional position along the strip (range -0.5 to +0.5).<BR>
   * BEWARE! The components are not Cartesian.<BR>
   * BEWARE! Neither coordinate may correspond to either local x or local y.<BR>
   * BEWARE! This involves ONLY strip-related measurements, not CSC wires! <BR>
   */
  virtual LocalPoint localPosition(const MeasurementPoint&) const;

  /** Strip in which a given LocalPoint lies. This is a float which
   * represents the fractional strip position within the detector.<BR>
   * Returns zero if the LocalPoint falls at the extreme low edge of the
   * detector or BELOW, and float(nstrips) if it falls at the extreme high
   * edge or ABOVE.
   */
  virtual float strip(const LocalPoint&) const;

  /**
   * Angle between strip and local y axis (measured clockwise from y axis)
   */
  float stripAngle(float strip) const;

  /**
   * Channel number corresponding to a strip or a LocalPoint.
   * Sometimes more than one strip is OR'ed into one channel.
   */
  virtual int channel(int strip) const = 0;
  virtual int channel(const LocalPoint& lp) const = 0;

  friend std::ostream & operator<<(std::ostream &, const OffsetRadialStripTopology &);

 private:
  /**
   * Transform from coordinates w.r.t. strip plane symmetry axes to
   * local coordinates
   */
  LocalPoint toLocal(float xprime, float yprime) const;
  /**
   * Transform from local coordinates to coordinates w.r.t. strip plane
   * symmetry axes
   */
  LocalPoint toPrime(const LocalPoint&) const;

  float theStripOffset; // fraction of a strip offset from sym about y
  float theCosOff; // cosine of angular offset
  float theSinOff; // sine of angular offset
};

#endif

