#ifndef SimVertex_H
#define SimVertex_H

#include "SimDataFormats/Vertex/interface/CoreSimVertex.h"
class SimVertex : public CoreSimVertex
{
 public:
  
  typedef CoreSimVertex Core;
  /// constructor
    SimVertex();
    SimVertex(const Hep3Vector & v, float tof);
    /// full constructor (position, time, index of parent in final vector)
    SimVertex(const Hep3Vector & v, float tof, int it);
    /// constructor from transient
    SimVertex(const CoreSimVertex & t, int it);
    /// index of the parent in the Event SimTrack container (-1 if no parent)
    int parentIndex() const { return  itrack; }
    bool noParent() const { return  itrack==-1; }
private: 
    int itrack;
};

#include <iosfwd>
std::ostream & operator <<(std::ostream & o , const SimVertex& v);
 

#endif
