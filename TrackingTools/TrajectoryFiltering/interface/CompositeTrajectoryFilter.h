#ifndef CompositeTrajectoryFilter_H
#define CompositeTrajectoryFilter_H

#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilter.h"
#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilterFactory.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TempTrajectory.h"

/** A TrajectoryFilter that stops reconstruction if P_t drops
 *  below some value at some confidence level.
 *  The CkfTrajectoryBuilder uses this class to
 *  implement the minimal P_t cut.
 */

class CompositeTrajectoryFilter : public TrajectoryFilter {
public:

  explicit CompositeTrajectoryFilter(){filters.clear();}

  explicit CompositeTrajectoryFilter( const edm::ParameterSet & pset)
  {
    //look for VPSet of filters
    std::vector<edm::ParameterSet> vpset=pset.getParameter<std::vector<edm::ParameterSet> >("filters");
    for (uint i=0;i!= vpset.size();i++)
      {filters.push_back(TrajectoryFilterFactory::get()->create(vpset[i].getParameter<std::string>("ComponentType"),
								vpset[i]));}
  }
  
  ~CompositeTrajectoryFilter() {}

  virtual bool qualityFilter( const Trajectory& traj) const { return QF<Trajectory>(traj);}
  virtual bool qualityFilter( const TempTrajectory& traj) const { return QF<TempTrajectory>(traj);}
 
  virtual bool toBeContinued( Trajectory& traj) const { return TBC<Trajectory>(traj);}
  virtual bool toBeContinued( TempTrajectory& traj) const { return TBC<TempTrajectory>(traj);}
  
  virtual std::string name() const {return "CompositeTrajectoryFilter";}

protected:
  template <class T> bool TBC( T& traj)const{
    uint i=0;
    uint n=filters.size();
    for (;i<n;i++){ if (!filters[i]->toBeContinued(traj)) return false; }
    return true;}

  template <class T> bool QF(const T& traj)const{
    uint i=0;
    uint n=filters.size();
    for (;i<n;i++){ if (!filters[i]->qualityFilter(traj)) return false; }
    return true;}
 protected:
  std::vector< const TrajectoryFilter *> filters;
  
};

#endif
