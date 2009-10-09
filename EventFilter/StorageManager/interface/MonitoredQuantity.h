// $Id: MonitoredQuantity.h,v 1.3 2009/07/20 13:06:10 mommsen Exp $
/// @file: MonitoredQuantity.h 

#ifndef StorageManager_MonitoredQuantity_h
#define StorageManager_MonitoredQuantity_h

#include "boost/shared_ptr.hpp"
#include "boost/thread/mutex.hpp"
#include <vector>
#include <math.h>

#include "EventFilter/StorageManager/interface/Utils.h"


namespace stor
{

  /**
   * This class keeps track of statistics for a set of sample values 
   * and provides timing information on the samples.
   *
   * $Author: mommsen $
   * $Revision: 1.3 $
   * $Date: 2009/07/20 13:06:10 $
   */

  class MonitoredQuantity
  {
    
  public:
    class Stats;

    enum DataSetType { FULL = 0,      // the full data set (all samples)
                       RECENT = 1 };  // recent data only

    static utils::duration_t ExpectedCalculationInterval();

    explicit MonitoredQuantity(double timeWindowForRecentResults = 10.0);

    /**
     * Adds the specified doubled valued sample value to the monitor instance.
     */
    void addSample(const double value = 1.0);

    /**
     * Adds the specified integer valued sample value to the monitor instance.
     */
    void addSample(const int value = 1);

    /**
     * Adds the specified uint32_t valued sample value to the monitor instance.
     */
    void addSample(const uint32_t value = 1);

    /**
     * Adds the specified uint64_t valued sample value to the monitor instance.
     */
    void addSample(const uint64_t value = 1);

    /**
     * Forces a calculation of the statistics for the monitored quantity.
     * The frequency of the updates to the statistics is driven by how
     * often this method is called.  It is expected that this method
     * will be called once per interval specified by
     * EXPECTED_CALCULATION_INTERVAL.
     */
    void calculateStatistics(utils::time_point_t currentTime = 
                             utils::getCurrentTime());

    /**
     * Resets the monitor (zeroes out all counters and restarts the
     * time interval).
     */
    void reset();

    /**
     * Enables the monitor (and resets the statistics to provide a
     * fresh start).
     */
    void enable();

    /**
     * Disables the monitor.
     */
    void disable();

    /**
     * Specifies a new time interval to be used when calculating
     * "recent" statistics.
     */
    void setNewTimeWindowForRecentResults(utils::duration_t interval);

    /**
     * Tests whether the monitor is currently enabled.
     */
//     bool isEnabled() const {return _enabled;}

    /**
     * Returns the length of the time window that has been specified
     * for recent results.  (This may be different than the actual
     * length of the recent time window which is affected by the
     * interval of calls to the calculateStatistics() method.  Use
     * a getDuration(RECENT) call to determine the actual recent
     * time window.)
     */
    utils::duration_t getTimeWindowForRecentResults() const {
      return _intervalForRecentStats;
    }

    /**
       Write all our collected statistics into the given Stats struct.
     */
    void getStats(Stats& stats) const;

  private:

    // Prevent copying of the MonitoredQuantity
    MonitoredQuantity(MonitoredQuantity const&);
    MonitoredQuantity& operator=(MonitoredQuantity const&);

    // Helper functions.
    void _reset_accumulators();
    void _reset_results();

    utils::time_point_t _lastCalculationTime;
    long long _workingSampleCount;
    double _workingValueSum;
    double _workingValueSumOfSquares;
    double _workingValueMin;
    double _workingValueMax;
    double _workingLastSampleValue;

    mutable boost::mutex _accumulationMutex;

    int _binCount;
    int _workingBinId;
    std::vector<long long> _binSampleCount;
    std::vector<double> _binValueSum;
    std::vector<double> _binValueSumOfSquares;
    std::vector<double> _binValueMin;
    std::vector<double> _binValueMax;
    std::vector<utils::duration_t> _binDuration;

    long long _fullSampleCount;
    double _fullSampleRate;
    double _fullValueSum;
    double _fullValueSumOfSquares;
    double _fullValueAverage;
    double _fullValueRMS;
    double _fullValueMin;
    double _fullValueMax;
    double _fullValueRate;
    utils::duration_t _fullDuration;

    long long _recentSampleCount;
    double _recentSampleRate;
    double _recentValueSum;
    double _recentValueSumOfSquares;
    double _recentValueAverage;
    double _recentValueRMS;
    double _recentValueMin;
    double _recentValueMax;
    double _recentValueRate;
    utils::duration_t _recentDuration;
    double _lastLatchedSampleValue;
    double _lastLatchedValueRate;

    mutable boost::mutex _resultsMutex;

    bool _enabled;
    utils::duration_t _intervalForRecentStats;  // seconds
  };

  struct MonitoredQuantity::Stats
  {
    long long fullSampleCount;
    double fullSampleRate;
    double fullValueSum;
    double fullValueSumOfSquares;
    double fullValueAverage;
    double fullValueRMS;
    double fullValueMin;
    double fullValueMax;
    double fullValueRate;
    double fullSampleLatency;
    utils::duration_t fullDuration;

    long long recentSampleCount;
    double recentSampleRate;
    double recentValueSum;
    double recentValueSumOfSquares;
    double recentValueAverage;
    double recentValueRMS;
    double recentValueMin;
    double recentValueMax;
    double recentValueRate;
    double recentSampleLatency;
    utils::duration_t recentDuration;
    std::vector<long long> recentBinnedSampleCounts;
    std::vector<double> recentBinnedValueSums;
    std::vector<utils::duration_t> recentBinnedDurations;

    double lastSampleValue;
    double lastValueRate;
    bool   enabled;

    long long getSampleCount(DataSetType t = FULL) const { return t == RECENT ? recentSampleCount : fullSampleCount; }
    double getValueSum(DataSetType t = FULL) const { return t == RECENT ? recentValueSum : fullValueSum; }
    double getValueAverage(DataSetType t = FULL) const { return t == RECENT ? recentValueAverage : fullValueAverage; }
    double getValueRate(DataSetType t = FULL) const { return t== RECENT ? recentValueRate : fullValueRate; }
    double getValueRMS(DataSetType t = FULL) const { return t == RECENT ? recentValueRMS : fullValueRMS; }
    double getValueMin(DataSetType t = FULL) const { return t == RECENT ? recentValueMin : fullValueMin; }
    double getValueMax(DataSetType t = FULL) const { return t == RECENT ? recentValueMax : fullValueMax; }
    utils::duration_t getDuration(DataSetType t = FULL) const { return t == RECENT ? recentDuration : fullDuration; }
    double getSampleRate(DataSetType t = FULL) const { return t == RECENT ? recentSampleRate : fullSampleRate; }
    double getSampleLatency(DataSetType t = FULL) const { double v=getSampleRate(t); return v  ? 1e6/v : INFINITY;}
    double getLastSampleValue() const { return lastSampleValue; }
    double getLastValueRate() const { return lastValueRate; }
    bool   isEnabled() const { return enabled; }
  };

} // namespace stor

#endif // StorageManager_MonitoredQuantity_h



/// emacs configuration
/// Local Variables: -
/// mode: c++ -
/// c-basic-offset: 2 -
/// indent-tabs-mode: nil -
/// End: -
