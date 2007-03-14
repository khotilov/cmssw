#ifndef MinL3Algorithm_H
#define MinL3Algorithm_H

/** \class MinL3Algorithm
 *  Implementation of the L3 Collaboration algorithm to solve a system Ax = B
 *  by minimization of |Ax-B| using an iterative linear approach
 *  This class is specific for the ECAL calibration
 *
 * 13.03.2007: R.Ofierzynski
 *  - implemented event weighting
 *
 * $Date: 2007/03/13 14:33:42 $
 * $Revision: 1.4 $
 * \author R.Ofierzynski, CERN
 */

#include <vector>
#include <iostream>

using namespace std;


class MinL3Algorithm
{
public:
  /// Default constructor
  /// kweight_ = event weight, squareMode_ = side length of the cluster square
  MinL3Algorithm(float kweight_ = 0., int squareMode_ = 5, int mineta_ = 1, int maxeta_ = 85, int minphi_ = 1, int maxphi_ = 20);

  /// Destructor
  ~MinL3Algorithm();

  /// method doing the full calibration running nIter number of times, recalibrating the event matrix after each iteration with the new solution
  /// returns the vector of calibration coefficients built from all iteration solutions
  /// >> to be used also as recipe on how to use the calibration methods one-by-one <<
  vector<float> iterate(const vector<vector<float> >& eventMatrix, const vector<int>& VmaxCeta, const vector<int>& VmaxCphi, const vector<float>& energyVector, const int& nIter, const bool& normalizeFlag = false);


  /// add event to the calculation of the calibration vector
  void addEvent(const vector<float>& eventSquare, const int& maxCeta, const int& maxCphi, const float& energy);

  /// recalibrate before next iteration: give previous solution vector as argument
  vector<float> recalibrateEvent(const vector<float>& eventSquare, const int& maxCeta, const int& maxCphi, const vector<float>& recalibrateVector); 

  /// get the solution at the end of the calibration
  vector<float> getSolution(bool resetsolution=true);

  /// reset for new iteration
  void resetSolution(); 

  /// method to translate from square indices to region indices
  int indexSqr2Reg(const int& sqrIndex, const int& maxCeta, const int& maxCphi);


private:

  float kweight;
  int squareMode;
  int mineta, maxeta, minphi, maxphi;
  int countEvents;
  int Nchannels, Nxtals;
  vector<float> wsum;
  vector<float> Ewsum;

};

#endif // MinL3Algorithm_H
