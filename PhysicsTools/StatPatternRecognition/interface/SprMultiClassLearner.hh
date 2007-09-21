// File and Version Information:
//      $Id: SprMultiClassLearner.hh,v 1.4 2007/05/14 18:08:08 narsky Exp $
//
// Description:
//      Class SprMultiClassLearner :
//          Implements a multi class learning algorithm described
//          in Allwein, Schapire and Singer, 
//          J. of Machine Learning Research 2000,
//          "Reducing multiclass to binary"
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Ilya Narsky                     Original author
//
// Copyright Information:
//      Copyright (C) 2005              California Institute of Technology
//
//------------------------------------------------------------------------
 
#ifndef _SprMultiClassLearner_HH
#define _SprMultiClassLearner_HH

#include "PhysicsTools/StatPatternRecognition/interface/SprAbsMultiClassLearner.hh"
#include "PhysicsTools/StatPatternRecognition/interface/SprTrainedMultiClassLearner.hh"
#include "PhysicsTools/StatPatternRecognition/src/SprMatrix.hh"

#include <vector>
#include <utility>
#include <iostream>
#include <cassert>

class SprAbsFilter;
class SprAbsClassifier;
class SprAbsTrainedClassifier;


class SprMultiClassLearner : public SprAbsMultiClassLearner
{
public:
  enum MultiClassMode { User, OneVsAll, OneVsOne };

  virtual ~SprMultiClassLearner();

  SprMultiClassLearner(SprAbsFilter* data, 
		       SprAbsClassifier* c,
		       const std::vector<int>& classes,
		       const SprMatrix& indicator,
		       MultiClassMode mode=User);

  /*
    Trains classifier on data. Returns true on success, false otherwise.
  */
  bool train(int verbose=0);

  /*
    Reset this classifier to untrained state.
  */
  bool reset();

  /*
    Replace training data.
  */
  bool setData(SprAbsFilter* data);

  /*
    Prints results of training.
  */
  void print(std::ostream& os) const;

  /*
    Make a trained classifier.
  */
  SprTrainedMultiClassLearner* makeTrained() const;

  /*
    Set parameters.
  */
  void setTrained(const SprMatrix& indicator, 
		  const std::vector<int>& classes,
		  const std::vector<
		  std::pair<const SprAbsTrainedClassifier*,bool> >& trained);

  // Show indicator matrix.
  void printIndicatorMatrix(std::ostream& os) const;

private:
  bool setClasses();
  void destroy();

  MultiClassMode mode_;
  SprMatrix indicator_;
  std::vector<int> mapper_;
  SprAbsClassifier* trainable_;
  std::vector<std::pair<const SprAbsTrainedClassifier*,bool> > trained_;
};

#endif
