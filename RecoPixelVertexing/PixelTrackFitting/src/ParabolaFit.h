#ifndef ParabolaFit_H
#define ParabolaFit_H

#include <vector>
using namespace std;

/** The parabola fit 
    y = parA + parB * x + parC * x*x
    see: R.L. Gluckstern, NIM 24 (1963) 381 
 */

class ParabolaFit {
public:
  struct Result { double parA, parB, parC; 
                  double varAA, varBB, varCC, varAB, varAC, varBC; } ;
  ParabolaFit() : doErr(true), hasValues(false), hasErrors(false) { }

  void addPoint(double x, double y, double weight=1.);

  void skipErrorCalculationByDefault() { doErr = false; }

  const Result & result(bool doErrors) const;

  double parA() const  { if(!hasValues) result(doErr); return theResult.parA; }
  double parB() const  { if(!hasValues) result(doErr); return theResult.parB; }
  double parC() const  { if(!hasValues) result(doErr); return theResult.parC; }
  double varAA() const { if(!hasErrors) result(true); return theResult.varAA; }
  double varBB() const { if(!hasErrors) result(true); return theResult.varBB; }
  double varCC() const { if(!hasErrors) result(true); return theResult.varCC; }
  double varAB() const { if(!hasErrors) result(true); return theResult.varAB; }
  double varAC() const { if(!hasErrors) result(true); return theResult.varAC; }
  double varBC() const { if(!hasErrors) result(true); return theResult.varBC; }

  double chi2() const;
  int    dof() const;

private:
  struct Column { double r1; double r2; double r3; };
  double det(const Column & c1, const Column & c2, const Column & c3) const;
  double fun(double x) const;

private:
  struct Point { double x; double y; double w; }; 
  vector<Point> points;
  bool doErr;
  mutable Result theResult;  mutable bool hasValues, hasErrors;
};

#endif
