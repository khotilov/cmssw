#ifndef TMom_H
#define TMom_H

#include<vector>

class TMom
{

 private:

  int nevt;
  double mean;
  double mean2;
  double mean3;
  double sum;
  double sum2;
  double sum3;
  double rms;
  double M3;
  double peak;
  double min;
  double max;
  int bing[101];
  std::vector<double> _cutLow;
  std::vector<double> _cutHigh;
  std::vector<double> _ampl;

  void init(double,double);
  void init(std::vector<double>,std::vector<double>);

 public:


  int _dimCut;


  // Default Constructor, mainly for Root
  TMom();

  // Default Constructor
  TMom(double, double);

  // Default Constructor
  TMom(std::vector<double>,std::vector<double>);

  // Destructor: Does nothing
  virtual ~TMom();

  void setCut(double, double);
  void setCut(std::vector<double>,std::vector<double>);
  void addEntry(double val);
  void addEntry(double val, std::vector<double> valcut);
  double getMean();
  double getMean2();
  double getMean3();
  int    getNevt();
  double getRMS();
  double getM3();
  double getMin();
  double getMax();
  std::vector<double> getPeak();

  //  ClassDef(TMom,1)
};

#endif
