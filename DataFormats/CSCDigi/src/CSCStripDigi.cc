/** \file
 * 
 *  $Date: 2006/03/03 13:50:52 $
 *  $Revision: 1.6 $
 *
 * \author M.Schmitt, Northwestern
 */
#include <DataFormats/CSCDigi/interface/CSCStripDigi.h>
#include <iostream>
#include <bitset>
#include <vector>

using namespace std;

// Constructors
CSCStripDigi::CSCStripDigi (int istrip, vector<int> vADCCounts){
  strip = istrip;
  ADCCounts = vADCCounts;
}

CSCStripDigi::CSCStripDigi (){
  vector<int> ZeroCounts(8,0);
  strip = 0;
  ADCCounts = ZeroCounts;
}

// Copy constructor
CSCStripDigi::CSCStripDigi(const CSCStripDigi& digi) {
 ADCCounts = digi.getADCCounts();
 strip = digi.getStrip();
}

// Assignment
CSCStripDigi& 
CSCStripDigi::operator=(const CSCStripDigi& digi){
  ADCCounts = digi.getADCCounts();
  strip = digi.getStrip();
  return *this;
}

// Comparison
bool
CSCStripDigi::operator == (const CSCStripDigi& digi) const {
  if ( getStrip() != digi.getStrip() ) return false;
  if ( getADCCounts().size() != digi.getADCCounts().size() ) return false;
  if ( getADCCounts() != digi.getADCCounts() ) return false;
  return true;
}

// Getters
int CSCStripDigi::getStrip() const { return strip; }
std::vector<int> CSCStripDigi::getADCCounts() const { return ADCCounts; }

// Setters
void CSCStripDigi::setStrip(int istrip) {
  strip = istrip;
}
void CSCStripDigi::setADCCounts(vector<int>vADCCounts) {
  bool badVal = false;
  for (int i=0; i<(int)vADCCounts.size(); i++) {
    if (vADCCounts[i] < 1) badVal = true;
  }
  if ( !badVal ) {
    ADCCounts = vADCCounts;
  } else {
    vector<int> ZeroCounts(8,0);
    ADCCounts = ZeroCounts;
  }
}

// Debug
void
CSCStripDigi::print() const {
  cout << "CSC Strip: " << getStrip() << "  ADC Counts: ";
  for (int i=0; i<(int)getADCCounts().size(); i++) {cout << getADCCounts()[i] << " ";}
  cout << "\n";
}




