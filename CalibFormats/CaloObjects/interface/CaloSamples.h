#ifndef CALOSAMPLES_H
#define CALOSAMPLES_H 1

#include <ostream>
#include "DataFormats/DetId/interface/DetId.h"

/** \class CaloSamples
    
Class which represents the charge/voltage measurements of an event/channel
with the ADC decoding performed.

$Date: 2006/01/05 16:26:05 $
$Revision: 1.3 $
*/
class CaloSamples {
public:
  CaloSamples();
  explicit CaloSamples(const DetId& id, int size);
  
  /// get the (generic) id
  DetId id() const { return id_; }

  /// get the size
  int size() const { return size_; }
  /// mutable operator to access samples
  double& operator[](int i) { return data_[i]; }
  /// const operator to access samples
  double operator[](int i) const { return data_[i]; }

  /// access presample information
  int presamples() const { return presamples_; }
  /// set presample information
  void setPresamples(int pre);

  /// multiply each item by this value
  CaloSamples& scale(double value);
  /// scale all samples
  CaloSamples& operator*=(double value) { return scale(value); }

  /// add a value to all samples
  CaloSamples& operator+=(double value);

  void setDetId( DetId detId ) { id_ = detId ; }

  void setSize( unsigned int size ) { size_ = size ; }

  bool isBlank() const ; // are the samples blank (zero?)

  void setBlank() ; // keep id, presamples, size but zero out data

  static const int MAXSAMPLES=10;
private:
  DetId id_;
  double data_[MAXSAMPLES]; // 
  int size_, presamples_;
};

std::ostream& operator<<(std::ostream& s, const CaloSamples& samps);

#endif
