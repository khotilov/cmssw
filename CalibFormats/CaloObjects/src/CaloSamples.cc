#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include <math.h>

CaloSamples::CaloSamples() : id_(), size_(0), presamples_(0) { setBlank() ; }

CaloSamples::CaloSamples(const DetId& id, int size) :
   id_          ( id   ) , 
   size_        ( size ) , 
   presamples_  ( 0    ) 
{
   setBlank() ;
}

void
CaloSamples::setPresamples( int pre ) 
{
   presamples_ = pre ;
}

CaloSamples& 
CaloSamples::scale( double value )
{
   for (int i=0; i<MAXSAMPLES; i++) data_[i]*=value;
   return (*this);
}

CaloSamples& 
CaloSamples::operator+=(double value) 
{  
   for (int i=0; i<MAXSAMPLES; i++) data_[i]+=value;
   return (*this);
}

bool 
CaloSamples::isBlank() const // are the samples blank (zero?)
{
   for( int i ( 0 ) ; i != MAXSAMPLES ; ++i )
   {
      if( 1.e-6 < fabs( data_[i] ) ) return false ;
   }
   return true ;
}

void 
CaloSamples::setBlank() // keep id, presamples, size but zero out data
{
   std::fill( data_ , data_ + MAXSAMPLES, (double)0.0 ) ;
}

std::ostream& operator<<(std::ostream& s, const CaloSamples& samples) {
  s << "DetId=" << samples.id().rawId();
  s << ", "<<  samples.size() << "samples" << std::endl;
  for (int i=0; i<samples.size(); i++)
    s << i << ":" << samples[i] << std::endl;
  return s;
}
