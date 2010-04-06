#ifndef Grid3D_h
#define Grid3D_h

/** \class Grid3D
 *
 *  Implementation of a 3D regular grid.
 *
*  \author T. Todorov
 */

#include "DataFormats/GeometryVector/interface/Basic3DVector.h"
#include "DataFormats/Math/interface/SSEVec.h"
#include "MagneticField/Interpolation/src/Grid1D.h"
#include <vector>

class Grid3D {
public:

  typedef Basic3DVector<float>  ReturnType;
  // typedef double   Scalar;
  typedef mathSSE::Vec3F ValueType;
  typedef float   Scalar;
 
  Grid3D() {}

  Grid3D( const Grid1D& ga, const Grid1D& gb, const Grid1D& gc,
	  std::vector<ValueType>& data) : 
    grida_(ga), gridb_(gb), gridc_(gc) {
     data_.swap(data);
     stride1_ = gridb_.nodes() * gridc_.nodes();
     stride2_ = gridc_.nodes();
  }

  Grid3D( const Grid1D& ga, const Grid1D& gb, const Grid1D& gc,
	  std::vector<ReturnType>& data);


  int index(int i, int j, int k) const {return i*stride1_ + j*stride2_ + k;}
  int stride1() const { return stride1_;}
  int stride2() const { return stride2_;}
  int stride3() const { return 1;}
  const ValueType& operator()(int i) const {
    return data_[i];
  }

  ReturnType operator()(int i, int j, int k) const {
    const ValueType& result = (*this)(index(i,j,k));
    return ReturnType(result);
  }

  const Grid1D& grida() const {return grida_;}
  const Grid1D& gridb() const {return gridb_;}
  const Grid1D& gridc() const {return gridc_;}

  const std::vector<ValueType>& data() const {return data_;}

  void dump() const;

private:

  Grid1D grida_;
  Grid1D gridb_;
  Grid1D gridc_;

  std::vector<ValueType> data_;

  int stride1_;
  int stride2_;


};

#endif
