#ifndef DataFormat_Math_SSERot_H
#define DataFormat_Math_SSERot_H


#include "DataFormats/Math/interface/SSEVec.h"

namespace mathSSE {

  template<typename T>
  struct OldRot { 
    T R11, R12, R13;
    T R21, R22, R23;
    T R31, R32, R33;
  }  __attribute__ ((aligned (16)));
  

  template<typename T>
  struct Rot3 {
    Vec3<T>  axis[3];

    Rot3() {
      axis[0].arr[0]=1;
      axis[1].arr[1]=1;
      axis[2].arr[2]=1;
    }
    
    Rot3( Vec3<T> x,  Vec3<T> y,  Vec3<T> z) {
      axis[0] =x;
      axis[1] =y;
      axis[2] =z;
    }

    Rot3( T xx, T xy, T xz, T yx, T yy, T yz, T zx, T zy, T zz) {
      axis[0].set(xx,xy,xz);
      axis[1].set(yx,yy,yz);
      axis[2].set(zx,zy,zz);
    }

    Rot3 transpose() const {
      return Rot3( axis[0].arr[0], axis[1].arr[0], axis[2].arr[0],
		   axis[0].arr[1], axis[1].arr[1], axis[2].arr[1],
		   axis[0].arr[2], axis[1].arr[2], axis[2].arr[2]
		   );
    }

    // toLocal...
    Vec3<T> rotate(Vec3<T> v) const {
      return transpose().rotateBack(v);
    }


    // toGlobal...
    Vec3<T> rotateBack(Vec3<T> v) const {
      return v.get1(0)*axis[0] +  v.get1(1)*axis[1] + v.get1(2)*axis[2];
    }

    Rot3 rotate(Rot3 const& r) const {
      Rot3 tr = transpose();
      return Rot3(tr.rotateBack(r.axis[0]),tr.rotateBack(r.axis[1]),tr.rotateBack(r.axis[2]));
    }

    Rot3 rotateBack(Rot3 const& r) const {
      return Rot3(rotateBack(r.axis[0]),rotateBack(r.axis[1]),rotateBack(r.axis[2]));
    }


  };

  typedef Rot3<float> Rot3F;

  typedef Rot3<double> Rot3D;
  
}

template<typename T>
inline  mathSSE::Rot3<T>  operator *(mathSSE::Rot3<T> const & rh, mathSSE::Rot3<T> const & lh) {
  //   return Rot3(lh.rotateBack(rh.axis[0]),lh.rotateBack(rh.axis[1]),lh.rotateBack(rh.axis[2]));
  return lh.rotateBack(rh);
}

#include <iosfwd>
std::ostream & operator<<(std::ostream & out, mathSSE::Rot3F const & v);
std::ostream & operator<<(std::ostream & out, mathSSE::Rot3D const & v);


#endif //  DataFormat_Math_SSERot_H
