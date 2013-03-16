#ifndef GeometryVector_oldBasic3DVector_h
#define GeometryVector_oldBasic3DVector_h
#if ( defined(IN_DICTBUILD) || defined(__CINT__) )  && !defined(__REFLEX__)
#define __REFLEX__
#endif
#include "DataFormats/GeometryVector/interface/Basic2DVector.h"
#include "DataFormats/GeometryVector/interface/Theta.h"
#include "DataFormats/GeometryVector/interface/Phi.h"
#include "DataFormats/GeometryVector/interface/PreciseFloatType.h"
#include "DataFormats/GeometryVector/interface/CoordinateSets.h"
#ifndef __REFLEX__ 
#include "DataFormats/Math/interface/SIMDVec.h"
#endif
#include <iosfwd>
#include <cmath>

namespace detailsBasic3DVector {
  inline float __attribute__((always_inline)) __attribute__ ((pure))
  eta(float x, float y, float z) { float t(z/std::sqrt(x*x+y*y)); return ::asinhf(t);} 
  inline double __attribute__((always_inline)) __attribute__ ((pure))
  eta(double x, double y, double z) { double t(z/std::sqrt(x*x+y*y)); return ::asinh(t);} 
  inline long double __attribute__((always_inline)) __attribute__ ((pure))
  eta(long double x, long double y, long double z) { long double t(z/std::sqrt(x*x+y*y)); return ::asinhl(t);} 
}


template < typename T> 
class Basic3DVector {
public:
  typedef  Basic3DVector<T> MathVector;


  typedef T                                   ScalarType;
  typedef Geom::Cylindrical2Cartesian<T>      Cylindrical;
  typedef Geom::Spherical2Cartesian<T>        Spherical;
  typedef Spherical                           Polar; // synonym

#ifdef __REFLEX__
  Basic3DVector() { for(int i=0;i!=4; ++i) theV[i]=0;}
  Basic3DVector( const T& x, const T& y, const T& z, const T& w=0) 
  { theV[0]=x; theV[1]=y; theV[2]=z; theV[3]=w;}
 Basic3DVector( const Basic3DVector & p) 
  { theV[0]=p.x(); theV[1]=p.y(); theV[2]=p.z();theV[3]=p.w();}
  template <class OtherPoint> 
  explicit Basic3DVector( const OtherPoint& p) 
  { theV[0]=p.x(); theV[1]=p.y(); theV[2]=p.z();theV[3]=0;}

#else
  /** default constructor uses default constructor of T to initialize the 
   *  components. For built-in floating-point types this means initialization 
   * to zero??? (force init to 0)
   */
  Basic3DVector() : theV{0,0, 0, 0} {}

  /// construct from cartesian coordinates
  Basic3DVector( const T& x, const T& y, const T& z, const T& w=0) : 
    theV{x,y,z,w} {}


  /// Copy constructor from same type. Should not be needed but for gcc bug 12685
  Basic3DVector( const Basic3DVector & p) : 
    theV[0]{p.x(),p.y(),p.z(),p.w()} {}


  /// Copy constructor and implicit conversion from Basic3DVector of different precision
  template <class U>
  Basic3DVector( const Basic3DVector<U> & p) : 
    theV[0]{p.x(),p.y(),p.z(),p.w()} {}

  /// constructor from 2D vector (X and Y from 2D vector, z set to zero)
  Basic3DVector( const Basic2DVector<T> & p) : 
    theV{p.x()), p.y(),0,0} {}

  /** Explicit constructor from other (possibly unrelated) vector classes 
   *  The only constraint on the argument type is that it has methods
   *  x(), y() and z(), and that these methods return a type convertible to T.
   *  Examples of use are
   *   <BR> construction from a Basic3DVector with different precision
   *   <BR> construction from a Hep3Vector
   *   <BR> construction from a coordinate system converter 
   */
  template <class OtherPoint> 
  explicit Basic3DVector( const OtherPoint& p) : 
    theV{p.x(),p.y(),p.z(),0} {}
#endif


#if  defined(USE_EXTVECT)
  template<typename U>
  Basic3DVector(Vec4<U> const& iv) :
    theV{iv[0], iv[1], iv[2], 0} {}
#elif  defined(USE_SSEVECT)
  // constructor from Vec4
  template<typename U>
  Basic3DVector(mathSSE::Vec4<U> const& iv) :
    theV{iv.arr[0], iv.arr[1], iv.arr[2], 0} {}
#endif  


  /** Deprecated construct from polar coordinates, use 
   *  <BR> Basic3DVector<T>( Basic3DVector<T>::Polar( theta, phi, r))
   *  instead. 
   */
  template <typename U>
  Basic3DVector( const Geom::Theta<U>& theta, 
		 const Geom::Phi<U>& phi, const T& r) {
    Polar p( theta.value(), phi.value(), r);
    theV[0] = p.x(); theV[1] = p.y(); theV[2] = p.z();
  }

  T operator[](int i) const { return theV[i];}
  T & operator[](int i) { return theV[i];}

  /// Cartesian x coordinate
  T x() const { return theV[0];}

  /// Cartesian y coordinate
  T y() const { return theV[1];}

  /// Cartesian z coordinate
  T z() const { return theV[2];}

  T w() const { return theV[3];}


  Basic2DVector<T> xy() const { return  Basic2DVector<T>(theV[0],theV[1]);}


  // equality
  bool operator==(const Basic3DVector& rh) const {
    return x()==rh.x() && y()==rh.y() && z()==rh.z();
  }

  /// The vector magnitude squared. Equivalent to vec.dot(vec)
  T mag2() const { return  x()*x() + y()*y()+z()*z();}

  /// The vector magnitude. Equivalent to sqrt(vec.mag2())
  T mag() const  { return std::sqrt( mag2());}

  /// Squared magnitude of transverse component 
  T perp2() const { return x()*x() + y()*y();}

  /// Magnitude of transverse component 
  T perp() const { return std::sqrt( perp2());}

  /// Another name for perp()
  T transverse() const { return perp();}

  /** Azimuthal angle. The value is returned in radians, in the range (-pi,pi].
   *  Same precision as the system atan2(x,y) function.
   *  The return type is Geom::Phi<T>, see it's documentation.
   */ 
  T barePhi() const {return std::atan2(y(),x());}
  Geom::Phi<T> phi() const {return Geom::Phi<T>(barePhi());}

  /** Polar angle. The value is returned in radians, in the range [0,pi]
   *  Same precision as the system atan2(x,y) function.
   *  The return type is Geom::Phi<T>, see it's documentation.
   */ 
  T bareTheta() const {return std::atan2(perp(),z());}
  Geom::Theta<T> theta() const {return Geom::Theta<T>(std::atan2(perp(),z()));}

  /** Pseudorapidity. 
   *  Does not check for zero transverse component; in this case the behavior 
   *  is as for divide-by zero, i.e. system-dependent.
   */
  // T eta() const { return -log( tan( theta()/2.));} 
  T eta() const { return detailsBasic3DVector::eta(x(),y(),z());} // correct 

  /** Unit vector parallel to this.
   *  If mag() is zero, a zero vector is returned.
   */
  Basic3DVector unit() const {
    T my_mag = mag2();
    if (my_mag==0) return *this;
    my_mag = T(1)/std::sqrt(my_mag);
    return *this * my_mag;
  }

  /** Operator += with a Basic3DVector of possibly different precision.
   */
  template <class U> 
  Basic3DVector& operator+= ( const Basic3DVector<U>& p) {
    theV[0] += p.x();
    theV[1] += p.y();
    theV[2] += p.z();
    theV[3] += p.w();
    return *this;
  } 

  /** Operator -= with a Basic3DVector of possibly different precision.
   */
  template <class U> 
  Basic3DVector& operator-= ( const Basic3DVector<U>& p) {
    theV[0] -= p.x();
    theV[1] -= p.y();
    theV[2] -= p.z();
    theV[3] -= p.w();
    return *this;
  } 

  /// Unary minus, returns a vector with components (-x(),-y(),-z())
  Basic3DVector operator-() const { return Basic3DVector(-x(),-y(),-z());}

  /// Scaling by a scalar value (multiplication)
  Basic3DVector& operator*= ( T t) {
    theV[0] *= t;
    theV[1] *= t;
    theV[2] *= t;
    theV[3] *= t;;
    return *this;
  } 

  /// Scaling by a scalar value (division)
  Basic3DVector& operator/= ( T t) {
    t = T(1)/t;
    theV[0] *= t;
    theV[1] *= t;   
    theV[2] *= t;
    theV[3] *= t;;
    return *this;
  } 

  /// Scalar product, or "dot" product, with a vector of same type.
  T dot( const Basic3DVector& v) const { 
    return x()*v.x() + y()*v.y() + z()*v.z();
  }

  /** Scalar (or dot) product with a vector of different precision.
   *  The product is computed without loss of precision. The type
   *  of the returned scalar is the more precise of the scalar types 
   *  of the two vectors.
   */
  template <class U> 
  typename PreciseFloatType<T,U>::Type dot( const Basic3DVector<U>& v) const { 
    return x()*v.x() + y()*v.y() + z()*v.z();
  }

  /// Vector product, or "cross" product, with a vector of same type.
  Basic3DVector cross( const Basic3DVector& v) const {
    return Basic3DVector( y()*v.z() - v.y()*z(), 
			  z()*v.x() - v.z()*x(), 
			  x()*v.y() - v.x()*y());
  }


  /** Vector (or cross) product with a vector of different precision.
   *  The product is computed without loss of precision. The type
   *  of the returned vector is the more precise of the types 
   *  of the two vectors.   
   */
  template <class U> 
  Basic3DVector<typename PreciseFloatType<T,U>::Type> 
  cross( const Basic3DVector<U>& v) const {
    return Basic3DVector<typename PreciseFloatType<T,U>::Type>( y()*v.z() - v.y()*z(), 
								z()*v.x() - v.z()*x(), 
								x()*v.y() - v.x()*y());
  }

private:
  T theV[4];

}  
#ifndef __CINT__
__attribute__ ((aligned (16)))
#endif
;


namespace geometryDetails {
  std::ostream & print3D(std::ostream& s, double x, double y, double z);
}

/// simple text output to standard streams
template <class T>
inline std::ostream & operator<<( std::ostream& s, const Basic3DVector<T>& v) {
  return geometryDetails::print3D(s, v.x(),v.y(), v.z());
}


/// vector sum and subtraction of vectors of possibly different precision
template <class T, class U>
inline Basic3DVector<typename PreciseFloatType<T,U>::Type>
operator+( const Basic3DVector<T>& a, const Basic3DVector<U>& b) {
  typedef Basic3DVector<typename PreciseFloatType<T,U>::Type> RT;
  return RT(a.x()+b.x(), a.y()+b.y(), a.z()+b.z(), a.w()+b.w());
}

template <class T, class U>
inline Basic3DVector<typename PreciseFloatType<T,U>::Type>
operator-( const Basic3DVector<T>& a, const Basic3DVector<U>& b) {
  typedef Basic3DVector<typename PreciseFloatType<T,U>::Type> RT;
  return RT(a.x()-b.x(), a.y()-b.y(), a.z()-b.z(), a.w()-b.w());
}

/// scalar product of vectors of same precision
template <class T>
inline T operator*( const Basic3DVector<T>& v1, const Basic3DVector<T>& v2) {
  return v1.dot(v2);
}

/// scalar product of vectors of different precision
template <class T, class U>
inline typename PreciseFloatType<T,U>::Type operator*( const Basic3DVector<T>& v1, 
						       const Basic3DVector<U>& v2) {
  return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
}

/** Multiplication by scalar, does not change the precision of the vector.
 *  The return type is the same as the type of the vector argument.
 */
template <class T>
inline Basic3DVector<T> operator*( const Basic3DVector<T>& v, T t) {
  return Basic3DVector<T>(v.x()*t, v.y()*t, v.z()*t, v.w()*t);
}

/// Same as operator*( Vector, Scalar)
template <class T>
inline Basic3DVector<T> operator*(T t, const Basic3DVector<T>& v) {
  return Basic3DVector<T>(v.x()*t, v.y()*t, v.z()*t, v.w()*t);
}

template <class T, typename S>
inline Basic3DVector<T> operator*(S t,  const Basic3DVector<T>& v) {
  return static_cast<T>(t)*v;
}

template <class T, typename S>
inline Basic3DVector<T> operator*(const Basic3DVector<T>& v, S t) {
  return static_cast<T>(t)*v;
}


/** Division by scalar, does not change the precision of the vector.
 *  The return type is the same as the type of the vector argument.
 */
template <class T, typename S>
inline Basic3DVector<T> operator/( const Basic3DVector<T>& v, S s) {
  T t = T(1)/s;
  return v*t;
}


typedef Basic3DVector<float> Basic3DVectorF;
typedef Basic3DVector<double> Basic3DVectorD;
typedef Basic3DVector<long double> Basic3DVectorLD;


#endif // GeometryVector_Basic3DVector_h


