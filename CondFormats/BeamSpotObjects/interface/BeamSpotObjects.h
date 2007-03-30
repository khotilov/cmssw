#ifndef BEAMSPOTOBJECTS_H
#define BEAMSPOTOBJECTS_H
/** \class BeamSpotObjects
 *
 * Reconstructed beam spot object. It provides position, error, and
 * width of the beam position.
 *
 * \author Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)
 *
 * \version $Id: BeamSpotObjects.h,v 1.2 2007/01/22 04:48:40 yumiceva Exp $
 *
 */

#include <math.h>

class BeamSpotObjects {
	
  public:

	/// default constructor
	BeamSpotObjects(){}
	
	virtual ~BeamSpotObjects(){}

	/// set XYZ position
	void SetPosition( double x, double y, double z) { 
		position_[0] = x;
		position_[1] = y;
		position_[2] = z;
	};
	/// set sigma Z, RMS bunch length
	void SetSigmaZ(double val) { sigmaZ_ = val; }
	/// set dxdz slope, crossing angle
	void Setdxdz(double val) { dxdz_ = val; }
	/// set dydz slope, crossing angle in XZ
	void Setdydz(double val) { dydz_ = val; }
	/// set average transverse beam width in YZ
	void SetBeamWidth(double val) { beamwidth_ = val; }
	/// set i,j element of the full covariance matrix 7x7
	void SetCovariance(int i, int j, double val) {
		covariance_[i][j] = val;
	}

	/// get X beam position
	double GetX() const { return position_[0]; }
	/// get Y beam position
	double GetY() const { return position_[1]; }
	/// get Z beam position
	double GetZ() const { return position_[2]; }
	/// get sigma Z, RMS bunch length
	double GetSigmaZ() const { return sigmaZ_; }
	/// get average transverse beam width
	double GetBeamWidth() const { return beamwidth_; }
	/// get dxdz slope, crossing angle in XZ
	double Getdxdz() const { return dxdz_; }
	/// get dydz slope, crossing angle in YZ
	double Getdydz() const { return dydz_; }
	/// get i,j element of the full covariance matrix 7x7
	double GetCovariance(int i, int j) const { return covariance_[i][j]; }
	/// get X beam position Error
	double GetXError() const { return sqrt(covariance_[0][0]); }
	/// get Y beam position Error
	double GetYError() const { return sqrt(covariance_[1][1]); }
	/// get Z beam position Error
	double GetZError() const { return sqrt(covariance_[2][2]); }
	/// get sigma Z, RMS bunch length Error
	double GetSigmaZError() const { return sqrt(covariance_[3][3]); }
	/// get average transverse beam width
	double GetBeamWidthError() const { return sqrt(covariance_[6][6]); }
	/// get dxdz slope, crossing angle in XZ Error
	double GetdxdzError() const { return sqrt(covariance_[4][4]); }
	/// get dydz slope, crossing angle in YZ Error
	double GetdydzError() const { return sqrt(covariance_[5][5]); }
	
  private:

	double position_[3];
	double sigmaZ_;
	double beamwidth_;
	double dxdz_;
	double dydz_;
	double covariance_[7][7];
	
};
#endif
