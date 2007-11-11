///////////////////////////////////////////////////////////////////////////////
// File: ClusterFP420.cc
// Date: 12.2006
// Description: ClusterFP420 for FP420
// Modifications: 
///////////////////////////////////////////////////////////////////////////////
#include "DataFormats/FP420Cluster/interface/ClusterFP420.h"
#include "DataFormats/FP420Digi/interface/HDigiFP420.h"
#include <iostream>
#include <cmath>
using namespace std;

//#define mydigidebug10
//#define mydigidebug11

//static float const_ps1s2[3] =  {0.050,.0125,0.005};// pitch, sigma_1channel, sigma_2channels
static float const_ps1s2[3] =  {0.050,.0139,0.0045};// pitch, sigma_1channel, sigma_2channels - Narrow
static float constW_ps1s2[3] = {0.400,.0920,0.0280};// pitch, sigma_1channel, sigma_2channels - Wide

// sense of zside here is X or Y type planes. Now we are working with X only, i.e. zside=2
ClusterFP420::ClusterFP420( unsigned int detid, unsigned int zside, const HDigiFP420Range& range, 
			    float & cog ,float & err ) :
  detId_(detid), zside_(zside), firstStrip_(range.first->strip())
{
  // For the range of strips in cluster assign adc(,its numbers i->strip()) calculate cog... 
  // strip   :    0-400 or 0-250
#ifdef mydigidebug11
   std::cout << "===================================== firstStrip = " << firstStrip_ << std::endl;
   std::cout << "==range.first->strip() = " << range.first->strip() << std::endl;
   std::cout << "==range.second->strip() = " << range.second->strip() << std::endl;
#endif

  amplitudes_.reserve( range.second - range.first);


  int striprange = 0;
  float sumx = 0.;
  float sumxx = 0.;
  float sumy = 0.;
  float sumyy = 0.;
  int suma = 0;

//  int lastStrip = -1;
  for (HDigiFP420Iter i=range.first; i!=range.second; i++) {
    striprange++;
#ifdef mydigidebug11
   std::cout << " current striprange = " << striprange << std::endl;
#endif
   /*
    /// check if digis consecutive: put amplitude=0 for 
      if (i!=ibeg && (difNarr(zside_,i,i-1) > 1 || difWide(zside_,i,i-1) > 1)   ){
    if (lastStrip>0 && i->strip() != lastStrip + 1) {
                 for (int j=0; j < i->strip()-(lastStrip+1); j++) {
		   amplitudes_.push_back( 0);
		 }
    }
    lastStrip = i->strip();
*/
    short amp = i->adc();       // FIXME: gain correction here

#ifdef mydigidebug11
   std::cout << " current strip = " << i->strip() << "  amp = " << amp << std::endl;
#endif

    amplitudes_.push_back(amp);
    if(zside_ == 1) {
      sumx += i->stripH()*amp;
      sumy += i->stripHW()*amp;
      suma += amp;
      sumxx += (i->stripH()) * (i->stripH()) * amp;
      sumyy += (i->stripHW()) * (i->stripHW()) * amp;
    }
    else if(zside_ == 2) {
      sumx += i->stripV()*amp;
      sumy += i->stripVW()*amp;
      suma += amp;
      sumxx += (i->stripV()) * (i->stripV()) * amp;
      sumyy += (i->stripVW()) * (i->stripVW()) * amp;
    }
    else {
      std::cout << " ClusterFP420: wrong zside_ = " << zside_ << std::endl;
    }

#ifdef mydigidebug11
   std::cout << " current sumx = " << sumx << std::endl;
   std::cout << " current sumy = " << sumy << std::endl;
   std::cout << " current suma = " << suma << std::endl;
   std::cout << " current barycenter = " << (sumx / static_cast<float>(suma) )  << std::endl;
   std::cout << " current barycenterW= " << (sumy / static_cast<float>(suma) )  << std::endl;
#endif
  } //for


  if(suma != 0) {
    barycenter_ = sumx / static_cast<float>(suma) ;
    barycerror_ = sumxx / static_cast<float>(suma) ;
    barycerror_ = fabs(barycerror_ - barycenter_*barycenter_) ;
#ifdef mydigidebug11
    std::cout << "barycerror_ = " << barycerror_ << "barycenter_ = " << barycenter_ << std::endl;
#endif
    barycenterW_ = sumy / static_cast<float>(suma) ;
    barycerrorW_ = sumyy / static_cast<float>(suma) ;
    barycerrorW_ = fabs(barycerrorW_ - barycenterW_*barycenterW_) ;
#ifdef mydigidebug11
    std::cout << "barycerrorW_ = " << barycerrorW_ << "barycenterW_ = " << barycenterW_ << std::endl;
#endif
  }
  else{
    barycenter_ = 1000000. ;
    barycerror_ = 1000000. ;
    barycenterW_ = 1000000. ;
    barycerrorW_ = 1000000. ;
  }

  /** The barycenter of the cluster, not corrected for Lorentz shift;
   *  it can means that should not be used as position estimate for tracking.
   */
  cog = barycenter_;// cog for Narrow pixels only

#ifdef mydigidebug11
   std::cout << "AT end: barycenter_ = " << barycenter_ << std::endl;
   std::cout << "AT end:  striprange = " << striprange << std::endl;
#endif





   /*

  float sumx0 = 0.;
  float sumxx = 0.;
  lastStrip = -1;
  for (HDigiFP420Iter i=range.first; i!=range.second; i++) {
#ifdef mydigidebug11
   std::cout << " current striprange = " << striprange << std::endl;
#endif
    /// check if digis consecutive
    if (lastStrip>0 && i->strip() != lastStrip + 1) {
                 for (int j=0; j < i->strip()-(lastStrip+1); j++) {
		   amplitudes_.push_back( 0);
		 }
    }
    lastStrip = i->strip();

    short amp = i->adc();       // FIXME: gain correction here

#ifdef mydigidebug11
   std::cout << " current strip = " << i->strip() << "  amp = " << amp << std::endl;
#endif

    sumx0 += (i->strip()-cog)*amp;
    sumxx += (i->strip()-cog) * (i->strip()-cog) * amp;


#ifdef mydigidebug11
   std::cout << " 2 current sumx0 = " << sumx0 << std::endl;
   std::cout << " 2 current sumxx = " << sumxx << std::endl;
#endif
  } //for


  if(suma != 0) {
    sumx0 = sumx0 / static_cast<float>(suma) ;
    sumxx = sumxx / static_cast<float>(suma);
    
    //barycerror_ = fabs(sumxx - sumx0*sumx0) ;

    //barycerror_ = (sumxx - sumx0*sumx0) ;
    //barycerror_ *= barycerror_ ;

      barycerror_ = sumxx ;

  }
  else{
    barycerror_ = 1000000. ;
  }

*/

#ifdef mydigidebug10
   std::cout << "pitchcommon = " << const_ps1s2[0] << " sigma1= " << const_ps1s2[1]  << " sigma2= " << const_ps1s2[2]  << std::endl;
#endif

   //
  if(barycerror_ == 0.0) {
    barycerror_ = const_ps1s2[1]/const_ps1s2[0];// 
  }
  else{
    barycerror_ = const_ps1s2[2]/const_ps1s2[0];//  
  }
    barycerror_ *= barycerror_;
   //
  if(barycerrorW_ == 0.0) {
    barycerrorW_ = constW_ps1s2[1]/constW_ps1s2[0];// 
  }
  else{
    barycerrorW_ = constW_ps1s2[2]/constW_ps1s2[0];// 
  }
    barycerrorW_ *= barycerrorW_;
   //

#ifdef mydigidebug11
   std::cout << "barycerror_ = " << barycerror_ << "barycerrorW_ = " << barycerrorW_ << std::endl;
#endif
         
   // change by hands:

	// number of station
	int  mysn0 = 3;
	
	// number of planes 
	int  mypn0 = 5;







  // unpack from detId_:
 	int  sScale = 2*mypn0;
	//	int  zScale=2;
	int  sector = (detId_-1)/sScale + 1 ;
	//	int  zmodule = (detId_ - (sector - 1)*sScale - 1) /zScale + 1 ;
	//////	int  zside = detId_ - (sector - 1)*sScale - (zmodule - 1)*zScale ;
	float a = 0.00001;


	if(mysn0 == 2) {
	  if(sector==2) {
	    a = 0.0026+((0.0075-0.0026)/7.)*(mypn0-2); // 8 m 
	      }
	}
	else if(mysn0 == 3) {
	  if(sector==2) {
	    a = 0.0012+((0.0036-0.0012)/7.)*(mypn0-2); // 4 m 
	      }
	  else if(sector==3) {
	    a = 0.0026+((0.0075-0.0026)/7.)*(mypn0-2); // 8 m 
	      }
	}
	else if(mysn0 == 4) {
	  if(sector==2) {
	    a = 0.0009+((0.0024-0.0009)/7.)*(mypn0-2); // 2.7 m 
	      }
	  else if(sector==3) {
	    a = 0.0018+((0.0050-0.0018)/7.)*(mypn0-2); // 5.4 m 
	      }
	  else if(sector==4) {
	    a = 0.0026+((0.0075-0.0026)/7.)*(mypn0-2); // 8.1 m 
	      }
	}

	barycerror_+=a*a;
	barycerrorW_+=a*a;

    /*

  if(detId_ < 21) {
    float a = 0.0001*(int((detId_-1)/2.)+1)/pitchall;
    barycerror_+=a*a;
  }
  else if(detId_ < 41) {
    float a = 0.0001*(int((detId_-21)/2.)+1)/pitchall;
           a +=0.0036; // 2.5 m
    //          a +=0.0052; // 4 m
    //  a +=0.0131;// 8. m
    barycerror_+=a*a;
  }
  else if(detId_ < 61) {
    float a = 0.0001*(int((detId_-41)/2.)+1)/pitchall;
             a +=0.0069;// 5 m  0.0059
    //          a +=0.0101;// 8. m
    //  a +=0.0241;// 16. m
    barycerror_+=a*a;
  }
  else if(detId_ < 81) {
    float a = 0.0001*(int((detId_-61)/2.)+1)/pitchall;
          a +=0.0131;// 7.5 m   0.0111
    //       a +=0.0151;// 12. m
    //  a +=0.0301;// 24. m
    barycerror_+=a*a;
  }
*/
#ifdef mydigidebug11
   std::cout << "AT end: barycerror_ = " << barycerror_ << std::endl;
#endif

  barycerror_ = sqrt(  barycerror_ );
  err = barycerror_;

  barycerrorW_ = sqrt(  barycerrorW_ );


#ifdef mydigidebug11
   std::cout << "AT end: err = " << err<< "   detId_= " << detId_ << std::endl;
#endif

}

