
#ifndef EcalSimAlgos_EcalCoder_h
#define EcalSimAlgos_EcalCoder_h 1

#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstantsMC.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EcalCorrelatedNoiseMatrix.h"

template<typename M> class CorrelatedNoisifier ;
class EcalMGPASample;
class EcalDataFrame;
class DetId;

#include<vector>

/* \class EEDigitizerTraits
 * \brief Converts CaloDataFrame in CaloTimeSample and vice versa.
 *
 */
class EcalCoder
{
   public:
      
      enum { NBITS         =   12 , // number of available bits
	     MAXADC        = 4095 , // 2^12 -1,  adc max range
	     ADCGAINSWITCH = 4079 , // adc gain switch
	     NGAINS        =    3   // number of electronic gains
      };

      /// ctor
      EcalCoder( bool                                 addNoise    , 
		 CorrelatedNoisifier<EcalCorrMatrix>* ebCorrNoise0 ,
		 CorrelatedNoisifier<EcalCorrMatrix>* eeCorrNoise0 = 0 ,
		 CorrelatedNoisifier<EcalCorrMatrix>* ebCorrNoise1 = 0 ,
		 CorrelatedNoisifier<EcalCorrMatrix>* eeCorrNoise1 = 0 ,
		 CorrelatedNoisifier<EcalCorrMatrix>* ebCorrNoise2 = 0 ,
		 CorrelatedNoisifier<EcalCorrMatrix>* eeCorrNoise2 = 0   ) ; // make EE version optional for tb compatibility
      /// dtor
      virtual ~EcalCoder() ;

      /// can be fetched every event from the EventSetup
      void setPedestals( const EcalPedestals* pedestals ) ;

      void setGainRatios( const EcalGainRatios* gainRatios ) ;

      void setFullScaleEnergy( double EBscale ,
			       double EEscale   ) ;

      void setIntercalibConstants( const EcalIntercalibConstantsMC* ical ) ; 
 

      /// from CaloSamples to EcalDataFrame
      virtual void analogToDigital( const CaloSamples& clf , 
				    EcalDataFrame&     df    ) const;
 
   private:

      /// limit on the energy scale due to the electronics range
      double fullScaleEnergy( const DetId & did ) const ;

      /// produce the pulse-shape
      void encode( const CaloSamples& caloSamples , 
		   EcalDataFrame&     df            ) const ;

//      double decode( const EcalMGPASample& sample , 
//		     const DetId&          detId    ) const ;

      /// not yet implemented
      //      void noisify( const EcalIntercalibConstantsMC* values ,
      //		    int                              size     ) const ;

      void findPedestal( const DetId& detId    , 
			 int          gainId   , 
			 double&      pedestal ,
			 double&      width      ) const ;
    
      void findGains( const DetId& detId, 
		      double       theGains[] ) const ;

      void findIntercalibConstant( const DetId& detId ,
				   double&      icalconst ) const ;
   
      const EcalPedestals* m_peds ;
      
      const EcalGainRatios* m_gainRatios ; // the electronics gains

      const EcalIntercalibConstantsMC* m_intercals ; //record specific for simulation of gain variation in MC

      double m_maxEneEB ; // max attainable energy in the ecal barrel
      double m_maxEneEE ; // max attainable energy in the ecal endcap
      
      bool m_addNoise ;   // whether add noise to the pedestals and the gains

      const CorrelatedNoisifier<EcalCorrMatrix>* m_ebCorrNoise[3] ;
      const CorrelatedNoisifier<EcalCorrMatrix>* m_eeCorrNoise[3] ;
};


#endif
