#ifndef DCCRAWDATADEFINITIONS_
#define DCCRAWDATADEFINITIONS_

enum globalFieds{

  B_MASK               =  1,
  HEADERLENGTH         =  9,
  HEADERSIZE           = 72,
  EMPTYEVENTSIZE       = 32,
	
  PHYSICTRIGGER        = 1,
  CALIBRATIONTRIGGER   = 2,
  TESTTRIGGER          = 3,
  TECHNICALTRIGGER     = 4,
      
  CH_ENABLED           = 0,
  CH_DISABLED          = 1,
  CH_TIMEOUT           = 2,
  CH_SUPPRESS          = 7,

  SRP_NREAD            = 0,
  SRP_NUMBFLAGS        = 68,
  SRP_BLOCKLENGTH      = 6,
  SRP_EB_NUMBFLAGS     = 68,
  
  BOEVALUE             = 0x5, 
  ERROR_EMPTYEVENT     = 0x1, 		
  TOWERH_SIZE          = 8, 
  TRAILER_SIZE         = 8,
  TCC_EB_NUMBTTS       = 68,
  TCCID_SMID_SHIFT_EB  = 27,
  
  //ARRAY SIZES
  NUMB_SM              = 54,
  NUMB_FE              = 68,
  NUMB_TCC             = 108,
  NUMB_XTAL            = 5,
  NUMB_STRIP           = 5
  
  
    
  
  

};



enum headerFields{ 
          
  H_FEDID_B            = 8,
  H_FEDID_MASK         = 0xFFF,
  
  H_BX_B               = 20,
  H_BX_MASK            = 0xFFF,
      
  H_L1_B               = 32,
  H_L1_MASK            = 0xFFFFFF,

  H_TTYPE_B            = 56,
  H_TTYPE_MASK         = 0xF,    

  H_EVLENGTH_MASK      = 0xFFFFFF,
      
  H_ERRORS_B           = 24,
  H_ERRORS_MASK        = 0xFF,

  H_RNUMB_B            = 32,
  H_RNUMB_MASK         = 0xFFFFFF,

  H_RTYPE_MASK         = 0xFFFFFFFF,

  H_SR_B               = 32,
  H_ZS_B               = 33,
  H_TZS_B              = 34,
        
  H_SRCHSTATUS_B       = 36,
  H_CHSTATUS_MASK      = 0xF,

  H_TCC1CHSTATUS_B     = 40, 
  H_TCC2CHSTATUS_B     = 44,
  H_TCC3CHSTATUS_B     = 48,
  H_TCC4CHSTATUS_B     = 52
     

};		

enum towerFields{ 
       
  TOWER_ID_MASK        = 0x7F,
  
  TOWER_NSAMP_MASK     = 0x7F,
  TOWER_NSAMP_B        = 8,  
      
  TOWER_BX_MASK        = 0xFFF,
  TOWER_BX_B           = 16,     
 
  TOWER_L1_MASK        = 0xFFF,
  TOWER_L1_B           = 32,
      
  TOWER_ADC_MASK       = 0xFFF,
  TOWER_DIGI_MASK      = 0x3FFF,
      
  TOWER_STRIPID_MASK   = 0x7,
      
  TOWER_XTALID_MASK    = 0x7,
  TOWER_XTALID_B       = 4,


  TOWER_LENGTH_MASK    = 0x1FF,
  TOWER_LENGTH_B       = 48

};	


enum tccFields{
  
   TCC_ID_MASK         = 0xFF,
 
   TCC_BX_MASK         = 0xFFF,
   TCC_BX_B            = 16,

   TCC_L1_MASK         = 0xFFF, 
   TCC_L1_B            = 32,  

   TCC_TT_MASK         = 0x7F,
   TCC_TT_B            = 48,

   TCC_TS_MASK         = 0xF,
   TCC_TS_B            = 55

};



enum srpFields{
  
   SRP_ID_MASK         = 0xFF,
 
   SRP_BX_MASK         = 0xFFF,
   SRP_BX_B            = 16,

   SRP_L1_MASK         = 0xFFF, 
   SRP_L1_B            = 32,  

   SRP_NFLAGS_MASK     = 0x7F,
   SRP_NFLAGS_B        = 48,
	
   SRP_SRFLAG_MASK     = 0x3

};


#endif	
