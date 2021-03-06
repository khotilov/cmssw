//************************************************************************
//* First test of CMKIN++ with the Pythia interface and HepMC output     *
//* into a POOL file                                                     *
//*                                                                      *
//* Based on testWriteHepMC.cpp by Vincenzo Innocente and Filip Moortgat *
//************************************************************************
//* author: Hector Naves Sordo -- Universidad de Oviedo / IFCA           *
//* email:  Hector.Naves.Sordo@cern.ch                                   *
//* First Release: 20/01/05                                              *
//*                                                                      *
//*      (Thanks to Anabel Diaz Rodriguez for her invaluable             *
//*         help and inspiration on make this code works.)               *
//************************************************************************

#include "Utilities/Configuration/interface/Architecture.h"
#include <iostream>
#include <memory>
#include <iomanip>
#include "GeneratorInterface/CMKINInterface/interface/Kikeys.h"
#include "GeneratorInterface/CMKINInterface/interface/Mcparam.h"
#include "GeneratorInterface/CMKINInterface/interface/Kigeneral.h"
#include "GeneratorInterface/CMKINInterface/interface/Kipythia.h"
#include "GeneratorInterface/Pythia6Interface/interface/Hepevt.h"
#include "GeneratorInterface/Pythia6Interface/interface/HepEventPythia6Reader.h"
#include "GeneratorInterface/CMKINInterface/interface/Kiutil.h"
#include "GeneratorInterface/CMKINInterface/interface/Kiffread.h"
#include "GeneratorInterface/CMKINInterface/interface/Kicwn.h"
#include "GeneratorInterface/CMKINInterface/interface/Kiuser.h"
#include "GeneratorInterface/CMKINInterface/interface/Util.h"
#include "GeneratorInterface/CMKINInterface/interface/Pystat.h"
#include "GeneratorInterface/CMKINInterface/interface/Pretauola.h"

// ******************************************************************
// From testWriteHepMC.cpp by Vincenzo Innocente and Filip Moortgat (CERN)
//
#include <fstream>
#include "CLHEP/HepMC/ConvertHEPEVT.h"
#include "CLHEP/HepMC/GenEvent.h"
#include "CLHEP/HepMC/WriteHepMC.h"

// best to include ParticleDataTableConfig after GenEvent
#include "CLHEP/HepMC/ParticleDataTableConfig.h"

#include "CARF/GeneratorApplication/interface/GeneratorFactory.h"

#include "Utilities/GenUtil/interface/CMSTime.h"

#include<vector>
#include<set>
#include<algorithm>
#include "CLHEP/Random/RandFlat.h"
#include<iostream>
//
// ******************************************************************

int isel;


using namespace std;


namespace {

  class MyFactory : public GeneratorFactory {
    HepMC::ConvertHEPEVT conv;

  public:

    MyFactory() 
    {

      // CMKIN initialisation (data cards, etc.)
      int i_flag = 0;
      ki_init(i_flag);


      // KSEL-dependent selection initialisation
      int flag = -1;
      bool DUMMY;
      ki_evt_sel(flag, DUMMY);

      // PYTHIA initialisation
      ki_pyt_init();


      // tauola usage
      Kitaumod myKitaumod;

      int itau = -1;
      if (myKitaumod.jakk1() >= 0) {
	PRETAUOLA(itau);
      }


      Kiextra myKiextra;
      Kiffgen myKiffgen;

      if (myKiextra.nNfil() == -99) {
	myKiextra.nNfil() = myKiffgen.nNsel();
      }

      //Initialise selected event counter
      isel = 0;
      //      int ifile = 0;


      //************************
      //From testWriteHepMC.cpp
      //************************
      //........................................Particle Data Table
      const char infile[] = "pythia.tbl";
      // open input file
      std::ifstream pdfile( infile );
      if( !pdfile ) {
	std::cerr << "cannot open " << infile << std::endl;
	exit(-1);
      }
      // construct empty PDT
      HepMCConfig::ParticleDataTable dataTable( "Pythia Table" );
      {
	// Construct table builder
	HepMCConfig::TableBuilder  tb(dataTable);
	// read the input - put as many here as you want
	if( !addPythiaParticles( pdfile, tb ) )
	  { std::cout << "error reading pythia table file " << std::endl; }
      }   // the tb destructor fills dataTable

    } //end MyFactory()






    virtual HepMC::GenEvent *  operator()() {

      int ihepc = 1;
      PYHEPC(&ihepc);
      return conv.getGenEventfromHEPEVT();
	  
    }


    virtual void run() {
//      int nev=1000;
      Kiffgen myKiffgen;
      Kiextra myKiextra;


      //      int ifile = 0;


      int nev = myKiffgen.nTrig();
      std::cout << "start a new run "
                << " with " << nev << " events" << std::endl;
      {
        newRun();
      }


      //**********************
      //*  Event generation  *
      //**********************

      //*************************
      // There is a Conflict with Hepevt common definition
      // in CLHEP/HepMC/ConvertHEPEVT.h.
      //
      // Use instead  hepevt_.nhep
      //*************************
      //      Hepevt myHepevt;
      Mcparam myMcparam;
      Kicwn myKicwn;
      Kitaumod myKitaumod;
      int itau;


      for (int iev=0; iev<nev;iev++) {
      
	hepevt_.nhep = 0;
	myMcparam.nvrmcp() = 0;
	PYEVNT();  // Generate Pythia event

	if (myKitaumod.jakk1() >= 0) {
	  itau = 0;
	  PRETAUOLA(itau);
	}

	int flag = 0;
	bool I_WANT;
	ki_evt_sel(flag, I_WANT);

	if (I_WANT) {
	  isel++;

	  // Show the first 10 events 
	  //	  if (isel <= 10) {
	  //	    int ilist = 1;
	  //	    PYLIST(&ilist);
	  //	  }

	  if (hepevt_.nhep == 0) {
	    kiu_fill_hepevt();
	  }
	  
	  if (myKicwn.idcwn() > 0 ) {
	    newGenerator();
	  }

	  else {
	    kir_anal();
	  }

	} // end if




        if ((isel % 10) == 1 || isel == myKiffgen.nNsel() ||
            (isel % myKiextra.nNfil()) == 0) {
          cout << "**KIMAIN** Nb. events triggered " << iev
               << ", selected =" << isel << endl;
	}

	if (isel == myKiffgen.nNsel()) break;


      } //end for

      //***********************
      //*  Termination Phase  *
      //***********************
      
      bool DUMMY;
      int flag = 1;
      ki_evt_sel(flag, DUMMY);  // KSEL-dependent selection termination

      int istat = 1;
      PYSTAT(istat);                // Write-out PYTHIA run summary

      if (myKitaumod.jakk1() > 0) {
	int itau = 1;
	PRETAUOLA(itau);
      }

      kir_end();

      int i_flag = 0;
      ki_end(i_flag);


    }

  };

  MyFactory bha;

}






















// -----------------------------------------------------------------------------
//  First steps toward CMKIN++
//  Kikeys.cc -> CMKIN++ specification of data cards variables 
// -----------------------------------------------------------------------------
//  Date: 2004/11/10
//  Revision: 0.002 (pre alpha!)
// -----------------------------------------------------------------------------
//  Author: H�ctor Naves Sordo (Email: naves@mail.cern.ch)
//          (Original kikeys.inc by Lucas Taylor 19/12/95)
// -----------------------------------------------------------------------------

#ifndef _WIN32
# define CFREAD     cfread_
# define KI_FFGEN   ki_ffgen_
# define KI_EXTRA   ki_extra_
# define KI_TOPREX  ki_toprex_
# define KI_TXPAR   ki_txpar_
# define KI_FCNC    ki_fcnc_
# define KI_TXDEC   ki_txdec_
# define KI_TAUMOD  ki_taumod_
#endif

//***********************
//* CFREAD common block *
//***********************

extern "C" {
  void *getaddr(void *);

  extern struct {
    double ffspace[50000];
  } CFREAD;
}

struct Cfread::_cfread * Cfread::__cfread = 0;

int Cfread::nDummy = 0;
double Cfread::dDummy = 0.0;

Cfread::Cfread(void)
{
}

Cfread::~Cfread(void)
{
}

void Cfread::init(void)
{
  __cfread = static_cast<struct Cfread::_cfread *>(getaddr(&CFREAD.ffspace[0]));
}

double &Cfread::ffspace(int n)
{
  if (__cfread == 0) init();
  if ((n<1)||(n>_nffword)) {
    dDummy = -999.0;
    return dDummy;
  }
  return Cfread::__cfread->ffspace[n-1];
}



//*************************
//* KI_FFGEN common block *
//*************************

extern "C" {
  void *getaddr(void *);

  extern struct {
    int nTrig;
    int nNsel;
    double energyCMS;
    int ksel;
    int icut[100];
    double rcut[100];
    double cut1;
    double cut2;
    int ikine;
    double pkine[10];
    int ihwg;
    int nrnhwg[2];
  } KI_FFGEN;
}

struct Kiffgen::_kiffgen * Kiffgen::__kiffgen = 0;

int Kiffgen::nDummy = 0;
double Kiffgen::dDummy = 0.0;

Kiffgen::Kiffgen(void)
{
}

Kiffgen::~Kiffgen(void)
{
}

void Kiffgen::init(void)
{
  __kiffgen = static_cast<struct  Kiffgen::_kiffgen *>(getaddr(&KI_FFGEN.nTrig));
}


int &Kiffgen::nTrig(void)
{
  if (__kiffgen == 0) init();
  return __kiffgen->nTrig; 
}

int &Kiffgen::nNsel(void)
{
  if (__kiffgen == 0) init();
  return __kiffgen->nNsel;
}

double &Kiffgen::energyCMS(void)
{
  if (__kiffgen == 0) init();
  return __kiffgen->energyCMS;
}

int &Kiffgen::ksel(void)
{
  if (__kiffgen == 0) init();
  return __kiffgen->ksel;
}

int &Kiffgen::icut(int n)
{
  if (__kiffgen == 0) init();
  if ((n<1)||(n>_ncut)) {
    nDummy = -999;
    return nDummy;
  }
  return __kiffgen->icut[n-1];
}
      
double &Kiffgen::rcut(int n)
{
  if (__kiffgen == 0) init();
  if ((n<1)||(n>_ncut)) {
    dDummy = -999.0;
    return dDummy;
  }
  return __kiffgen->rcut[n-1];
}
  
double &Kiffgen::cut1(void)
{
  if (__kiffgen == 0) init();
  return __kiffgen->cut1;
}

double &Kiffgen::cut2(void)
{
  if (__kiffgen == 0) init();
  return __kiffgen->cut2;
}

int &Kiffgen::ikine(void)
{
  if (__kiffgen == 0) init();
  return __kiffgen->ikine;
}

double &Kiffgen::pkine(int n)
{
  if (__kiffgen == 0) init();
  if ((n<1)||(n>_lenPkine)) {
    dDummy = -999.0;
    return dDummy;
  }
  return __kiffgen->pkine[n-1];
}

int &Kiffgen::ihwg(void)
{
  if (__kiffgen == 0) init();
  return __kiffgen->ihwg;
}

int &Kiffgen::nrnhwg(int n)
{
  if (__kiffgen == 0) init();
  if ((n<1)||(n>_lenNrnhwg)) {
    nDummy = -999;
    return nDummy;
  }
  return __kiffgen->nrnhwg[n-1];
}




//*************************
//* KI_EXTRA common block *
//*************************

extern "C" {
  void *getaddr(void *);

  extern struct {
    int nNfil;
  } KI_EXTRA;
}

struct Kiextra::_kiextra *Kiextra::__kiextra = 0;

int Kiextra::nDummy = 0;
double Kiextra::dDummy = 0.0;
                                                                                
Kiextra::Kiextra(void)
{
}
                                                                                
Kiextra::~Kiextra(void)
{
}

void Kiextra::init(void)
{
  __kiextra = static_cast<struct Kiextra::_kiextra *>(getaddr(&KI_EXTRA.nNfil));
}

int &Kiextra::nNfil(void)
{
  if (__kiextra == 0) init();
  return __kiextra->nNfil;
}



//*************************
//* KI_TXPAR common block *
//*************************

extern "C" {
  void *getaddr(void *);

  extern struct {
    int ipars[200];
    double rpars[200];
  } KI_TXPAR;
}

struct Kitxpar::_kitxpar *Kitxpar::__kitxpar = 0;

int Kitxpar::nDummy = 0;
double Kitxpar::dDummy = 0.0;

Kitxpar::Kitxpar(void)
{
}

Kitxpar::~Kitxpar(void)
{
}

void Kitxpar::init(void)
{
  __kitxpar = static_cast<struct Kitxpar::_kitxpar *>(getaddr(&KI_TXPAR.ipars[0]));
}

int &Kitxpar::ipars(int n)
{
  if (__kitxpar == 0) init();
  if ((n<1)||(n>_lenPars)) {
    nDummy = -999;
    return nDummy;
  }
  return __kitxpar->ipars[n-1];
}

double &Kitxpar::rpars(int n)
{
  if (__kitxpar == 0) init();
  if ((n<1)||(n>_lenPars)) {
    dDummy = -999.0;
    return dDummy;
  }
  return __kitxpar->rpars[n-1];
}



//************************
//* KI_FCNC common block *
//************************

extern "C" {
  void *getaddr(void *);
                                                                                       
  extern struct {
    double gfcncs[2][8];
  } KI_FCNC;
}

struct Kifcnc::_kifcnc *Kifcnc::__kifcnc = 0;

int Kifcnc::nDummy = 0;
double Kifcnc::dDummy = 0.0;
                                                                                
Kifcnc::Kifcnc(void)
{
}
                                                                                
Kifcnc::~Kifcnc(void)
{
}
                                                                                
void Kifcnc::init(void)
{
  __kifcnc = static_cast<struct Kifcnc::_kifcnc *>(getaddr(&KI_FCNC.gfcncs[0][0]));
}

double &Kifcnc::gfcncs(int i, int j)
{
  if (__kifcnc == 0) init();
  if ( ((i<1)||(i>_lenGfcncs)) ||
       ((j<1)||(j>_depthGfcncs)) )
    {
      dDummy = -999.0;
      return dDummy;
    }
  else
    return Kifcnc::__kifcnc->gfcncs[i-1][j-1];

}



//*************************
//* KI_TXDEC common block *
//*************************

extern "C" {
  void *getaddr(void *);

  extern struct {
    int midx[50][4];
  } KI_TXDEC;
}

struct Kitxdec::_kitxdec *Kitxdec::__kitxdec = 0;


int Kitxdec::nDummy = 0;
double Kitxdec::dDummy = 0.0;
                                                                                
Kitxdec::Kitxdec(void)
{
}
                                                                                
Kitxdec::~Kitxdec(void)
{
}
                                                                                
void Kitxdec::init(void)
{
  __kitxdec = static_cast<struct Kitxdec::_kitxdec *>(getaddr(&KI_TXDEC.midx[0][0]));
}

int &Kitxdec::midx(int i, int j)
{
  if (__kitxdec == 0) init();
  if ( ((i<1)||(i>_lenMidx)) ||
       ((j<1)||(j>_depthMidx)) )
    {
      nDummy = -999;
      return nDummy;
    }
  else
    return Kitxdec::__kitxdec->midx[i-1][j-1];
}



//**************************
//* KI_TAUMOD common block *
//**************************


extern "C" {
  void *getaddr(void *);

  extern struct {
    int jakk1;
    int jakk2;
    int mdtau;
  } KI_TAUMOD;
}

struct Kitaumod::_kitaumod *Kitaumod::__kitaumod = 0;

int Kitaumod::nDummy = 0;
double Kitaumod::dDummy = 0.0;
                                                                                
Kitaumod::Kitaumod(void)
{
}
                                                                                
Kitaumod::~Kitaumod(void)
{
}
                                                                                
void Kitaumod::init(void)
{
  __kitaumod = static_cast<struct Kitaumod::_kitaumod *>(getaddr(&KI_TAUMOD.jakk1));
}

int &Kitaumod::jakk1(void)
{
  if (__kitaumod == 0) init();
  return __kitaumod->jakk1;

}

int &Kitaumod::jakk2(void)
{
  if (__kitaumod == 0) init();
  return __kitaumod->jakk2;
}

int &Kitaumod::mdtau(void)
{
  if (__kitaumod == 0) init();
  return __kitaumod->mdtau;
}



extern "C" void *getaddr(void *);

void *getaddr(void *item)
{
	return item;
}




// -----------------------------------------------------------------------------
//  First steps toward CMKIN++
//  Mcparam.cc -> CMKIN++ specification of data cards variables 
// -----------------------------------------------------------------------------
//  Date: 2004/11/15
//  Revision: 0.002 (pre alpha!)
// -----------------------------------------------------------------------------
//  Author: H�ctor Naves Sordo (Email: naves@mail.cern.ch)
//          (Original kikeys.inc by Lucas Taylor 19/12/95)
// -----------------------------------------------------------------------------

#ifndef _WIN32
# define MC_PARAM     mc_param_
#endif


//*************************
//* MC_PARAM common block *
//*************************

extern "C" {
  void *getaddr(void *);

  extern struct {
    int ivercms;
    int irnmcp;
    int ievmcp;
    double wgtmcp;
    double xsecn;
    int ifilter;
    int nvrmcp;
    double varmcp[200];
  } MC_PARAM;
}

struct Mcparam::_mcparam * Mcparam::__mcparam = 0;

int Mcparam::nDummy = 0;
double Mcparam::dDummy = 0.0;

Mcparam::Mcparam(void)
{
}

Mcparam::~Mcparam(void)
{
}

void Mcparam::init(void)
{
  __mcparam = static_cast<struct  Mcparam::_mcparam *>(getaddr(&MC_PARAM.ivercms));
}


int &Mcparam::ivercms(void)
{
  if (__mcparam == 0) init();
  return __mcparam->ivercms; 
}

int &Mcparam::irnmcp(void)
{
  if (__mcparam == 0) init();
  return __mcparam->irnmcp;
}

int &Mcparam::ievmcp(void)
{
  if (__mcparam == 0) init();
  return __mcparam->ievmcp;
}

double &Mcparam::wgtmcp(void)
{
  if (__mcparam == 0) init();
  return __mcparam->wgtmcp;
}

double &Mcparam::xsecn(void)
{
  if (__mcparam == 0) init();
  return __mcparam->xsecn;
}

int &Mcparam::ifilter(void)
{
  if (__mcparam == 0) init();
  return __mcparam->ifilter;
}

int &Mcparam::nvrmcp(void)
{
  if (__mcparam == 0) init();
  return __mcparam->nvrmcp;
}
      
double &Mcparam::varmcp(int n)
{
  if (__mcparam == 0) init();
  if ((n<1)||(n>_nmxmcp)) {
    dDummy = -999.0;
    return dDummy;
  }
  return __mcparam->varmcp[n-1];
}
  
                                                                                


// -----------------------------------------------------------------------------
//  First steps toward CMKIN++
//  Kicwn.cc -> CMKIN++ specification of data cards variables 
// -----------------------------------------------------------------------------
//  Date: 2004/11/15
//  Revision: 0.002 (pre alpha!)
// -----------------------------------------------------------------------------
//  Author: H�ctor Naves Sordo (Email: naves@mail.cern.ch)
//          (Original kikeys.inc by Lucas Taylor 19/12/95)
// -----------------------------------------------------------------------------





#ifndef _WIN32
# define PAWC      pawc_
# define KICWN     kicwn_
#endif


//*************************
//*   PAWC common block   *
//*************************

extern "C" {
  void *getaddr(void *);

  extern struct {
    double hpawc[2000000];
  } PAWC;
}

struct Pawc::_pawc * Pawc::__pawc = 0;

int Pawc::nDummy = 0;
double Pawc::dDummy = 0.0;

Pawc::Pawc(void)
{
}

Pawc::~Pawc(void)
{
}

void Pawc::init(void)
{
  __pawc = static_cast<struct  Pawc::_pawc *>(getaddr(&PAWC.hpawc[0]));
}

double &Pawc::hpawc(int n)
{
  if (__pawc == 0) init();
  if ((n<1)||(n>_nwpawc)) {
    dDummy = -999.0;
    return dDummy;
  }
  return __pawc->hpawc[n-1];
}





//*************************
//*   KICWN common block   *
//*************************
                                                                                      
extern "C" {
  void *getaddr(void *);
                                                                                      
  extern struct {
    int idcwn;
  } KICWN;
}
                                                                                      
struct Kicwn::_kicwn * Kicwn::__kicwn = 0;
                                                                                      
int Kicwn::nDummy = 0;
double Kicwn::dDummy = 0.0;
                                                                                      
Kicwn::Kicwn(void)
{
}
                                                                                      
Kicwn::~Kicwn(void)
{
}
                                                                                      
void Kicwn::init(void)
{
  __kicwn = static_cast<struct  Kicwn::_kicwn *>(getaddr(&KICWN.idcwn));
}

int &Kicwn::idcwn(void)
{
  if (__kicwn == 0) init();
  return __kicwn->idcwn; 
}

                                                                                


