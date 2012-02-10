#include "../interface/TestProposal.h"
#include "../interface/DebugProposal.h"
#include "../interface/VerticalInterpPdf.h"
#include "../interface/VerticalInterpHistPdf.h"
#include "../interface/AsymPow.h"
#include "../interface/CombDataSetFactory.h"
#include "../interface/TH1Keys.h"
#include "../interface/RooSimultaneousOpt.h"
#include "../interface/SimpleCacheSentry.h"
#include "../interface/th1fmorph.h"
#include "../interface/HZZ4LRooPdfs.h"
#include "../interface/HZZ2L2QRooPdfs.h"
#include "../interface/HGGRooPdfs.h"
#if ROOT_VERSION_CODE < ROOT_VERSION(5,29,99)
#include "../interface/FlexibleInterpVar.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class TestProposal+;
#pragma link C++ class DebugProposal+;
#pragma link C++ class VerticalInterpPdf+;
#pragma link C++ class VerticalInterpHistPdf+;
#pragma link C++ class AsymPow+;
#pragma link C++ class CombDataSetFactory+;
#pragma link C++ class TH1Keys+;
#pragma link C++ class RooSimultaneousOpt+;
#pragma link C++ class SimpleCacheSentry+;
#pragma link C++ function th1fmorph;
#pragma link C++ class RooqqZZPdf+;
#pragma link C++ class RooggZZPdf+;
#pragma link C++ class RooqqZZPdf_v2+;
#pragma link C++ class RooggZZPdf_v2+;
#pragma link C++ class RooBetaFunc_v2+;
#pragma link C++ class Roo4lMasses2D+;
#pragma link C++ class Roo4lMasses2D_Bkg+;
#pragma link C++ class Roo4lMasses2D_BkgGGZZ+;
#pragma link C++ class RooFourMuMassShapePdf2+;
#pragma link C++ class RooFourEMassShapePdf2+;
#pragma link C++ class RooTwoETwoMuMassShapePdf2+;
#pragma link C++ class RooFourMuMassRes+;
#pragma link C++ class RooFourEMassRes+;
#pragma link C++ class RooTwoETwoMuMassRes+;
#pragma link C++ class RooRelBW1+;
#pragma link C++ class RooRelBWUF+;
#pragma link C++ class RooRelBWUF_SM4+;
#pragma link C++ class RooRelBWUFParam+;
#pragma link C++ class RooRelBW+;
#pragma link C++ class RooDoubleCB+;
#pragma link C++ class RooCB+;
#pragma link C++ class RooFermi+;
#pragma link C++ class Triangle+;
#pragma link C++ class RooLevelledExp+;
#pragma link C++ class RooPower+;
#if ROOT_VERSION_CODE < ROOT_VERSION(5,29,99)
#pragma link C++ class RooStats::HistFactory::FlexibleInterpVar+;
#endif

#endif
