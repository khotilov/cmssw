// Original Author: Gero Flucke
// last change    : $Date: 2010/10/26 21:34:25 $
// by             : $Author: flucke $
#ifndef COMPAREMILLEPEDE_H
#define COMPAREMILLEPEDE_H

class PlotMillePede;

class GFHistManager;
class CompareMillePede

{
 public:
  CompareMillePede(const char *fileName1, const char *fileName2, Int_t iter1 = 2, Int_t iter2 = 2,
		   Int_t hieraLevel = 0);// iter=1/2: singlerun/merged; heiraLev: -1 ignore, 0 lowest level, etc.
  virtual ~CompareMillePede();

  void DrawPedeParam(Option_t *option = "", unsigned int nNonRigidParam = 12);//"add": keep old canvas, "free1/2": if free param in file 1/2 
  void DrawPedeParamVsLocation(Option_t *option = "", unsigned int nNonRigidParam = 12);//"add": keep old canvas, "free1/2": if free param in file 1/2 

  void DrawParam(Option_t *option="");//"add": keep old canvas, "free1/2": if free param in file 1/2 
  void DrawParamVsLocation(Option_t *option="");//"add": keep old canvas, "free1/2": if free param in file 1/2 
  void DrawParamDeltaMis(Option_t *option="");//"add": keep old canvas, "free1/2": if free param in file 1/2 
  void DrawParamDeltaMisVsLoc(Option_t *option="");//"add": keep old canvas, "free1/2": if free param in file 1/2
  void DrawNumHits(Option_t *opt="");//"add": keep old canvas

  void DrawAbsPos(Option_t *opt="");//"start": at start (else end), "add": keep old canvas

  bool IsConsistent(); // check correct order of alignables, if false draw some hists
  TString DeltaPar(UInt_t iPar) const; // par_2 - par_1
  TString DeltaParBySigma(UInt_t iPar, const PlotMillePede *sigmaSource) const;
  TString DeltaMisPar(UInt_t iPar) const; // abs(misalignment_2)-abs(misalignment_1)
  TString DeltaMisParBySigma(UInt_t iPar, const PlotMillePede *sigmaSource) const;
  TString DeltaPos(UInt_t iPos) const;

  void AddIsFreeSel(TString &sel, const TString &option, UInt_t iPar) const;

  void SetSubDetId(Int_t subDetId); // 1-6 are TPB, TPE, TIB, TID, TOB, TEC, -1 means: take all
  void AddSubDetId(Int_t subDetId); // 1-6 are TPB, TPE, TIB, TID, TOB, TEC
  void SetAlignableTypeId(Int_t alignableTypeId);//detunit=1,det=2,rod=3,etc. from AlignableObjectIdType (-1: all)
  void SetHieraLevel(Int_t hieraLevel); // select hierarchical level (-1: all)
  void AddAdditionalSel(const char *selection);// special select; StripDoubleOr1D,StripRphi,StripStereo
  void AddAdditionalSel(const TString &xyzrPhiNhit, Float_t min, Float_t max); // x,y,z,r,phi,Nhit
  //  const TString GetAdditionalSel () const { return fAdditionalSel;}
  void ClearAdditionalSel ();


  TString TitleAdd() const;

  PlotMillePede* GetPlotMillePede1() {return fPlotMp1;}
  PlotMillePede* GetPlotMillePede2() {return fPlotMp2;}
  GFHistManager* GetHistManager() { return fHistManager;}

  static const unsigned int kNpar; // number of parameters we have...

 private: 
  CompareMillePede() : fPlotMp1(0), fPlotMp2(0), fHistManager(0) {}

  Int_t PrepareAdd(bool addPlots);

  PlotMillePede *fPlotMp1;
  PlotMillePede *fPlotMp2;

  GFHistManager *fHistManager;
};

#endif
