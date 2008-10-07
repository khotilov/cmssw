//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Nov 18 13:55:12 2007 by ROOT version 5.14/00f
// from TTree t1/analysis tree
// found on file: analysis_my.root
//////////////////////////////////////////////////////////

#ifndef JPTRootAnalysisEta_h
#define JPTRootAnalysisEta_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class JPTRootAnalysisEta {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Double_t        EtaGen1;
   Double_t        PhiGen1;
   Double_t        EtaRaw1;
   Double_t        PhiRaw1;
   Double_t        EtGen1;
   Double_t        EtRaw1;
   Double_t        EtZSP1;
   Double_t        EtJPT1;
   Double_t        DRMAXgjet1;
   Double_t        EtaGen2;
   Double_t        PhiGen2;
   Double_t        EtaRaw2;
   Double_t        PhiRaw2;
   Double_t        EtGen2;
   Double_t        EtRaw2;
   Double_t        EtZSP2;
   Double_t        EtJPT2;
   Double_t        DRMAXgjet2;

   // List of branches
   TBranch        *b_EtaGen1;   //!
   TBranch        *b_PhiGen1;   //!
   TBranch        *b_EtaRaw1;   //!
   TBranch        *b_PhiRaw1;   //!
   TBranch        *b_EtGen1;   //!
   TBranch        *b_EtRaw1;   //!
   TBranch        *b_EtZSP1;   //!
   TBranch        *b_EtJPT1;   //!
   TBranch        *b_DRMAXgjet1;   //!
   TBranch        *b_EtaGen2;   //!
   TBranch        *b_PhiGen2;   //!
   TBranch        *b_EtaRaw2;   //!
   TBranch        *b_PhiRaw2;   //!
   TBranch        *b_EtGen2;   //!
   TBranch        *b_EtRaw2;   //!
   TBranch        *b_EtZSP2;   //!
   TBranch        *b_EtJPT2;   //!
   TBranch        *b_DRMAXgjet2;   //!

   JPTRootAnalysisEta(TTree *tree=0);
   virtual ~JPTRootAnalysisEta();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef JPTRootAnalysisEta_cxx
JPTRootAnalysisEta::JPTRootAnalysisEta(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("analysis.root");
    //     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("analysisIterTrk_2.root");
    if (!f) {
      f = new TFile("analysis.root");
      //        f = new TFile("analysisIterTrk_2.root");
      }
    tree = (TTree*)gDirectory->Get("t1");
   }

/*
  if (tree == 0) {
    TChain * chain = new TChain("t1","");
    
     // eff
    //    chain->Add("/afs/cern.ch/user/a/anikiten/scratch0/CMSSW_1_6_9/src/JetMETCorrections/JetPlusTrack/test/RESULT/analysisIterTrkOutLeakMuonsZSP152noNu_*.root");
    chain->Add("/afs/cern.ch/user/a/anikiten/scratch0/CMSSW_1_6_9/src/JetMETCorrections/JetPlusTrack/test/RESULT/analysisIterTrk_*_resp.root");
    tree = chain;
*/

   Init(tree);
}

JPTRootAnalysisEta::~JPTRootAnalysisEta()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JPTRootAnalysisEta::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t JPTRootAnalysisEta::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void JPTRootAnalysisEta::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EtaGen1", &EtaGen1, &b_EtaGen1);
   fChain->SetBranchAddress("PhiGen1", &PhiGen1, &b_PhiGen1);
   fChain->SetBranchAddress("EtaRaw1", &EtaRaw1, &b_EtaRaw1);
   fChain->SetBranchAddress("PhiRaw1", &PhiRaw1, &b_PhiRaw1);
   fChain->SetBranchAddress("EtGen1", &EtGen1, &b_EtGen1);
   fChain->SetBranchAddress("EtRaw1", &EtRaw1, &b_EtRaw1);
   fChain->SetBranchAddress("EtZSP1", &EtZSP1, &b_EtZSP1);
   fChain->SetBranchAddress("EtJPT1", &EtJPT1, &b_EtJPT1);
   fChain->SetBranchAddress("DRMAXgjet1", &DRMAXgjet1, &b_DRMAXgjet1);
   fChain->SetBranchAddress("EtaGen2", &EtaGen2, &b_EtaGen2);
   fChain->SetBranchAddress("PhiGen2", &PhiGen2, &b_PhiGen2);
   fChain->SetBranchAddress("EtaRaw2", &EtaRaw2, &b_EtaRaw2);
   fChain->SetBranchAddress("PhiRaw2", &PhiRaw2, &b_PhiRaw2);
   fChain->SetBranchAddress("EtGen2", &EtGen2, &b_EtGen2);
   fChain->SetBranchAddress("EtRaw2", &EtRaw2, &b_EtRaw2);
   fChain->SetBranchAddress("EtZSP2", &EtZSP2, &b_EtZSP2);
   fChain->SetBranchAddress("EtJPT2", &EtJPT2, &b_EtJPT2);
   fChain->SetBranchAddress("DRMAXgjet2", &DRMAXgjet2, &b_DRMAXgjet2);
   Notify();
}

Bool_t JPTRootAnalysisEta::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void JPTRootAnalysisEta::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t JPTRootAnalysisEta::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef JPTRootAnalysisEta_cxx
