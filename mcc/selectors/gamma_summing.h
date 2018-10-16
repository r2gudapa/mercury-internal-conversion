//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 25 13:18:27 2016 by ROOT version 5.34/24
// from TTree FragmentTree/FragmentTree
// found on file: fragment07844_000.root
//////////////////////////////////////////////////////////

#ifndef gamma_summing_h
#define gamma_summing_h

#include "TChain.h"
#include "TFile.h"

#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"

// Header file for the classes stored in the TTree if any.
#include "TGriffin.h"
//#include "TSceptar.h"
// #include "TZeroDegree.h"
// #include "TLaBr.h"
// #include "TTAC.h"
#include "TPaces.h"
#include "TGRSISelector.h"
#include "TGraphErrors.h"

// Fixed size dimensions of array or collections stored in the TTree if any.

class gamma_summing : public TGRSISelector {

 public :
   TGriffin    * fGrif;
  // TSceptar    * fScep;
  // TZeroDegree * fZDS;
  // TLaBr       * fLaBr;
  // TTAC        * fTAC;
  TPaces      * fPaces;

   gamma_summing(TTree * /*tree*/ =0) : TGRSISelector(), fGrif(0), fPaces(0) /*,fZDS(0), fLaBr(0), fTAC(0), fScep(0)*/ {
      SetOutputPrefix("EfficiencyEvent");
   }
   virtual ~gamma_summing() { }
   virtual Int_t   Version() const { return 2; }
   void CreateHistograms();
   void FillHistograms();
   void InitializeBranches(TTree *tree);

   ClassDef(gamma_summing,2);
};

#endif

#ifdef gamma_summing_cxx
void gamma_summing::InitializeBranches(TTree* tree) {
   if (!tree) return;
   tree->SetBranchAddress("TGriffin", &fGrif);
   //tree->SetBranchAddress("TSceptar", &fScep);
   // tree->SetBranchAddress("TZeroDegree", &fZDS);
   // tree->SetBranchAddress("TLaBr", &fLaBr);
   // tree->SetBranchAddress("TTAC", &fTAC);
   tree->SetBranchAddress("TPaces", &fPaces);

}

#endif // #ifdef gamma_summing_cxx
