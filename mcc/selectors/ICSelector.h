//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 25 13:18:27 2016 by ROOT version 5.34/24
// from TTree FragmentTree/FragmentTree
// found on file: fragment07844_000.root
//////////////////////////////////////////////////////////

#ifndef ICSelector_h
#define ICSelector_h

#include "TChain.h"
#include "TFile.h"

#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"

// Header file for the classes stored in the TTree if any.
#include "TGriffin.h"
#include "TSceptar.h"
#include "TPaces.h"
#include "TZeroDegree.h"
#include "TGRSISelector.h"

std::vector<std::pair<double, double>> PACESEnergyCorrectionA();
std::vector<std::pair<double, double>> PACESEnergyCorrectionB();
std::vector<std::pair<double, double>> PACESEnergyCorrectionC();
std::vector<std::pair<double, double>> PACESEnergyCorrectionD();
std::vector<std::pair<double, double>> PACESEnergyCorrectionE();

// TODO: add px, py, pbins, pDetectors
int gx=0, gy=1000, gbins=1000; //griffin energy bins 
int gDetectors = 16; //number of germanium detectors
int gCry = 4; //number of germanium crystals
int px=0, py=1000, pbins=1000;
int pDetectors = 5;

class ICSelector : public TGRSISelector { //Must be same name as .C and .h

 public :
   TGriffin * fGrif; //Pointers to spot that events will be
   TPaces * fPaces;

   std::vector<std::pair<double, double>> fPACESEnergyCorrectionA;
   std::vector<std::pair<double, double>> fPACESEnergyCorrectionB;
   std::vector<std::pair<double, double>> fPACESEnergyCorrectionC;
   std::vector<std::pair<double, double>> fPACESEnergyCorrectionD;
   std::vector<std::pair<double, double>> fPACESEnergyCorrectionE;
   
   std::vector<TGriffin> fLastGrif;
   std::vector<TPaces> fLastPaces;

   ICSelector(TTree * /*tree*/ =0) : TGRSISelector(), fGrif(0), fPaces(0) {

	 SetOutputPrefix("ICSelector"); //Changes prefix of output file

  	 fPACESEnergyCorrectionA = PACESEnergyCorrectionA();
         fPACESEnergyCorrectionB = PACESEnergyCorrectionB();
         fPACESEnergyCorrectionC = PACESEnergyCorrectionC();
         fPACESEnergyCorrectionD = PACESEnergyCorrectionD();
         fPACESEnergyCorrectionE = PACESEnergyCorrectionE();
    }

//These functions are expected to exist
   virtual ~ICSelector() { }
   virtual Int_t   Version() const { return 2; }
   void CreateHistograms();
   void FillHistograms();
   void InitializeBranches(TTree *tree);

   ClassDef(ICSelector,2); //Makes ROOT happier
};

#endif

#ifdef ICSelector_cxx
void ICSelector::InitializeBranches(TTree* tree)
{
   if(!tree) return;
   if(tree->SetBranchAddress("TGriffin", &fGrif) == TTree::kMissingBranch) {
		fGrif = new TGriffin;
	}
   if(tree->SetBranchAddress("TPaces", &fPaces) == TTree::kMissingBranch) {
		fPaces = new TPaces;
	}
}

#endif // #ifdef ICSelector_cxx

//Map of residuals to correct PACES non-linearities [Efine = Ecoarse - Eres]
//We produce a vector containing <Ecoarse, Eres> for each detector
std::vector<std::pair<double, double>> PACESEnergyCorrectionA()
{
   std::vector<std::pair<double, double>> result;
   const int dp0 = 11; double ex0[dp0]; double ey0[dp0];	//poor res.
   //PACES0	**corrected energies
   ex0[0]=70.97; ex0[1]=72.92; ex0[2]=131.29; ex0[3]=141.49; ex0[4]=174.16; ex0[5]=196.10; ex0[6]=328.46; ex0[7]=339.95; ex0[8]=358.60; ex0[9]=397.88; ex0[10]=504.83;
   ey0[0]=2.08; ey0[1]=2.11; ey0[2]=-1.21; ey0[3]=-1.61; ey0[4]=-1.24; ey0[5]=-1.17; ey0[6]=-0.24; ey0[7]=-0.25; ey0[8]=-0.10; ey0[9]=0.92; ey0[10]=0.73;

   for(int i = 0; i < int(dp0); ++i) {
      result.push_back(std::make_pair(ex0[i], ey0[i]));
   }
   return result;
}
std::vector<std::pair<double, double>> PACESEnergyCorrectionB()
{
   std::vector<std::pair<double, double>> result;
   const int dp1 = 16; double ex1[dp1]; double ey1[dp1];
   //PACES1	**corrected energies
   ex1[0]=32.85; ex1[1]=44.21; ex1[2]=47.08; ex1[3]=69.22; ex1[4]=71.18; ex1[5]=132.43; ex1[6]=142.92; ex1[7]=175.07; ex1[8]=197.12; ex1[9]=245.19; ex1[10]=307.35; ex1[11]=328.71; 
   ex1[12]=339.44; ex1[13]=359.06; ex1[14]=397.48; ex1[15]=504.23;
   ey1[0]=-0.05; ey1[1]=0.03; ey1[2]=0.14; ey1[3]=0.32; ey1[4]=0.36; ey1[5]=-0.07; ey1[6]=-0.18; ey1[7]=-0.33; ey1[8]=-0.15; ey1[9]=-0.37; ey1[10]=0.05; ey1[11]=0.01;
   ey1[12]=-0.76; ey1[13]=0.36; ey1[14]=0.52; ey1[15]=0.13;

   for(int i = 0; i < int(dp1); ++i) {
      result.push_back(std::make_pair(ex1[i], ey1[i]));
   }
   return result;
}
std::vector<std::pair<double, double>> PACESEnergyCorrectionC()
{
   std::vector<std::pair<double, double>> result;
   const int dp2 = 11; double ex2[dp2]; double ey2[dp2];	//poor res.
   //PACES2	**corrected energies
   ex2[0]=71.16; ex2[1]=73.05; ex2[2]=130.91; ex2[3]=141.36; ex2[4]=174.06; ex2[5]=196.43; ex2[6]=328.46; ex2[7]=339.24; ex2[8]=359.37; ex2[9]=397.77; ex2[10]=504.83;
   ey2[0]=2.27; ey2[1]=2.23; ey2[2]=-1.59; ey2[3]=-1.74; ey2[4]=-1.34; ey2[5]=-0.84; ey2[6]=-0.24; ey2[7]=-0.96; ey2[8]=0.67; ey2[9]=0.81; ey2[10]=0.73;

   for(int i = 0; i < int(dp2); ++i) {
      result.push_back(std::make_pair(ex2[i], ey2[i]));
   }
   return result;
}
std::vector<std::pair<double, double>> PACESEnergyCorrectionD()
{
   std::vector<std::pair<double, double>> result;
   const int dp3 = 16; double ex3[dp3]; double ey3[dp3];
   //PACES3	**corrected energies
   ex3[0]=30.77; ex3[1]=42.83; ex3[2]=45.25; ex3[3]=67.80; ex3[4]=69.87; ex3[5]=131.95; ex3[6]=143.09; ex3[7]=175.72; ex3[8]=199.65; ex3[9]=249.38; ex3[10]=312.52; ex3[11]=334.54; 
   ex3[12]=344.68; ex3[13]=358.51; ex3[14]=388.78; ex3[15]=498.17;
   ey3[0]=-2.13; ey3[1]=-1.35; ey3[2]=-1.69; ey3[3]=-1.09; ey3[4]=-0.95; ey3[5]=-0.55; ey3[6]=-0.01; ey3[7]=0.32; ey3[8]=2.38; ey3[9]=3.82; ey3[10]=5.23; ey3[11]=5.84;
   ey3[12]=4.48; ey3[13]=-0.19; ey3[14]=-8.18; ey3[15]=-5.93;
   
   for(int i = 0; i < int(dp3); ++i) {
      result.push_back(std::make_pair(ex3[i], ey3[i]));
   }
   return result;
}
std::vector<std::pair<double, double>> PACESEnergyCorrectionE()
{
   std::vector<std::pair<double, double>> result;
   const int dp4 = 16; double ex4[dp4]; double ey4[dp4];	//poor res.
   //PACES4	**corrected energies
   ex4[0]=34.22; ex4[1]=45.72; ex4[2]=47.76; ex4[3]=72.63; ex4[4]=74.59; ex4[5]=136.26; ex4[6]=144.82; ex4[7]=167.33; ex4[8]=191.10; ex4[9]=240.29; ex4[10]=304.98; ex4[11]=326.92; 
   ex4[12]=338.30; ex4[13]=358.74; ex4[14]=399.13; ex4[15]=510.70;
   ey4[0]=1.32; ey4[1]=1.54; ey4[2]=0.82; ey4[3]=3.74; ey4[4]=3.78; ey4[5]=3.76; ey4[6]=1.72; ey4[7]=-8.07; ey4[8]=-6.17; ey4[9]=-5.27; ey4[10]=-2.32; ey4[11]=-1.78;
   ey4[12]=-1.90; ey4[13]=0.04; ey4[14]=2.17; ey4[15]=6.60;

   for(int i = 0; i < int(dp4); ++i) {
      result.push_back(std::make_pair(ex4[i], ey4[i]));
   }
   return result;
}
