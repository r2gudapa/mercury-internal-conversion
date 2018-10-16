#define gamma_summing_cxx
// The class definition in gamma_summing.h has been generated automatically
#include "gamma_summing.h"

//DEFS
int maxMult = 6;
double thresh = 0;
int rchan = 64;
double bthresh = 200;		//applied to gainmatched beta energy
bool bhit, abhit;
double gePromptLo = 40.;
double gePromptHi = 100.;
double geRandomLo = 200.;
double geRandomHi = 1000.;
double ggPrompt = 400.;
double ggRandomLo = 400.;
double ggRandomHi = 2000.;
double tlow=-40.,thigh=40.;
double bgalo=60.,bgahi=120.;	//background timing gates (either side of prompt peak), 120/175 suggested
double bgblo=-120.,bgbhi=-60.;	//background timing gates (either side of prompt peak), -175/-120 suggested
double gbin=5000,gx=0,gy=5000;
double sbin=1000,sx=0,sy=5000;
double gbPrompt1=-10,gbPrompt2=40;
double ggbBG1=200,ggbBG2=2000;
double dist[16] = {113.7, 114.8, 112.5, 115.5, 110., 110., 110., 110., 110., 110., 110., 110., 110., 113.3, 113.5, 110.};

// previously: < 80
bool GGCoincidenceCondition(TGriffinHit* hit_one, TGriffinHit* hit_two){
   return TMath::Abs(hit_one->GetTime() - hit_two->GetTime()) < ggPrompt; 
}

// previously: < ggbBG2 and > ggbBG1
bool GGBGCondition(TGriffinHit* hit_one, TGriffinHit* hit_two){
   return TMath::Abs(hit_one->GetTime() - hit_two->GetTime()) < ggRandomHi && TMath::Abs(hit_one->GetTime() - hit_two->GetTime()) > ggRandomLo; 
}

bool GECoincidenceCondition(TPacesHit* hit_paces, TGriffinHit* hit_griffin){	
	double ptime = hit_paces->GetTime();
	double gtime = hit_griffin->GetTime();
	if ((gtime - ptime) < 0) { // hit1 is griffin - LO
		return TMath::Abs(gtime - ptime) < PromptLo;
	} else {
		return TMath::Abs(gtime - ptime) < PromptHi;
	}
}

bool GEBGCondition(TPacesHit* hit_paces, TGriffinHit* hit_griffin){		
   return (TMath::Abs(hit_paces->GetTime() - hit_griffin->GetTime()) < RandomHi) && 
		(TMath::Abs(hit_paces->GetTime() - hit_griffin->GetTime()) > RandomLo); 
}

//bool CoincidenceCondition(TGriffinHit* hit_one, TSceptarHit* hit_two){	
//   return ((hit_one->GetTimeStamp() - hit_two->GetTimeStamp()) < gbPrompt2) && ((hit_one->GetTimeStamp() - hit_two->GetTimeStamp()) > gbPrompt1); 
//}

void gamma_summing::CreateHistograms() {

	//Scintillator
//	fH1["sE"] = new TH1D("sE","#beta Singles",sbin,sx,sy);
	//HPGe - Single Hits
//	fH1["gE"] = new TH1D("gE","#gamma Singles; Energy (keV); Counts per keV",gbin,gx,gy);
//	fH2["gEchan"] = new TH2D("gEchan","#gamma Singles; Channel Number; Energy (keV)",64,0,64,gbin,gx,gy);
//	fH2["gQchan"] = new TH2D("gQchan","#gamma Singles; Channel Number; Charge",64,0,64,8000,0,8000);
//	fH1["gb"] = new TH1D("gb","#beta-#gamma coincidence",gbin,gx,gy);
//	fH1["gb_dt"] = new TH1D("gb_dt","#beta-#gamma coincidence time difference; Time difference (10 ns); Counts per 10 ns",200,-100,100);

	//HPGe - Addback Hits
	fH1["aE"] = new TH1D("aE","#gamma Addback; Energy (keV); Counts per keV",gbin,gx,gy);
	fH1["aE_pileup"] = new TH1D("aE_pileup","#gamma Addback (pileup); Energy (keV); Counts per keV",gbin,gx,gy);
	//fH1["ae"] = new TH1D("ae","e^{-}-#gamma coincidence",gbin,gx,gy);
	//fH1["ae_dt"] = new TH1D("ae", "e^{-}-#gamma coincidence time difference; Time difference (10 ns); Counts per 10ns",200,-100,100);
//	fH1["ab"] = new TH1D("ab","#beta-#gamma coincidence",gbin,gx,gy);
//	fH1["ab_dt"] = new TH1D("ab_dt","#beta-#gamma coincidence time difference; Time difference (10 ns); Counts per 10 ns",200,-100,100);
	
	fH1["aa_dt"] = new TH1D("aa_dt","#gamma-#gamma addback coincidence time difference; Time difference (ns); Counts per ns",4400,-2200,2200);
	fH2["aa"] = new TH2D("aa","#gamma-#gamma coincidence; Energy (keV); Energy (kev)",gbin,gx,gy,gbin,gx,gy);
	fH2["aa_bg"] = new TH2D("aa_bg","#gamma-#gamma coincidence (TR); Energy (keV); Energy (kev)",gbin,gx,gy,gbin,gx,gy);

	fH1["ae_dt"] = new TH1D("ae_dt", "e^{-}-#gamma addback coincidence time difference; Time difference (ns); Counts per ns",4400,-2200,2200);
	fH2["ae"] = new TH2D("ae","e^{-}-#gamma coincidence; Energy (keV); Energy (keV);"gbin,gx,gy,gbin,gx,gy);
	fH2["ae_bg"] = new TH2D("ae_bg","e^{-}-#gamma coincidence (TR); Energy (keV); Energy (keV)",gbin,gx,gy,gbin,gx,gy);

	for(auto it : fH1) {
		GetOutputList()->Add(it.second);
	}
	for(auto it : fH2) {
		GetOutputList()->Add(it.second);
	}
	for(auto it : fHSparse) {
		GetOutputList()->Add(it.second);
	}

}


void gamma_summing::FillHistograms() {

	const int n = 11;
double x[n] = {71.55, 256.89, 523.79, 1121.12, 1435.82, 1460.82, 1472.36, 1519.30, 1553.77, 1590.85, 2205.84};
double dx[n] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
double y[64][n] = 
	{ 
		{1.08, 0.5, 0.74, 0.85, 0.94, 0.4, 0.36, 0.63, 0.52, 0.54, 0.41},
		{-0.23, 0.23, 0.73, 1.35, -0.55, -0.51, -0.84, -0.55, -0.45, -0.47, 0.8},
		{0.87, 0.52, 0.76, 0.65, 1.21, 0.95, 0.11, 1.07, 1.07, 0.92, 0.94},
		{1.83, 0.05, 0.49, 0.7, 0.74, 0.43, 1.37, 0.19, 0.05, -0.29, 0.22},
		{0.89, 0.34, 0.86, 0.96, 1.05, 1.14, 1.93, 0.75, 0.9, 0.72, 1.03},
		{1.33, 0.36, 0.14, 0.75, 1.34, 0.93, 1.43, 0.86, 1.05, 0.73, -0.81},
		{0.78, 0.72, 0.59, 0.59, 0.3, 0.58, -0.04, 0.35, 0.48, -0.14, 0.54},
		{1.06, 1.36, 1.47, 0.72, 1.28, 1.61, 1.34, 0.95, 1.15, 1.35, -0.49},
		{0.88, 1.23, 1.06, 1.07, 0.78, 0.74, 1.44, 0.69, 0.85, 0.83, 0.49},
		{0.79, 0.23, 0.22, 0.91, 0.57, 1.51, 1.53, 0.5, 0.37, -0.15, 0.30},
		{0.86, 1.15, 0.75, 0.74, 0.9, 1.3, 0.49, 0.95, 0.75, 0.70, 0.46},
		{1.27, 0.28, 0.31, 0.96, 1.44, 1.19, 1.42, 1.34, 1.12, 0.85, 0.44},
		{0.62, 0.44, 0.89, 1.00, 0.68, 0.65, 1.19, 0.87, 0.87, 0.83, 0.76},
		{-0.54, -0.06, 0.49, 0.71, 0.35, -0.36, 0.15, -0.57, -0.48, -0.46, 0.08},
		{0.56, 0.38, 0.64, 0.34, 0.53, 0.35, 0.77, 0.16, 0.18, 0.04, 0.70},
		{0.51, 0.61, 0.45, 1.02, 1.46, 1.50, 0.47, 0.12, -0.08, -0.44, -0.34},
		{0.46, 0.27, 0.48, 0.66, 0.72, 0.71, 1.41, 0.52, 0.54, 0.51, 1.12},
		{-0.01, -0.05, 0.43, 1.07, -0.18, -0.05, -0.18, -0.2, -0.03, -0.02, 0.14},
    		{0.63, 0.53, 0.05, 0.54, 1.34, 0.53, 0.58, 0.53, 0.59, 0.55, 0.99},
		{1.07, -0.04, 0.21, 1.13, 0.34, 0.69, 0.79, -0.05, -0.09, -0.07, 0.84},
		{0.76, 0.89, 0.64, 0.77, 0.84, 0.91, 1.58, 1.07, 1.03, 0.92, 0.81},
		{1.57, -0.15, 0.16, 0.8, 0.29, 1.83, 1.26 ,0.51, 0.53, 0.49, 0.37},
		{0.65, 1.01, 0.66, 0.73, 0.49, 1.34, -0.1, 0.66, 0.69, 0.61, 0.30},
		{1.49, 0.01, 0.16, 0.81, 1.53, 1.7, 1.72, 1.55, 1.6, 1.5, 0.49},
		{0.64, 0.8, 0.64, 0.77, 0.67, 0.67, 1.24, 0.76, 0.83, 0.81, 0.77},
		{0.23, 0.24, 0.61, 1.28, 0.01, -0.03, 0.42, -0.08, -0.05, -0.13, 0.71},
		{1.25, 1.36, 1.01, 0.74, 0.84, 0.87, 0.91, 0.7, 0.77, 0.78, 0.25},
		{0.25, 1.05, -0.16, 0.96, 1.23, 1.25, 1.24, 0.49, 0.31, 0.2, 1.52},
		{1.13, 0.97, 0.65, 0.73, 0.26, 0.88, 0.68, 0.45, 0.45, 0.44, 0.83},
		{0.96, -0.26, -0.11, 0.75, 0.69, 0.65, -0.53, 0.23, 0.19, -0.06, -0.36},
		{0.52, 0.65, 0.85, 0.52, 0.97, 1.20, 1.59, 0.37, 0.37, 0.48, 1.09},
		{0.58, 0.79, 1.02, 0.53, 1.05, 1.27, 1.79, 1.38, 1.26, 1.33, -0.32},
		{1.03, 0.73, 0.55, 1.06, 0.52, 0.89, 1.28, 0.63, 0.75, 0.69, 0.40},
		{-0.03, -0.05, 0.44, 0.63, 0.36, 0.72, 0.28, 0.43, 0.31, 0.52, 0.56},
		{0.43, 0.45, 0.62, 0.63, 0.86, 0.46, 0.84, 0.45, 0.41, 0.31, 0.26},
		{1.44, 0.04, 0.29, 1.12, -0.35, -0.46, -0.43, -0.25, -0.18, -0.16, 0.83},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		{0.8, 0.78, 0.84, 0.7, 0.95, 0.71, 1.07, 0.54, 0.7, 0.56, 0.53},
		{0.31, 0.14, 0.69, 0.96, 1.13, 1.74, 0.95, 1.2, 1.24, 1.14, -0.2,},
		{0.47, 0.72, 0.67, 0.76, 0.58, 0.78, 0.56, 0.51, 0.61, 0.47, 0.62},
		{1.62, 1.88, -0.43, 0.61, 1.25, 1.24, 0.21, 0.96, 1.08, 0.48, 0.19},
		{0.39, 0.6, 0.59, 0.42, 0.66, 0.18, -0.56, 0.36, 0.22, 0.12, 0.33},
		{0.36, 0.41, 1.03, 1.4, 0.92, 0.23, 0.83, 0.14, 0.25, 0.45, 1.25},
		{0.56, 0.52, 0.48, 0.45, 0.91, 0.54, 0.61, 0.37, 0.38, 0.28, 0.12},
		{1.05, 0.24, 0.49, 1.07, 2.02, 1.8, 1.11, 1.56, 1.51, 1.32, 0.58},
		{0.71, 0.71, 0.70, 0.75, 0.56, 0.51, 1.57, 0.47, 0.52, 0.48, 0.59},
		{-0.15, 0.19, 0.52, 1.28, 0.43, 0.2, -0.23, 0.4, 0.36, 0.45, 1.15},
		{0.93, 0.88, 0.99, 1.06, 0.91, 1.00, 0.13, 1.08, 0.99, 0.91, 0.59},
		{1.88, 0.65, 0.23, 1.01, 0.35, 0.43, 0.7, -0.24, -0.43, -0.56, -0.06},
		{1.54, 1.02, 0.3, 0.18, 0.68, 0.97, 0.69, 0.49, 0.63, 0.65, 0.51},
		{0.07, 0.23, 0.76, 0.82, -0.22, 0.00, -0.18 -0.47, -0.43, -0.54, 0.61},
		{0.79, 0.54, 0.5, 0.64, 0.88, 0.85, 0.52, 0.64, 0.55, 0.54, 0.29},
		{1.4, 1.69, 0.3, 0.41, 1.12, 0.58, 1.19, -0.02, -0.29, -0.62, -0.07},
		{0.49, 0.65, 0.86, 0.93, 0.74, 0.74, 0.64, 0.69, 0.7, 0.73, 0.86},
		{-0.17, 0.16, 0.4, 0.76, -0.55, -0.58, -0.8, -0.73, -0.59, -0.58, -0.22},
		{1.06, 1.18, 1.15, 0.72, 0.83, 1.06, 1.18, 0.9, 0.9, 0.85, 0.78},
		{-0.38, 0.00, 0.47, 1.03, 0.57, 0.24, 0.07, -0.42, -0.87, -1.32, -0.58},
		{0.68, 0.87, 0.45, 0.47, 0.62, 0.35, 0.3, 0.3, 0.32, 0.38, 0.1},
		{0.91, 0.32, 0.31, 0.67, -0.87, -0.87, -1.57, -1.02, -1.05, -1.02, -0.72},
		{0.97, 1.01, 0.47, 0.69, 0.82, 0.84, 0.37, 0.67, 0.8, 0.69, 0.71},
		{0.82, 0.23, 0.29, 0.95, -0.15, -0.02, 0.62, -0.11, -0.06, -0.07, 0.8}
	};	
	
double  dy[64][n] =
	{
		{0.01, 0.01, 0.01, 0.01, 0.21, 0.34, 0.54, 0.02, 0.02, 0.04, 0.24},
		{0.01, 0.01, 0.01, 0.02, 0.09, 0.19, 0.33, 0.02, 0.02, 0.04, 0.18},
		{0.01, 0.01, 0.01, 0.01, 0.16, 0.17, 0.7, 0.02, 0.01, 0.03, 0.03},
		{0.02, 0.01, 0.02, 0.01, 1.77, 0.12, 0.31, 0.04, 0.02, 0.03, 0.37},
		{0.03, 0.01, 0.01, 0.02, 0.14, 0.1, 0.34, 0.02, 0.02, 0.03, 0.15},
		{0.01, 0.01, 0.01, 0.03, 0.12, 0.12, 0.59, 0.02, 0.06, 0.04, 0.17},
		{0.01, 0.01, 0.00, 0.02, 0.06, 0.31, 0.37, 0.02, 0.03, 0.02, 0.18},
		{0.01, 0.02, 0.02, 0.03, 0.17, 0.12, 0.27, 0.02, 0.03, 0.04, 0.16},
		{0.01, 0.02, 0.01, 0.02, 0.15, 0.14, 1.41, 0.02, 0.03, 0.07, 0.21},
		{0.01, 0.03, 0.01, 0.03, 0.21, 0.17, 0.38, 0.02, 0.01, 0.04, 0.17},
		{0.01, 0.01, 0.01, 0.01, 0.12, 0.14, 0.44, 0.02, 0.01, 0.02, 0.18},
		{0.01, 0.01, 0.04, 0.01, 0.08, 0.1, 0.58, 0.02, 0.01, 0.05, 0.17},
		{0.01, 0.01, 0.01, 0.02, 0.23, 0.12, 0.82, 0.04, 0.01, 0.02, 0.17},
		{0.02, 0.01, 0.02, 0.02, 0.26, 0.32, 0.35, 0.02, 0.04, 0.08, 0.21},
		{0.01, 0.01, 0.02, 0.01, 0.2, 0.09, 0.31, 0.03, 0.01, 0.02, 0.31},
		{0.02, 0.01, 0.01, 0.01, 0.49, 0.01, 0.49, 0.02, 0.02, 0.03, 0.44},
		{0.01, 0.01, 0.01, 0.02, 0.29, 0.09, 0.46, 0.02, 0.01, 0.02, 0.41},
		{0.01, 0.01, 0.01, 0.02, 0.08, 0.24, 0.26, 0.02, 0.01, 0.02, 0.29},
		{0.01, 0.01, 0.01, 0.02, 0.05, 0.19, 0.47, 0.02, 0.01, 0.03, 0.22},
		{0.01, 0.01, 0.01, 0.01, 0.14, 0.26, 0.19, 0.02, 0.01, 0.02, 0.19},
		{0.02, 0.01, 0.01, 0.02, 0.11, 0.16, 0.6, 0.06, 0.02, 0.04, 0.24},
		{0.02, 0.01, 0.01, 0.02, 0.09, 0.26, 0.35 ,0.04, 0.01, 0.04, 0.25},
		{0.02, 0.01, 0.01, 0.01, 0.13, 0.35, 0.76, 0.02, 0.02, 0.03, 0.33},
		{0.02, 0.01, 0.02, 0.03, 0.35, 0.24, 0.19, 0.05, 0.03, 0.07, 0.43}, 
		{0.01, 0.01, 0.02, 0.01, 0.06, 0.1, 0.34, 0.04, 0.01, 0.02, 0.17},
		{0.01, 0.01, 0.01, 0.01, 0.05, 0.21, 0.33, 0.05, 0.01, 0.06, 0.25},
		{0.01, 0.02, 0.01, 0.01, 0.08, 0.16, 0.23, 0.02, 0.01, 0.02, 0.18},
		{0.01, 0.01, 0.01, 0.01, 0.02, 0.18, 0.79, 0.04, 0.02, 0.1, 0.23}, 
		{0.01, 0.01, 0.02, 0.02, 0.1, 0.21, 0.24, 0.02, 0.01, 0.06, 0.49},
		{0.02, 0.01, 0.01, 0.02, 0.12, 0.24, 0.39, 0.03, 0.02, 0.05, 0.53}, 
		{0.03, 0.01, 0.02, 0.04, 0.16, 0.24, 0.21, 0.02, 0.01, 0.03, 0.22},
		{0.04, 0.01, 0.01, 0.02, 0.17, 0.21, 0.6, 0.03, 0.01, 0.03, 0.19}, 
		{0.01, 0.01, 0.01, 0.02, 0.07, 0.7, 0.28, 0.03, 0.01, 0.05, 0.14},
		{0.01, 0.01, 0.01, 0.02, 0.18, 0.14, 1.61, 0.03, 0.02, 0.04, 0.29},
		{0.03, 0.01, 0.01, 0.01, 0.11, 0.12, 0.23, 0.07, 0.01, 0.02, 0.2},
		{0.02, 0.01, 0.01, 0.01, 0.09, 0.13, 0.42, 0.02, 0.01, 0., 0.19},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		{0.01, 0.01, 0.02, 0.03, 0.16, 0.13, 0.35, 0.02, 0.02, 0.06, 0.32},
		{0.02, 0.01, 0.01, 0.02, 0.09, 0.28, 0.27, 0.02, 0.02, 0.03, 0.16},
		{0.02, 0.01, 0.01, 0.01, 0.09, 0.11, 0.55, 0.02, 0.02, 0.03, 0.21},
		{0.02, 0.01, 0.01, 0.01, 0.34, 0.17, 0.64, 0.04, 0.05, 0.05, 0.22},
		{0.01, 0.01, 0.01, 0.01, 0.38, 0.15, 0.37, 0.03, 0.02, 0.03, 0.17},
		{0.02, 0.01, 0.02, 0.02, 0.18, 0.15, 0.3, 0.03, 0.02, 0.03, 0.33},
		{0.01, 0.01, 0.01, 0.01, 0.75, 0.02, 0.59, 0.07, 0.01, 0.02, 0.2},
		{0.02, 0.04, 0.02, 0.03, 0.32, 0.19, 0.36, 0.02, 0.01, 0.05, 0.18},
		{0.01, 0.01, 0.01, 0.01, 0.08, 0.12, 0.19, 0.02, 0.01, 0.02, 0.29},
		{0.01, 0.02, 0.01, 0.02, 0.21, 0.37, 0.25, 0.03, 0.02, 0.03, 0.32},
		{0.01, 0.01, 0.01, 0.03, 0.1, 0.34, 0.32, 0.02, 0.01, 0.02, 0.17},
		{0.02, 0.01, 0.01, 0.03, 0.06, 0.35, 0.6, 0.03, 0.04, 0.06, 0.21},
		{0.03, 0.01, 0.01, 0.01, 0.11, 0.11, 0.19, 0.02, 0.01, 0.03, 0.31},
		{0.03, 0.01, 0.01, 0.03, 0.06, 0.22, 0.25, 0.03, 0.02, 0.03, 0.3},
		{0.01, 0.01, 0.02, 0.04, 0.02, 0.3, 0.57, 0.03, 0.01, 0.04, 0.66},
		{0.02, 0.01, 0.01, 0.01, 0.1, 0.17, 0.25, 0.02, 0.01, 0.02, 0.18},
		{0.02, 0.02, 0.01, 0.02, 0.23, 0.16, 0.26, 0.04, 0.01, 0.02, 0.19},
		{0.02, 0.01, 0.01, 0.03, 0.09, 0.16, 1.02, 0.03, 0.02, 0.03, 0.52},
		{0.01, 0.01, 0.01, 0.01, 0.07, 0.24, 0.47, 0.02, 0.03, 0.03, 0.27},
		{0.02, 0.01, 0.01, 0.04, 0.05, 0.18, 0.31, 0.05, 0.0, 0.02, 0.16},
		{0.01, 0.01, 0.01, 0.01, 0.07, 0.18, 0.21, 0.05, 0.02, 0.03, 0.39},
		{0.01, 0.03, 0.01, 0.02, 0.09, 0.13, 0.11, 0.05, 0.01, 0.03, 0.21},
		{0.01, 0.01, 0.01, 0.01, 0.06, 0.25, 0.27, 0.02, 0.4, 0.03, 0.13},
		{0.01, 0.01, 0.03, 0.01, 0.06, 0.25, 0.31, 0.05, 0.01, 0.03, 0.5}
	};

	//things that we need::
	//improved energy resolution
		//very noticeable high-E tails - yes this is pileup. Only accept single hits.
		//non-linearity correction:: uncorrected appears *better* ?? maybe due to extrapolation above ~2MeV with non-zero gradient (a TGraph 'feature')
			//abandon NL correction for now.
	//channel 61 looks very bad, get rid of this? Yes. Channels also 23, 27, 49 have poor resolution.

	//should check a subrun from each run to make sure response doesn't change.
//	if(fScep){
//	for(auto i = 0; i < fScep->GetMultiplicity(); ++i){
//	   double se = fScep->GetSceptarHit(i)->GetEnergy();
//	   fH1["sE"]->Fill(se);
//	}
//	}
	//-----------------------------------------------------------------------------------------
	//SINGLE HITS------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------------
	for(auto i = 0; i < fGrif->GetMultiplicity(); ++i){
		bhit = false;
		double k1 = fGrif->GetGriffinHit(i)->GetKValue();
		int ch1 = fGrif->GetGriffinHit(i)->GetArrayNumber()-1;
		
		if(k1<700) continue;
		//if(ch1==61) continue;

		double e1 = fGrif->GetGriffinHit(i)->GetEnergy();
		double q1 = fGrif->GetGriffinHit(i)->GetCharge();
		double t1 = fGrif->GetGriffinHit(i)->GetTimeStamp();
		double cfd1 = fGrif->GetGriffinHit(i)->GetTime();
	
		//g
		//fH2["gEchan"]->Fill(ch1,e1);
		//fH2["gQchan"]->Fill(ch1,q1);
		//fH1["gE"]->Fill(e1);

		//gb
//		if(fScep){
//		for(auto j = 0; j < fScep->GetMultiplicity(); ++j){
//		   double se = fScep->GetSceptarHit(j)->GetEnergy();
//		   double st = fScep->GetSceptarHit(j)->GetTimeStamp();
//		   if(bhit) continue;
//		   if(se < bthresh) continue;
//
//		   fH1["gb_dt"]->Fill(t1-st);
  //         	   if(CoincidenceCondition(fGrif->GetGriffinHit(i),fScep->GetSceptarHit(j))){
//		      fH1["gb"]->Fill(e1);
//		      bhit = true;
//		   }//end beta-gamma coincidence
//		}//end Sceptar Event
//		}
	}//end gamma1

	//-----------------------------------------------------------------------------------------
	//ADDBACK HIT------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------------
	for(auto i = 0; i < fGrif->GetAddbackMultiplicity(); ++i){
		abhit = false;
		double k1 = fGrif->GetAddbackHit(i)->GetKValue();
		int ch1 = fGrif->GetAddbackHit(i)->GetArrayNumber()-1;
	
		double e1 = fGrif->GetAddbackHit(i)->GetEnergy();
		double t1 = fGrif->GetAddbackHit(i)->GetTimeStamp();
		double cfd1 = fGrif->GetAddbackHit(i)->GetTime();
		double d1 = fGrif->GetAddbackHit(i)->GetDetector();
		double c1 = fGrif->GetAddbackHit(i)->GetCrystal();
		TVector3 v1 = (TGriffin::GetPosition(d1,c1,dist[int(d1-1)]));

		if(k1<700) fH1["aE_pileup"]->Fill(e1);
		//if(ch1==61) continue;	

		//a
		fH1["aE"]->Fill(e1);

		//aa - for summing corrections only
		for(auto j = 0; j < fGrif->GetAddbackMultiplicity(); ++j){
		   if(i==j) continue;
		   double k2 = fGrif->GetAddbackHit(j)->GetKValue();
		   int ch2 = fGrif->GetAddbackHit(j)->GetArrayNumber()-1;

		   if(k2<700) continue;
		   //if(ch2==61) continue;

		   double e2 = fGrif->GetAddbackHit(j)->GetEnergy();
		   double cfd2 = fGrif->GetAddbackHit(j)->GetTime();
		   double d2 = fGrif->GetAddbackHit(j)->GetDetector();
		   double c2 = fGrif->GetAddbackHit(j)->GetCrystal();
		   TVector3 v2 = (TGriffin::GetPosition(d2,c2,dist[int(d2-1)]));
		   double angle = v1.Angle(v2)* 180. / TMath::Pi();

		   if(TMath::Abs(angle-180.)<0.001){
		      fH1["aa_dt"]->Fill(cfd2-cfd1);
           	      if(CoincidenceCondition(fGrif->GetAddbackHit(i),fGrif->GetAddbackHit(j))){
		         fH2["aa"]->Fill(e1,e2);
		      }
		      if(BGCondition(fGrif->GetAddbackHit(i),fGrif->GetAddbackHit(j))){
		         fH2["aa_bg"]->Fill(e1,e2);
		      }
		   }
		}//end gamma2 coincidence

		if(fPaces) {
			for(auto j = 0; j < fPaces->GetMultiplicity(); ++j) {
				double pe = fPaces->GetPacesHit(j)->GetEnergy();
				double st
			}
		}

		//ab
//		if(fScep){
//		for(auto j = 0; j < fScep->GetMultiplicity(); ++j){
//		   double se = fScep->GetSceptarHit(j)->GetEnergy();
//		   double st = fScep->GetSceptarHit(j)->GetTimeStamp();
//		   if(abhit) continue;
//		   if(se < bthresh) continue;
//
//		   fH1["ab_dt"]->Fill(t1-st);
  //         	   if(CoincidenceCondition(fGrif->GetAddbackHit(i),fScep->GetSceptarHit(j))){
//		      fH1["ab"]->Fill(e1);
//		      abhit = true;
//		   }//end beta-gamma coincidence
//		}//end Sceptar Event
//		}
	}//end gamma1
}
