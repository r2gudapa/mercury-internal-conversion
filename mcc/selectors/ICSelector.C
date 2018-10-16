#define ICSelector_cxx
// The class definition in ICSelector.h has been generated automatically
// selector to create prompt and time random coincidence matrices between PACES and GRIFFIN hits
#include "ICSelector.h"

std::vector<double> xtmp, ytmp;

double PromptLo = 40.;
double PromptHi = 100.;
double RandomLo = 200.;
double RandomHi = 1000.;
double ggPrompt = 400.;
double ggRandomLo = 400.;
double ggRandomHi = 2000.;
double eSparseMin[3] = {0,0,0};
double eSparseMax[3] = {300,300,300};
const int eSparseBins[3] = {600,600,600};

bool CoincidenceCondition(TPacesHit* hit_paces, TGriffinHit* hit_griffin){	
	double ptime = hit_paces->GetTime();
	double gtime = hit_griffin->GetTime();
	if ((gtime - ptime) < 0) { // hit1 is griffin - LO
		return TMath::Abs(gtime - ptime) < PromptLo;
	} else {
		return TMath::Abs(gtime - ptime) < PromptHi;
	}
}

bool BGCondition(TPacesHit* hit_paces, TGriffinHit* hit_griffin){		
   return (TMath::Abs(hit_paces->GetTime() - hit_griffin->GetTime()) < RandomHi) && 
		(TMath::Abs(hit_paces->GetTime() - hit_griffin->GetTime()) > RandomLo); 
}


bool GGCoincidenceCondition(TGriffinHit* hit_one, TGriffinHit* hit_two) {
	return TMath::Abs(hit_one->GetTime() - hit_two->GetTime()) < ggPrompt;
}

bool GGBGCondition(TGriffinHit* hit_one, TGriffinHit* hit_two) {
	return (TMath::Abs(hit_one->GetTime() - hit_two->GetTime()) < ggRandomHi) &&
		(TMath::Abs(hit_one->GetTime() - hit_two->GetTime()) > ggRandomLo);
}

double corr = 0.;

void ICSelector::CreateHistograms() {
	//Define Histograms

        fH1["gE"] = new TH1D("gE","GRIFFIN singles; Energy (keV); Counts",gbins,gx,gy);
	fH1["pEt"] = new TH1D("pEt", "PACES singles (total); Energy (keV); Counts",pbins,px,py);
	fH1["pEp"] = new TH1D("pEp", "PACES singles (partial); Energy (keV); Counts",pbins,px,py);
	fH1["geCFD"] = new TH1D("geCFD","#gamma-PACES time difference; [ns]; Counts per ns",6000,-3000,3000);
	
	fH2["gChan"] = new TH2D("gChan", "GRIFFIN Energy vs. Channel; Channel Number; Energy (keV)",64,0,64,gbins,gx,gy);
	fH2["pChan"] = new TH2D("pChan", "PACES Energy vs. Channel; Channel Number; Energy (keV)", 5,0,5,pbins,px,py);

	// PACES - GRIFFIN COINCIDENCES ------------------------------------------------------------------------------------------------
	// all detectors
	fH2["gePromptTotal"] =
		new TH2D("gePromptTotal",
			 Form(" #gamma-e^{-}, |#Deltat_{#gamma-e^{-}}| < %.1f || %.1f: All Detectors",PromptLo,PromptHi),gbins,gx,gy,pbins,px,py);

	fH2["geBGtotal"] = 
		new TH2D("geBGtotal", 
			 Form("#gamma-e^{-}, #Deltat_{#gamma-e^{-}} = %.1f - %.1f: All Detectors",RandomHi,RandomLo),
				gbins,gx,gy,pbins,px,py);
	// detector 1
	fH2["gePromptDet1"] =
		new TH2D("gePromptDet1",
			 Form(" #gamma-e^{-}, |#Deltat_{#gamma-e^{-}}| < %.1f || %.1f: Detector 1",PromptLo,PromptHi),gbins,gx,gy,pbins,px,py);
	
	fH2["geBGdet1"] = 
		new TH2D("geBGdet1", 
			 Form("#gamma-e^{-}, #Deltat_{#gamma-e^{-}} = %.1f - %.1f: Detector 1",RandomHi,RandomLo),
				gbins,gx,gy,pbins,px,py);
	// detector 3
	fH2["gePromptDet3"] =
		new TH2D("gePromptDet3",
			 Form(" #gamma-e^{-}, |#Deltat_{#gamma-e^{-}}| < %.1f || %.1f: Detector 3",PromptLo,PromptHi),gbins,gx,gy,pbins,px,py);

	fH2["geBGdet3"] = 
		new TH2D("geBGdet3", 
			 Form("#gamma-e^{-}, #Deltat_{#gamma-e^{-}} = %.1f - %.1f: Detector 3",RandomHi,RandomLo),
				gbins,gx,gy,pbins,px,py);
	// detector 4
	fH2["gePromptDet4"] =
		new TH2D("gePromptDet4",
			 Form(" #gamma-e^{-}, |#Deltat_{#gamma-e^{-}}| < %.1f || %.1f: Detector 4",PromptLo,PromptHi),gbins,gx,gy,pbins,px,py);

	fH2["geBGdet4"] = 
		new TH2D("geBGdet4", 
			 Form("#gamma-e^{-}, #Deltat_{#gamma-e^{-}} = %.1f - %.1f: Detector 4",RandomHi,RandomLo),
				gbins,gx,gy,pbins,px,py);
	// all "good" detectors
	fH2["gePromptPartial"] =
		new TH2D("gePromptPartial",
			 Form(" #gamma-e^{-}, |#Deltat_{#gamma-e^{-}}| < %.1f || %.1f: Detectors 1, 3, 4",PromptLo,PromptHi),gbins,gx,gy,pbins,px,py);

	fH2["geBGpartial"] = 
		new TH2D("geBGpartial", 
			 Form("#gamma-e^{-}, #Deltat_{#gamma-e^{-}} = %.1f - %.1f: Detectors 1, 3, 4",RandomHi,RandomLo),
				gbins,gx,gy,pbins,px,py);
	// -----------------------------------------------------------------------------------------------------------------------------	

	fH2["ggPrompt"] =
		new TH2D("ggPrompt", 
			  Form(" #gamma-#gamma, |#Deltat_{#gamma-#gamma}| < %.1f", ggPrompt), gbins, gx, gy, gbins, gx, gy);

	fH2["ggBG"] = 
		new TH2D("ggBG",
			  Form("#gamma-#gamma, #Deltat_{#gamma-#gamma} = %.1f - %.1f", ggRandomHi, ggRandomLo), gbins, gx, gy, gbins, gx, gy);

//	fHSparse["ggeSparse"] = new THnSparseF("ggeSparse","THnSparse #gamma-#gamma-e^{-} 3D matrix; GRIFFIN-E1 (keV); GRIFFIN-E2 (keV); PACES-E (keV)",3,eSparseBins,eSparseMin,eSparseMax);

	//Send histograms to Output list to be added and written.
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

void ICSelector::FillHistograms() {
	if(fPaces) {
		for(auto p = 0; p < fPaces->GetMultiplicity(); ++p) {
			auto phit = fPaces->GetPacesHit(p);
		 	double pdep = phit->GetEnergy();
			int pdet = phit->GetDetector();
			double ptime = phit->GetTime();
		
			// non-linearity correction
			switch (pdet) {
				case 0: 
					for (int p = 0; p < int(fPACESEnergyCorrectionA.size()); p++) {
						xtmp.push_back(fPACESEnergyCorrectionA[p].first); ytmp.push_back (fPACESEnergyCorrectionA[p].second);
					}
					break;
				case 1:	
					for (int p = 0; p < int(fPACESEnergyCorrectionB.size()); p++) {
						xtmp.push_back(fPACESEnergyCorrectionB[p].first); ytmp.push_back (fPACESEnergyCorrectionB[p].second);
					}
					break;
				case 2:
					for (int p = 0; p < int(fPACESEnergyCorrectionC.size()); p++) {
						xtmp.push_back(fPACESEnergyCorrectionC[p].first); ytmp.push_back (fPACESEnergyCorrectionC[p].second);
					}
					break;
				case 3:	
					for (int p = 0; p < int(fPACESEnergyCorrectionD.size()); p++) {
						xtmp.push_back(fPACESEnergyCorrectionD[p].first); ytmp.push_back (fPACESEnergyCorrectionD[p].second);
					}
					break;
				case 4:
					for (int p = 0; p < int(fPACESEnergyCorrectionE.size()); p++) {
						xtmp.push_back(fPACESEnergyCorrectionE[p].first); ytmp.push_back (fPACESEnergyCorrectionE[p].second);
					}
					break;
					
			}
			TGraph *tmpgraph = new TGraph(xtmp.size(), &(xtmp[0]), &(ytmp[0]));
			corr = tmpgraph->Eval(pdep);
			double pdepc = pdep - corr;
			
			fH1["pEt"]->Fill(pdepc);
			if ((pdet == 1) || (pdet == 3) || (pdet == 4)) {
				// fill "good" paces detector sum
				fH1["pEp"]->Fill(pdepc);
			}

			delete tmpgraph; xtmp.clear(); ytmp.clear();
			fH2["pChan"]->Fill(pdet, pdepc);
		
			if(fGrif) {
				for(auto g = 0; g < fGrif->GetMultiplicity(); ++g) {
					auto ghit = fGrif->GetGriffinHit(g);
		 			double gdep = ghit->GetEnergy();
					double gtime = ghit->GetTime();
					int gdet = ghit->GetArrayNumber() - 1;
	
					fH1["geCFD"]->Fill((gtime-ptime));

					if(CoincidenceCondition(phit, ghit)) {
						fH2["gePromptTotal"]->Fill(gdep, pdepc);

						if((pdet == 1) || (pdet == 3) || (pdet == 4)) {
							
							fH2["gePromptPartial"]->Fill(gdep, pdepc);

							if(pdet == 1) {
								fH2["gePromptDet1"]->Fill(gdep, pdepc);

							} else if (pdet == 3) {
								fH2["gePromptDet3"]->Fill(gdep, pdepc);

							} else {
								fH2["gePromptDet4"]->Fill(gdep, pdepc);
							
							}
						}		
					}

					if(BGCondition(phit, ghit)) {
						fH2["geBGtotal"]->Fill(gdep, pdepc);
						
						if((pdet == 1) || (pdet == 3) || (pdet == 4)) {
							
							fH2["geBGpartial"]->Fill(gdep, pdepc);

							if(pdet == 1) {
								fH2["geBGdet1"]->Fill(gdep, pdepc);

							} else if (pdet == 3) {
								fH2["geBGdet3"]->Fill(gdep, pdepc);

							} else {
								fH2["geBGdet4"]->Fill(gdep, pdepc);
							
							}
						}		
					}
 
					for(auto g2 = 0; g2 < fGrif->GetMultiplicity(); ++g2) {
						auto ghit2 = fGrif->GetGriffinHit(g2);
		 				double gdep2 = ghit2->GetEnergy();
						double gtime2 = ghit2->GetTime();
						int gdet2 = ghit2->GetArrayNumber() - 1;

						double coords[3] = {gdep, gdep2, pdepc};
						/*if(GGCoincidenceCondition(ghit, ghit2)) {
							if(CoincidenceCondition(phit, ghit2) && CoincidenceCondition(phit, ghit)) {
								fHSparse["ggeSparse"]->Fill(coords, 1);
							}
						}*/
					}
				}	 
			}
		} 
	} 
	if(fGrif) {
		for(auto g1 = 0; g1 < fGrif->GetMultiplicity(); g1++) {
			auto ghit1 = fGrif->GetGriffinHit(g1);
		 	double gdep1 = ghit1->GetEnergy();
			double gtime1 = ghit1->GetTime();
			int gdet1 = ghit1->GetArrayNumber() - 1;

			fH1["gE"]->Fill(gdep1);
			fH2["gChan"]->Fill(gdet1, gdep1);

			for(auto g2 = 0; g2 < fGrif->GetMultiplicity(); g2++) {
				auto ghit2 = fGrif->GetGriffinHit(g2);
				double gdep2 = ghit2->GetEnergy();

				if(GGCoincidenceCondition(ghit1, ghit2)) {
					fH2["ggPrompt"]->Fill(gdep1, gdep2);
				}
				if(GGBGCondition(ghit1, ghit2)) {
					fH2["ggBG"]->Fill(gdep1, gdep2);
				}
			}
		}
	}
}		



