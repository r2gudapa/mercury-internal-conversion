void calibrate() {
	//'good' base calibration = run 9984
	//const char* cpath = "/data3/mbowry/paces/radium/data/ra228/calfiles/ra228_HPGe_online_run09984.cal";			//old with channels 37 and 61 swapped
	//const char* cpath = "/data3/mbowry/paces/radium/data/ra228/calfiles/ra228_HPGe_online_run09984_updated.cal";		//fixed
	const char* cpath = "/data3/rgudapati/run1989.cal";	//digitizer type set to GRF16, correct address for ch.37,61

	//all ra-228 BEAM ON runs + subruns----------------------------------------------------------------------
	const char* fpath = "/data3/rgudapati/Hg198/";
	const int nr = 4;
	int run[nr]={9189,9190,9191,9192};
	int sub[nr]={50,44,29,44};
	//-----------------------------------------------------------------------------------------------
	//eu-152 calibration run
	/*const char* fpath = "/data3/mbowry/paces/radium/data/eu152_r1103/";
	const int nr = 1;
	int run[nr]={10076};
	int sub[nr]={20};*/
	//-----------------------------------------------------------------------------------------------
	//ba-133 calibration run
	/*const char* fpath = "/data3/mbowry/paces/radium/data/ba133_r794/";
	const int nr = 1;
	int run[nr]={10092};
	int sub[nr]={12};*/
	//-----------------------------------------------------------------------------------------------
	//co-60 calibration run
	/*const char* fpath = "/data3/mbowry/paces/radium/data/co60_r1105/";
	const int nr = 1;
	int run[nr]={10097};
	int sub[nr]={10};*/
	//-----------------------------------------------------------------------------------------------
	//co-56 calibration run
	/*const char* fpath = "/data3/mbowry/paces/radium/data/co56_r1122/";
	const int nr = 1;
	int run[nr]={10094};
	int sub[nr]={8};*/
	//-----------------------------------------------------------------------------------------------

	for(int i=0;i<int(nr);i++){
		for(int j=0;j<sub[i];j++){
   		   TFile *rf = new TFile(Form("%sanalysis%05i_%03i.root",fpath,run[i],j),"read");
 		   TChannel::ReadCalFile(Form("%s",cpath));
		   TChannel::WriteToRoot();
		   rf->Close();
		}
	}
}
