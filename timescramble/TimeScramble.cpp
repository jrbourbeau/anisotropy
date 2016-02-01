
// 11/2015
// James Bourbeau
// Time scrambling code for CR background

// Standard C++ libraries
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iterator>

// ROOT libraries
#include <TMath.h>
#include <TStopwatch.h>
#include <TChain.h>
#include <TH1.h>
#include <TRandom.h>

// Boost headers
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// Custom header files
#include <SimpleDST.h>

// healpix-cxx headers
#include <healpix-cxx/cxxsupport/fitshandle.h>
#include <healpix-cxx/healpix/healpix_map.h>
#include <healpix-cxx/healpix/healpix_map_fitsio.h>
#include <healpix-cxx/cxxsupport/pointing.h>

// astro headers
#include <astro/astro.h>
#include <astro/Astro_Time.h>
#include <astro/Astro_Coords.h>
#include <astro/Astro_Detector.h>

// photospilne headers
#include <photospline/splinetable.h>
#include <photospline/bspline.h>

using namespace std;
namespace po = boost::program_options;

const double deg_to_rad = TMath::DegToRad();
const Double_t hours_to_days = 1 / 24.;
const Double_t minute = hours_to_days / 60.;
const Double_t second = minute / 60.;
const Double_t millisecond = 1e-3 * second;
const Double_t microsecond = 1e-6 * second;

void TimeScramble(po::variables_map variables_map, string inFile);
//void TimeScramble(po::variables_map variables_map, string inFile, string datapath);
bool FilterCut(po::variables_map variables_map, SimpleDST dst);
int ITenergyCut(po::variables_map variables_map, SimpleDST dst, vector<float> ebins);
int ITcompCut(po::variables_map variables_map, SimpleDST dst, vector<string> cbins);
int ITs125Cut(po::variables_map variables_map, SimpleDST dst, vector<float> sbins);
int ICenergyCut(po::variables_map variables_map, SimpleDST dst, struct splinetable table, double zenith, vector<float> ebins);

int main(int argc, char* argv[]){

	po::options_description description("\nAllowed options");
	description.add_options()
	// Options used for all configurations
	("help,h", "Produce help message")
	("version,v", "Display the version number")
	("overwrite,o", "Overwrite existing maps and .fit files")
	// NOTE: boost version can't parse multiple multitoken options.
	// Multiple infiles workaround by reading infiles from text file
	//("infiles", po::value< vector<string> >()->multitoken(),"")
	//("datapath", po::value<string>(), "Path to directory containing root files to be analyzed")
	("batchfile", po::value<string>(), "Text file with input root filenames")
	("batch_idx", po::value<string>(), "Line number to read in txt file")
	("outbase", po::value<string>(), "Base name for outfile")
	("outdir", po::value<string>(), "Path to directory for outfile")
	("config", po::value<string>(), "Detector configuration")
	("integrationtime", po::value<string>(), "Integration time in hours")
	("method", po::value<string>(), "Sidereal, Anti, Solar, Extended")
	// IceCube specific options
	("spline", po::value<string>(), "File containing spline tables")
	// Icetop specific options
	("filter", po::value<string>(), "Filter for IceTop data")
	("comp", po::value< vector<string> >()->multitoken(), "Comp bins")
	("sbins", po::value< vector<string> >()->multitoken(), "S125 bins")
	("emin", po::value<string>(), "Minimum reconstructed energy")
	// Options for either detector
	("ebins", po::value< vector<string> >()->multitoken(), "Energy bins")
	;

	po::variables_map variables_map;
	po::store(po::parse_command_line(argc, argv, description), variables_map);
	po::notify(variables_map);

	if (variables_map.count("help")) {
		cout << description << "\n";
		return 1;
	}
	if (variables_map.count("version")) {
		cout << "\nVerision 2.0 of TimeScramble.\n" << endl;
		return 1;
	}
	
	// Check for all necessary parameters
	vector<string> option_keys;
	option_keys.push_back("outbase");
	option_keys.push_back("config");
	option_keys.push_back("integrationtime");
	option_keys.push_back("method");
	size_t nKeys = option_keys.size();
	for (size_t i = 0; i < nKeys; ++i) {
		if (not variables_map.count(option_keys[i])) {
			cerr << "\nUsage: " << option_keys[i] << " parameter not defined.\n";
			return 1;
		}
		if (variables_map["method"].as<string>()!= "Sidereal" && variables_map["method"].as<string>()!= "Anti" 
			&& variables_map["method"].as<string>()!= "Solar" && variables_map["method"].as<string>()!= "Extended") {
			cerr << "\n Invalid method parameter '" << variables_map["method"].as<string>() << "' entered...\n" << endl;
			return 1;
		}
	}

	// Read in filelist from batch_idx element of batchfile
	ifstream batch_file(variables_map["batchfile"].as<string>().c_str());
	if(!batch_file.is_open()){	// Check that batchfile was able to be opened
		cerr << "Couldn't open batchfile..." << endl;
		return 1;
	}
	string fileList;
	vector<string> inFiles;	// inFile is vector containing names of input files to be analyzed
	if (variables_map.count("batch_idx")){ 
		int batch_idx = atoi(variables_map["batch_idx"].as<string>().c_str());
		for (int i = 0; i < batch_idx; ++i){
			getline(batch_file, fileList);
			//inFiles.push_back(fileList);
			//cout << "fileList = " << fileList << endl;
		}
		getline(batch_file, fileList);
		// Convert fileList string to vector
		istringstream iss(fileList);
		copy(istream_iterator<string>(iss),
			istream_iterator<string>(),
			back_inserter(inFiles));
	}
	else {
		while(getline(batch_file, fileList)){
			inFiles.push_back(fileList);			
		}
	}
	batch_file.close();
	
	cout << "\nInput file(s) are: " << endl;
	for (size_t i = 0; i < inFiles.size(); ++i) cout << inFiles[i] << endl;
	
	if (variables_map.count("ebins")) {
		vector<string> test = variables_map["ebins"].as< vector<string> >();
		cout << "\nEbin values:" << endl;
		for (size_t i = 0; i < test.size(); ++i) cout << " " << test[i];
		cout << endl;
	}

	if (variables_map.count("sbins")) {
		vector<string> test = variables_map["sbins"].as< vector<string> >();
		cout << "\nS125 bin values:" << endl;
		for (size_t i = 0; i < test.size(); ++i) cout << " " << test[i];
		cout << endl;
	}
	
	if (variables_map.count("comp")) {
		vector<string> test = variables_map["comp"].as< vector<string> >();
		cout << "Comp bin values:" << endl;
		for (size_t i = 0; i < test.size(); ++i) cout << " " << test[i];
		cout << endl;
	}

	/*string datapath;
	if (variables_map.count("datapath")) {
		datapath = variables_map["datapath"].as<string>();
	}*/

	for (size_t i = 0; i < inFiles.size(); ++i) {
		//cout << "\nInput file is: " << datapath+inFiles[i] << endl;
		cout << "\nInput file is: " << inFiles[i] << endl;
		//TimeScramble(variables_map, inFiles[i], datapath);
		TimeScramble(variables_map, inFiles[i]);
	}
	
	return 0;
}

void TimeScramble(po::variables_map variables_map, string inFile) {
//void TimeScramble(po::variables_map variables_map, string inFile, string datapath) {
	
	/*TStopwatch timer;
	timer.Start();*/

	// Read input parameters
	string outbase = variables_map["outbase"].as<string>();
	string outdir;
	if (variables_map.count("outdir")){
		outdir = variables_map["outdir"].as<string>();
		boost::filesystem::path p(outdir.c_str());
		if (!boost::filesystem::is_directory(p)){
			cout << "\nCannot find the output directory '" << outdir << "'...will create this directory now..." << endl;
			boost::filesystem::create_directory(p);
		}
	}
	else outdir = "";
	string config = variables_map["config"].as<string>();
	string method = variables_map["method"].as<string>();
	const Int_t nInt = atoi(variables_map["integrationtime"].as<string>().c_str());
	string input_file = inFile;
	string detector = config.substr(0,2);

	// Read in spline tables if provided
	struct splinetable table;
	if (variables_map.count("spline")) {
		string splineFile = variables_map["spline"].as<string>();
		readsplinefitstable(splineFile.c_str(), &table);
	}

	// Energy binning setup
	unsigned int nMaps = 1;
	vector<float> ebins, sbins;
	vector<string> ebinstr, sbinstr, cbins;
	if (variables_map.count("ebins")) {
		ebinstr = variables_map["ebins"].as<vector <string> >();
		for (size_t n = 0; n < ebinstr.size(); ++n) ebins.push_back(atof(ebinstr[n].c_str()));
		nMaps = ebins.size() - 1;
	}
	if (variables_map.count("sbins")) {
		sbinstr = variables_map["sbins"].as<vector <string> >();
		for (size_t n = 0; n < sbinstr.size(); ++n) sbins.push_back(atof(sbinstr[n].c_str()));
		nMaps = sbins.size() - 1;
	}
	if (variables_map.count("comp")) {
		cbins = variables_map["comp"].as<vector <string> >();
		nMaps = cbins.size();
	}

	Astro::IceCubeDetector ice;
	Astro::LocalCoord local;
	Astro::EqCoord eqApparent;

	int NSide = 64;

	// Allow for a map for each energy/composition bin
	cout << "\nNumber of maps: " << nMaps << endl;
	vector< Healpix_Map<double> > LocalMapInt(nMaps);
	vector< Healpix_Map<double> > LocalMap(nMaps);
	vector< Healpix_Map<double> > DataMapInt(nMaps);
	vector< Healpix_Map<double> > DataMap(nMaps);
	vector< Healpix_Map<double> > BGMap(nMaps);

	for (unsigned int i = 0; i < nMaps; ++i) {
		LocalMapInt[i].SetNside(NSide, RING);
		LocalMapInt[i].fill(0.);
		LocalMap[i].SetNside(NSide, RING);
		LocalMap[i].fill(0.);
		DataMapInt[i].SetNside(NSide, RING);
		DataMapInt[i].fill(0.);
		DataMap[i].SetNside(NSide, RING);
		DataMap[i].fill(0.);
		BGMap[i].SetNside(NSide, RING);
		BGMap[i].fill(0.);
	}

	pointing sphereDir;
	int pixelID;

	stringstream sstr;
	sstr.str("");
	
	// Determine date from the input_file name
	//
	// IMPORTANT: FILES MUST END WITH "yyyy-mm-dd.root"!!!!
	//
	size_t idx0 = input_file.find_last_of('.') - 10;
	string year   = input_file.substr(idx0, 4);
	string month  = input_file.substr(idx0+5, 2);
	string day    = input_file.substr(idx0+8, 2);
	string inbase = input_file.substr(0, idx0);
	string yyyymmdd = year+"-"+month+"-"+day;

	int yyyy = atoi(year.c_str());
	int mm = atoi(month.c_str());
	int dd = atoi(day.c_str());

	// Get info from previous and next day files
	Astro::Time time(yyyy, mm, dd, 0, 0, 0);
	double mjd_input_file = time.GetMJD();

	time.SetTime(mjd_input_file - 1);	// Previous day
	string date_prev = time.GetCalendarString();
	string year_prev = date_prev.substr(date_prev.find('/')-4, 4);
	string month_prev = date_prev.substr(date_prev.find('/')+1, 2);
	string day_prev = date_prev.substr(date_prev.find('/')+4, 2);
	string yyyymmdd_prev = year_prev + "-" + month_prev + "-" + day_prev;
	
	string file_prev = inbase + yyyymmdd_prev + ".root";
	//if ((detector=="IT") && (year_prev!=year) && (detector!="IT81-III")) {
	//	idx0 = file_prev.find('/'+year+'/') + 1;
	//	file_prev.replace(idx0, 4, year_prev);
	//}
	ifstream ifile_prev(file_prev.c_str());
	//ifstream ifile_prev((datapath+file_prev).c_str());
	if(!ifile_prev.is_open()){	// Check that ifile_prev was able to be opened
		cerr << "Couldn't open ifile_prev" << endl;
	}

	time.SetTime(mjd_input_file + 1);	// Next day
	string date_post = time.GetCalendarString();
	string year_post = date_post.substr(date_post.find('/')-4, 4);
	string month_post = date_post.substr(date_post.find('/')+1, 2);
	string day_post = date_post.substr(date_post.find('/')+4, 2);
	string yyyymmdd_post = year_post + "-" + month_post + "-" + day_post;

	string file_post = inbase + yyyymmdd_post + ".root";
//	if ((detector=="IT") && (year_post!=year) && (detector!="IT81-III")) {
//		idx0 = file_post.find('/'+year+'/') + 1;
//		file_post.replace(idx0, 4, year_post);
//	}
	ifstream ifile_post(file_post.c_str());
	//ifstream ifile_post((datapath+file_post).c_str());
	if(!ifile_post.is_open()){	// Check that ifile_post was able to be opened
		cout << "ifile_post = " << file_post.c_str() << endl;
		//cout << "ifile_post = " << (datapath+file_post).c_str() << endl;
		cerr << "Couldn't open ifile_post" << endl;
	}

	const char* masterTree;
	if (detector == "IC") masterTree = "CutDST";
	if (detector == "IT") masterTree = "master_tree";

	// Initialize the chain and read data
	TChain *cutDST = new TChain(masterTree);
	if (ifile_prev) {
		cutDST->Add(file_prev.c_str());
		//cutDST->Add((datapath+file_prev).c_str());
		cout << "Open TChain: " << yyyymmdd_prev.c_str() << endl;
	}
	cutDST->Add(inFile.c_str());
	//cutDST->Add((datapath+inFile).c_str());
	cout << "Open TChain: " << yyyymmdd.c_str() << endl;
	if (ifile_post) {
		cutDST->Add(file_post.c_str());
		//cutDST->Add((datapath+file_post).c_str());
		cout << "Open TChain: " << yyyymmdd_post.c_str() << endl;
	}
	
	SimpleDST dst(cutDST, config);

	//cout << "Number of chained files: " << cutDST->GetNtrees() << endl;

	Long64_t nEntries = cutDST->GetEntries();	// Sum of number of entriesfor each TTree in the TChain 
	vector<Long64_t> nEvents(nMaps, 0);
	vector<Long64_t> nUsedEvents(nMaps, 0);

	const int nBGResample = 20;
	const double alpha = 1. / nBGResample;
	const double pi = TMath::Pi();

	Double_t startMJD, mjd2, mjd1=0;
	double zenith, azimuth, theta, phi, rndMJD;
	bool pass_cuts;

	// Integration time
	//const Double_t dt = nInt * hours_to_days;
	Double_t dt = nInt * hours_to_days;
	cout << "Integration time = " << nInt << " (hours) " << dt << " (day)\n";
	cout << "Reading " << nEntries << " entries...\n";

	cutDST->GetEntry(0);
	mjd1 = mjd_input_file;
	mjd2 = mjd1 + 1;
	startMJD = mjd1;
	
	// Setup histograms for storing time information
	vector<TH1D*> histMJD(nMaps);
	const char* histName;
	for (unsigned int i = 0; i < nMaps; ++i) {
		sstr.str("");
		sstr << "histMJD_" << i;
		histName = sstr.str().c_str();
		histMJD[i] = new TH1D(histName, ";modified julian day;events", Int_t((mjd2 - mjd1) / (10. * second)), mjd1, mjd2);
	}

	// Track the local coordinates
	vector< vector<Double_t> > LocCoord_theta(nMaps);
	vector< vector<Double_t> > LocCoord_phi(nMaps);
	
	double percent_completed_indicator = 0;
	int progress_bar_width = 70;

	for (Long64_t jentry = 0; jentry < nEntries; ++jentry) {
		
		cutDST->GetEntry(jentry);

		// Print progress bar to the console
		if (double(jentry)/nEntries >= percent_completed_indicator) {
			if (percent_completed_indicator == 0) cout << "\n" << "Progress:" << endl;
			cout << "[";
			int pos = progress_bar_width*percent_completed_indicator;
			for (int i = 0; i < progress_bar_width; ++i) {
				if (i < pos) cout << "=";
				else if (i == pos) cout << ">";
				else cout << " ";
			}
			cout << "] " << int(percent_completed_indicator * 100.0) << " %\r" << endl;
			percent_completed_indicator += 0.1;
		}

		if (detector == "IC") {
			zenith = dst.LLHZenithDeg * deg_to_rad;
			azimuth = dst.LLHAzimuthDeg * deg_to_rad;
		}
		if (detector == "IT") {
			zenith = dst.Zenith;
			azimuth = dst.Azimuth;
		}

		// Basic time check
		if (dst.ModJulDay < mjd1) continue;	
			// ISSUE:
			// This should skip all events from the previous file
			// If so, why did we load the previous file to begin with?
		
		//
		// Additional checks on data
		//
		pass_cuts = true;
		int mapIdx = 0;

		// Reconstruction cuts
		if (!dst.isReco || zenith != zenith || azimuth != azimuth) pass_cuts = false;

		// IceTop filter cut
		if (variables_map.count("filter")) {
			bool temp = FilterCut(variables_map, dst);
			if (not temp)
			pass_cuts = false;
		}
		
		// Energy cuts for IceTop and IceCube
		if (detector == "IT" && variables_map.count("ebins")) mapIdx = ITenergyCut(variables_map, dst, ebins);
		if (variables_map.count("spline"))	mapIdx = ICenergyCut(variables_map, dst, table, zenith, ebins);
		if (variables_map.count("sbins")) mapIdx = ITs125Cut(variables_map, dst, sbins);
		if (variables_map.count("comp")) mapIdx = ITcompCut(variables_map, dst, cbins);

		if (mapIdx == -1) pass_cuts = false;

		// If event passes cuts and within integration time, add to skymaps
		if (pass_cuts && dst.ModJulDay <= (startMJD+dt)) {

			// Store local coordinates
			++nEvents[mapIdx];
			LocCoord_theta[mapIdx].push_back(zenith);
			LocCoord_phi[mapIdx].push_back(azimuth);
			sphereDir.theta = zenith;
			sphereDir.phi = azimuth;
			pixelID = LocalMapInt[mapIdx].ang2pix(sphereDir);
			LocalMapInt[mapIdx][pixelID] += 1.0;

			// Calculate equatorial coordinates
			time.SetTime(dst.ModJulDay);
			local.SetLocalRad(zenith, azimuth);
			if (method == "Anti") eqApparent = ice.LocalToEquatorial_FromAntiSid(local, time);
			//if (method == "Extended") eqApparent = ice.LocalToEquatorial_FromExtendedSid(local, time);
			if (method == "Sidereal") eqApparent = ice.LocalToEquatorial(local, time);
			if (method == "Solar") eqApparent = ice.LocalToEquatorial_FromSolar(local, time);
			
			// Write to map
			sphereDir.theta = pi/2. - eqApparent.GetDecRad();
			sphereDir.phi = eqApparent.GetRaRad();

			// Solar coordinates need a 180 deg flip in phi (definition difference)
			if (method == "Solar") sphereDir.phi -= pi;
			while (sphereDir.phi < 0) sphereDir.phi += 2.*pi;
			pixelID = DataMapInt[mapIdx].ang2pix(sphereDir);
			DataMapInt[mapIdx][pixelID] += 1.0;

			// Store time
			histMJD[mapIdx]->Fill(dst.ModJulDay);
		}


		if ((dst.ModJulDay > (startMJD + dt)) || (dst.ModJulDay > mjd2) || (jentry + 1 == nEntries)) {

			for (unsigned int mEntry = 0; mEntry < nMaps; ++mEntry) {

				if (variables_map.count("ebins")) {
					cout << "Working on energy bin " << ebinstr[mEntry] << "-" << ebinstr[mEntry+1] << "GeV..." << endl;
				}
				if (variables_map.count("sbins")) {
					cout << "Working on s125 bin " << sbinstr[mEntry] << " to " << sbinstr[mEntry+1] << "s125..." << endl;
				}
				nUsedEvents[mEntry] += (nEvents[mEntry]);

				// Scramble the time
				cout << "  Scrambling time for (" << nBGResample << "x " << nEvents[mEntry] << " events)..." << endl;
				gRandom->SetSeed(0);

				for (Long64_t iEntry = 0; iEntry < (Long64_t)(nEvents[mEntry]); iEntry++) {

					// Get local coordinates
					theta = LocCoord_theta[mEntry][iEntry];
					phi = LocCoord_phi[mEntry][iEntry];
					local.SetLocalRad(theta, phi);

					for (int k=0; k<nBGResample; ++k) {

						// Generate new equatorial coordinates
						rndMJD = histMJD[mEntry]->GetRandom();
						time.SetTime(rndMJD);
						if (method == "Anti") eqApparent = ice.LocalToEquatorial_FromAntiSid(local, time);
						//if (method == "ext") eqApparent = ice.LocalToEquatorial_FromExtendedSid(local, time);
						if (method == "Sidereal") eqApparent = ice.LocalToEquatorial(local, time);
						if (method == "Solar") eqApparent = ice.LocalToEquatorial_FromSolar(local, time);

						// Write to map
						sphereDir.theta = (pi/2. - eqApparent.GetDecRad());
						sphereDir.phi = eqApparent.GetRaRad();
						if (method == "Solar") sphereDir.phi -= pi;
						while (sphereDir.phi < 0) sphereDir.phi += 2.*pi;

						pixelID = DataMapInt[mEntry].ang2pix(sphereDir);
						BGMap[mEntry][pixelID] += 1.0;
					}
				}

				// Update the data map for this time interval
				for (int i = 0; i < DataMap[mEntry].Npix(); ++i) {
					DataMap[mEntry][i] += DataMapInt[mEntry][i];
					LocalMap[mEntry][i] += LocalMapInt[mEntry][i];
				}
			}	// End of mEntry loop
			
			//cout << "jentry : " << jentry << endl;
			jentry = jentry - 1;
			startMJD += dt;
			
			if (startMJD + second >= mjd2) break;
			else {
				cout << "new startMJD :" << setprecision(12) << startMJD << endl;
				for (unsigned kEntry = 0; kEntry < nMaps; ++kEntry) {
					LocCoord_phi[kEntry].erase(LocCoord_phi[kEntry].begin(), LocCoord_phi[kEntry].end());
					LocCoord_theta[kEntry].erase(LocCoord_theta[kEntry].begin(), LocCoord_theta[kEntry].end());
					histMJD[kEntry]->Reset();
					DataMapInt[kEntry].fill(0);
					LocalMapInt[kEntry].fill(0);
					nEvents[kEntry] = 0;
				}
			}
		}
	} // End of jentry loop
  
	for (unsigned int m=0; m<nMaps; ++m) {

		// Scale background map
		for (int i=0; i<BGMap[m].Npix(); ++i) BGMap[m][i] *= alpha;

		cout << "Read " << nEntries << " events" << "\n" << "Used " << nUsedEvents[m] << " events" << endl;

		// Save BG, Data, and Local maps in one file
		arr<std::string> colname(3);	// Using arr instead of vector because that's what prepare_Healpix_fitsmap requires
		colname[0] = "data map";
		colname[1] = "background map";
		colname[2] = "local map";

		sstr.str("");
		sstr << outdir + "/" + outbase;
		if (variables_map.count("ebins")) sstr << "_" << ebinstr[m] << "-" << ebinstr[m+1] << "GeV";
		if (variables_map.count("sbins")) sstr << "_" << sbinstr[m] << "to" << sbinstr[m+1] << "s125";
		if (variables_map.count("comp")) sstr << "_" << cbins[m];
		sstr << "_" << yyyymmdd << ".fits";
		// If overwrite option given, and file with same name exists, remove the existing file
		boost::filesystem::path outfile_path((sstr.str()).c_str());
		if (variables_map.count("overwrite") && boost::filesystem::exists(outfile_path)){
			boost::filesystem::remove(outfile_path);
		}
		
		fitshandle fitsOut;
		fitsOut.create(sstr.str().c_str());

		fitsOut.add_comment("Maps: data, bg, local");
		prepare_Healpix_fitsmap(fitsOut, DataMap[m], FITSUTIL<double>::DTYPE, colname);
		fitsOut.write_column(1, DataMap[m].Map());
		fitsOut.write_column(2, BGMap[m].Map());
		fitsOut.write_column(3, LocalMap[m].Map());
		fitsOut.close();
	}
	// Clean up
	delete cutDST;

	// To avoid messages like "Warning in <TROOT::Append>: Replacing existing TH1: histMJD_0 (Potential memory leak)"	
	// This might lead to some errors...not sure yet. 
	for (unsigned int i=0; i < nMaps; ++i) delete histMJD[i];

	/*timer.Stop();
	printf("RT=%7.3f s, Cpu=%7.3f s\n",timer.RealTime(),timer.CpuTime());*/

}

bool FilterCut(po::variables_map variables_map, SimpleDST dst) {

	string filter = variables_map["filter"].as<string>();
	string config = variables_map["config"].as<string>();
	if (filter=="STA3" && dst.isSTA3) return true;
	if (config=="IT59" || config=="IT73" || config=="IT81") {
		if ((filter=="STA8" && dst.isSTA8) || (filter=="NotSTA8" && dst.isSTA3 && !dst.isSTA8)) return true;
	}
	if (config=="IT81-II" || config=="IT81-III") {
		if ((filter=="STA8" && dst.nStations>=8) || (filter=="NotSTA8" && dst.isSTA3 && dst.nStations<8)) return true;
	}
	return false;
}

int ICenergyCut(po::variables_map variables_map, SimpleDST dst, splinetable t, double zenith, vector<float> ebins) {

	// Setup basic parameters
	double x = cos(zenith);
	double y = log10(dst.NChannels);

	// Boundary check (energy cut tables go to 0.3 in cos(zenith))
	if (x < 0.3) return -1;

	// Catch additional outliers
	double coords[2] = {x, y};
	int centers[t.ndim];
	if (tablesearchcenters(&t, coords, centers) != 0) {
		cout << "Variables outside of table boundaries" << endl;
		cout << "x: " << x << " y: " << y << endl;
		return -1;
	}

	// Calculate reconstructed energy
	double median = ndsplineeval(&t, coords, centers, 0);
	// Make sure we're in the energy bin range
	if ((median < ebins[0]) || (median > ebins.back())) return -1;

	// Get energy bin
	int ebin = 0;
	while (median > ebins[ebin+1]) ebin += 1;
	
	return ebin;
}

int ITenergyCut(po::variables_map variables_map, SimpleDST dst, vector<float> ebins) {

	// Get most likely energy value
	//
	// Original line commented out. dst.pEnergy and dst.fEnergy don't exist.
	//
	//original: double llhEnergy = (dst.pLLH >= dst.fLLH) ? dst.pEnergy : dst.fEnergy;
	//double llhEnergy = (dst.pLLH >= dst.fLLH) ? dst.pLLH : dst.fLLH;
	double llhEnergy = dst.llhEnergy;
	double logEnergy = log10(llhEnergy);

	// Make sure we're in the energy bin range
	if ((logEnergy < ebins[0]) || (logEnergy > ebins.back())) return -1;

	// Get energy bin
	int ebin = 0;
	while (logEnergy > ebins[ebin+1]) ebin += 1;

	return ebin;
}

int ITcompCut(po::variables_map variables_map, SimpleDST dst, vector<string> cbins) {

	// Get most likely composition
	Text_t llh_comp = dst.comp;
	float llh_energy = log10(dst.llhEnergy);

	// Read in minimum energy if provided
	float emin = -1.;
	if (variables_map.count("emin")) emin = atof(variables_map["emin"].as<string>().c_str());
	if ((emin == -1) || (llh_energy < emin)) return -1;

	// Get composition bin
	int cbin = -1;
	for (unsigned int i=0; i < cbins.size(); ++i) {
		if (llh_comp == cbins[i]) cbin = i;
	}

	return cbin;
}

int ITs125Cut(po::variables_map variables_map, SimpleDST dst, vector<float> sbins) {

	// Get desired s125 value
	double s125 = (dst.nStations >= 5) ? dst.s125 : dst.ss125;
	double logS125 = log10(s125);

	// Make sure we're in the bin range
	if ((logS125 < sbins[0]) || (logS125 > sbins.back())) return -1;

	// Get s125 bin
	int sbin = 0;
	while (logS125 > sbins[sbin+1]) sbin += 1;

	return sbin;
}

