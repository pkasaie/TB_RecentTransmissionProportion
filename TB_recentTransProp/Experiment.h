
#include <cstdlib>
#include <cstring>
#include <climits>
#include <iostream>
#include <string>
#include <iterator>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <math.h>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


using namespace std;

class Experiment{

public:
	Experiment(int simtime );


	void readExperimentInfo(string filename);

	vector<double> testShortcut(gsl_rng *rng, vector<vector<int>> vvDNAS, double coverage, int duration, int startTime);//return sampleSize and #secondary cases; do not provide info on the number of clusteres
	vector<vector<double>>	runExperiments(gsl_rng *rng, vector<vector<int>> vvDNAS, int numBootstrap);



	int T;//simulation time steps (years)
	 
	vector<int> vDuration;//durations
	vector<double> vCoverage;//coverage



	


};