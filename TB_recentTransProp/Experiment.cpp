#include <iterator>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <math.h>
#include <vector>
#include <list>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Experiment.h"
#include "methods.h"

using namespace std;

Experiment::Experiment(int simtime){
	T = simtime;
	
	
}


void Experiment::readExperimentInfo(string filename){
	ifstream file;
	file.open(filename);
	string line;
	
	while (getline(file, line)){
		if (line == "") continue;
		istringstream iss(line);
		int num1;	double num2;

		while (iss >> num1 >> num2){
			vDuration.push_back(num1);
			vCoverage.push_back(num2);
		}
	}
};

vector<double> Experiment::testShortcut(gsl_rng *rng, vector<vector<int>> vvDNAS, double coverage, int duration, int startTime){
	vector<double> vResults;
	
		//collect sample
		list<int> lSample;
		for (int t = startTime; t < startTime + duration; t++){
			for (int o = 1; o < vvDNAS[t].size(); o++){//start from element1 (element 0 is the time)
				if (gsl_rng_uniform(rng) < coverage)
					lSample.push_back(vvDNAS[t][o]);
			}
		}
		//analyse sample using n-1
		int sampleSize = lSample.size();
		lSample.sort();
		lSample.unique();
		int numSecondary = sampleSize - lSample.size();

		vResults.push_back(sampleSize);
		vResults.push_back(numSecondary);
		vResults.push_back((double)numSecondary / sampleSize);

	
	return vResults;
};

vector<vector<double>> Experiment::runExperiments(gsl_rng *rng, vector<vector<int>> vvDNAS, int numBootstrap ){
	//rows represent each experimental point
	//in each row we have the mean value of bootstraped experiments
	
	int numTests = vDuration.size();
	vector<vector<double>> vvTestsMeanResults;
	for (int i = 0; i < numTests; i++){
		double p = vCoverage[i];
		int d = vDuration[i];
		vector<vector<double>> vvBootsResults;//a set of boosts results: B *  O(3outputs)
		for (int b = 0; b < numBootstrap; b++){
			//generate new start time
			int s = gsl_rng_uniform_int(rng, T-d );
			vvBootsResults.push_back(testShortcut(rng, vvDNAS, p, d, s));
			//summarize results
			//vMeanRes
		}
		vector<double> vMean = methods::meanVec2DColumns(vvBootsResults);//mean of numBoots tests with p and d (changing starttime)
		vvTestsMeanResults.push_back(vMean);
	}
	return vvTestsMeanResults;
};


