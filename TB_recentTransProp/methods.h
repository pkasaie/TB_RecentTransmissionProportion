//Methods.h
#include <string>
#include <vector>
#include <list>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>   //for time function

#include "Person.h"
#include "Community.h"
#include "params.h"



using namespace std;

class methods{
 
public:

	static void printVector3D(vector<vector < vector<int>>> vvvOutputs);
	static void saveAgeDist(string s, bool b, vector<vector<vector<int>>> vvvAgeDists);//static methods can be used outside the scope of community
	static void saveSingleOutputAcrossReps(string s, bool b, vector<vector < vector<int>>> vvvOutputs, int outputId);
	static void SaveSimulationInfo(string s, bool b, vector<int> vRngSeeds);
	static void SaveSimulationInfo(string s, bool b);


	static int sumVec1D(vector<int> vec);
	static vector<int> sumVector2DAcrossRows(vector<vector<int>> vec);
	static vector<int> sumVector2DAcrossColumn(vector<vector<int>> vec);

	static vector<vector<double>> MeanVector3DAcrossReps(vector<vector < vector<int>>> vvvOutputs);


	static void saveVector1D(string s, bool b, vector<double> vec);
	static void saveVector1D(string s, bool b, vector<int> vec);

	static void saveVector2D(string s, bool b, vector<vector<int>> vec);
	static void saveVector2D(string s, bool b, vector<vector<double>> vec);
	static void saveVector2D(string filename, bool b, vector<vector<double>> vec, string firstline);
	static void saveVector2D(string filename, bool b, vector<vector<int>> vec, string firstline);
		
	static void saveVector3DCumulative(string s, bool b, vector<vector<vector<int>>> vvvDNAStrains);
	static void saveVector3DCumulative(string filename, bool b, vector<vector<vector<int>>> vvvDNAStrains,string firstline);

	static double meanVec1D(vector<int> vec);

	static double meanVec1D(vector<double> vec);
	static double varVec1D(vector<double> vec);

	static vector<double>  meanVec2DColumns(vector<vector<int>> vec);
	static vector<double>  meanVec2DColumns(vector<vector<double>> vec);
	static vector<double>  meanVec2DColumns(vector<vector<double>> vec,int start);
	static vector<double>  meanVec2DColumns(vector<vector<int>> vec, int start);
	static double meanVec1DStartEnd(vector<int> vec, int start, int end);

	
	static void runSimulationParallel( bool isCollectingDNA, bool isCountAllStrains);
	static double runSimulation(int reps, double lambda,int goal);//run the simulation for rep times,return mean results
	
	
	static void runSimulation(gsl_rng * rng, int tid,  double lambda, double exRate, bool isCollectingDNA, bool isCountAllStrains);
	static void runSimulation(gsl_rng * rng, int tid, double lambda, double exRate, bool isCollectingDNA, bool isCountAllStrains,bool isAnalyzingDist);
	static void runSimulation(gsl_rng * rng, int tid, double lambda, double exRate, double spRate, bool isCollectingDNA, bool isCountAllStrains);

	//runs simulation and all it's replications on a single thread
};