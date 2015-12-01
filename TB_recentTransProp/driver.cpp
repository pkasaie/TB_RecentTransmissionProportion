 //TB_RecentTransmissionProportion
//Developed: Parastu Kasaie 

// driver.cpp
// driver for dengue epidemic model

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
#include "methods.h"
#include <omp.h>
 
//we have both migration (all people) and exchange rate (for LLTB and REC) people; migrationRate is global and Fixed


#define  scenarios 1

using namespace std;
int main(int argc, char *argv[]) {

	cout << "MODEL RUN:" << endl;
	
	bool isAnalyzingDist = false;
	bool isCollectDNA = false;
	bool isCountAllStrains = true;
	REPLICATIONS =1 ;
	MAXRUNTIME =6000;
	TRANSIENT =1200;
	POPULATIONSIZE = 100000;
	MIGRATE = 0;

	//list of scenarios:
	
	//vector<double> vExRate = { .025, .1, .5, .025, .1, .5, .025, .1, .5, .025, .1, .5, .025, .1, .5, .025, .1, .5, .025, .1, .5, .025, .1, .5 };
	vector<double> vExRate = { .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1 };
	//vector<double> vLambda = { 1.66, 1.8, 1.93, 2.06, 2.15, 2.24, 2.31, 2.39, 1.44, 1.55, 1.705, 1.82, 1.95, 2.04, 2.13, 2.2, 1.24, 1.35, 1.463, 1.58, 1.7, 1.82, 1.92, 2, 0.585, 0.625, 0.665, 0.708, 0.755, 0.805, 0.86, 0.92, 0.37, 0.39, 0.41, 0.432, 0.448, 0.465, 0.49, 0.515};
	//vector<double> vSPRate = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10, 10, 10};
	vector<double> vLambda = { 1.65 	};
	vector<double> vSPRate = { .1};

	
	gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus2);
	double  lam, exR,spR; lam = exR =spR= 0;
	int tid, x;//thread's id

	//each thread reads one scenario, runs the simulation for REP times, save the results.   
	/*
#pragma omp parallel shared(vLambda,vExRate,vSPRate ) private(x,tid,rng, lam,exR)
	{
	#pragma omp for  ordered  nowait //----------schedule tasks in order
	for (x = 0; x < scenarios; x++)
	{
	#pragma omp critical
	{
	 */
	 x=0;
	 cout << x << endl;
				tid = omp_get_thread_num();
				exR = vExRate[x];
				lam = vLambda[x];
				spR = vSPRate[x]/12000.0;
				cout << "Reading thread" << tid << " L: " << lam<< "slowProgR: " << spR << endl;
				rng = gsl_rng_alloc(gsl_rng_taus2);
				gsl_rng_set(rng, time(0)*tid);
	 //	}
		methods::runSimulation(rng, tid, lam, exR,spR, isCollectDNA, isCountAllStrains);
		 
/*	}
#pragma omp barrier
		printf("finished %d\n", tid);

	} */
	 
	/*
//	cout << "plrease enter replications, MaxTime, transient, popsize: " << endl;
	//cin >> REPLICATIONS >> MAXRUNTIME >> TRANSIENT >> POPULATIONSIZE;
	REPLICATIONS = 1;
	MAXRUNTIME = 6000;
	TRANSIENT = 1200;
	POPULATIONSIZE = 100000;
	MIGRATE = 0;
	double lam, exR;
	lam = 1.463;
	exR = .025;
//	cout << "plrease enter lambda, ExRate: " << endl;
//	cin >> lam >> exR;
	gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus2);
	methods::runSimulation(rng, 0, lam, exR, false, true);
	*/
	cout << "click to end:" << endl;
	cin >> x;
	return 0;



}; 


