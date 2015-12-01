// Community.h
// Manages relationship between people, locations, and mosquitoes
#ifndef __COMMUNITY_H
#define __COMMUNITY_H


#include <string>
#include <vector>
#include <list>

using namespace std;

class Person;


class Community {
public:

	void restart(gsl_rng *rng);
	vector<int> runModel_incVec(gsl_rng *rng);//runs the model and retunr the vector of incidences
	double runModel_incMean(gsl_rng *rng);//runs the model and resurns the mean incidence after the transient period
	void runModel_MeanOutputs(gsl_rng *rng, int rep, bool isCollectDNA, bool isCountAllStrains);
	void  runModel_FullOutputs(gsl_rng *rng, int rep, bool isCollectDNA, bool isCountAllStrains);
	void  runModel_FullOutputsDist(gsl_rng *rng, int rep, bool isCollectDNA, bool isCountAllStrains, bool isAnalyzingDist);
	void clusterDistAnalysis();


	Community(int tid, double lambda, double exRate);
	Community(int tid, double lambda, double exRate, double spRate);
	virtual ~Community();

	void createPopulation(gsl_rng *rng, int N);//create a hypo.population of size N
	void initialInfection(gsl_rng *rng, double initialStates[]);//set initial healthstatus with unique strains
	void tick(gsl_rng *rng, bool isCollectDNA);//execute model's events at each step


	//populaiton dynamics..................................................................................
	void contactPoisson(gsl_rng *rng);
	void progress(gsl_rng *rng, bool isCollectDNA);
	void deathBirth(gsl_rng *rng);//people die due to old ag or natural mortality;new-borns are added to the population
	void migrationSwaping(gsl_rng *rng);//swap each Migout with a migIn from our population
	void exchange(gsl_rng *rng);//swap each Migout with a migIn from our population




	//outputs:----------------------------------------------------------------------------
	vector<int> returnAgeDist();
	double returnVectorAverage(vector<int> vec);
	void saveRunOutputs(int r);//saves this run's ourputs
	void saveDieaseDuration(bool b);
	void printDiseaseVector();
	//vector<int>  generateUniqueIDs(gsl_rng *rng, int MaxBound, int n);
	
	int  countTotalCirculatingStrains();


	//helpers--------------
	vector<int> returnDiseaseVector();//check all, returns current disease vector
	//vector<vector<int>>  returnOutputs();

	
	vector<vector<int>> returnOutputs(){ return vvOut; }
	vector<vector<int>> returnDNAStrains(){ return vvDNAStrains; }
	vector<int> returnIncidence(){ return vvOut[4]; }

private:
	double _LAMBDA;
	double _EXCHANGERATE;
	double _SPRATE;

	vector<Person*> _vPopulation; //vector of population
	//vector <vector<Person*>> _vPopulationAgeCohort; //vector of population's age groups and everyone in each group
	int _nMonth;           //global variable across community to keep track of time-steps
	int _nLastStrain;//global variable across community to keep track of strains
	int _nNextID; // unique ID to assign to the next Person allocated
	int _nThreadId;

	int numContacts;//number of contacts successful or not
	int numSusExposures;
	int numEnterELTB;//num of successful contacts leading to infection and moving to ELTB
	int	numIncidence;
	int numFastProg;
	int numSlowProg;
	int numRelapses;
	int numImmATB;//number of immigrants to population who are ATB
	int numMigrated;
	int numFatality;
	int numMortality;
	int numOldageMortality;
	int numNewborn;
	int numATBfromREC;
	int numExchanged;

	vector<int> vTotalcirculatingStrains; //to keep track of the total #of revoving Strains in the population

	vector<int> vPeopleDisDuration;//when someone dies or recover (not mogration)


	vector<int> vDNAStrains;//1D contains DNAstrains at one point of time
	vector<vector<int>> vvDNAStrains;//2D contains DNAstrains at multiple points of time in 1 relication

	vector<vector<int>> vvOut;//2D containing all summary output vectors at one rep
	vector<vector<int>> vvAgeDists;


};
#endif
