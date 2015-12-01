//Methods.cpp
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

using namespace std;
 

void methods::saveAgeDist(string s, bool b, vector<vector<vector<int>>> vvvAgeDists){
	//compute average size of each group at each point of time across all replications
	//out:  T   * #ageGroups
	ofstream outFile;
	if (b){
		outFile.open(s);
	}
	else outFile.open(s, ios_base::app);
	int R = vvvAgeDists.size();
	//sum across replications
	vector<vector<int>> vAgeDistsSums;
	vAgeDistsSums = vvvAgeDists[0];

	for (int i = 0; i < vvvAgeDists[0].size(); i++){
		for (int j = 0; j < vvvAgeDists[0][0].size(); j++){
			for (int r = 1; r < R; r++){		//cout<<t<<" "<<o<<" "<<r<< endl;	
				vAgeDistsSums[i][j] = vAgeDistsSums[i][j] + vvvAgeDists[r][i][j];
			}
		}
	}
	for (int i = 0; i < vvvAgeDists[0].size(); i++){
		for (int j = 0; j < vvvAgeDists[0][0].size(); j++){
			outFile << (double)vAgeDistsSums[i][j] / R << " ";
		}
		outFile << endl;
	}
	outFile.close();

};
void methods::saveSingleOutputAcrossReps(string s, bool b, vector<vector < vector<int>>> vvvOutputs, int outputId){
	//writes out a specific output across time and for each replication
	//R * T

	int R = vvvOutputs.size();
	int T = vvvOutputs[0].size();

	ofstream f;
	if (b) {//newFile
		f.open(s);
	}
	else
		f.open(s, ios_base::app);

	//sum across replications

	for (int r = 0; r < R; r++){
		for (int t = 0; t < T; t++){
			f << vvvOutputs[r][t][outputId] << " ";
		}
		f << endl;
	}

	f.close();
};
void methods::SaveSimulationInfo(string s, bool b, vector<int> vRngSeeds){
/*	ofstream f;
	if (b) {//newFile
		f.open(s);
	}
	else		f.open(s, ios_base::app);
	time_t t = time(0);   // get time now
	//struct tm * now = localtime(&t);
	//f << "SystemTime: " << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' << now->tm_mday << endl;
	f << "Relications-R " << REPLICATIONS << endl;
	f << "MaxRunTime-T " << MAXRUNTIME << endl;

	f << "Lambda: " <<   endl;
	f << "MigrationRate: " <<  endl;
	f << "--------------------------------------------" << endl;
	f << "PopSize N " << POPULATIONSIZE << endl;
	f << "InitialStates - " << INITIALSTATES[0] << " " << INITIALSTATES[1] << " " << INITIALSTATES[2] << " " << INITIALSTATES[3] << " " << endl;
	f << "CumFastProgRate: " << CUMFASTPROG << endl;
	f << "SlowProgRate: " << SLOWPROGRATE << endl;
	f << "RecoveryRate: " << RECOVERYRATE << endl;
	f << "RelapseRate: " << RELAPSE << endl;
	f << "FatalityRate: " << FATALITY << endl;
	f << "MaximumInfectiousness: " << MAXINFECTIOUSNESS << endl;
	f << "MaxDuration: " << MAXINFDURATION << endl;
	f << "ImmunityEL: " << IMMUNITYELTB << endl;
	f << "ImmunityLL: " << IMMUNITYLLTB << endl;
	f << "ImmunityATB: " << IMMUNITYATB << endl;
	f << "ImmunityREC: " << IMMUNITYREC << endl;

	f << "--------------------------------------------" << endl;
	f << "RNG SEEDS: " << endl;
	for (int i = 0; i < vRngSeeds.size(); i++)
		f << vRngSeeds[i] << endl;
	f.close();*/

};
void methods::SaveSimulationInfo(string s, bool b ){
	/*ofstream f;
	if (b) {//newFile
		f.open(s);
	}
	else		f.open(s, ios_base::app);
	time_t t = time(0);   // get time now
	//struct tm * now = localtime(&t);
	//f << "SystemTime: " << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' << now->tm_mday << endl;
	f << "Relications-R " << REPLICATIONS << endl;
	f << "MaxRunTime-T " << MAXRUNTIME << endl;

	f << "Lambda: " <<  endl;
	f << "MigrationRate: " <<   endl;
	f << "--------------------------------------------" << endl;
	f << "PopSize N " << POPULATIONSIZE << endl;
	f << "InitialStates - " << INITIALSTATES[0] << " " << INITIALSTATES[1] << " " << INITIALSTATES[2] << " " << INITIALSTATES[3] << " " << endl;
	f << "CumFastProgRate: " << CUMFASTPROG << endl;
	f << "SlowProgRate: " << SLOWPROGRATE << endl;
	f << "RecoveryRate: " << RECOVERYRATE << endl;
	f << "RelapseRate: " << RELAPSE << endl;
	f << "FatalityRate: " << FATALITY << endl;
	f << "MaximumInfectiousness: " << MAXINFECTIOUSNESS << endl;
	f << "MaxDuration: " << MAXINFDURATION << endl;
	f << "ImmunityEL: " << IMMUNITYELTB << endl;
	f << "ImmunityLL: " << IMMUNITYLLTB << endl;
	f << "ImmunityATB: " << IMMUNITYATB << endl;
	f << "ImmunityREC: " << IMMUNITYREC << endl;

	f << "--------------------------------------------" << endl;
	f << "RNG SEEDS: " << endl;
	 
	f.close();*/

};

void methods::saveVector1D(string s, bool b,  vector<double> vec){
	ofstream f;
	if (b) {//newFile
		f.open(s);
	}
	else
		f.open(s, ios_base::app);
	int I = vec.size();
	for (int i = 0; i < I; i++){
		
			f << vec[i] << " ";

		
	}
	f << endl;
	f.close();
};
void methods::saveVector1D(string s, bool b, vector<int> vec){
	ofstream f;
	if (b) {//newFile
		f.open(s);
	}
	else
		f.open(s, ios_base::app);
	int I = vec.size();
	for (int i = 0; i < I; i++){

		f << vec[i] << " ";


	}
	f << endl;
	f.close();
};

void methods::saveVector2D(string s, bool b, vector<vector<double>> vec){
	ofstream f;
	if (b) {//newFile
		f.open(s);
	}
	else
		f.open(s, ios_base::app);
	int I = vec.size();
	for (int i = 0; i < I; i++){
		for (int j = 0; j < vec[i].size(); j++){
			f << vec[i][j] << " ";

		}
		f << endl;
	}
	f.close();
};
void methods::saveVector2D(string s, bool b, vector<vector<int>> vec){
	ofstream f;
	if (b) {//newFile
		f.open(s);
	}
	else
		f.open(s, ios_base::app);
	int I = vec.size();
	for (int i = 0; i < I; i++){
		for (int j = 0; j < vec[i].size(); j++){
			f << vec[i][j] << " ";

		}
		f << endl;
	}
	f.close();
};
void methods::saveVector2D(string filename, bool b, vector<vector<double>> vec, string firstline){
	ofstream f;
	if (b) {//newFile
		f.open(filename);
	}
	else
		f.open(filename, ios_base::app);
	f << firstline << endl << endl;

	int I = vec.size();
	for (int i = 0; i < I; i++){
		for (int j = 0; j < vec[i].size(); j++){
			f << vec[i][j] << " ";

		}
		f << endl;
	}
	f.close();
};
void methods::saveVector2D(string filename, bool b, vector<vector<int>> vec, string firstline){
	ofstream f;
	if (b) {//newFile
		f.open(filename);
	}
	else
		f.open(filename, ios_base::app);
	f << firstline << endl << endl;

	int I = vec.size();
	for (int i = 0; i < I; i++){
		for (int j = 0; j < vec[i].size(); j++){
			f << vec[i][j] << " ";

		}
		f << endl;
	}
	f.close();
};
void methods::saveVector3DCumulative(string s, bool b, vector<vector<vector<int>>> vvvDNAStrains){
	ofstream f;
	if (b) {//newFile
		f.open(s);
	}
	else
		f.open(s, ios_base::app);
	int R = vvvDNAStrains.size();
	int T = vvvDNAStrains[0].size();

	for (int r = 0; r < R; r++){
		for (int t = 0; t < T; t++){
			for (int s = 0; s < vvvDNAStrains[r][t].size(); s++){
				f << vvvDNAStrains[r][t][s] << " ";
			}
			f << endl;
		}

	}
	f.close();

};
void methods::saveVector3DCumulative(string filename, bool b, vector<vector<vector<int>>> vvvDNAStrains, string firstline){
	ofstream f;
	if (b) {//newFile
		f.open(filename);
	}
	else
		f.open(filename, ios_base::app);
	f << firstline << endl;
	f << endl;
	int R = vvvDNAStrains.size();
	int T = vvvDNAStrains[0].size();

	for (int r = 0; r < R; r++){
		for (int t = 0; t < T; t++){
			for (int s = 0; s < vvvDNAStrains[r][t].size(); s++){
				f << vvvDNAStrains[r][t][s] << " ";
			}
			f << endl;
		}

	}
	f.close();

};



int methods::sumVec1D(vector<int> vec){
	int sum = 0;
	int I = vec.size();
	for (int i = 0; i < I; i++){
		sum = sum + vec[i];
	}
	return sum;
};
vector<int> sumVector2DAcrossRows(vector<vector<int>> vec){
	vector<int> vSum = vec[0];
	int I = vec.size();
	int J = vec[0].size();
	for (int i = 1; i < I; i++){
		for (int j = 0; j < J; j++){
			vSum[j] = vSum[j] + vec[i][j];
		}
	}
	return vSum;
};
vector<int> sumVector2DAcrossColumn(vector<vector<int>> vec){
	vector<int> vSum;
	int sum = 0;
	int I = vec.size();
	int J = vec[0].size();
	for (int i = 0; i < I; i++){
		for (int j = 0; j < J; j++){
			sum = sum + vec[i][j];
		}
		vSum.push_back(sum);
		sum = 0;
	}
	return vSum;
};


double methods::varVec1D(vector<double> vec){

	double sum = 0;
	int I = vec.size();
	for (int i = 0; i < I; i++){
		sum = sum + vec[i];
	}
	double mean= (double) sum / I;
	sum = 0;
	 
	for (int i = 0; i < I; i++){
		sum = sum + (vec[i] - mean)*(vec[i] - mean);
	}

	return (double)sum / (I- 1);

};
double methods::meanVec1D(vector<int> vec){
	int sum = 0;
	int I = vec.size();
	for (int i = 0; i < I; i++){
		sum = sum + vec[i];
	}
	return (double)sum / I;
};
double methods::meanVec1D(vector<double> vec){
	double sum = 0;
	int I = vec.size();
	for (int i = 0; i < I; i++){
		sum = sum + vec[i];
	}
	return (double)sum / I;
};
double methods::meanVec1DStartEnd(vector<int> vec, int start, int end){
	int sum = 0;
	for (int i = start; i < end; i++){
		sum = sum + vec[i];
	}
	return (double)sum / (end - start);
};
vector<double> methods::meanVec2DColumns(vector<vector<int>> vec){
	vector<int> sum;
	int rows = vec.size();
	int cols = vec[0].size();
	sum = vec[0];
	for (int i = 1; i < rows; i++){
		for (int j = 0; j<cols; j++){
			sum[j] = sum[j] + vec[i][j];
		}
	}
	vector<double> mean;
	for (int j = 0; j < cols; j++){
		mean[j] = (double)sum[j] / rows;
	}
	return mean;
};
vector<double> methods::meanVec2DColumns(vector<vector<double>> vec){
	vector<double> sum;
	int rows = vec.size();
	int cols = vec[0].size();
	sum = vec[0];
	for (int i = 1; i < rows; i++){
		for (int j = 0; j<cols; j++){
			sum[j] = sum[j] + vec[i][j];
		}
	}
	vector<double> mean;
	for (int j = 0; j < cols; j++){
		mean.push_back((double)sum[j] / rows);
	}
	return mean;
};
vector<double> methods::meanVec2DColumns(vector<vector<double>> vec,int start){
	vector<double> sum;
	int rows = vec.size();
	int cols = vec[0].size();
	sum = vec[0];
	for (int i = start; i < rows; i++){
		for (int j = 0; j<cols; j++){
			sum[j] = sum[j] + vec[i][j];
		}
	}
	vector<double> mean;
	int denom = rows - start;
	for (int j = 0; j < cols; j++){
		mean.push_back((double)sum[j] / denom);
	}
	return mean;
};
vector<double> methods::meanVec2DColumns(vector<vector<int>> vec, int start){
	vector<int> sum;
	int rows = vec.size();
	int cols = vec[0].size();
	sum = vec[start];
	for (int i = start; i < rows; i++){
		for (int j = 0; j<cols; j++){
			sum[j] = sum[j] + vec[i][j];
		}
	}
	vector<double> mean;
	int denom = rows - start;
	for (int j = 0; j < cols; j++){
		mean.push_back((double)sum[j] / denom);
	}
	return mean;
};
vector<vector<double>> methods::MeanVector3DAcrossReps(vector<vector < vector<int>>> vvvOutputs){
	//VECTOR SHould not include the transient period
	//computs average mean value  across time anf reports it for each replications
	//out:  R * O(num outputs)
	int R = vvvOutputs.size();
	int T = vvvOutputs[0].size();
	int O = vvvOutputs[0][0].size();//#outcomes

	//sum across time, add all to the first elements
	vector<vector<int>> vvSums = vvvOutputs[0];
	for (int r = 1; r < R; r++){
		for (int t = 0; t < T; t++){
			for (int o = 0; o < O; o++){
				vvSums[t][o] = vvSums[t][o] + vvvOutputs[r][t][o];

			}
		}
	}
	vector<vector<double>> vvMean;
	vector<double> temp;
	for (int t = 0; t < T; t++){
		for (int o = 0; o < O; o++){
			temp.push_back((double)vvSums[t][o] / R);
		}
		vvMean.push_back(temp);
		temp.clear();
	}
	return vvMean;
};



void methods::printVector3D(vector<vector<vector<int>>> vvvDNAStrains){
	int R = vvvDNAStrains.size();
	int T = vvvDNAStrains[0].size();

	for (int r = 0; r < R; r++){
		cout << "i:   " << r << "---------------" << endl;
		for (int t = 0; t < T; t++){
			cout << "j=" << t << "--------------" << endl;
			for (int s = 0; s < vvvDNAStrains[r][t].size(); s++){
				cout << vvvDNAStrains[r][t][s] << " ";
			}
			cout << endl;
		}

	}


};

/*
double methods::runSimuCalib(int reps, double lambda,int goal){
	int REPL = reps;
	LAMBDA = lambda;

	//start parallelization------------
	int iCPU = omp_get_num_procs();//determne #cores
	omp_set_num_threads(iCPU);//set #of threads
	

	int tid;//thread's id
	double incMean;//mean incidence from a rep
	vector<double> vIncMeans;//shared vector of mean incidence; filled by all threads

	vector < vector < int>> vvInc;
	vector < vector < vector<int>>> vvvInc;
	gsl_rng * rng;//

	Community * com;
	int repsPerThread = TOTALREPS / iCPU;
	cout << "number of threads" << iCPU <<" each running: " <<repsPerThread<<" simulations"<< endl;
	int r;//counter
	vector<int> vtemp;
#pragma omp parallel shared(repsPerThread,vvvInc,vIncMeans ) private(tid,rng,com,r,incMean,vvInc,vtemp)
	{
		for (int j = 0; j < iCPU; j++){
#pragma omp single nowait 
			{

				tid = omp_get_thread_num();
				//printf("Hello World from thread = %d\n", tid);

				rng = gsl_rng_alloc(gsl_rng_taus2);
				gsl_rng_set(rng, time(0)*tid);
				com = new Community(tid);

				for (r = 0; r < repsPerThread; r++){
					printf("newRep from thread %d\n", tid);
					if (r>0)
						com->restart(rng);
					else
						com->createPopulation(rng, POPULATIONSIZE);

					com->initialInfection(rng, INITIALSTATES);

					//incMean = com->runModel_incMean(rng);
					//vIncMeans.push_back(d);
					vtemp = com->runModel_incVec(rng);
					vvInc.push_back(vtemp);
					incMean = methods::meanVec1DStartEnd(vtemp, TRANSIENT/12, MAXRUNTIME/12);
					vIncMeans.push_back(incMean);
				}
				delete com;
			//	printf("byebye from thread = %d\n", tid);
			}
		}
		vvvInc.push_back(vvInc);
	}

	// end of parallel section 
	double inc = methods::meanVec1D(vIncMeans);
	
	ostringstream filename;
	ostringstream firstline;

	filename << "calibIncidences-G-" << goal <<".txt";
	firstline << "Lambda:  " << lambda; 
	saveVector3DCumulative(filename.str(), false, vvvInc,firstline.str());



	return inc;



};

void methods::runSimulationParallel(bool isCollectingDNA,bool isCountAllStrains){

	  
	 int iCPU = omp_get_num_procs();//determne #cores
	 omp_set_num_threads(iCPU);//set #of threads


	 int tid;//thread's id
	 double incMean;//mean incidence from a rep
	 vector<double> vIncMeans;//shared vector of mean incidence; filled by all threads
	 
	 vector < vector<int>>  vvCumOutputs;//cumulative outputs across reps
	 vector<vector < vector<int>>> vvvOutputs;//3D contains DNAstrains at multiple points of time in multiple relications
	 vector<vector < vector<int>>> vvvDNAStrains;//3D contains DNAstrains at multiple points of time in multiple relications


	 gsl_rng * rng;//
	 Community * com;
	 int repsPerThread = REPLICATIONS / iCPU;
	 cout << "number of threads" << iCPU << " each running: " << repsPerThread << " simulations" << endl;
	 int r;//counter
	 vector<int> vtemp;


#pragma omp parallel shared(repsPerThread,vvvOutputs,vvvDNAStrains ) private(tid,rng,com,r)
	 {
		 for (int j = 0; j < iCPU; j++){
#pragma omp single nowait 
			 {

				 tid = omp_get_thread_num();
				 //printf("Hello World from thread = %d\n", tid);

				 rng = gsl_rng_alloc(gsl_rng_taus2);
				 gsl_rng_set(rng, time(0)*tid);
				 com = new Community(tid);

				 for (r = 0; r < repsPerThread; r++){
					 printf("newRep from thread %d\n", tid);
					 if (r>0)
						 com->restart(rng);
					 else
						 com->createPopulation(rng, POPULATIONSIZE);
					 //infect
					 com->initialInfection(rng, INITIALSTATES);
					 //run
					 com->runModelFull(rng, isCollectingDNA, isCountAllStrains);
					 //saveoutputs
					 vvvOutputs.push_back(com->returnOutputs());
					 vvvDNAStrains.push_back(com->returnDNAStrains());
					 
					 //record cumulativeoutputs
					 vvCumOutputs.push_back(com->returnCumOutputs());
					}
				 delete com;
				 //	printf("byebye from thread = %d\n", tid);
			 }
		 }
		
	 }

  ostringstream ext;
  ext << "-L-" << LAMBDA << "-M-" << MIGRATE;
  //sLambda << LAMBDA;
	 methods::SaveSimulationInfo("outSimulationSetupInfo" + ext.str() + ".txt", true, vtemp);
	 //single outputs
//	 methods::saveSingleOutputAcrossReps("outRepBase-Incidence-" + ext.str() + ".txt", true, vvvOutputs, 4);
	// methods::saveSingleOutputAcrossReps("outRepBase-FP-" + ext.str() + ".txt", true, vvvOutputs, 5);
	 //summary all outputs
	 methods::saveVector2D("outSummary" + ext.str() + ".txt", true, methods::MeanVector3DAcrossReps(vvvOutputs));

	 //DNAs
	 if (isCollectingDNA) methods::saveVector3DCumulative("DNASamples" + ext.str() + ".txt", true, vvvDNAStrains);

	 //save cumulativeoutputs
	 ostringstream est;
	 est << LAMBDA << " " << MIGRATE;
	 saveVector2D("cumOutputs.txt", false, vvCumOutputs, est.str());
	 


 };
 */

void methods::runSimulation(gsl_rng * rng, int tid, double lambda, double exRate, bool isCollectingDNA, bool isCountAllStrains, bool isAnalyzingDist){
	//create community with lambda, migrate/run simulation for rep times
	//*we don't use a 3D output vector anymore. instead just save the output means at the end of each run. 
	//*Can add the lines if wish to save mean outputs across time
	/*
	vector<vector < vector<int>>> vvvOutputs;//3D contains DNAstrains at multiple points of time in multiple relications
	//vector<double> vMeanofMeans;//vector of meanoutput values across all replications of this model

	Community * com = new Community(tid , lambda, exRate);
	
	for (int r = 0; r < REPLICATIONS; r++)
	{
		cout << "running rep " << r << " by " << tid << endl;
		if (r>0)
			com->restart(rng);
		else
			com->createPopulation(rng, POPULATIONSIZE);

		com->initialInfection(rng, INITIALSTATES);
		//com->runModel_MeanOutputs(rng, r, isCollectingDNA, isCountAllStrains);

		com->runModel_FullOutputs(rng, r, isCollectingDNA, isCountAllStrains, isAnalyzingDist);
		vvvOutputs.push_back(com->returnOutputs());
		
		//write DNA samples at the end of last file
		ostringstream ext;
		ext << "-L-" << lambda << "-ExR-" << exRate << "-SP-" << 0.001;
		if (isCollectingDNA) methods::saveVector2D("DNASamples" + ext.str() + ".txt", false, com->returnDNAStrains());
		//write mean outputs from this rep
		//vector<double> vTemp = com->returnMeanOutputs();
		//saveVector1D("meanOutputs.txt", false, vTemp);
		
		//record means to comput overall means
		//int I = vTemp.size();
		//if (r == 0) vMeanofMeans = vTemp;
		//else{
		//	for (int i = 0; i < I; i++){
		//		vMeanofMeans[i] += vTemp[i];
		//	}
		//}
	}

	//int I = vMeanofMeans.size();
	//for (int i = 0; i < I; i++){
	//	vMeanofMeans[i] = vMeanofMeans[i] / REPLICATIONS;
	//}
	//methods::saveVector1D("meanOfMeans.txt", false, vMeanofMeans);

	
	ostringstream ext;
	ext << "-L-" << lambda << "-ExR-" << exRate << "-SP-" << 0.001;
	methods::saveVector2D("outSummary" + ext.str() + ".txt", true, methods::MeanVector3DAcrossReps(vvvOutputs));
 
	//cout << "bybye " << tid << endl;
	delete com;
	*/
};

void methods::runSimulation(gsl_rng * rng, int tid, double lambda, double exRate,double spRate, bool isCollectingDNA, bool isCountAllStrains){
	//create community with lambda, migrate/run simulation for rep times
	//*we don't use a 3D output vector anymore. instead just save the output means at the end of each run. 
	//*Can add the lines if wish to save mean outputs across time

	vector<vector < vector<int>>> vvvOutputs;//3D contains DNAstrains at multiple points of time in multiple relications
	//vector<double> vMeanofMeans;//vector of meanoutput values across all replications of this model

	Community * com = new Community(tid, lambda, exRate,spRate);

	for (int r = 0; r < REPLICATIONS; r++)
	{
		cout << "running rep " << r << " by " << tid << endl;
		if (r>0)
			com->restart(rng);
		else
			com->createPopulation(rng, POPULATIONSIZE);

		com->initialInfection(rng, INITIALSTATES);
		
		//com->runModel_MeanOutputs(rng, r, isCollectingDNA, isCountAllStrains);
			
		com->runModel_FullOutputs(rng, r, isCollectingDNA, isCountAllStrains);
		vvvOutputs.push_back(com->returnOutputs());
		
		//write DNA samples at the end of last file
		ostringstream ext;
		double spr = floor(spRate * 12000000) / 1000;
		ext << "-L-" << lambda << "-ExR-" << exRate << "-SP-" << spr;
		if (isCollectingDNA) methods::saveVector2D("DNASamples" + ext.str() + ".txt", false, com->returnDNAStrains());
		
		//record means to comput overall means
		//int I = vTemp.size();
		//if (r == 0) vMeanofMeans = vTemp;
		//else{
		//	for (int i = 0; i < I; i++){
		//		vMeanofMeans[i] += vTemp[i];
		//	}
		//}
	}

	//int I = vMeanofMeans.size();
	//for (int i = 0; i < I; i++){
	//	vMeanofMeans[i] = vMeanofMeans[i] / REPLICATIONS;
	//}
	//methods::saveVector1D("meanOfMeans.txt", false, vMeanofMeans);

	

ostringstream ext;
	double spr = floor(spRate * 12000000) / 1000;
	ext << "-L-" << lambda << "-ExR-" << exRate << "-SP-" << spr;
	methods::saveVector2D("outSummary" + ext.str() + ".txt", true, methods::MeanVector3DAcrossReps(vvvOutputs));
	
	//cout << "bybye " << tid << endl;

	delete com;
};