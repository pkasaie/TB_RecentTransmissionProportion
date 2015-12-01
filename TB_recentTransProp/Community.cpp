// Community.cpp

#include <cstdlib>
#include <cstring>
#include <climits>
#include <iostream>
#include <string>
#include <iterator>
#include <fstream>
#include <sstream>
#include <algorithm> //std::sort
#include <math.h>
#include <list>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>   //for time function
#include "Person.h"
#include "Community.h"
#include "params.h"
#include "methods.h"

using namespace std;

int initialPopSize;


Community::Community(int tid, double lambda, double exRate) {
	_LAMBDA = lambda;
	_EXCHANGERATE = exRate;

	_nThreadId = tid;
	_nNextID = 1;
	_nMonth = 1;
	_nLastStrain = 0;


	numContacts = 0; //number of contacts
	numSusExposures = 0;
	numEnterELTB = 0;//number of ELTB infections
	numIncidence = 0;//number of ATB infections
	numFastProg = 0;
	numSlowProg = 0;
	numATBfromREC = 0;//including those rec people who reinfect+relapse
	numRelapses = 0;
	numImmATB = 0;
	numMigrated = 0;
	numFatality = 0;
	numMortality = 0;
	numOldageMortality = 0;
	numNewborn = 0;
	numExchanged = 0;

	cout << "community constructed by " << tid << " L-" << lambda << "-ExR-" << exRate << endl;
};
Community::Community(int tid, double lambda, double exRate, double spRate) {
	_LAMBDA = lambda;
	_EXCHANGERATE = exRate;
	_SPRATE = spRate;

	_nThreadId = tid;
	_nNextID = 1;
	_nMonth = 1;
	_nLastStrain = 0;


	numContacts = 0; //number of contacts
	numSusExposures = 0;
	numEnterELTB = 0;//number of ELTB infections
	numIncidence = 0;//number of ATB infections
	numFastProg = 0;
	numSlowProg = 0;
	numATBfromREC = 0;//including those rec people who reinfect+relapse
	numRelapses = 0;
	numImmATB = 0;
	numMigrated = 0;
	numFatality = 0;
	numMortality = 0;
	numOldageMortality = 0;
	numNewborn = 0;
	numExchanged = 0;

	cout << "community constructed by " << tid << " L-" << lambda << "-ExR-" << exRate << endl;
};

Community::~Community(){
	int N = _vPopulation.size();
	for (int i = 0; i < N; i++){//killing and destroying people
		Person * p = _vPopulation.back();
		delete p;
		_vPopulation.pop_back();
	}
};

void Community::restart(gsl_rng *rng) {
	//_labmda stay the same
	//_migrate stay the same
	//__nThreadId syas the same

	_nNextID = 1;
	_nMonth = 1;
	_nLastStrain = 0;



	numContacts = 0; //number of contacts
	numSusExposures = 0;
	numEnterELTB = 0;//number of ELTB infections
	numIncidence = 0;//number of ATB infections
	numFastProg = 0;
	numSlowProg = 0;
	numATBfromREC = 0;//including those rec people who reinfect+relapse
	numRelapses = 0;
	numImmATB = 0;
	numMigrated = 0;
	numFatality = 0;
	numMortality = 0;
	numOldageMortality = 0;
	numNewborn = 0;
	numExchanged = 0;

	//--------------population------------------------
	int N = _vPopulation.size();
	if (N < POPULATIONSIZE){
		for (int i = 0; i < POPULATIONSIZE - N; i++){
			_vPopulation.push_back(new Person(i));

		}
	}
	else {
		for (int i = 0; i < N - POPULATIONSIZE; i++){
			delete _vPopulation.back();
			_vPopulation.pop_back();
		}
	}
	//restart population

	int age = 0;
	for (int i = 0; i < POPULATIONSIZE; i++){
		age = gsl_rng_uniform_int(rng, MAXAGE);
		_vPopulation[i]->reborn(age, _nNextID);
		_nNextID++;
	}
	//------vectors-------------------------
	
	vvOut.clear();
	vvAgeDists.clear();
	vTotalcirculatingStrains.clear();
	vPeopleDisDuration.clear();
	vDNAStrains.clear();
	vvDNAStrains.clear();
	//printf("population reborn \n");

	//cout << "community restart by " << _threadid<< " L-" << _LAMBDA << " M-" << _MIGRATE << endl;
};

void Community::createPopulation(gsl_rng *rng, int N){
	for (int i = 0; i < N; i++){
		int age = gsl_rng_uniform_int(rng, MAXAGE);
		_vPopulation.push_back(new Person(age, _nNextID));
		_nNextID++;

	}
	initialPopSize = N;
	//	printf("population created \n");
};
void Community::initialInfection(gsl_rng *rng, double initialStates[]){
	int N = _vPopulation.size();
	double pEL = initialStates[0];
	double pLL = initialStates[1];
	double pATB = initialStates[2];
	double pREC = initialStates[3];
	double pSUS = 1 - pEL - pLL - pATB - pREC;
	for (int i = 0; i < N; i++){
		Person *p = _vPopulation[i];
		double r = gsl_rng_uniform(rng);
		if (r > pSUS){//is infected?
			if (r < pSUS + pLL){//---------------LLTB
				p->enterLLTB(rng, _nMonth, _SPRATE);
				_nLastStrain++;	p->setStrain(_nLastStrain);
			}
			else {
				if (r < pSUS + pLL + pEL){//---------------ELTB
					int t = _nMonth - gsl_rng_uniform_int(rng, 60);//randomly assign the previous time of entering ELTB
					p->enterELTB(rng, t);
					_nLastStrain++;	p->setStrain(_nLastStrain);
					numEnterELTB++;
				}
				else{//---------------REC
					if (r < pSUS + pLL + pEL + pREC){
						p->enterREC(rng, _nMonth);
						_nLastStrain++;	p->setStrain(_nLastStrain);
					}
					else{//----------------------ATB
						int t = _nMonth - gsl_rng_uniform_int(rng, 10);//randomly assign the previous time of entering ATB
						p->enterATB(rng, t);
						_nLastStrain++;	p->setStrain(_nLastStrain);
						numIncidence++;
					}

				}
			}

		}
	}
	//cout<<"initial infection successful."<<endl;
};

//PopulatonDynamics------------------------------------------------------------------------------
void Community::contactPoisson(gsl_rng *rng){
	int C = 0; //num contacts
	double r = 0;
	double ptran = 0;
	int N = _vPopulation.size();
	vector<Person*> vATB;
	for (int i = 0; i < N; i++){
		Person *p = _vPopulation[i];
		if (p->getState() == 3) vATB.push_back(p);
	}
	int A = vATB.size();
	for (int i = 0; i < A; i++){
		Person *p = vATB[i];
		//generate number of contacts
		C = gsl_ran_poisson(rng, _LAMBDA); //vec.push_back(C);
		double pInf = p->getInfectiousness();
		int pStrain = p->getStrain();
		for (int c = 0; c < C; c++){

			Person *q = _vPopulation[gsl_rng_uniform_int(rng, N)];
			ptran = pInf*(1 - q->getImmunity());
			numContacts++;

			if (gsl_rng_uniform(rng) < ptran){
				if (q->getState() == 0) numSusExposures++;
				q->setStrain(pStrain);
				q->enterELTB(rng, _nMonth);
				numEnterELTB++;

			}
		}
	}

};
void Community::progress(gsl_rng *rng, bool isCollectDNA){//disease progress at the end of month
	//check the timing of scheduled events and execute them if it's the time
	//!!! _nMonth >= eventTime!!!
	list<int> temp;//list of dead people
	int N = _vPopulation.size();
	for (int i = 0; i < N; i++){
		Person *p = _vPopulation[i];
		int pState = p->getState();
		switch (pState){
		case(0) : break;//susceptible
		case (1) : {//ELTB-check: fastProg or proceeding to LLTB
					   if (_nMonth >= p->getTimeOfFastProg()){
						   numFastProg++;
						   numIncidence++;
						   if (p->_bRecovered) numATBfromREC++;
						   p->enterATB(rng, _nMonth);
					   }
					   else if ((_nMonth - p->getELTBTime()) >= 60){
						   p->enterLLTB(rng, _nMonth, _SPRATE);
					   }
					   break;
		};
		case (2) : {//LLTB-check: slow prog
					   if (_nMonth >= p->getTimeOfSlowProg()){
						   numSlowProg++;
						   numIncidence++;
						   if (p->_bRecovered) numATBfromREC++;
						   p->enterATB(rng, _nMonth);
					   }
					   break;
		};
		case(3) : {	//ATB-check: fatality or recovery
					  if (_nMonth >= p->getTimeOfFatality()) {//case-fatlity
						  numFatality++;
						  vPeopleDisDuration.push_back(_nMonth - p->getATBTime());
						  temp.push_back(i);
					  }
					  else if (_nMonth >= p->getTimeOfRecovery()){//recovery
						  //DNA Sampling--------if we're collecting samples
						  if (isCollectDNA)	vDNAStrains.push_back(p->getStrain());
						  p->enterREC(rng, _nMonth);
						  vPeopleDisDuration.push_back(_nMonth - p->getATBTime());
					  }
					  else
						  p->updateInfectiousness(_nMonth);

					  break;
		};
		case (4) : {//REC check: relapse
					   if (_nMonth >= p->getTimeOfRelapse()) {//relapse
						   numRelapses++;
						   numIncidence++;
						   numATBfromREC++;
						   p->enterATB(rng, _nMonth);
					   }
					   break;
		};
		default: cout << "meh";
			break;
		};
	};

	temp.sort();//start removing people from the end of population to save the order
	while (!temp.empty()){
		int id = temp.back();
		_vPopulation.erase(_vPopulation.begin() + id);
		temp.pop_back();
	}
};
void  Community::deathBirth(gsl_rng *rng){
	//list<int> tempID;//list of dead people
	vector<int> vTempID;
	int N = _vPopulation.size();
	double r;
	for (int i = 0; i < N; i++){
		Person *p = _vPopulation[i];
		int myState = p->getState();
		//without age:
		/*
		r = gsl_rng_uniform(rng);
		if (r < NATURALMORTALITY){
		numMortality++;
		if (p->getState() == 3) 	{ vPeopleDisDuration.push_back(_nMonth - p->getATBTime()); }
		temp.push_back(i);
		}
		else p->addAge();
		}
		*/
		//with age
		if ((p->getAge() >= MAXAGE)){
			numOldageMortality++;
			if (myState == 3) 	{ vPeopleDisDuration.push_back(_nMonth - p->getATBTime()); }
			//tempID.push_back(i);
			vTempID.push_back(i);
		}
		else {
			r = gsl_rng_uniform(rng);
			double d;
			if (p->getAge() == 0) d = AGESPECIFICMORTALITY[0];
			else d = (double)AGESPECIFICMORTALITY[p->getAge() / 5 + 1] / 5;
			if (r < d){   // NATURALMORTALITY){
				numMortality++;
				if (myState == 3) 	{ vPeopleDisDuration.push_back(_nMonth - p->getATBTime()); }
				//tempID.push_back(i);
				vTempID.push_back(i);
			}
			else p->addAge();
		}
	}

	/*
	//initial algorithm: erase dead people:add newborns
	tempID.sort();//start removing people from the end of population to save the order
	while (!tempID.empty()){
	int id = tempID.back();
	delete _vPopulation[id];
	_vPopulation.erase(_vPopulation.begin() + id);
	tempID.pop_back();
	}

	//adding new-borns:
	//int newborns = gsl_ran_poisson(rng, BIRTHMEAN);

	int newborns = gsl_ran_poisson(rng, initialPopSize - _vPopulation.size());
	numNewborn = numNewborn + newborns;
	for (int i = 0; i < newborns; i++){
	_vPopulation.push_back(new Person());//with age 0
	}*/
	//problem: erase is a bottlneck:we shold minimize it
	//lets decide about number of newborns yet; replace ded people with them; and only erase the rest of dedpeople if there is any remaining

	int numDead = vTempID.size();
	int numBirth = gsl_ran_poisson(rng, initialPopSize - _vPopulation.size() + numDead);

	if (numBirth < numDead){
		for (int i = 0; i < numBirth; i++){//reborn dead person
			_vPopulation[vTempID[i]]->reborn(_nNextID);
			_nNextID++;
		}
		for (int i = numBirth; i < numDead; i++){//erase the rest of them
			int id = vTempID.back();//we start erasing from the end to front
			vTempID.pop_back();
			delete _vPopulation[id];
			_vPopulation.erase(_vPopulation.begin() + id);
		}
	}
	else{
		for (int i = 0; i < numDead; i++){//reborn all dead person
			_vPopulation[vTempID[i]]->reborn(_nNextID);
			_nNextID++;
		}
		for (int i = numDead; i < numBirth; i++){//add new people
			_vPopulation.push_back(new Person(_nNextID));//with age 0
			_nNextID++;
			numNewborn++;
		}

	}




};
void Community::migrationSwaping(gsl_rng *rng){
	double r = 0;
	//MIGRATION OUT
	int N = _vPopulation.size();
	for (int i = 0; i < N; i++){
		r = gsl_rng_uniform(rng);
		if (r < MIGRATE){
			numMigrated++;
			//select random person to leave
			int id1 = gsl_rng_uniform_int(rng, N);
			Person *p = _vPopulation[id1];
			//select a random person to replace him
			int id2 = gsl_rng_uniform_int(rng, N);
			Person *q = _vPopulation[id2];
			//swap
			*p = *q; //keeping the original copy of Person p, and coping all attributes from person q
			p->setId(_nNextID);
			_nNextID++;
			//p->setBithTime(_nMonth);
			int pState = p->getState();
			if ((pState > 0)){// && (p->getState() < 4)){
				_nLastStrain++;  p->setStrain(_nLastStrain);
			}

			//p is leaving population
			//for this new person who's just entered: we don't restart the disease state to keep the timings in order
			//but we need to update the stats as if some new sick person entered our population
			if (pState == 3){
				numImmATB++;
			}

			//cout<<_vPopulation[id1]->getID()<<endl;
			//cout<<_vPopulation[id2]->getID()<<endl;
		}
	}
};
void Community::exchange(gsl_rng *rng){
	//exchange LLTB people for new strains
	int N = _vPopulation.size();
	for (int i = 0; i < N; i++){
		Person *p = _vPopulation[i];
		if (p->getState() == 2) {//Only LLTB
			if (gsl_rng_uniform(rng) < _EXCHANGERATE){
				numExchanged++;
				_nLastStrain++;  p->setStrain(_nLastStrain);

			}
		}
	}
	//cout << numExchanged << endl;
};

//outputs----------------------------------------------------------------------------------------------------
vector<int> Community::returnAgeDist(){
	vector<int> ageVec;
	for (int i = 0; i < 100; i++){ ageVec.push_back(0); }

	for (int i = 0; i < _vPopulation.size(); i++){
		ageVec[_vPopulation[i]->getAge()]++;
	}
	return ageVec;


}
vector<int> Community::returnDiseaseVector(){	//rearrenge disease vectors
	vector<int> vDS;
	for (int i = 0; i < 5; i++) vDS.push_back(0);
	for (int i = 0; i < _vPopulation.size(); i++){
		int s = _vPopulation[i]->getState();
		vDS[s]++;
	}
	return vDS;
};
void Community::saveRunOutputs(int r){
	ofstream outFile;
	stringstream filename;
	filename << "Rep-" << r << ".txt";
	outFile.open(filename.str());
	//outFile<<"time "<<" "<<"Incidence "<<"FastProg "<<"SlowProg "<<"natMort "<<"Fatality "<<endl;}
	for (int i = 0; i < vvOut.size(); i++){
		for (int j = 0; j < vvOut[0].size(); j++){
			outFile << vvOut[i][j] << " ";
		}
		outFile << endl;
	}
	outFile.close();
};
void Community::saveDieaseDuration(bool b){
	ofstream outFile;
	if (b) {//newFile
		outFile.open("outDD.txt");
	}
	//outFile<<"time "<<"Incidence "<<"FastProg "<<"SlowProg "<<"natMort "<<"Fatality "<<endl;}
	else
		outFile.open("outDD.txt", ios_base::app);
	for (int i = 0; i < vPeopleDisDuration.size(); i++){
		outFile << vPeopleDisDuration[i] << endl;
	}
	outFile.close();
};
void Community::printDiseaseVector(){
	vector<int> vDS;
	vDS = returnDiseaseVector();
	for (int i = 0; i < vDS.size(); i++){
		cout << vDS[i] << " ";
	}
	cout << endl;
};
double Community::returnVectorAverage(vector<int> vec){
	int sum = 0;
	if (vec.size()>0){
		for (int i = 0; i < vec.size(); i++){
			sum = sum + vec[i];
		}
		return (double)sum / vec.size();
	}
	else return 0;
}

int  Community::countTotalCirculatingStrains(){//counts # of circulating strains in the whole population EL LL ATB
	vector<int> vStrains;
	//list<int> strains;
	int strain = 0; int j = 0;
	int I = _vPopulation.size();
	for (int i = 0; i < I; i++){
		//if ((_vPopulation[i]->getState() >0)){///only for E,LL,ATB infected people not REC
		vStrains.push_back(_vPopulation[i]->getStrain());
		//}
	}

	//// sort(vStrains.begin(), vStrains.end());     //(12 26 32 33 45 53 71 80)
	// unique(vStrains.begin(), vStrains.end());

	std::sort(vStrains.begin(), vStrains.end());
	vStrains.erase(std::unique(vStrains.begin(), vStrains.end()), vStrains.end());

	return vStrains.size();

};
//-------------------------------------------------------------------------------
void Community::tick(gsl_rng *rng, bool isCollectDNA){


	contactPoisson(rng);//contactHomo(rng);
	progress(rng, isCollectDNA);

	if (_nMonth % 12 == 0) {
		//	cout << numIncidence << endl;
		deathBirth(rng);
		migrationSwaping(rng);


		//save DNA strains
		if (isCollectDNA){
			vDNAStrains.insert(vDNAStrains.begin(), _nMonth / 12);
			vvDNAStrains.push_back(vDNAStrains);
			vDNAStrains.clear();
		};


		//updating outputs
		vector<int> temp; vvOut.push_back(temp);
		vvOut.back().push_back(_nMonth);//year//0
		vvOut.back().push_back(numContacts);	numContacts = 0;//1
		vvOut.back().push_back(numSusExposures);	numSusExposures = 0;//2
		vvOut.back().push_back(numEnterELTB);		numEnterELTB = 0;//3
		vvOut.back().push_back(numIncidence);	numIncidence = 0;//4
		vvOut.back().push_back(numFastProg);	numFastProg = 0;//5
		vvOut.back().push_back(numSlowProg);	numSlowProg = 0;//6
		vvOut.back().push_back(numRelapses);	numRelapses = 0;
		vvOut.back().push_back(numImmATB);	numImmATB = 0;
		vvOut.back().push_back(numMigrated); numMigrated = 0;
		vvOut.back().push_back(numFatality);	numFatality = 0;
		vvOut.back().push_back(numMortality);	numMortality = 0;
		vvOut.back().push_back(numOldageMortality);	numOldageMortality = 0;
		vvOut.back().push_back(numNewborn);	numNewborn = 0;

		temp = returnDiseaseVector();
		for (int i = 0; i < 5; i++){
			vvOut.back().push_back(temp[i]);
		}
		vvOut.back().push_back(_vPopulation.size());
		//disease duration
		vvOut.back().push_back(returnVectorAverage(vPeopleDisDuration));
		vPeopleDisDuration.clear();
		//
		vvOut.back().push_back(numATBfromREC);	numATBfromREC = 0;//those who r sick and have been previously treated


		//if ((_nMonth % 300) == 0) vvAgeDists.push_back(returnAgeDist());

		//DNA srains circulating
		if ((_nMonth % 300) == 0)			vvOut.back().push_back(countTotalCirculatingStrains()); 		else vvOut.back().push_back(0);


	}

	_nMonth++;
};
//------------------------------------------------------------------------------------------------------

vector<int> Community::runModel_incVec(gsl_rng *rng){//CALIBRATION: runs the model and return the vector of incidences
	vector<int> vInc;
	for (int t = 0; t < MAXRUNTIME; t++){
		contactPoisson(rng);//contactHomo(rng);
		progress(rng, false);
		if (_nMonth % 12 == 0) {
			//	cout << numIncidence << endl;
			deathBirth(rng);
			migrationSwaping(rng);

			vInc.push_back(numIncidence);	numIncidence = 0;//4
		};



		_nMonth++;
	}
	return vInc;

};
double Community::runModel_incMean(gsl_rng *rng){//runs the model and resurns the mean incidence after the transient period
	vector<int> vInc;
	for (int t = 0; t < MAXRUNTIME; t++){
		contactPoisson(rng);//contactHomo(rng);
		progress(rng, false);
		if (_nMonth % 12 == 0) {
			//	cout << numIncidence << endl;
			deathBirth(rng);
			migrationSwaping(rng);


			vInc.push_back(numIncidence);	numIncidence = 0;//4
		};



		_nMonth++;
	}
	int sum = 0;

	for (int i = TRANSIENT / 12; i < MAXRUNTIME / 12; i++){
		sum = sum + vInc[i];
	}

	return (double)sum / (MAXRUNTIME - TRANSIENT);

};
void Community::runModel_FullOutputsDist(gsl_rng *rng, int curRep, bool isCollectDNA, bool isCountAllStrains, bool isAnalyzingDist){
	//runs the model, records all outputs and means of outputs
	int numOutputs = 24;
	vector<double> vMeans;
	for (int i = 0; i < numOutputs; i++) vMeans.push_back(0);//initialize

	for (int t = 0; t < MAXRUNTIME; t++){
		contactPoisson(rng);//contactHomo(rng);
		progress(rng, isCollectDNA);

		if (_nMonth % 12 == 0) {
			deathBirth(rng);
			//migrationSwaping(rng);
			exchange(rng);

			if (isCollectDNA){//save DNA strains
				vDNAStrains.insert(vDNAStrains.begin(), _nMonth / 12);//add the year to the begining of line-Canbe SKIPPED
				vvDNAStrains.push_back(vDNAStrains);//save DNA strain types
				vDNAStrains.clear();
			};

			//updating outputs
			vector<int> temp; vvOut.push_back(temp);
			vvOut.back().push_back(_nMonth);//year//0
			vvOut.back().push_back(numContacts);	numContacts = 0;//1
			vvOut.back().push_back(numSusExposures);	numSusExposures = 0;//2
			vvOut.back().push_back(numEnterELTB);		numEnterELTB = 0;//3
			vvOut.back().push_back(numIncidence);	numIncidence = 0;//4//4
			vvOut.back().push_back(numFastProg);	numFastProg = 0;//5//5
			vvOut.back().push_back(numSlowProg); numSlowProg = 0;//6	//6
			vvOut.back().push_back(numRelapses);	numRelapses = 0;//7
			vvOut.back().push_back(numImmATB);	numImmATB = 0;//8
			vvOut.back().push_back(numMigrated); numMigrated = 0;//9
			vvOut.back().push_back(numFatality);	numFatality = 0;//10
			vvOut.back().push_back(numMortality);	numMortality = 0;//11
			vvOut.back().push_back(numOldageMortality);	numOldageMortality = 0;//12
			vvOut.back().push_back(numNewborn);	numNewborn = 0; //13

			temp = returnDiseaseVector();
			for (int i = 0; i < 5; i++){//14 -18
				vvOut.back().push_back(temp[i]);
			}
			vvOut.back().push_back(_vPopulation.size());//19
			vvOut.back().push_back(returnVectorAverage(vPeopleDisDuration));		vPeopleDisDuration.clear();//20
			vvOut.back().push_back(numATBfromREC);	numATBfromREC = 0;//21//those who r sick and have been previously treated

			//if ((_nMonth % 300) == 0) vvAgeDists.push_back(returnAgeDist());
			//DNA srains circulating

			if (isCountAllStrains){
				int nCirculatingStrains = 0;
				if ((_nMonth % 600) == 0)
					nCirculatingStrains = countTotalCirculatingStrains();
				vvOut.back().push_back(nCirculatingStrains);//22
				//	cout << nCirculatingStrains << endl;
			}

			vvOut.back().push_back(numExchanged); numExchanged = 0;//23

			//add this to last 
			if (_nMonth > TRANSIENT){
				for (int i = 0; i < numOutputs; i++) vMeans[i] += vvOut.back()[i];//initialize
			}

			//end of year
		}
		if (isAnalyzingDist){
			if (((_nMonth % 240) == 0) || (_nMonth == 1)){ clusterDistAnalysis(); }
		}

		_nMonth++;
	}

	//compute means:
	int var = (MAXRUNTIME - TRANSIENT) / 12;
	//1-for all outputs
	for (int i = 0; i < numOutputs - 1; i++) vMeans[i] /= var;//initialize
	//2-for TB strains that are collected every 600 steps
	var = (MAXRUNTIME - TRANSIENT) / 600;	vMeans[numOutputs - 1] /= var;

	vMeans.push_back(_LAMBDA);
	vMeans.push_back(_EXCHANGERATE);
	vMeans.push_back(curRep);
	vMeans.push_back(gsl_rng_get(rng));
	methods::saveVector1D("meanOutputs.txt", false, vMeans);
};
void Community::clusterDistAnalysis(){
	vector<int> vSample;

	//collect a sample of all ATB people
	int I = _vPopulation.size();
	for (int i = 0; i < I; i++){
		if (_vPopulation[i]->getState() == 3){
			vSample.push_back(_vPopulation[i]->getStrain());
		}
	}

	//sort:
	std::sort(vSample.begin(), vSample.end());

	//grouping the clusters:-------------------
	vector<int> vClustSizes;//size of clusters 1 if unique, >1 if clusters
	int temp = 111111111111111;
	I = vSample.size();
	for (int i = 0; i < I; i++){
		if (vSample[i] == temp){//same cluster-not valid for fir element
			vClustSizes.back()++;
		}
		else{//new cluster
			temp = vSample[i];
			vClustSizes.push_back(1);
		}
	}

	std::sort(vClustSizes.begin(), vClustSizes.end());
	vector<int> vClustSizeFreq;//size of clusters 1 if unique, >1 if clusters
	int S = vClustSizes.back();//max cluster size
	for (int i = 0; i < S; i++){ vClustSizeFreq.push_back(0); }

	I = vClustSizes.size();
	for (int i = 0; i < I; i++){
		vClustSizeFreq[vClustSizes[i] - 1]++;
	}
	//save different file for difffernt scenarios	
	stringstream filename;
	filename << "clusterDist-" << _LAMBDA << "-" << _EXCHANGERATE << ".txt";
	methods::saveVector1D(filename.str(), false, vClustSizeFreq);

};


void Community::runModel_FullOutputs(gsl_rng *rng, int curRep, bool isCollectDNA, bool isCountAllStrains){
	//runs the model, records all outputs and means of outputs
	
	//cout << "FULL apprach" << MAXRUNTIME << " " << TRANSIENT << endl; 
	
	vvOut.clear();
	int numOutputs = 24;
	vector<double> vMeanOutputs;
	for (int i = 0; i < numOutputs; i++) vMeanOutputs.push_back(0);//initialize
	int myvar1 = 0; int myvar2 = 0;
	for (int t = 0; t < MAXRUNTIME; t++)
	{
		contactPoisson(rng);//contactHomo(rng);
		progress(rng, isCollectDNA);
		if (_nMonth % 12 == 0)
		{
			deathBirth(rng);
			exchange(rng);
			if (isCollectDNA){//save DNA strains
				vDNAStrains.insert(vDNAStrains.begin(), _nMonth / 12);//add the year to the begining of line-Canbe SKIPPED
				vvDNAStrains.push_back(vDNAStrains);//save DNA strain types
				vDNAStrains.clear();
			};
			
				vector<int> vTemp;
				vTemp.push_back(_nMonth);//year//0
				vTemp.push_back(numContacts);	numContacts = 0;//1
				vTemp.push_back(numSusExposures);	numSusExposures = 0;//2
				vTemp.push_back(numEnterELTB);		numEnterELTB = 0;//3
				vTemp.push_back(numIncidence);	numIncidence = 0;//4//4
				vTemp.push_back(numFastProg);	numFastProg = 0;//5//5
				vTemp.push_back(numSlowProg); numSlowProg = 0;//6	//6
				vTemp.push_back(numRelapses);	numRelapses = 0;//7
				vTemp.push_back(numImmATB);	numImmATB = 0;//8
				vTemp.push_back(numMigrated); numMigrated = 0;//9
				vTemp.push_back(numFatality);	numFatality = 0;//10
				vTemp.push_back(numMortality);	numMortality = 0;//11
				vTemp.push_back(numOldageMortality);	numOldageMortality = 0;//12
				vTemp.push_back(numNewborn);	numNewborn = 0; //13
				vector<int> vDV = returnDiseaseVector();
				for (int i = 0; i < 5; i++){//14 -18
					vTemp.push_back(vDV[i]);
				}
				vTemp.push_back(_vPopulation.size());//19
				vTemp.push_back(returnVectorAverage(vPeopleDisDuration));		vPeopleDisDuration.clear();//20
				vTemp.push_back(numATBfromREC);	numATBfromREC = 0;//21//those who r sick and have been previously treated
				//DNA srains circulating
				vTemp.push_back(numExchanged); numExchanged = 0;//22
				int nCirculatingStrains = 0;
				if (isCountAllStrains){
					if (_nMonth % 600==0){
						cout << _nMonth<<endl;
						myvar1++;
						nCirculatingStrains = countTotalCirculatingStrains();
					}
				}
				vTemp.push_back(nCirculatingStrains);//23
				vvOut.push_back(vTemp);
				//add this to last
				
				int I = vTemp.size();
				if (_nMonth > TRANSIENT)
				{
					myvar2++;
				for (int i = 0; i < I; i++) vMeanOutputs[i] += vTemp[i];
			//	cout << _nMonth << " inc " << vMeanOutputs[4] << " " << vTemp[4] << endl;
				}
			
		}
		_nMonth++;
	}
	cout << "myvar " << myvar1<<" "<<myvar2;
	//compute means:
	int var = (MAXRUNTIME - TRANSIENT) / 12;
	//cout << _nMonth++<<" "<< "first line:  " << vMeanOutputs[4] << "var " << var << endl;
	//1-for all outputs
	for (int i = 0; i < numOutputs - 1; i++) vMeanOutputs[i] /= var;//initialize
	//2-for TB strains that are collected every 600 steps
	var = (MAXRUNTIME - TRANSIENT) / 600;	vMeanOutputs[numOutputs - 1] /= var;

	vMeanOutputs.push_back(_LAMBDA);
	vMeanOutputs.push_back(_EXCHANGERATE);
	double spr = floor(_SPRATE * 12000000) / 1000;	vMeanOutputs.push_back(spr);

	vMeanOutputs.push_back(curRep);
	vMeanOutputs.push_back(gsl_rng_get(rng));
	methods::saveVector1D("meanOutputs.txt", false, vMeanOutputs);

};

void Community::runModel_MeanOutputs(gsl_rng *rng, int curRep, bool isCollectDNA, bool isCountAllStrains){
	//runs the model, only records output means during the steady state
	//cout << "MEAN apprach"<<MAXRUNTIME<<" "<<TRANSIENT<<endl;

	int numOutputs = 24;
	vector<double> vMeanOutputs;
	for (int i = 0; i < numOutputs; i++) vMeanOutputs.push_back(0);//initialize

	for (int t = 0; t < MAXRUNTIME; t++)
	{
		contactPoisson(rng);//contactHomo(rng);
		progress(rng, isCollectDNA);
		if (_nMonth % 12 == 0)
		{
			deathBirth(rng);
			exchange(rng);
			if (isCollectDNA){//save DNA strains
				vDNAStrains.insert(vDNAStrains.begin(), _nMonth / 12);//add the year to the begining of line-Canbe SKIPPED
				vvDNAStrains.push_back(vDNAStrains);//save DNA strain types
				vDNAStrains.clear();
			};
			
				vector<int> vTemp;
				vTemp.push_back(_nMonth);//year//0
				vTemp.push_back(numContacts);	numContacts = 0;//1
				vTemp.push_back(numSusExposures);	numSusExposures = 0;//2
				vTemp.push_back(numEnterELTB);		numEnterELTB = 0;//3
				vTemp.push_back(numIncidence);	numIncidence = 0;//4//4
				vTemp.push_back(numFastProg);	numFastProg = 0;//5//5
				vTemp.push_back(numSlowProg); numSlowProg = 0;//6	//6
				vTemp.push_back(numRelapses);	numRelapses = 0;//7
				vTemp.push_back(numImmATB);	numImmATB = 0;//8
				vTemp.push_back(numMigrated); numMigrated = 0;//9
				vTemp.push_back(numFatality);	numFatality = 0;//10
				vTemp.push_back(numMortality);	numMortality = 0;//11
				vTemp.push_back(numOldageMortality);	numOldageMortality = 0;//12
				vTemp.push_back(numNewborn);	numNewborn = 0; //13
				vector<int> vDV = returnDiseaseVector();
				for (int i = 0; i < 5; i++){//14 -18
					vTemp.push_back(vDV[i]);
				}
				vTemp.push_back(_vPopulation.size());//19
				vTemp.push_back(returnVectorAverage(vPeopleDisDuration));		vPeopleDisDuration.clear();//20
				vTemp.push_back(numATBfromREC);	numATBfromREC = 0;//21//those who r sick and have been previously treated
				//DNA srains circulating
				vTemp.push_back(numExchanged); numExchanged = 0;//22
				int nCirculatingStrains = 0;
				if (isCountAllStrains){
					if (_nMonth % 600){
						nCirculatingStrains = countTotalCirculatingStrains();
					}
				}
				vTemp.push_back(nCirculatingStrains);//23
				
				//add this to last
				if (_nMonth > TRANSIENT){
				int I = vTemp.size();
				for (int i = 0; i < I; i++) vMeanOutputs[i] += vTemp[i];
				//cout << _nMonth << " " << vMeanOutputs[4] << endl;
			}
		}
		_nMonth++;
	}

	//compute means:
	int var = (MAXRUNTIME - TRANSIENT) / 12;
	//cout << "first line:  " << vMeanOutputs[4] << "var " << var << endl;
	//1-for all outputs
	for (int i = 0; i < numOutputs - 1; i++) vMeanOutputs[i] /= var;//initialize
	//2-for TB strains that are collected every 600 steps
	var = (MAXRUNTIME - TRANSIENT) / 600;	vMeanOutputs[numOutputs - 1] /= var;

	vMeanOutputs.push_back(_LAMBDA);
	vMeanOutputs.push_back(_EXCHANGERATE);
	double spr = floor(_SPRATE * 12000000) / 1000;	vMeanOutputs.push_back(spr);

	vMeanOutputs.push_back(curRep);
	vMeanOutputs.push_back(gsl_rng_get(rng));
	methods::saveVector1D("meanOutputs.txt", false, vMeanOutputs);

};

