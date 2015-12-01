// Person.cpp
#include <cstdlib>
#include <cstring>
#include <climits>
#include <iostream>
#include <string>

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Person.h"
#include "params.h"




using namespace std;


Person::Person(int id) {
	_bRecovered = false;
	_nID = id;
	//_nBirthTime=time;    
	_nAge = 0;
	_dInfectiousness = 0;  	_dImmunity = 0;
	_nState = 0;
	_nStrain = 0;

	_nELTBTime = 1000000;  _nLLTBTime = 1000000;	_nATBTime = 1000000;	_nRECTime = 1000000;
	_tOfSlowProg = 1000000;	_tOfFastProg = 1000000; 	_tOfRecovery = 1000000; _tOfFatality = 1000000; _tOfRelapse = 1000000;
	//printf("person %d born  \n", _nID);

};

Person::Person(int age, int id) {
	_bRecovered = false;
	_nID = id;
	//_nBirthTime=time;    
	_nAge = age;
	_dInfectiousness = 0;  	_dImmunity = 0;
	_nState = 0;
	_nStrain = 0;

	_nELTBTime = 1000000;  _nLLTBTime = 1000000;	_nATBTime = 1000000;	_nRECTime = 1000000;
	_tOfSlowProg = 1000000;	_tOfFastProg = 1000000; 	_tOfRecovery = 1000000; _tOfFatality = 1000000; _tOfRelapse = 1000000;
	//printf("person %d born  \n", _nID);
};

void Person::reborn(int id) {
	_bRecovered = false;
	_nID = id;
	//_nBirthTime=time;    
	_nAge = 0;
	_dInfectiousness = 0;  	_dImmunity = 0;
	_nState = 0;
	_nStrain = 0;

	_nELTBTime = 1000000;  _nLLTBTime = 1000000;	_nATBTime = 1000000;	_nRECTime = 1000000;
	_tOfSlowProg = 1000000;	_tOfFastProg = 1000000; 	_tOfRecovery = 1000000; _tOfFatality = 1000000; _tOfRelapse = 1000000;
	//cout<<_nID<<" constrcuted";
};
void Person::reborn(int age, int id) {
	_bRecovered = false;
	_nID = id;
	//_nBirthTime=time;    
	_nAge = age;
	_dInfectiousness = 0;  	_dImmunity = 0;
	_nState = 0;
	_nStrain = 0;

	_nELTBTime = 1000000;  _nLLTBTime = 1000000;	_nATBTime = 1000000;	_nRECTime = 1000000;
	_tOfSlowProg = 1000000;	_tOfFastProg = 1000000; 	_tOfRecovery = 1000000; _tOfFatality = 1000000; _tOfRelapse = 1000000;
	//printf("person %d reborn  \n", _nID);
};

Person::~Person() {
	//	cout<<_nID<<" deconstructed";
};


void Person::enterELTB(gsl_rng *rng, int time){
	_nState = 1;
	_dImmunity = IMMUNITYELTB;
	_dInfectiousness = 0;
	_nELTBTime = time;
	//_tOfFastProgression
	_tOfFastProg = 1000000;
	double r = gsl_rng_uniform(rng);
	if (r < FASTPROGRATECOEF[4] * CUMFASTPROG){//yearly rate
		for (int j = 0; j < 5; j++){
			if (r < FASTPROGRATECOEF[j] * CUMFASTPROG) _tOfFastProg = j * 12 + gsl_rng_uniform_int(rng, 12) + time;
			break;
		}
	}

	/*	_tOfFastProg=0;
	int Y=gsl_ran_geometric(rng,FASTPROGRATE);
	_tOfFastProg= Y+time;//(Y-1)*12+gsl_rng_uniform_int(rng,12)+ time;//time of progression
	*/
};
void Person::enterLLTB(gsl_rng *rng, int time,double spRate){
	_nState = 2;
	_dImmunity = IMMUNITYLLTB;
	_nLLTBTime = time;
	//_tOfSlowProg
	_tOfSlowProg = 1000000;
	int Y = gsl_ran_geometric(rng, spRate);//generates the number of years to slowProg
	_tOfSlowProg = Y + time;   // (Y-1)*12+gsl_rng_uniform_int(rng,12)+ time;//time of progression
};

void Person::enterATB(gsl_rng *rng, int time){
	_nState = 3;
	_dImmunity = IMMUNITYATB;
	_nATBTime = time;
	//timeto recovery	
	_tOfRecovery = 10000000;
	int Y = gsl_ran_geometric(rng, RECOVERYRATE);
	_tOfRecovery = Y + time;   //(Y-1)*12+gsl_rng_uniform_int(rng,12)+time;//timeof recovery
	//cout<<_tOfRecovery<<endl;
	//time to fatality
	_tOfFatality = 10000000;
	Y = gsl_ran_geometric(rng, FATALITY);
	_tOfFatality = Y + time;    //(Y-1)*12+gsl_rng_uniform_int(rng,12)+time;//timeof recovery
};
void Person::enterREC(gsl_rng *rng, int time){
	_bRecovered = true;
	_nState = 4;
	//_nStrain=0;
	_dImmunity = IMMUNITYREC;
	_nRECTime = time;
	//time to relapse
	_tOfRelapse = 10000000;
	if (RELAPSE > 0){
		int Y = gsl_ran_geometric(rng, RELAPSE);
		if (Y < 24){
			_tOfRelapse = Y + time;
		}
	}
};

void Person::updateInfectiousness(int t){

	int  d = (t)-(_nATBTime);
	if (d < MAXINFDURATION)
		_dInfectiousness = (MAXINFECTIOUSNESS* d / MAXINFDURATION);
	else
		_dInfectiousness = MAXINFECTIOUSNESS;
};
