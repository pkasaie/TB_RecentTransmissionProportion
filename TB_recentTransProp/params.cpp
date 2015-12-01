 #include <string>
#include <vector>
#include "params.h"

//double LAMBDA ;// = 1.45;// [] = { 1.45, 2, 3, 4, 1, .5 };
//double SLOWPROGRATE;// = .001 / 12;// annual risk of reactivation per ferebee(1970)


double MIGRATE;// = 0.1;// 0.10;

int REPLICATIONS;// 0;//replications
int MAXRUNTIME ;// = 1200;// 00;				//simulation period
int TRANSIENT;
//int TRANSIENTPERIOD = 100;//years
//int DNASamplingStartMonth = TRANSIENTPERIOD*12+1;//month

int POPULATIONSIZE;// = 100000;//population size 
double INITIALSTATES[] = { 0.045, 0.46, 0.002, 0.06 };




//inputs:



//mortality for age groups{AGELT1	AGE1-4	AGE5-9	AGE10-14	AGE15-19	AGE20-24	AGE25-29	AGE30-34	AGE35-39	AGE40-44	AGE45-49	AGE50-54	AGE55-59	AGE60-64	AGE65-69	AGE70-74	AGE75-79	AGE80-84	AGE85-89 90++}
//http://apps.who.int/gho/data/view.main.60740?lang=en
double AGESPECIFICMORTALITY[20] = { 0.0472, 0.0148, 0.00618, 0.00453, 0.00722, 0.01006, 0.01058, 0.01334, 0.01748, 0.02227, 0.03097, 0.04654, 0.06698, 0.11444, 0.16777, 0.26212, 0.33928, 0.43998, 0.5684, 1 };
int MAXAGE = 90;// 75;



//person
//we try to keep the rates at montly values if purely yearly average disease duration will be 5-6 months
double CUMFASTPROG = 0.138;//adults//cumulative value in 5 years  
double FASTPROGRATECOEF[5] = { 0.604, 0.85, 0.93, 0.98, 1 };//computed based on vynnycky
//double FASTPROGRATE=.03/12;
double RECOVERYRATE = .9 / 12;//if put to yearly rate, average disease duration will be 5-6 months!!!
double FATALITY = 0.12 / 12;
double RELAPSE = 0.06 / 12;// 0.06 / 12;


double MAXINFECTIOUSNESS = 1;//maxInfectiousness
int MAXINFDURATION = 9;//month

double IMMUNITYELTB = 0.5;
double IMMUNITYLLTB = 0.5;
double IMMUNITYATB = 1;
double IMMUNITYREC = 0;






