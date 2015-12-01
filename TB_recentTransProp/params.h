#include <string>
#include <vector>
using namespace std;

extern vector<vector<double>> vRECEIPE;

//extern double LAMBDA;
extern double MIGRATE;

extern int REPLICATIONS;
extern int MAXRUNTIME;
extern int TRANSIENT;
extern int POPULATIONSIZE;
extern double INITIALSTATES[];

//community

//extern double SLOWPROGRATE;// = .001 / 12;// annual risk of reactivation per ferebee(1970)


extern int MAXAGE;

extern double AGESPECIFICMORTALITY[];


//person
extern double CUMFASTPROG; // age-specific cumulative value in 5 years  
extern double FASTPROGRATECOEF[];// [5] = { 0.604, 0.85, 0.93, 0.98, 1 };//computed based on vynnycky

extern double RECOVERYRATE;// = .95 / 12;//if put to yearly rate, average disease duration will be 5-6 months!!!
extern double FATALITY;// = 0.12 / 12;

extern double MAXINFECTIOUSNESS;// = 1;//maxInfectiousness
extern int MAXINFDURATION;// = 9;//month

extern double IMMUNITYELTB;// = 0.5; 
extern double IMMUNITYLLTB;// = 0.5;
extern double IMMUNITYATB;//= 1; 
extern double IMMUNITYREC;// = 0;

extern double RELAPSE;

