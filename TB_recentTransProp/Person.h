 // Person.h
// modeling a single individual.

#ifndef __PERSON_H
#define __PERSON_H





class Person {
public:
	Person(int id);
	Person(int age, int id);
	~Person();
	void reborn(int id);
	void reborn(int age, int id);

	inline int getID() { return _nID; }	                  inline void setId(int n){ _nID = n; }
	//	int getBithTime(){return _nBirthTime;} void setBithTime(int t){ _nBirthTime=t;}
	inline int getAge() { return _nAge; }			inline void setAge(int n) { _nAge = n; }	inline void addAge(){ _nAge++; }
	inline int getState() { return _nState; }		inline void setState(int n) { _nState = n; }
	//int getGender(){return _nGender;}	void setGender(int g) {_nGender=g;}
	inline double getImmunity(){ return _dImmunity; }					inline void setImmunity(double d){ _dImmunity = d; }
	inline double getInfectiousness() { return _dInfectiousness; }		inline void setInfectiousness(double d){ _dInfectiousness = d; }

	inline int getELTBTime() { return _nELTBTime; }		inline void setELTBTime(int n) { _nELTBTime = n; }
	inline int getLLTBTime() { return _nLLTBTime; }		inline void setLLTBTime(int n) { _nLLTBTime = n; }
	inline int getATBTime() { return _nATBTime; }		inline  void setATBTime(int n) { _nATBTime = n; }
	inline int getRECTime() { return _nRECTime; }		inline  void setRECTime(int n) { _nRECTime = n; }
	inline int getTimeOfFastProg(){ return _tOfFastProg; }
	inline int getTimeOfSlowProg(){ return _tOfSlowProg; }
	inline int getTimeOfFatality(){ return _tOfFatality; }
	inline int getTimeOfRecovery(){ return _tOfRecovery; }
	inline int getTimeOfRelapse(){ return _tOfRelapse; }



	inline int getStrain() { return _nStrain; }			 inline void setStrain(int s) { _nStrain = s; }



	void enterELTB(gsl_rng *rng, int time);
	void enterLLTB(gsl_rng *rng, int time);
	void enterLLTB(gsl_rng *rng, int time,double spRate);
	void enterATB(gsl_rng *rng, int time);
	void enterREC(gsl_rng *rng, int time);
	void updateInfectiousness(int t);


	bool _bRecovered;


protected:
	int _nID;            // unique identifier
	int _nAge;

	int _nState;
	int _nStrain;
    
	double _dInfectiousness;  //0 or >0
	double _dImmunity;

	int _nELTBTime;
    int _nLLTBTime;
    int _nATBTime;
    int _nRECTime;
	int _tOfSlowProg;
	int _tOfFastProg;
	int _tOfRecovery;
	int _tOfFatality;
	int _tOfRelapse;


};

#endif
