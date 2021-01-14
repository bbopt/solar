#ifndef _CLOCK_H_
#define _CLOCK_H_

class Clock
{
private:
	int _numberOfIncrements; 
	int _currentTime;		//minutes
	int _sizeOfIncrements;	//minutes
	int _incrementsCounter;

public:
	Clock(int, int, int);
	~Clock(){}

	//-------------------------------------------------------------------
	void fTimeIncrement();

	//set methods
	void set_currentTime(int currentTime)           { _currentTime = currentTime; }
	void set_sizeOfIncrements(int sizeOfIncrements) { _sizeOfIncrements = sizeOfIncrements; }
	void set_incrementsCounter(int incrementsCounter) { _incrementsCounter = incrementsCounter; }
	void fResetTime();


	//get methods
	int  get_currentTime() const            { return _currentTime;}
	int  get_incrementsCounter() const      { return _incrementsCounter; }
	int get_sizeOfIncrements() const		{ return _sizeOfIncrements; }
	int get_numberOfIncrements() const		{ return _numberOfIncrements; }
};

#endif
