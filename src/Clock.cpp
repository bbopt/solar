#include "Clock.hpp"

Clock::Clock(int numberOfIncrement, int currentTime, int sizeOfIncrements)
:_numberOfIncrements(numberOfIncrement),
_currentTime(currentTime),
_sizeOfIncrements(sizeOfIncrements),
_incrementsCounter(0)
{}

void Clock::fTimeIncrement()
{
	_currentTime += _sizeOfIncrements;
	++_incrementsCounter;
}

void Clock::fResetTime()
{
	_currentTime = 0;
	_incrementsCounter = 0;
}
