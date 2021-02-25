/*-------------------------------------------------------------------------------*/
/*  SOLAR - The solar thermal power plant simulator                              */
/*  https://github.com/bbopt/solar                                               */
/*                                                                               */
/*  Miguel Diago, Sebastien Le Digabel, Mathieu Lemyre-Garneau, Bastien Talgorn  */
/*                                                                               */
/*  Polytechnique Montreal / GERAD                                               */
/*  sebastien.le-digabel@polymtl.ca                                              */
/*                                                                               */
/*  This program is free software: you can redistribute it and/or modify it      */
/*  under the terms of the GNU Lesser General Public License as published by     */
/*  the Free Software Foundation, either version 3 of the License, or (at your   */
/*  option) any later version.                                                   */
/*                                                                               */
/*  This program is distributed in the hope that it will be useful, but WITHOUT  */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  */
/*  for more details.                                                            */
/*                                                                               */
/*  You should have received a copy of the GNU Lesser General Public License     */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.         */
/*                                                                               */
/*-------------------------------------------------------------------------------*/
#ifndef __TIME_MANAGER_H__
#define __TIME_MANAGER_H__

class Time_Manager {
  
private:

  int _numberOfIncrements; 
  int _currentTime;       // minutes
  int _sizeOfIncrements;  // minutes
  int _incrementsCounter;

public:

  Time_Manager ( int numberOfIncrements, int currentTime, int sizeOfIncrements ) :
    _numberOfIncrements ( numberOfIncrements ) ,
    _currentTime        ( currentTime        ) ,
    _sizeOfIncrements   ( sizeOfIncrements   ) ,
    _incrementsCounter  ( 0                  )   {}

  Time_Manager ( const Time_Manager & t ) :
    _numberOfIncrements ( t._numberOfIncrements ) ,
    _currentTime        ( t._currentTime        ) ,
    _sizeOfIncrements   ( t._sizeOfIncrements   ) ,
    _incrementsCounter  ( t._incrementsCounter  )   {}
  
  ~Time_Manager ( void ) {}

  void fTimeIncrement ( void );
  void fResetTime     ( void );
  
  // set methods:
  void set_currentTime       ( int currentTime       ) { _currentTime       = currentTime;       }
  void set_sizeOfIncrements  ( int sizeOfIncrements  ) { _sizeOfIncrements  = sizeOfIncrements;  }
  void set_incrementsCounter ( int incrementsCounter ) { _incrementsCounter = incrementsCounter; }

  //get methods:
  int get_currentTime        ( void ) const { return _currentTime;        }
  int get_incrementsCounter  ( void ) const { return _incrementsCounter;  }
  int get_sizeOfIncrements   ( void ) const { return _sizeOfIncrements;   }
  int get_numberOfIncrements ( void ) const { return _numberOfIncrements; }
};

#endif
