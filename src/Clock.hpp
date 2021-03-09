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
#ifndef __CLOCK__
#define __CLOCK__

#include <ctime>

#ifdef _MSC_VER
#pragma warning(disable:4275)
#pragma warning(disable:4251)
#endif


class Clock {
        
private:
        
  time_t              _real_t0;          ///< Wall clock time measurement.
  clock_t             _CPU_t0;           ///< CPU time measurement.
  static const double _D_CLOCKS_PER_SEC; ///< System constant for CPU time measurement.
        
public:
        
  /// Constructor.
  Clock ( void ) : _CPU_t0 ( clock() ) { time (&_real_t0); }
        
  /// Copy constructor.
  /**
     \param c The copied object -- \b IN.
  */
  Clock ( const Clock & c ) : _real_t0 ( c._real_t0 ) , _CPU_t0 ( c._CPU_t0 ) {}
        
  /// Affectation operator.
  /**
     \param  c The right-hand side object -- \b IN.
     \return \c *this as the result of the affectation.
  */
  Clock & operator = ( const Clock & c )
  {
    _real_t0 = c._real_t0;
    _CPU_t0  = c._CPU_t0;
    return *this;
  }
        
  /// Destructor.
  virtual ~Clock ( void ) {}
        
  /// Reset the clock.
  void reset ( void )
  {
    time ( &_real_t0 );
    _CPU_t0 = clock();
  }
        
  /// Get wall clock time.
  /**
     \return The wall clock time.
  */
  int get_real_time ( void ) const;
        
  /// Get the CPU time.
  /**
     \return The CPU time.
  */
  double get_CPU_time ( void ) const
  {
    return ( static_cast<double>(clock()) - _CPU_t0 ) / _D_CLOCKS_PER_SEC;
  }
};

#endif
