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
#ifndef __CENTRAL_RECEIVER_H__
#define __CENTRAL_RECEIVER_H__

#include "MoltenSalt.hpp"
#include "Constants.hpp"
#include <vector>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <cmath>

class CentralReceiver {

private:
  
  // fluid conditions at inlet and outlet:
  MoltenSalt * _input;
  MoltenSalt * _output;

  // input design parameters:
  double _apertureHeight;
  double _apertureWidth;
  double _insulationThickness;
  double _tubesInsideDiameter;
  double _tubesOutsideDiameter;
  int    _numberOfTubes;

  // consequential attributes:
  int    _numberOfPasses;
  double _receiverSurfaceArea;
  double _receiverEfficiency;
  
  double computeEmissionLosses   ( double ) const;
  double computeConvectionLosses ( double ) const;
  double computeConductionLosses ( double ) const;
  double fSolveForT              ( double, double, double, double, double) const;

  // simulation data:
  std::vector<double> _losses;
  std::vector<double> _efficiency;
  std::vector<double> _surfaceTemperature;
  std::vector<double> _msRate;

  double computeReflectionLosses ( double Q_in ) const {
    return Q_in * RECEIVER_SURF_REFLECTIVITY * ((_apertureHeight*_apertureWidth)/_receiverSurfaceArea);
  }
  
public:
  CentralReceiver  ( MoltenSalt*, MoltenSalt*, double, double, double, double, double, int );
  ~CentralReceiver ( void ) {}
 
  double computeEnergyToFluid ( double Q_in );

  //function shall provide an option to calculate if the output temperature
  //is imposed or if the mass flow is imposed.
  
  //for imposed temperature the mass flow will vary but varying the mass flow
  //effectively changes the convection transfer... Or not, because the
  // temperature being fixed means that the temperature distribution should be
  //the same wether the flow is fast or slow.
  
  //For imposted flow the output temperature will vary and the higher it is,
  //the less heat transfer will go on. The temperature distribution must be
  //determined along the tubes so that the steady-state transfer may be determined.
  //that is, the amount of energy absorbed to the fluid isn't the same in the transient
  //regime as in the steady-state regime. Though the difference shouldn't be much...
  
  double computeYieldPressure ( void ) const {
    return ( _tubesOutsideDiameter - _tubesInsideDiameter ) *
      SS316_YIELD_PRESSURE / (0.5*(_tubesInsideDiameter + _tubesOutsideDiameter));
  }
  
  double computePressureInTubes ( void ) const;
 
};

#endif
