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
#ifndef __HELPFUNCTIONS_H__
#define __HELPFUNCTIONS_H__

#include "Scenario.hpp"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>

void display_usage    ( std::ostream & out );
void display_info     ( std::ostream & out , const std::string & version );
void display_help     ( std::ostream & out , const std::vector<Problem> & , const std::string & pb_id );
void display_help     ( std::ostream & out , const std::vector<Problem> & );
void display_problems ( std::ostream & out , const std::vector<Problem> & );

void print_maxNrg_H1     ( std::ostream & out );
void print_minSurf_H1    ( std::ostream & out );
void print_minCost_C1    ( std::ostream & out );
void print_minCost_C2    ( std::ostream & out );
void print_maxComp_HTF1  ( std::ostream & out );
void print_minCost_TS    ( std::ostream & out );
void print_maxEff_RE     ( std::ostream & out );
void print_maxHF_minCost ( std::ostream & out );
void print_maxNrg_minPar ( std::ostream & out );


#endif
