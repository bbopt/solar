#ifndef _HELPFUNCTIONS_H_
#define _HELPFUNCTIONS_H_

#include "Scenario.hpp"
#include <algorithm>
#include <cstring>
#include <iostream>
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
