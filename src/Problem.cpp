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
#include "Problem.hpp"

/*---------------------------------------------------------------------*/
/*  is a specific output of the problem stochastic or deterministic ?  */
/*  output_index is equal to -1 by default: in this case the function  */
/*  returns true if the problem is stochastic                          */
/*---------------------------------------------------------------------*/
bool Problem::is_stochastic ( int output_index ) const {

  // SOLAR1:
  if ( _pb_id == "MAXNRG_H1" &&
       ( output_index == -1 || output_index == 0 ) )
    return true;

  // SOLAR2:
  if ( _pb_id == "MINSURF_H1" &&
       ( output_index == -1 || output_index == 2 || output_index == 7 || output_index == 8 || output_index == 9 ) )
    return true;

  // SOLAR3:  
  if ( _pb_id == "MINCOST_C1" &&
       ( output_index == -1 || output_index == 2 || output_index == 6 || output_index == 7 || output_index == 8 || output_index == 13 ) )
    return true;
  
  // SOLAR4:
  if ( _pb_id == "MINCOST_C2" &&
       ( output_index == -1 || output_index == 2 || output_index == 6 || output_index == 7 || output_index == 8 || output_index == 9 || output_index == 13 ) )
    return true; 
 
  // SOLAR5: Deterministic
  if ( _pb_id == "MAXCOMP_HTF1" )
    return false;

  // SOLAR6: Deterministic
  if ( _pb_id == "MINCOST_TS" )
    return false;

  // SOLAR7:
  if ( _pb_id == "MAXEFF_RE" &&
       ( output_index == -1 || output_index == 0 || output_index == 2 || output_index == 6 ) )
    return true;  
 
  // SOLAR8:
  if ( _pb_id == "MAXHF_MINCOST" &&
       ( output_index == -1 || output_index == 6 || output_index == 9 || output_index == 10 ) )
    return true;
 
  // SOLAR9:
  if ( _pb_id == "MAXNRG_MINPAR"  &&
       ( output_index == -1 || output_index ==  0 || output_index ==  1 || output_index ==  3 ||
	 output_index ==  8 || output_index ==  9 || output_index == 10 || output_index == 11 ||
	 output_index == 15 ) )
    return true;

  // SOLAR10: Deterministic
  if ( _pb_id == "MINCOST_UNCONSTRAINED" )
    return false;
  
  return false;
}

/*---------------------------------------------------------*/
/*              create the list of problems                */
/*---------------------------------------------------------*/
void create_problems ( std::vector<Problem> & problems ) {

  // pb_id, f, n, m
  problems.push_back ( Problem ( "MAXNRG_H1"            , "total solar energy on the receiver"   ,  1, 1,  9,  5 ) ); // #1
  problems.push_back ( Problem ( "MINSURF_H1"           , "total heliostats field surface"       ,  2, 1, 14, 12 ) ); // #2
  problems.push_back ( Problem ( "MINCOST_C1"           , "total investment cost"                ,  3, 1, 20, 13 ) ); // #3
  problems.push_back ( Problem ( "MINCOST_C2"           , "total investment cost"                ,  4, 1, 29, 16 ) ); // #4
  problems.push_back ( Problem ( "MAXCOMP_HTF1"         , "compliance to a demand profile"       ,  5, 1, 20, 12 ) ); // #5
  problems.push_back ( Problem ( "MINCOST_TS"           , "cost of storage"                      ,  6, 1,  5,  6 ) ); // #6
  problems.push_back ( Problem ( "MAXEFF_RE"            , "receiver efficiency"                  ,  7, 1,  7,  6 ) ); // #7
  problems.push_back ( Problem ( "MAXHF_MINCOST"        , "heliostat field performance and cost" ,  8, 2, 13,  9 ) ); // #8
  problems.push_back ( Problem ( "MAXNRG_MINPAR"        , "power and losses"                     ,  9, 2, 29, 17 ) ); // #9
  problems.push_back ( Problem ( "MINCOST_UNCONSTRAINED", "cost of storage + penalties"          , 10, 1,  5,  0 ) ); // #10
}

/*---------------------------------------------------------*/
/*              find a problem from its pb_id              */
/*---------------------------------------------------------*/
const Problem * find_problem ( const std::vector<Problem> & problems , const std::string & pb_id ) {

  std::string id = toupper(pb_id);
  int  i , nb_pb = static_cast<int>(problems.size());

  // 1. id may be an integer in [1;problems.size()]:
  // -----------------------------------------------
  if  ( string_to_int (id,i) && i > 0 && i <= nb_pb )
    return &problems[i-1];

  // 2. id may be with format "solarx":
  // ----------------------------------
  if ( id.length() > 5 ) {
    std::string s = id.substr(0,5);  
    if ( s == "SOLAR" ) {
      std::string s_index = id.substr (5,id.length());
      if ( !string_to_int(s_index,i) )
	return NULL;
      if ( i > 0 && i <= nb_pb )
	return &problems[i-1];
    }
  }
  
  // 3. we search for id == problem.pb_id:
  // -------------------------------------
  for ( i = 0 ; i < nb_pb ; ++i )
    if ( problems[i].get_pb_id() == id )
      return &problems[i];
  
  return NULL;
}
