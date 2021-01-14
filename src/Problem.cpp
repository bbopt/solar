#include "Problem.hpp"

/*---------------------------------------------------------*/
/*              create the list of problems                */
/*---------------------------------------------------------*/
void create_problems ( std::vector<Problem> & problems ) {

  // pb_id, f, n, m
  problems.push_back ( Problem ( "MAXNRG_H1"    , "total solar energy on the receiver"   , 1,  9,  5 ) ); // #1
  problems.push_back ( Problem ( "MINSURF_H1"   , "total heliostats field surface"       , 2, 14, 13 ) ); // #2
  problems.push_back ( Problem ( "MINCOST_C1"   , "total investment cost"                , 3, 20, 13 ) ); // #3
  problems.push_back ( Problem ( "MINCOST_C2"   , "total investment cost"                , 4, 29, 16 ) ); // #4
  problems.push_back ( Problem ( "MAXCOMP_HTF1" , "compliance to a demand profile"       , 5, 20, 12 ) ); // #5
  problems.push_back ( Problem ( "MINCOST_TS"   , "cost of storage"                      , 6,  5,  6 ) ); // #6
  problems.push_back ( Problem ( "MAXEFF_RE"    , "receiver efficiency"                  , 7,  7,  6 ) ); // #7
  problems.push_back ( Problem ( "MAXHF_MINCOST", "heliostat field performance and cost" , 8, 13,  9 ) ); // #8
  problems.push_back ( Problem ( "MAXNRG_MINPAR", "power and losses"                     , 9, 29, 17 ) ); // #9
 
}

/*---------------------------------------------------------*/
/*              find a problem from its pb_id              */
/*---------------------------------------------------------*/
const Problem * find_problem ( const std::vector<Problem> & problems , const std::string & pb_id ) {

  std::string id = toupper(pb_id);
  size_t i;

  // 1. id may be an integer in [1;problems.size()]:
  // -----------------------------------------------
  try {
    i = std::stoi(id);
  }
  catch ( ... ) {
    i = 0;
  }
  if ( i > 0 && i <= problems.size() )
    return &problems[i-1];

  // 2. id may be with format "solarx":
  // ----------------------------------
  if ( id.length() > 5 ) {
    std::string s = id.substr (0,5);  
    if ( s == "SOLAR" ) {
      std::string s_index = id.substr (5,id.length());
      try {
	i = std::stoi(s_index);
      }
      catch ( ... ) {
	i = 0;
      }
      if ( i > 0 && i <= problems.size() )
	return &problems[i-1];
    }
  }
  
  // 3. we search for id == problem.pb_id:
  // -------------------------------------
  for ( size_t i = 0 ; i < problems.size() ; ++i )
    if ( problems[i].get_pb_id() == id )
      return &problems[i];
  
  return NULL;
}

/*------------------------------------------*/
/*                  toupper                 */
/*------------------------------------------*/
std::string toupper ( std::string s ) {
  for ( size_t i = 0 ; i < s.size() ; ++i )
    s[i] = std::toupper(s[i]);
  return s;
}
