#ifndef _PROBLEM_H_
#define _PROBLEM_H_

#include <vector>
#include <iostream>

class Problem {

private:

  std::string _pb_id;
  std::string _f_description;
  int         _index , _n , _m;

public:

  Problem ( const std::string & pb_id, const std::string & f, int index , int n, int m )
    : _pb_id(pb_id), _f_description(f), _index(index) , _n(n), _m(m) {}
  
  Problem ( const Problem & p ) : _pb_id         ( p._pb_id         ) ,
				  _f_description ( p._f_description ) ,
				  _index         ( p._index         ) ,
				  _n             ( p._n             ) ,
				  _m             ( p._m             )   {}
  ~Problem ( void ) {}

  const std::string & get_pb_id         ( void ) const { return _pb_id; }
  const std::string & get_f_description ( void ) const { return _f_description; }

  int get_index ( void ) const { return _index; }
  int get_n     ( void ) const { return _n; }
  int get_m     ( void ) const { return _m; }
  
};


// creates the list of problems:
void create_problems ( std::vector<Problem> & problems );

// find a problem from its pb_id:
const Problem * find_problem ( const std::vector<Problem> & problems , const std::string & pb_id );

// put a sting in upper cases:
std::string toupper ( std::string );

#endif
