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
#include "sampling.hpp"

/*-----------------------------------------------------------------*/
/*  the sampling: evaluates one instance with random or LH points  */
/*-----------------------------------------------------------------*/
void sampling ( void ) {

  // Example of use:
  //   ./solar10 > solar10.txt  &
  //   disown %1
  
  // list of problems:
  std::vector<Problem> problems;
  create_problems ( problems );

  // number of sampling points:
  size_t p = 30;
  
  // set the instance:
  int instance = 6;

  int n = 0 , m = 0;
  const Problem * pb = NULL;
  double        * lb = NULL;
  double        * ub = NULL;
  
  if ( instance == 1 ) {
    n  = 9;
    m  = 5;
    pb = find_problem ( problems, "1" );
    lb = new double[n];
    ub = new double[n];

    lb[0] = 1;  ub[0] = 40;
    lb[1] = 1;  ub[1] = 40;
    lb[2] = 20; ub[2] = 250;
    lb[3] = 1;  ub[3] = 30;
    lb[4] = 1;  ub[4] = 30;
    lb[5] = 1;  ub[5] = 10000; // artificial upper bound based on observations
    lb[6] = 1;  ub[6] = 89;
    lb[7] = 0;  ub[7] = 20;
    lb[8] = 1;  ub[8] = 20;
  }
  else if ( instance == 2 ) {
    n  = 14;
    m  = 12;
    pb = find_problem ( problems, "2" );
    lb = new double[n];
    ub = new double[n];
    lb[0] = 1;      ub[0] = 40;
    lb[1] = 1;      ub[1] = 40;
    lb[2] = 20;     ub[2] = 250;
    lb[3] = 1;      ub[3] = 30;
    lb[4] = 1;      ub[4] = 30;
    lb[5] = 1;      ub[5] = 10000; // artificial upper bound based on observations
    lb[6] = 1;      ub[6] = 89;
    lb[7] = 0;      ub[7] = 20;
    lb[8] = 1;      ub[8] = 20;
    lb[9] = 793;    ub[9] = 995 ;
    lb[10] = 1;     ub[10] = 9424;
    lb[11] = 0.01;  ub[11] = 5;
    lb[12] = 0.005; ub[12] = 0.1 ;
    lb[13] = 0.005; ub[13] = 0.1;   
  }
  else if ( instance == 3 ) {
    n  = 20;
    m  = 13;
    pb = find_problem ( problems, "3" );
    lb = new double[n];
    ub = new double[n];
    lb[0] = 1;      ub[0] = 40;
    lb[1] = 1;      ub[1] = 40;
    lb[2] = 20;     ub[2] = 250;
    lb[3] = 1;      ub[3] = 30;
    lb[4] = 1;      ub[4] = 30;
    lb[5] = 1;      ub[5] = 10000; // artificial upper bound based on observations
    lb[6] = 1;      ub[6] = 89;
    lb[7] = 0;      ub[7] = 20;
    lb[8] = 1;      ub[8] = 20;
    lb[9] = 793;    ub[9] = 995;
    lb[10] = 1;     ub[10] = 50;
    lb[11] = 1;     ub[11] = 30;
    lb[12] = 0.01;  ub[12] = 5;
    lb[13] = 0.01;  ub[13] = 5;
    lb[14] = 495;   ub[14] = 650;
    lb[15] = 1;     ub[15] = 9424;
    lb[16] = 0.01;  ub[16] = 5;
    lb[17] = 0.005; ub[17] = 0.1;
    lb[18] = 0.005; ub[18] = 0.1;
    lb[19] = 1;     ub[19] = 8;  
  }
  else if ( instance == 4 ) {
    n  = 29;
    m  = 16;
    pb = find_problem ( problems, "4" );
    lb = new double[n];
    ub = new double[n];
    lb[0] = 1;      ub[0] = 40;
    lb[1] = 1;      ub[1] = 40;
    lb[2] = 20;     ub[2] = 250;
    lb[3] = 1;      ub[3] = 30;
    lb[4] = 1;      ub[4] = 30;
    lb[5] = 1;      ub[5] = 10000; // artificial upper bound based on observations
    lb[6] = 1;      ub[6] = 89;
    lb[7] = 0;      ub[7] = 20;
    lb[8] = 1;      ub[8] = 20;
    lb[9] = 793;    ub[9] = 995;
    lb[10] = 1;     ub[10] = 50;
    lb[11] = 1;     ub[11] = 30;
    lb[12] = 0.01;  ub[12] = 5;
    lb[13] = 0.01;  ub[13] = 5;
    lb[14] = 495;   ub[14] = 650;
    lb[15] = 1;     ub[15] = 7853;
    lb[16] = 0.01;  ub[16] = 5;
    lb[17] = 0.005; ub[17] = 0.1;
    lb[18] = 0.006; ub[18] = 0.1;
    lb[19] = 0.007; ub[19] = 0.2;
    lb[20] = 0.5;   ub[20] = 10;
    lb[21] = 0.005; ub[21] = 0.1;
    lb[22] = 0.006; ub[22] = 0.1;
    lb[23] = 0.15;  ub[23] = 0.4;
    lb[24] = 2;     ub[24] = 20;     // artificial upper bound based on observations
    lb[25] = 1;     ub[25] = 100000; // artificial upper bound based on observations
    lb[26] = 1;     ub[26] = 10;
    lb[27] = 1;     ub[27] = 9;
    lb[28] = 1;     ub[28] = 8;
  }
  else if ( instance == 5 ) {
    n  = 20;
    m  = 12;
    pb = find_problem ( problems, "5" );
    lb = new double[n];
    ub = new double[n];
    lb[0] = 793;  ub[0] = 995;
    lb[1] = 1;    ub[1] = 30;
    lb[2] = 1;    ub[2] = 30;
    lb[3] = 0.01; ub[3] = 2;
    lb[4] = 0.01; ub[4] = 2;
    lb[5] = 495;  ub[5] = 650;
    lb[6] = 1;    ub[6] = 1884;
    lb[7] = 0.1;  ub[7] = 2;
    lb[8] = 0.005;  ub[8] = 0.1;
    lb[9] = 0.005;  ub[9] = 0.1;
    lb[10] = 0.006; ub[10] = 0.2;
    lb[11] = 0.5;   ub[11] = 10;
    lb[12] = 0.005; ub[12] = 0.1;
    lb[13] = 0.006; ub[13] = 0.1;
    lb[14] = 0.15;  ub[14] = 0.4;
    lb[15] = 2;     ub[15] = 20;     // artificial upper bound based on observations
    lb[16] = 1;     ub[16] = 100000; // artificial upper bound based on observations
    lb[17] = 1;     ub[17] = 10;
    lb[18] = 1;     ub[18] = 9;
    lb[19] = 1;     ub[19] = 8;
  }
  else if ( instance == 6 ) {
    n  = 5;
    m  = 6;
    pb = find_problem ( problems, "6" );
    lb = new double[n];
    ub = new double[n];
    lb[0] = 793;  ub[0] = 995;
    lb[1] = 2;    ub[1] = 50;
    lb[2] = 2;    ub[2] = 30;
    lb[3] = 0.01; ub[3] = 5;
    lb[4] = 0.01; ub[4] = 5;
  }
  else if ( instance == 7 ) {
    n  = 7;
    m  = 6;
    pb = find_problem ( problems, "7" );
    lb = new double[n];
    ub = new double[n];
    lb[0] = 1;      ub[0] = 30;
    lb[1] = 1;      ub[1] = 30;
    lb[2] = 793;    ub[2] = 995;
    lb[3] = 1;      ub[3] = 8567;
    lb[4] = 0.01;   ub[4] = 5;
    lb[5] = 0.005;  ub[5] = 0.1;
    lb[6] = 0.0055; ub[6] = 0.1;
  }
  else if ( instance == 8 ) {
    n  = 13;
    m  = 9;
    pb = find_problem ( problems, "8" );
    lb = new double[n];
    ub = new double[n];
    lb[0] = 1;      ub[0] = 40;
    lb[1] = 1;      ub[1] = 40;
    lb[2] = 20;     ub[2] = 250;
    lb[3] = 1;      ub[3] = 30;
    lb[4] = 1;      ub[4] = 30;
    lb[5] = 1;      ub[5] = 10000; // artificial upper bound based on observations
    lb[6] = 1;      ub[6] = 89;
    lb[7] = 0;      ub[7] = 20;
    lb[8] = 1;      ub[8] = 20;
    lb[9] = 1;      ub[9] = 7853;
    lb[10] = 0.01;  ub[10] = 5;
    lb[11] = 0.005; ub[11] = 0.1;
    lb[12] = 0.006; ub[12] = 0.1;    
  }
  else if ( instance == 9 ) {
    n  = 29;
    m  = 17;
    pb = find_problem ( problems, "9" );
    lb = new double[n];
    ub = new double[n];
    lb[0] = 1;      ub[0] = 40;
    lb[1] = 1;      ub[1] = 40;
    lb[2] = 20;     ub[2] = 250;
    lb[3] = 1;      ub[3] = 30;
    lb[4] = 1;      ub[4] = 30;
    lb[5] = 1;      ub[5] = 10000; // artificial upper bound based on observations
    lb[6] = 1;      ub[6] = 89;
    lb[7] = 0;      ub[7] = 20;
    lb[8] = 1;      ub[8] = 20;
    lb[9] = 793;    ub[9] = 995;
    lb[10] = 1;     ub[10] = 50;
    lb[11] = 1;     ub[11] = 30;
    lb[12] = 0.01;  ub[12] = 5;
    lb[13] = 0.01;  ub[13] = 5;
    lb[14] = 495;   ub[14] = 650;
    lb[15] = 1;     ub[15] = 7853;
    lb[16] = 0.01;  ub[16] = 5;
    lb[17] = 0.005; ub[17] = 0.1;
    lb[18] = 0.006; ub[18] = 0.1;
    lb[19] = 0.007; ub[19] = 0.2;
    lb[20] = 0.5;   ub[20] = 10;
    lb[21] = 0.005; ub[21] = 0.1;
    lb[22] = 0.006; ub[22] = 0.1;
    lb[23] = 0.15;  ub[23] = 0.4;
    lb[24] = 2;     ub[24] = 20;     // artificial upper bound based on observations
    lb[25] = 1;     ub[25] = 100000; // artificial upper bound based on observations
    lb[26] = 1;     ub[26] = 10;
    lb[27] = 1;     ub[27] = 9;
    lb[28] = 1;     ub[28] = 8; 
  }
  else if ( instance == 10 ) {
    n  = 5;
    m  = 0;
    pb = find_problem ( problems, "10" );
    lb = new double[n];
    ub = new double[n];
    lb[0] = 793;  ub[0] = 995;
    lb[1] = 2;    ub[1] = 50;
    lb[2] = 2;    ub[2] = 30;
    lb[3] = 0.01; ub[3] = 5;
    lb[4] = 0.01; ub[4] = 5;
  }
    

  std::vector<double *> points;
    
  // get_random_points ( p, n, lb, ub, points );
  get_LH_points ( p, n, lb, ub, points );

  delete [] lb;
  delete [] ub;

  // rounding of discrete variables:
  for ( size_t k = 0 ; k < p ; ++k ) {
    if ( instance == 1 ) {
      points[k][5] = myround(points[k][5]);
    }
    else if ( instance == 2 ) {
      points[k][ 5] = myround(points[k][ 5]);
      points[k][10] = myround(points[k][10]);
    }
    else if ( instance == 3 ) {
      points[k][ 5] = myround(points[k][ 5]);
      points[k][15] = myround(points[k][15]);
      points[k][19] = myround(points[k][19]);
    }
    else if ( instance == 4 ) {
      points[k][ 5] = myround(points[k][ 5]);
      points[k][15] = myround(points[k][15]);
      points[k][24] = myround(points[k][24]);
      points[k][25] = myround(points[k][25]);
      points[k][26] = myround(points[k][26]);
      points[k][27] = myround(points[k][27]);
      points[k][28] = myround(points[k][28]);
    }  
    else if ( instance == 5 ) {
      points[k][ 6] = myround(points[k][ 6]);
      points[k][15] = myround(points[k][15]);
      points[k][16] = myround(points[k][16]);
      points[k][17] = myround(points[k][17]);
      points[k][18] = myround(points[k][18]);
      points[k][19] = myround(points[k][19]);
    }
    else if ( instance == 7 ) {
      points[k][3] = myround(points[k][3]);
    }
    else if ( instance == 8 ) {
      points[k][5] = myround(points[k][5]);
      points[k][9] = myround(points[k][9]);
    }
    else if ( instance == 9 ) {
      points[k][ 5] = myround(points[k][ 5]);
      points[k][15] = myround(points[k][15]);
      points[k][24] = myround(points[k][24]);
      points[k][25] = myround(points[k][25]);
      points[k][26] = myround(points[k][26]);
      points[k][27] = myround(points[k][27]);
      points[k][28] = myround(points[k][28]);
    }
  }

  // For only displaying the LH points:
  // {
  //   for ( size_t k = 0 ; k < p ; ++k ) {
  //     for ( int i = 0 ; i < n ; ++i )
  //  	std::cout << points[k][i] << " ";
  //      std::cout << std::endl;
  //   }
  //   for ( size_t k = 0 ; k < p ; ++k )
  //     delete [] points[k];
  //   return;
  // }
  
  std::vector<double> all_f;
  std::vector<double> all_f1;
  std::vector<double> all_f2;
  std::vector<int>    all_i;
  std::vector<double> succ_f;
  std::vector<int>    succ_i;
  int    nb_evals = 0;
  double best_f   = 1E20;
  
  // display format:
  std::cout << "# x F(x) cnt_eval simu_completed time:\n\n";
  
  for ( size_t k = 0 ; k < p ; ++k ) {
  
    std::ostringstream output_stream;
    Evaluator evaluator ( *pb, output_stream );
    evaluator.set_x ( points[k] );

    bool        simulation_completed = false; // output flag
    bool        cnt_eval             = false; // output flag
    double      time;
    
    // the evaluation:
    {
      Clock clock;
      
      int         x_index      = 0;     // only one point
      int         seed         = 0;     // random seed
      double      fidelity     = 1.0;   // full fidelity
      int         replications = 1;     // one replication
      bool        verbose      = false; // no display
      std::string error_msg;                    // error message
      evaluator.eval_x ( x_index, seed, fidelity, replications, simulation_completed, cnt_eval, error_msg, verbose );

      time = clock.get_CPU_time();
    }
  
    // get the outputs: this puts the outputs into output_stream:
    evaluator.display_outputs();
    std::istringstream outputs ( output_stream.str() );

    double f, f1, f2;

    if ( instance == 8 || instance == 9 )
      outputs >> f1 >> f2;
    else
      outputs >> f;   
  
    if ( cnt_eval )
      ++nb_evals;
    
    bool feasible = false;
    if ( cnt_eval && simulation_completed ) {
      double c;
      feasible = true;
      for ( int j = 0 ; j < m ; ++j ) {
	outputs >> c;
	if ( c > 0.0 ) {	  
	  feasible = false;
	  break;
	}
      }
    }

    if ( feasible ) {

      // biobjective:
      if ( instance == 8 || instance == 9 ) {
	all_f1.push_back ( f1 );
	all_f2.push_back ( f2 );
	all_i.push_back  ( nb_evals );
      }

      // single-objective:
      else {
      	all_f.push_back ( f );
	all_i.push_back ( nb_evals );
	if ( f < best_f ) {
	  best_f = f;
	  succ_f.push_back ( f );
	  succ_i.push_back ( nb_evals );
	}
      }     
    }
    
    // display:
    std::cout << k << " ";
    for ( int i = 0 ; i < n ; ++i )
      std::cout << points[k][i] << " ";
    std::cout << output_stream.str() << " " << cnt_eval << " " << simulation_completed << " " << time << std::endl;
    
  }

  // display all feasible evals:
  std::cout << "\nall feasible evals:\n\n";
  for ( size_t k = 0 ; k < all_i.size(); ++k ) {
    std::cout << all_i[k] << " ";
    if ( instance == 8 || instance == 9 )
      std::cout << all_f1[k] << " " << all_f2[k] << std::endl;
    else
      std::cout << all_f[k] << std::endl;
  }
  
  // display successes:
  if ( instance != 8 && instance != 9 ) {
    std::cout << "\nsuccesses:\n\n";
    for ( size_t k = 0 ; k < succ_i.size(); ++k )
      std::cout << succ_i[k] << " " << succ_f[k] << std::endl;
  }
 
  // delete points:
  for ( size_t k = 0 ; k < p ; ++k )
    delete [] points[k];
}

void get_random_points ( size_t p, size_t n, const double * lb, const double * ub, std::vector<double *> & points ) {
  for ( size_t k = 0 ; k < p ; ++k ) {
    double * x = new double[n];
    for (size_t i = 0 ; i < n ; ++i ) {
      x[i] = RNG::rand(lb[i],ub[i]);
    }   
    points.push_back(x);
  }
}

void get_LH_points (size_t p, size_t n, const double * lb, const double * ub, std::vector<double *> & points ) {
    size_t         pm1 = p-1;
  Random_Pickup ** rps = new Random_Pickup *[n];
  for (size_t k = 0 ; k < p ; ++k ) {
    double * x = new double[n];
    for (size_t i = 0 ; i < n ; ++i ) {
      if ( k==0 )
	rps[i] = new Random_Pickup(static_cast<int>(p));
      x[i] = lb[i] +
	(ub[i]-lb[i]) *
	( rps[i]->pickup() + RNG::rand(0.0,1.0) ) / p;
      if ( k==pm1 )
	delete rps[i];
    }
    points.push_back(x);
  }
  delete [] rps;
}

  // // Test LHS:
  // {
  //   int p = 100, n = 2;
  
  //   std::vector<double *> points;

  //   double lb[] = {0,0};
  //   double ub[] = {100,100};
    
  //   get_LH_points ( p, n, lb, ub, points );

  //   for ( int k = 0 ; k < p ; ++k ) {
  //     for ( int i = 0 ; i < n ; ++i )
  // 	std::cout << points[k][i] << ",";
  //     std::cout << std::endl;
  //     delete [] points[k];
  //   }
     
  // }
