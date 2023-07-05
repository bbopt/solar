# SOLAR v0.5.7 (July 2023)
The **SOLAR** blackbox optimization benchmarking framework.

### Compilation
In the following, `$SOLAR_HOME` denotes the directory where **SOLAR** has been downloaded.
Go to the [$SOLAR_HOME/src](src) directory and type `make`. It will generate the binary
executable `solar` located in [$SOLAR_HOME/bin](bin).

For Windows users, it is possible to directly use a binary [executable](./bin/solar_WINDOWS.exe).

### Validation
Your **SOLAR** installation must be validated. For this, type `./bin/solar -check`.
This will execute several tests during 10 to 20 minutes.

If the validation fails, please send an email to nomad@gerad.ca with the full output.

### Execution
Type `./bin/solar` and you will be guided with the following help:

```
Run SOLAR (basic)   : solar pb_id x.txt (add -v for verbose mode)
Run SOLAR (advanced): solar pb_id x.txt -seed=S -fid=F -rep=R -v
     S: Random seed: integer >=0 or "diff"; Default=0
     F: Fidelity: real in ]0;1]; Default=1.0 (truth)
     R: Replications: integer >= 1 ; Default=1

Validation: solar -check (can take several minutes)
Help(1)   : solar -h
Help(2)   : solar -h pb_id
Info      : solar -i
```
To run a simulation, type `$SOLAR_HOME/bin/solar pb_id x.txt -seed=S -fid=F -rep=R -v (optional)`.

The different options are:

```
 pb_id: Problem ID (see list of problems below)

 x.txt: Input vector: Point at which the simulator is evaluated
        Values separated with spaces
        It is possible to specify several vectors: Use one line for each

    -v: Verbose option

     S: Random seed:
          Some SOLAR instances are stochastic. This parameter impacts the value of stochastic outputs
          The seed is a natural integer
          If SOLAR is run twice at the same point with the same seed, it will give the same outputs
          The default value is 0
          Use -seed=diff to let SOLAR use a different random seed each time
          The random number generator can be validated by running 'solar -check'

     F: Fidelity of the simulator
          Real value in ]0;1]
          Default: 1.0, which corresponds to the "true blackbox", or the "truth"
          Any value in ]0;1[ corresponds to a "static surrogate" of the truth
          The execution time increases with the fidelity
          A good default static surrogate is -prec=0.5

     R: Number of replications
          Integer >= 1, default=1
          Number of times that the simulator is run at the same point
          Each replication uses a different random seed dependent on the -seed option
          The mean value of stochastic outputs is displayed
          It is not possible to use R>1 with deterministic instances

Help for a problem: solar pb_id or solar -h pb_id
```

The list of instances is:

```
	#	          pb_id	                              obj.(f)	# of objectives(p)	   # of var.(n)	# of constr.(m)
	1	      MAXNRG_H1	   total solar energy on the receiver	              1	              9	              5
	2	     MINSURF_H1	       total heliostats field surface	              1	             14	             13
	3	     MINCOST_C1	                total investment cost	              1	             20	             13
	4	     MINCOST_C2	                total investment cost	              1	             29	             16
	5	   MAXCOMP_HTF1	       compliance to a demand profile	              1	             20	             12
	6	     MINCOST_TS	                      cost of storage	              1	              5	              6
	7	      MAXEFF_RE	                  receiver efficiency	              1	              7	              6
	8         MAXHF_MINCOST	 heliostat field performance and cost	              2	             13	              9
	9    	  MAXNRG_MINPAR	                     power and losses	              2	             29	             17
       10 MINCOST_UNCONSTRAINED	          cost of storage + penalties	              1	              5	              0
```
List of best know values for single-objective instances (one replication, full fidelity, default seed of zero):
```
	SOLAR1 	-902,503.692418
	SOLAR2 	841,839.671915
	SOLAR3 	70,813,885.0684
	SOLAR4 	108,197,236.146
	SOLAR5 	-28.8817193932
	SOLAR6 	43,954,935.1836
	SOLAR7 	-4,972.88703862
	SOLAR10	42.542949
```
`SOLAR10` best solutions found with the help of J. Larson and GOOMBAH in the [IBCDFO package](https://github.com/POptUS/IBCDFO).

### Example

The command `$SOLAR_HOME/bin/solar 1 ./tests/1_MAXNRG_H1/x0.txt` should display

`-122505.5978 -10881140.57 -1512631.39776 -134 -4.5 0`

which corresponds to the feasible point
*(8, 8, 150, 7, 7, 250, 45, 0.5, 5)*
of value *-122,505.5978*.

Other points and NOMAD parameters files can be found in the
[./tests](tests) directory.

It is also possible to modify the `main()` function in [./src/main.cpp](main.cpp) to call SOLAR from a code. A minimal example is provided.

### How to cite

```
@mastersthesis{MScMLG,
  Author = {Mathieu Lemyre Garneau},
  Title  = {{Modelling of a solar thermal power plant for benchmarking blackbox optimization solvers}},
  School = {Polytechnique Montr\'eal},
  Year   = {2015},
  Url    = {https://publications.polymtl.ca/1996/},
  Note   = {Text available at \url{https://publications.polymtl.ca/1996},
            code available at \url{https://github.com/bbopt/solar}}
}
```
(article to come)
