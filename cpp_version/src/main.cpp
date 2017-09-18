# define MAX_POP_SIZE 10000
# define VERBOSE false
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <ctime>
# include "rf.h"
#include <algorithm> 
using namespace std;

//
//  Each GENOTYPE is a member of the population, with
//  gene: a string of variables,
//  fitness: the fitness
//  upper: the variable upper bounds,
//  lower: the variable lower bounds,
//  rfitness: the relative fitness,
//  cfitness: the cumulative fitness.
//

std::vector<genotype> population, newpopulation, bestpopulation;
Forest* forest;
bool reload_data;
size_t pop_size, max_gens, n_vars, n_features, best_gen;
size_t tournament_size;
double lambda, averageFitness, bestAverageFitness, overallPredictionError, overallCorrelation;
std::vector<int> pop_series;
std::ostream *verbose_out, *trees_out;
char filename[20];
std::string base_path = "D:\\ERF_Project\\logs\\";
std::string input_file_path;
std::string case_weight_file_path;
int main (int argc, char* argv[]);
void crossover ( int &seed );
void elitist ( );
void saiyajin( size_t gen );
bool evaluate ( int &seed );
int i4_uniform_ab ( int a, int b, int &seed );
void initialize ( ifstream& input, int &seed );
void keep_the_best ( );
void mutate ( int &seed );
double r8_uniform_ab ( double a, double b, int &seed );
void report ( int generation );
void selector ( int &seed );
void timestamp ( );
void Xover ( int one, int two, int &seed );

//****************************************************************************80

int main (int argc, char* argv[])

//****************************************************************************80
//
//  Purpose:
//
//    MAIN supervises the genetic algorithm.
//
//  Discussion:
//
//    Each generation involves selecting the best 
//    members, performing crossover & mutation and then 
//    evaluating the resulting population, until the terminating 
//    condition is satisfied   
//
//    This is a simple genetic algorithm implementation where the 
//    evaluation function takes positive values only and the 
//    fitness of an individual is the same as the value of the 
//    objective function.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Zbigniew Michalewicz,
//    Genetic Algorithms + Data Structures = Evolution Programs,
//    Third Edition,
//    Springer, 1996,
//    ISBN: 3-540-60676-9,
//    LC: QA76.618.M53.
//
//  Parameters:
//
//    max_gens is the maximum number of generations.
//
//    n_vars is the number of problem variables.
//
//    PMUTATION is the probability of mutation.
//
//    pop_size is the population size. 
//
//    PXOVER is the probability of crossover.                          
//
{
	if (argc > 1) {
		std::vector<genotype> NullPointer;
		rf(argc, argv, NullPointer, true);
	}
	else{
	  string config_filename = "C:\\Users\\i-fen\\Documents\\ERF_Project\\forest_init_config.txt";
	  uint generation;
	  int seed;
	  double elapsed_secs;
	  ifstream input;
	  clock_t begin, end, time_check1, time_check2;
	  seed = 12345678;
	  input.open(config_filename.c_str());
	  reload_data = true;
	  while (!input.eof()) {
		  begin = clock();
		  initialize(input, seed);
		  timestamp();
		  evaluate(seed);
		  keep_the_best();
		  report(0);
		  saiyajin(0);
		  reload_data = false;
		  cout << "breed generations";
		  time_check1 = clock();
		  *verbose_out << endl << "Initilization time: " << double(time_check1 - begin) / CLOCKS_PER_SEC << " seconds" << endl;
		  for (generation = 0; generation < max_gens; ++generation) {
			  selector(seed);
			  crossover(seed);
			  mutate(seed);
			  if (evaluate(seed)) {
				  if (generation == 0) {
					  time_check2 = clock();
					  *verbose_out << "Evaluation time: " << double(time_check2 - time_check1) / CLOCKS_PER_SEC << "second" << endl << endl;
				  }
				  report(generation + 1);
				  elitist();
				  saiyajin(generation + 1);
				  cout << ".";
			  }
			  else {
				  break;
			  }
		  }
		  end = clock();
		  elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		  population = bestpopulation;
		  evaluate(seed);
		  population.clear();
		  newpopulation.clear();
		  bestpopulation.clear();
		  reload_data = true;
		  *verbose_out << "Best average fitness: " << bestAverageFitness << " at generation " << best_gen << endl;
		  *verbose_out << "Forest and confusion for the best population are generated." << endl;
		  *verbose_out << "Total elapsed time: " << elapsed_secs << " seconds" << endl << endl;
		  cout << endl;
	  }
	  //
	  //  Terminate.
	  //
	  timestamp();
	  input.close();
	  cout << "finished." << endl << "Have a good day!" << endl;
  }
  return 0;
}
//****************************************************************************80

void crossover ( int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    CROSSOVER selects two parents for the single point crossover.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int FIRST, is a count of the number of members chosen.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  const double a = 0.0;
  const double b = 1.0;
  uint mem = 0;
  int one = 0;
  int first = 0;
  double x;

  for ( mem = 0; mem < pop_size; ++mem )
  {
    x = r8_uniform_ab ( a, b, seed );

    if ( x < PXOVER )
    {
      ++first;

      if ( first % 2 == 0 )
      {
        Xover ( one, mem, seed );
      }
      else
      {
        one = mem;
      }

    }
  }
  return;
}
//****************************************************************************80

void elitist ( )

//****************************************************************************80
// 
//  Purpose:
//
//    ELITIST stores the best member of the previous generation.
//
//  Discussion:
//
//    The best member of the previous generation is stored as 
//    the last in the array. If the best member of the current 
//    generation is worse then the best member of the previous 
//    generation, the latter one would replace the worst member 
//    of the current population.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double BEST, the best fitness value.
//
//    Local, double WORST, the worst fitness value.
//
{
  uint i, best_mem = 0, worst_mem = 0;
  double best, worst;

  best = population[0].fitness;
  worst = population[0].fitness;

  for ( i = 0; i < pop_size - 1; ++i )
  {
    if ( population[i+1].fitness < population[i].fitness )
    {

      if ( best <= population[i].fitness )
      {
        best = population[i].fitness;
        best_mem = i;
      }

      if ( population[i+1].fitness <= worst )
      {
        worst = population[i+1].fitness;
        worst_mem = i + 1;
      }

    }
    else
    {

      if ( population[i].fitness <= worst )
      {
        worst = population[i].fitness;
        worst_mem = i;
      }

      if ( best <= population[i+1].fitness )
      {
        best = population[i+1].fitness;
        best_mem = i + 1;
      }

    }

  }
// 
//  If the best individual from the new population is better than 
//  the best individual from the previous population, then 
//  copy the best from the new population; else replace the 
//  worst individual from the current population with the 
//  best one from the previous generation                     
//
  if ( population[pop_size].fitness <= best )
  {
    for ( i = 0; i < n_vars; i++ )
    {
      population[pop_size].gene[i] = population[best_mem].gene[i];
    }
    population[pop_size].fitness = population[best_mem].fitness;
  }
  else
  {
    for ( i = 0; i < n_vars; i++ )
    {
      population[worst_mem].gene[i] = population[pop_size].gene[i];
    }
    population[worst_mem].fitness = population[pop_size].fitness;
  } 

  return;
}
void saiyajin(size_t gen)

//****************************************************************************80
// 
//  Purpose:
//
//    SAIYAJIN stores the best population through evolution
//
//  Discussion:
//
//    If the current population has a lower prediction error than the previous
//    one, then store it as the best generation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double BEST, the best fitness value.
//
//    Local, double WORST, the worst fitness value.
//
{
	if (averageFitness > bestAverageFitness) {
		for (size_t i = 0; i < pop_size; ++i) {
			bestpopulation[i] = population[i];
		}
		best_gen = gen;
		bestAverageFitness = averageFitness;
	}

	return;
}
//****************************************************************************80

bool evaluate ( int& seed)

//****************************************************************************80
// 
//  Purpose:
//
//    EVALUATE implements the user-defined valuation function
//
//  Discussion:
//
//    Each time this is changed, the code has to be recompiled.
//    The current function is:  x[1]^2-x[1]*x[2]+x[3]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
{
	char pop_size_char[4];
	char file_path_char[512];
	char case_weight_file_path_char[512];
	char seed_char[8];
	char output_file_path[512];
	sprintf(pop_size_char, "%d", pop_size);
	sprintf(seed_char, "%d", seed);
	strcpy(file_path_char, input_file_path.c_str());
	strcpy(case_weight_file_path_char, case_weight_file_path.c_str());
	strcpy(output_file_path, base_path.c_str());
	strcat(output_file_path, filename);
	char * args[] = { "ranger", // dummy argument
		//"--verbose",
		"--file",
		file_path_char,
		"--depvarname",
		"\"label\"",
		"--treetype",
		"1",
		"--ntree",
		pop_size_char,
		"--seed",
		seed_char,
		"--caseweights",
		case_weight_file_path_char,
		"--write",
		"--outprefix",
		output_file_path,
		"--nthreads",
	    "8"};
	int argc = sizeof(args) / sizeof(*args);

	forest = rf(argc, args, population, reload_data);
	// check if forest is successfully created
	double prediction_error, max_correlation;
	if (forest) {
		for (size_t member = 0; member < pop_size; ++member) {
			prediction_error = forest->getPredictionErrorOfTree(member);
			max_correlation = forest->getMaxCorrelation(member);
			population[member].fitness = max(0.0, 1 - prediction_error - (lambda * max_correlation));
		}
		overallPredictionError = forest->getOverallPredictionError();
		overallCorrelation = forest->getOverallCorrelation();
		delete forest;
	}
	else {
		std::cout <<  "*" ;
		*verbose_out << "Something is wrong. Probably insufficient memory." << endl;
		return false;
	}
	return true;
	
}
//****************************************************************************80

int i4_uniform_ab ( int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}
//****************************************************************************80

void initialize ( ifstream& input, int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    INITIALIZE initializes the genes within the variables bounds. 
//
//  Discussion:
//
//    It also initializes (to zero) all fitness values for each
//    member of the population. It reads upper and lower bounds 
//    of each variable from the input file `gadata.txt'. It 
//    randomly generates values between these bounds for each 
//    gene of each genotype in the population. The format of 
//    the input file `gadata.txt' is 
//
//      var1_lower_bound var1_upper bound
//      var2_lower_bound var2_upper bound ...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the input file.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  double lbound;
  double ubound;


  if ( !input )
  {
    cerr << "\n";
    cerr << "INITIALIZE - Fatal error!\n";
    cerr << "  Cannot open the input file!\n";
    exit ( 1 );
  }

  // Verbose output to logfile if non-verbose mode
  if (VERBOSE) {
	  verbose_out = &std::cout;
  }
  else {
	  std::ofstream* logfile = new std::ofstream();
	  std::ofstream* treesfile = new std::ofstream();
	  char report_filename[27];
	  char tree_filename[26];
	  std::string file_path, tree_file_path;
	  time_t t = time(0);
	  strftime(filename, sizeof(filename), "%Y%m%d%H%M%S", gmtime(&t));
	  strcpy(report_filename, filename);
	  strcpy(tree_filename, filename);
	  strcat(report_filename, "_report.txt");
	  strcat(tree_filename, "_trees.txt");
	  file_path = base_path + std::string(report_filename);
	  tree_file_path = base_path + std::string(tree_filename);
	  logfile->open(file_path);
	  treesfile->open(tree_file_path);
	  if (!logfile->good()) {
		  throw std::runtime_error("Could not write to logfile.");
	  }
	  verbose_out = logfile;
	  trees_out = treesfile;
	  cout << "Writing output to " << file_path << endl;
  }

//  Initialize hyper parameters
  input >> input_file_path >> case_weight_file_path;
  input >> n_vars >> n_features >> max_gens >> pop_size >> lambda;
  *verbose_out << "Input file: " << input_file_path << endl;
  *verbose_out << "Max generations: " << max_gens << endl;
  *verbose_out << "Population size: " << pop_size << endl;
  *verbose_out << "Lambda         : " << lambda << endl;
  *verbose_out << "N. of variables: " << n_vars << endl;
  *verbose_out << "N. of features : " << n_features << endl;
  *verbose_out << "Parameter ranges: " << endl;
  if (!VERBOSE) { // print message if verbose mode for debugging
	  cout << "Max generations: " << max_gens << endl;
	  cout << "Population size: " << pop_size << endl;
	  cout << "Lambda         : " << lambda << endl;
	  cout << "N. of variables: " << n_vars << endl;
	  cout << "N. of features : " << n_features << endl;
  }
  if (n_vars < 2)
  {
	  *verbose_out << "\n";
	  *verbose_out << "  The crossover modification will not be available,\n";
	  *verbose_out << "  since it requires 2 <= n_vars.\n";
  }

//  Initialize population
//  +1 to store the best individual
  population = std::vector<genotype>(pop_size + 1, genotype(n_vars));
  newpopulation = std::vector<genotype>(pop_size + 1, genotype(n_vars));
  bestpopulation = std::vector<genotype>(pop_size, genotype(n_vars));
// 
//  Initialize variables within the bounds 
//
  for ( size_t i = 0; i < n_vars; ++i )
  {
	  if (i < n_features) {
		  input >> lbound >> ubound;
		  *verbose_out << "   " << i + 1 << ": [" << lbound << ", " << ubound << "]" << endl;
	  }
	  else {
		  // bounds for select weights
		  lbound = 1;
		  ubound = 999;
	  }
      for ( size_t j = 0; j < pop_size; j++ ) {
		  population[j].lower[i] = lbound;
		  population[j].upper[i] = ubound;
          population[j].gene[i] = uint(r8_uniform_ab ( lbound, ubound, seed ));
      }
  }

  // initilize pop series
  for ( size_t i = 0; i < pop_size; ++i) {
	  pop_series.push_back(i);
  }
  // set tournament size
  tournament_size = int(ceil(pow(pop_size, 0.25)));
  // 
  bestAverageFitness = 0;
}
//****************************************************************************80

void keep_the_best ( )

//****************************************************************************80
// 
//  Purpose:
//
//    KEEP_THE_BEST keeps track of the best member of the population. 
//
//  Discussion:
//
//    Note that the last entry in the array Population holds a 
//    copy of the best individual.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int CUR_BEST, the index of the best individual.
//
{
  int cur_best;
  uint mem, i;

  cur_best = 0;

  for ( mem = 0; mem < pop_size; mem++ )
  {
    if ( population[pop_size].fitness < population[mem].fitness )
    {
      cur_best = mem;
      population[pop_size].fitness = population[mem].fitness;
    }
  }
// 
//  Once the best member in the population is found, copy the genes.
//
  for ( i = 0; i < n_vars; i++ )
  {
    population[pop_size].gene[i] = population[cur_best].gene[i];
  }

  return;
}
//****************************************************************************80

void mutate ( int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    MUTATE performs a random uniform mutation. 
//
//  Discussion:
//
//    A variable selected for mutation is replaced by a random value 
//    between the lower and upper bounds of this variable.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  const double a = 0.0;
  const double b = 1.0;
  uint i, j;
  double lbound;
  double ubound;
  double x;

  for ( i = 0; i < pop_size; i++ )
  {
    for ( j = 0; j < n_vars; j++ )
    {
      x = r8_uniform_ab ( a, b, seed );
      if ( x < PMUTATION )
      {
        lbound = population[i].lower[j];
        ubound = population[i].upper[j];  
        population[i].gene[j] = r8_uniform_ab ( lbound, ubound, seed );
      }
    }
  }

  return;
}
//****************************************************************************80

double r8_uniform_ab ( double a, double b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_AB returns a scaled pseudorandom R8.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_AB, a number strictly between A and B.
//
{
  int i4_huge = 2147483647;
  int k;
  double value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  value = ( double ) ( seed ) * 4.656612875E-10;

  value = a + ( b - a ) * value;

  return value;
}
//****************************************************************************80

void report ( int generation )

//****************************************************************************80
// 
//  Purpose:
//
//    REPORT reports progress of the simulation. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double avg, the average population fitness.
//
//    Local, best_val, the best population fitness.
//
//    Local, double square_sum, square of sum for std calc.
//
//    Local, double stddev, standard deviation of population fitness.
//
//    Local, double sum, the total population fitness.
//
//    Local, double sum_square, sum of squares for std calc.
//
{
  uint i;
  double avg;
  double best_val;
  double square_sum;
  double stddev;
  double sum;
  double sum_square;

 /* *verbose_out << "Fitness of each tree:" << endl;
  for (size_t i = 0; i < pop_size; ++i)
  {
	  *verbose_out << i << "(";
	  for (size_t j = 0; j < n_vars; ++j) {
		  *verbose_out << population[i].gene[j] << ',';
	  }
	  *verbose_out << ") " << population[i].fitness << endl;
  }*/

  if ( generation == 0 )
  {
    *verbose_out << "\n";
    *verbose_out << "  Generation       Best            Average       Standard         Overall         Overall\n";
    *verbose_out << "  number           value           fitness       deviation       OOB Error      Correlation\n";
    *verbose_out << "\n";
  }

  sum = 0.0;
  sum_square = 0.0;

  for ( i = 0; i < pop_size; i++ )
  {
    sum = sum + population[i].fitness;
    sum_square = sum_square + population[i].fitness * population[i].fitness;
  }

  avg = sum / ( double ) pop_size;
  square_sum = avg * avg * pop_size;
  stddev = sqrt ( ( sum_square - square_sum ) / ( pop_size - 1 ) );
  best_val = population[pop_size].fitness;

  *verbose_out << "  " << setw(8) << generation 
       << "  " << setw(14) << best_val 
       << "  " << setw(14) << avg 
	   << "  " << setw(14) << stddev
	   << "  " << setw(14) << overallPredictionError
	   << "  " << setw(14) << overallCorrelation
	   << "\n";

  // print trees
  *trees_out << "%" << generation << "th generation" << endl;
  for (size_t i = 0; i < pop_size; i++) {
	  for (size_t j = 0; j < n_vars; j++) {
		  *trees_out << population[i].gene[j];
		  *trees_out << ' ';
	  }
	  *trees_out << population[i].fitness;
	  *trees_out << endl;
  }
  verbose_out->flush();
  averageFitness = avg;
  return;
}
//****************************************************************************80

void selector ( int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    SELECTOR is the selection function.
//
//  Discussion:
//
//    Standard proportional selection for maximization problems incorporating 
//    the elitist model.  This makes sure that the best member always survives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  //const double a = 0.0;
  //const double b = 1.0;
  uint i, k;
  //int j;
  uint picked, winner = 0;
  uint mem;
  //double p;
  double sum;
//
//  Find the total fitness of the population.
//
  sum = 0.0;
  for ( mem = 0; mem < pop_size; mem++ )
  {
    sum = sum + population[mem].fitness;
  }
//
//  Calculate the relative fitness of each member.
//
  //for ( mem = 0; mem < pop_size; mem++ )
  //{
  //  population[mem].rfitness = population[mem].fitness / sum;
  //}

  // tournament selection
  std::srand(seed);
  for (mem = 0; mem < pop_size; mem++) {
	  std::random_shuffle(pop_series.begin(), pop_series.end());
	  for (k = 0; k < tournament_size; k++) {
		  picked = pop_series.at(k);
	      if (k == 0 || population[picked].fitness > population[winner].fitness) {
			  winner = picked;
		  }
	  }
	  newpopulation[mem] = population[winner];
  }
  

// 
////  Calculate the cumulative fitness.
////
//  population[0].cfitness = population[0].rfitness;
//  for ( mem = 1; mem < pop_size; mem++ )
//  {
//    population[mem].cfitness = population[mem-1].cfitness +       
//      population[mem].rfitness;
//  }
//// 
////  Select survivors using cumulative fitness. 
////
//  for ( i = 0; i < pop_size; i++ )
//  { 
//    p = r8_uniform_ab ( a, b, seed );
//    if ( p < population[0].cfitness )
//    {
//      newpopulation[i] = population[0];      
//    }
//    else
//    {
//      for ( j = 0; j < pop_size; j++ )
//      { 
//        if ( population[j].cfitness <= p && p < population[j+1].cfitness )
//        {
//          newpopulation[i] = population[j+1];
//		    break;
//        }
//      }
//    }
//  }
// 

//  Overwrite the old population with the new one.
//
  for ( i = 0; i < pop_size; i++ )
  {
    population[i] = newpopulation[i]; 
  }

  return;     
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2003
//
//  Author:
//
//    John Burkardt
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  *verbose_out << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

void Xover ( int one, int two, int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    XOVER performs crossover of the two selected parents. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int point, the crossover point.
//
//  Parameters:
//
//    Input, int ONE, TWO, the indices of the two parents.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  size_t i;
  double t;
// 
//  Select the crossover point.
//
//  point = i4_uniform_ab ( 0, n_vars - 1, seed );
//
//  Swap genes in positions 0 through POINT-1.
//
//  for ( i = 0; i < point; i++ )
//  {
//    t                       = population[one].gene[i];
//    population[one].gene[i] = population[two].gene[i];
//    population[two].gene[i] = t;
//  }

  // uniform crossover
  const double xover_rate = 0.5;
  for (i = 0; i < n_vars; i++)
  {
	  if (r8_uniform_ab(0.0, 1.0, seed) < xover_rate) {
		  t = population[one].gene[i];
		  population[one].gene[i] = population[two].gene[i];
		  population[two].gene[i] = t;
	  }
  }
  
  return;
}