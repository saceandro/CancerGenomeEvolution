#include <iostream>
#include <fstream>
#include <math.h>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
using namespace std;

void generate_binom(ofstream& f, double mu, unsigned int M, unsigned int n)
{
  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  for (int i=0; i<n; ++i)
    {
      unsigned int m = gsl_ran_binomial(r, mu, M);
      f << m << "\t" << M << endl;
    }
  
  gsl_rng_free (r);
}

int main(int argc, char** argv)
{
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 5)
    {
      cerr << "usage: ./binomial_generate (outfile) mu M n" << endl;
      exit(EXIT_FAILURE);
    }
  
  ofstream f (argv[1]);
  double mu;
  unsigned int M, n;
  mu = atof(argv[2]);
  M = atoi(argv[3]);
  n = atoi(argv[4]);
  
  generate_binom(f, mu, M, n);
  f.close();
  
  return 0;
}
