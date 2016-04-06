#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "./mt19937-64/mt64.h"
using namespace std;

typedef vector<double> Vdouble;
typedef vector<unsigned int> Vuint;

double calc_mu(Vuint& a, unsigned int total_cn, unsigned int max_subtype)
{
  unsigned int sum = 0;
  for (Vuint::iterator it = a.begin(); it != a.end(); ++it)
    sum += *it;
  return ((double) sum) / (max_subtype + 1) / total_cn;
}

void generate_binom(ofstream& f, Vdouble& kappa, unsigned int M, unsigned int n, unsigned int total_cn, unsigned int max_subtype)
{
  Vdouble kappa_cum (kappa.size(), 0);
  
  for (int r=1; r<=total_cn; ++r)
    {
      kappa_cum[r] = kappa[r];
    }

  //cerr << kappa_cum[0] << endl;
  //cerr << kappa_cum[1] << endl;
  for (int r=2; r<=total_cn; ++r)
    {
      kappa_cum[r] = kappa_cum[r-1] + kappa_cum[r];
      //cerr << kappa_cum[r] << endl;
    }

  init_genrand64(1);

  Vuint cns (max_subtype+1, 0);
  
  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  for (int i=0; i<n; ++i)
    {
      for (int i=0; i<=max_subtype; ++i)
        {
          double x = genrand64_real2();
          //cerr << x << endl;
          for (int r=1; r<=total_cn; ++r)
            {
              if (x < kappa_cum[r])
                {
                  cns[i] = r;
                  //cerr << cns[i] << " ";
                  break;
                }
            }
          //cerr << endl;
        }

      double mu = calc_mu(cns, total_cn, max_subtype);
      //cerr << "mu: " << mu << endl;
      
      unsigned int m = gsl_ran_binomial(r, mu, M);
      f << m << "\t" << M << endl;
    }
  
  gsl_rng_free (r);
}

int main(int argc, char** argv)
{
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 6)
    {
      //cerr << "usage: ./mixture_cn_generate (outfile) M n max_subtype (infile)" << endl;
      exit(EXIT_FAILURE);
    }
  
  ofstream f (argv[1]);
  unsigned int M, n, total_cn, max_subtype;
  Vdouble kappa;
  
  M = atoi(argv[2]);
  n = atoi(argv[3]);
  max_subtype = atoi(argv[4]);
  ifstream g (argv[5]);

  g >> total_cn;
  
  kappa.assign(total_cn, 0);

  for (int r=1; r<=total_cn; ++r)
    g >> kappa[r];
  
  generate_binom(f, kappa, M, n, total_cn, max_subtype);
  f.close();
  
  return 0;
}
