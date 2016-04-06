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

class state 
{
public:
  std::vector<int> total_cn;
  std::vector<int> variant_cn;

  // use default constructor
  // state (int, std::vector<int>, std::vector<int>, double);
};

typedef std::vector<state*> states;

class params
{
public:
  std::vector<double> pi;
  std::vector<std::vector<double> > kappa;

  // use default constructor
  // params (std::vector<double>, std::vector<std::vector<double> >);
};

typedef struct _hyperparams
{
  unsigned int MAX_SUBTYPE;
  unsigned int TOTAL_CN;
}
  hyperparams;

void write_params(std::ofstream& f, params& pa)
{
  f << endl << "pi" << endl;
  for (int c=1; c<pa.pi.size(); ++c)
    f << pa.pi[c] << "\t";
  f << endl;

  f << endl << "kappa" << endl;
  for (int c=1; c<pa.pi.size(); ++c)
    {
      for (int d=1; d<=c; ++d)
        f << pa.kappa[c][d] << "\t";
      f << endl;
    }
  return;
}

void init_state(state& st, hyperparams& hpa)
{
  st.total_cn.assign(hpa.MAX_SUBTYPE + 1, 0);
  st.variant_cn.assign(hpa.MAX_SUBTYPE + 1, 0);
}

void init_params(params& pa, hyperparams& hpa)
{
  pa.pi.assign(hpa.TOTAL_CN + 1, 0);
  pa.kappa.assign(hpa.TOTAL_CN + 1, Vdouble (hpa.TOTAL_CN + 1, 0));
}

double calc_mu(state& st, hyperparams& hpa)
{
  unsigned int denom = 0;
  unsigned int num = 0;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      denom += st.total_cn[i];
      num += st.variant_cn[i];
    }
    
  return ((double) num) / denom;
}

void generate_binom(ofstream& f, params& pa, unsigned int M, unsigned int n, hyperparams hpa)
{
  params pa_cum;
  init_params(pa_cum, hpa);
  
  for (int l=1; l<=hpa.TOTAL_CN; ++l)
    {
      pa_cum.pi[l] = pa.pi[l];
    }

  for (int l=2; l<=hpa.TOTAL_CN; ++l)
    {
      pa_cum.pi[l] = pa_cum.pi[l-1] + pa_cum.pi[l];
    }

  for (int l=1; l<=hpa.TOTAL_CN; ++l)
    {
      for (int r=1; r<=l; ++r)
        {
          pa_cum.kappa[l][r] = pa.kappa[l][r];
        }

      for (int r=2; r<=l; ++r)
        {
          pa_cum.kappa[l][r] = pa_cum.kappa[l][r-1] + pa_cum.kappa[l][r];
        }
    }

  init_genrand64(1);

  state st;
  init_state(st, hpa);
  
  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  for (int i=0; i<n; ++i)
    {
      for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
        {
          double x = genrand64_real2();
          for (int l=1; l<=hpa.TOTAL_CN; ++l)
            {
              if (x < pa_cum.pi[l])
                {
                  st.total_cn[i] = l;
                  break;
                }
            }
          
          double y = genrand64_real2();
          for (int r=1; r<=st.total_cn[i]; ++r)
            {
              if (y < pa_cum.kappa[st.total_cn[i]][r])
                {
                  st.variant_cn[i] = r;
                  break;
                }
            }
        }

      double mu = calc_mu(st, hpa);
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
      cerr << "usage: ./mixture_total_variant_generate (outfile) M n max_subtype (infile)" << endl;
      exit(EXIT_FAILURE);
    }
  
  ofstream f (argv[1]);
  unsigned int M, n;
  hyperparams hpa;
  params pa;
  
  M = atoi(argv[2]);
  n = atoi(argv[3]);
  hpa.MAX_SUBTYPE = atoi(argv[4]);
  ifstream g (argv[5]);

  g >> hpa.TOTAL_CN;

  init_params(pa, hpa);

  for (int l=1; l<=hpa.TOTAL_CN; ++l)
    g >> pa.pi[l];

  for (int l=1; l<=hpa.TOTAL_CN; ++l)
    for (int r=1; r<=l; ++r)
      g >> pa.kappa[l][r];
  
  generate_binom(f, pa, M, n, hpa);
  f.close();
  
  return 0;
}
