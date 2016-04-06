#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "../mt19937-64/mt64.h"
using namespace std;

typedef vector<double> Vdouble;
typedef vector<int> Vint;
typedef vector<Vdouble> VVdouble;

class state 
{
public:
  Vint total_cn;
  Vint variant_cn;

  // use default constructor
  // state (int, std::vector<int>, std::vector<int>, double);
};

typedef std::vector<state*> states;

class param
{
public:
  Vdouble pi;
  VVdouble kappa;

  param (Vdouble _pi, VVdouble _kappa) : pi(_pi), kappa(_kappa) {}
};

typedef vector<param*> params;

typedef struct _hyperparams
{
  Vdouble alpha;
  VVdouble beta;
  int MAX_SUBTYPE;
  int TOTAL_CN;
}
  hyperparams;

void write_params(std::ofstream& f, params& pa, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double sum = 0;
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        f << pa[i]->pi[l] << "\t";
      f << endl;

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (int r=1; r<=l; ++r)
            f << pa[i]->kappa[l][r] << "\t";
          f << endl;
        }
      f << endl;
    }
}

void init_state(state& st, hyperparams& hpa)
{
  st.total_cn.assign(hpa.MAX_SUBTYPE + 1, 0);
  st.variant_cn.assign(hpa.MAX_SUBTYPE + 1, 0);
}

void init_params(params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    pa.push_back(new param (Vdouble (hpa.TOTAL_CN + 1, 0), VVdouble (hpa.TOTAL_CN + 1, Vdouble (hpa.TOTAL_CN + 1, 0))));
}

void delete_params(params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    delete pa[i];
}

void init_hyperparams(hyperparams& hpa)
{
  hpa.alpha.assign(hpa.TOTAL_CN + 1, 0.05);
  hpa.beta.assign(hpa.TOTAL_CN + 1, Vdouble (hpa.TOTAL_CN + 1, 0.05));
}

double calc_mu(state& st, hyperparams& hpa)
{
  int denom = 0;
  int num = 0;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      denom += st.total_cn[i];
      num += st.variant_cn[i];
    }
    
  return ((double) num) / denom;
}

void generate_params(params& pa, hyperparams& hpa, gsl_rng* r)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      gsl_ran_dirichlet(r, hpa.TOTAL_CN, &hpa.alpha[1], &pa[i]->pi[1]);

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        gsl_ran_dirichlet(r, l, &hpa.beta[l][1], &pa[i]->kappa[l][1]);
    }
}

void generate_binom(ofstream& f, int M, int n, params& pa, hyperparams& hpa, int seed, gsl_rng* r)
{
  params pa_cum;
  init_params(pa_cum, hpa);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          pa_cum[i]->pi[l] = pa[i]->pi[l];
        }

      for (int l=2; l<=hpa.TOTAL_CN; ++l)
        {
          pa_cum[i]->pi[l] = pa_cum[i]->pi[l-1] + pa_cum[i]->pi[l];
        }

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (int r=1; r<=l; ++r)
            {
              pa_cum[i]->kappa[l][r] = pa[i]->kappa[l][r];
            }

          for (int r=2; r<=l; ++r)
            {
              pa_cum[i]->kappa[l][r] = pa_cum[i]->kappa[l][r-1] + pa_cum[i]->kappa[l][r];
            }
        }
    }

  init_genrand64(seed);

  state st;
  init_state(st, hpa);

  for (int k=0; k<n; ++k)
    {
      st.total_cn[0] = 2;
      st.variant_cn[0] = 0;
      
      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        {
          double x = genrand64_real2();
          for (int l=1; l<=hpa.TOTAL_CN; ++l)
            {
              if (x < pa_cum[i]->pi[l])
                {
                  st.total_cn[i] = l;
                  break;
                }
            }
          
          double y = genrand64_real2();
          for (int r=1; r<=st.total_cn[i]; ++r)
            {
              if (y < pa_cum[i]->kappa[st.total_cn[i]][r])
                {
                  st.variant_cn[i] = r;
                  break;
                }
            }
        }

      double mu = calc_mu(st, hpa);
      // cerr << "mu: " << mu << endl;
      
      unsigned int m = gsl_ran_binomial(r, mu, M);
      f << m << "\t" << M << endl;
    }
  
  gsl_rng_free (r);
}

int main(int argc, char** argv)
{
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  cerr << scientific;
  
  if (argc != 8)
    {
      cerr << "usage: ./lda_generate max_subtype total_cn M n seed (pi_kappa infile) (reads outfile)" << endl;
      exit(EXIT_FAILURE);
    }
  
  hyperparams hpa;
  int M, n, seed;
  
  hpa.MAX_SUBTYPE = atoi(argv[1]);
  hpa.TOTAL_CN = atoi(argv[2]);
  M = atoi(argv[3]);
  n = atoi(argv[4]);
  seed = atoi(argv[5]);

  ifstream f (argv[6]);

  params pa;
  init_params(pa, hpa);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          f >> pa[i]->pi[l];
        }

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (int r=1; r<=l; ++r)
            {
              f >> pa[i]->kappa[l][r];
            }
        }
    }

  ofstream g (argv[7]);
  
  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  generate_binom(g, M, n, pa, hpa, seed, r);

  delete_params(pa, hpa);
  
  f.close();
  g.close();
  
  return 0;
}
