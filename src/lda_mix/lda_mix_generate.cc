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
  double n;
  Vdouble pi;
  VVdouble kappa;

  param (double _n, Vdouble _pi, VVdouble _kappa) : n(_n), pi(_pi), kappa(_kappa) {}
};

typedef vector<param*> params;

typedef struct _hyperparams
{
  Vdouble gamma;
  Vdouble alpha;
  VVdouble beta;
  int MAX_SUBTYPE;
  int TOTAL_CN;
}
  hyperparams;

void write_params(std::ofstream& f, params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    f << pa[i]->n << "\t";
  f << endl << endl;
  
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
    pa.push_back(new param (0, Vdouble (hpa.TOTAL_CN + 1, 0), VVdouble (hpa.TOTAL_CN + 1, Vdouble (hpa.TOTAL_CN + 1, 0))));
}

void delete_params(params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    delete pa[i];
}

void init_hyperparams(hyperparams& hpa)
{
  hpa.gamma.assign(hpa.MAX_SUBTYPE + 1, 0.1);
  hpa.alpha.assign(hpa.TOTAL_CN + 1, 0.1);
  hpa.beta.assign(hpa.TOTAL_CN + 1, Vdouble (hpa.TOTAL_CN + 1, 0.1));
}

double calc_mu(state& st, params& pa, hyperparams& hpa)
{
  double denom = 0;
  double num = 0;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      denom += pa[i]->n * st.total_cn[i];
      num += pa[i]->n * st.variant_cn[i];
    }
    
  return num / denom;
}

void generate_params(params& pa, hyperparams& hpa, gsl_rng* r)
{
  Vdouble v (hpa.MAX_SUBTYPE + 1, 0);
  for (int i=0; i<1024; ++i) // for appropriet random number generation
    gsl_ran_dirichlet(r, hpa.MAX_SUBTYPE + 1, &hpa.gamma[0], &v[0]);
  gsl_ran_dirichlet(r, hpa.MAX_SUBTYPE + 1, &hpa.gamma[0], &v[0]);
  
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    pa[i]->n = v[i];
  
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
  for (int i=0; i<2048; ++i)
    genrand64_real2();
  
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

      double mu = calc_mu(st, pa, hpa);
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
      cerr << "usage: ./lda_mix_generate max_subtype total_cn M n seed (n_pi_kappa outfile) (reads outfile)" << endl;
      exit(EXIT_FAILURE);
    }
  
  hyperparams hpa;
  int M, n, seed;
  
  hpa.MAX_SUBTYPE = atoi(argv[1]);
  hpa.TOTAL_CN = atoi(argv[2]);
  M = atoi(argv[3]);
  n = atoi(argv[4]);
  seed = atoi(argv[5]);

  ofstream f (argv[6]);
  ofstream g (argv[7]);
  f << scientific;
  
  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  // T = gsl_rng_default; // default is mt19937
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);

  init_hyperparams(hpa);
  
  params pa;
  init_params(pa, hpa);
  
  generate_params(pa, hpa, r);
  write_params(f, pa, hpa);
  
  generate_binom(g, M, n, pa, hpa, seed, r);

  delete_params(pa, hpa);
  
  f.close();
  g.close();
  
  return 0;
}
