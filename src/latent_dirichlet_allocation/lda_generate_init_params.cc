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
  hpa.alpha.assign(hpa.TOTAL_CN + 1, 1.0);
  hpa.beta.assign(hpa.TOTAL_CN + 1, Vdouble (hpa.TOTAL_CN + 1, 1.0));
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

int main(int argc, char** argv)
{
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  cerr << scientific;
  
  if (argc != 4)
    {
      cerr << "usage: ./lda_generate_init_params max_subtype total_cn (pi_kappa outfile)" << endl;
      exit(EXIT_FAILURE);
    }
  
  hyperparams hpa;
  
  hpa.MAX_SUBTYPE = atoi(argv[1]);
  hpa.TOTAL_CN = atoi(argv[2]);

  ofstream f (argv[3]);
  f << scientific;
  
  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, time(NULL));

  init_hyperparams(hpa);
  
  params pa;
  init_params(pa, hpa);
  
  generate_params(pa, hpa, r);
  write_params(f, pa, hpa);
  
  delete_params(pa, hpa);
  
  f.close();

  return 0;
}
