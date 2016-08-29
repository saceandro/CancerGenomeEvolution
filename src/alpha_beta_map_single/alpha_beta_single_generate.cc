#include "setting.hh"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <xmmintrin.h>
#include "../../util/enumtree_wf.hh"
using namespace std;

extern void variant_fraction_partition(int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, VLog& vf, VLog& vf_numerator, VLog& vf_denominator, Log& partition);
extern void set_gegen(VVLog &gegen);
extern void set_gegen_integral(VLog &gegen_int, VLog &gegen_int_err);

void write_params(std::ofstream& f, params& pa, hyperparams& hpa)
{
  f << pa.pa[1]->u.eval() << endl;
}

void write_t_n(std::ofstream& f, subtype& st, hyperparams& hpa)
{
  f << st.t.eval() << endl;

  f << st.n.eval() << endl;
}

double calc_mu(subtype& st, hyperparams& hpa)
{
    
  return (st.x / Log(2)).eval();
}

bool strictly_greater_than(double i, double j)
{
  return (i > j);
}

int strictly_less_than(const void* i, const void* j)
{
  double* a = (double*) i;
  double* b = (double*) j;
  
  return (*a < *b);
}

void generate_params(params& pa, hyperparams& hpa, subtype& st, gsl_rng* rng)
{
  for (int i=0; i<1024; ++i) // for appropriate random number generation
    gsl_rng_uniform(rng);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    pa.pa[i]->u = Log(gsl_rng_uniform(rng));

  // for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     for (int j=0; j<tr[i].children.size(); ++j)
  //       pa.pa[i]->beta[j] = Log(gsl_rng_uniform(rng));
    // }

  st.t = pa.pa[1]->u;
  st.n = Log(1);
  // calc_n(pa, hpa, tr);

  // Vdouble v (hpa.MAX_SUBTYPE, 0);
  
  // gsl_ran_dirichlet(rng, hpa.MAX_SUBTYPE, &hpa.gamma[1], &v[0]);
  // for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     pa.pa[i]->xi = Log(v[i-1]);
  //   }
}

void generate_binom(ofstream& f, int M, int n, params& pa, hyperparams& hpa, subtype& st, int seed, gsl_rng* rng, VVLog& gegen, VLog& gegen_int)
{
  st.total_cn = 2;
  st.variant_cn = 1;

  VLog vf(FRACTIONS + 1, Log(0));
  VLog vf_numerator(FRACTIONS + 1, Log(0));
  VLog vf_denominator(FRACTIONS + 1, Log(0));
  Log partition = Log(0);

  variant_fraction_partition(0, 1, Log(1), st.t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf, vf_numerator, vf_denominator, partition);

  VLog vf_cum(FRACTIONS + 1, Log(0));
  for (int s=1; s<=FRACTIONS; ++s)
    {
      vf_cum[s] = vf[s];
      cerr << vf[s].eval() << "\t";
    }
  cerr << endl;
  
  for (int s=2; s<=FRACTIONS; ++s)
    {
      vf_cum[s] = vf_cum[s-1] + vf_cum[s];
    }
  
  for (int k=0; k<n; ++k)
    {
      Log z = Log(gsl_rng_uniform(rng));
      for (int s=1; s<=FRACTIONS; ++s)
        {
          if (z < vf_cum[s])
            {
              st.x = Log(((double) s) / FRACTIONS);
              break;
            }
        }
      
      double mu = calc_mu(st, hpa);
      // cerr << "mu: " << mu << endl;
      
      unsigned int m = gsl_ran_binomial(rng, mu, M);
      f << m << "\t" << M << "\t" << st.x.eval() << endl;
    }
}

int main(int argc, char** argv)
{
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  gsl_set_error_handler_off ();
  cerr << scientific;
  
  if (argc != 6)
    {
      cerr << "usage: ./alpha_beta_single_generate M n seed (u outfile) (reads outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  int M, n, seed, MAX_SUBTYPE, TOTAL_CN, MAX_TREE;
  
  MAX_SUBTYPE = 1;
  TOTAL_CN = 2;
  M = atoi(argv[1]);
  n = atoi(argv[2]);
  seed = atoi(argv[3]);

  // trees tr;
  // trees_cons(tr, MAX_SUBTYPE);

  // MAX_TREE = tr.size();

  subtype st (1, 2, 1, 0, Log(0), Log(0), Log(1), Log(0), Log(0), NULL, NULL, std::vector<subtype*>(0));
  ofstream f (argv[4]);
  ofstream g (argv[5]);

  f << scientific;
  g << scientific;
  
  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  // T = gsl_rng_default; // default is mt19937
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed);

  VVLog gegen;
  gegen.assign(FRACTIONS+1, VLog(GEGEN_MAX+1, Log(0)));

  set_gegen(gegen);

  VLog gegen_int (FRACTIONS+1, Log(0));
  VLog gegen_int_err (FRACTIONS+1, Log(0));

  set_gegen_integral(gegen_int, gegen_int_err);

  hyperparams hpa (MAX_SUBTYPE, TOTAL_CN, MAX_TREE);
   
  params pa (hpa);
  
  generate_params(pa, hpa, st, r);
  
  write_params(f, pa, hpa);
  
  generate_binom(g, M, n, pa, hpa, st, seed, r, gegen, gegen_int);

  gsl_rng_free (r);
  
  f.close();
  g.close();
  
  return 0;
}
