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
#include "../../util/enumtree_wf_n.hh"
using namespace std;

extern void variant_fraction_partition(int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, VLog& vf, VLog& vf_numerator, VLog& vf_denominator, Log& partition);
extern void set_gegen(VVLog &gegen);
extern void set_gegen_integral(VLog &gegen_int, VLog &gegen_int_err);

void write_params(std::ofstream& f, params& pa, hyperparams& hpa, subtypes& tr)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f << pa.pa[i]->u.eval() << "\t";
    }
  f << endl;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f << pa.pa[i]->n.eval() << "\t";
    }
  f << endl;
}

void read_params(std::ifstream& f, params& pa, hyperparams& hpa, subtypes& tr)
{
  double a;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f >> a;
      pa.pa[i]->u = Log(a);
    }

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f >> a;
      pa.pa[i]->n = Log(a);
    }
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f >> a;
      pa.pa[i]->xi = Log(a);
    }
}

void write_t_n(std::ofstream& f, subtypes& st, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    f << st[i].t.eval() << "\t";
  f << endl;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    f << st[i].n.eval() << "\t";
  f << endl << endl;
}

double calc_mu(subtypes& tr, hyperparams& hpa, int q)
{
  return (tr[q].n * tr[q].x / Log(2)).eval(); // corrected
}

// double calc_mu(subtypes& st, hyperparams& hpa)
// {
//   Log denom;
//   Log num;
//   for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
//     {
//       denom += st[i].n * Log(st[i].total_cn);
//       num += st[i].n * st[i].x * Log(st[i].variant_cn);
//     }
    
//   return (num / denom).eval();
// }

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

void generate_params(params& pa, hyperparams& hpa, subtypes& tr, gsl_rng* rng)
{
  for (int i=0; i<1024; ++i) // for appropriate random number generation
    gsl_rng_uniform(rng);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    pa.pa[i]->u = Log(gsl_rng_uniform(rng));

  Vdouble w (hpa.MAX_SUBTYPE + 1, 0);
  gsl_ran_dirichlet(rng, hpa.MAX_SUBTYPE + 1, &hpa.gamma[0], &w[0]);
  
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->n = Log(w[i]);
    }

  calc_t(pa, hpa, tr);
  calc_n(pa, hpa, tr);

  Vdouble v (hpa.MAX_SUBTYPE, 0);
  gsl_ran_dirichlet(rng, hpa.MAX_SUBTYPE, &hpa.gamma[1], &v[0]);
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->xi = Log(v[i-1]);
    }
}

void generate_binom(ofstream& f, int M, int n, params& pa, hyperparams& hpa, subtypes& tr, int seed, gsl_rng* rng, VVLog& gegen, VLog& gegen_int)
{
  VVLog vf (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0)));
  VVLog vf_numerator (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0)));
  VVLog vf_denominator (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0)));
  VLog partition (hpa.MAX_SUBTYPE + 1, Log(0));

  for (int q=1; q<=hpa.MAX_SUBTYPE; ++q)
    variant_fraction_partition(0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q], vf_numerator[q], vf_denominator[q], partition[q]);

  VVLog vf_cum(hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0)));

  for (int q=1; q<=hpa.MAX_SUBTYPE; ++q)
    {
      for (int s=1; s<=FRACTIONS; ++s)
        {
          vf_cum[q][s] = vf[q][s];
        }
      
      for (int s=2; s<=FRACTIONS; ++s)
        {
          vf_cum[q][s] = vf_cum[q][s-1] + vf_cum[q][s];
        }
    }
  
  params pa_cum (hpa);
  int q;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa_cum.pa[i]->xi = pa.pa[i]->xi;
    }

  for (int i=2; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa_cum.pa[i]->xi = pa_cum.pa[i-1]->xi + pa_cum.pa[i]->xi;
    }
  
  for (int k=0; k<n; ++k)
    {
      Log w = Log(gsl_rng_uniform(rng));
      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        {
          if (w < pa_cum.pa[i]->xi)
            {
              q = i;
              break;
            }
        }
      
      // for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
      //   {
      //     tr[i].total_cn = 2;
      //     tr[i].variant_cn = 1;
      //   }
      
      Log z = Log(gsl_rng_uniform(rng));
      for (int s=1; s<=FRACTIONS; ++s)
        {
          if (z < vf_cum[q][s])
            {
              tr[q].x = Log(((double) s) / FRACTIONS);
              break;
            }
        }
      
      double mu = calc_mu(tr, hpa, q);
      // cerr << "mu: " << mu << endl;
      
      unsigned int m = gsl_ran_binomial(rng, mu, M);
      f << m << "\t" << M << "\t" << q << "\t" << endl;
    }
}

int main(int argc, char** argv)
{
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  gsl_set_error_handler_off ();
  cerr << scientific;
  
  if (argc != 9)
    {
      cerr << "usage: ./generate_reads max_subtype total_cn M n seed (u_n_xi infile) (reads outfile) topology" << endl;
      exit(EXIT_FAILURE);
    }

  int M, n, seed, MAX_SUBTYPE, TOTAL_CN, MAX_TREE;
  
  MAX_SUBTYPE = atoi(argv[1]);
  TOTAL_CN = atoi(argv[2]);
  M = atoi(argv[3]);
  n = atoi(argv[4]);
  seed = atoi(argv[5]);

  trees tr;
  trees_cons(tr, MAX_SUBTYPE);
  MAX_TREE = tr.size();

  ifstream f (argv[6]);
  ofstream h (argv[7]);

  int a = atoi(argv[8]);
  
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

  read_params(f, pa, hpa, tr[a]);
  calc_t(pa, hpa, tr[a]);
  calc_n(pa, hpa, tr[a]);

  generate_binom(h, M, n, pa, hpa, tr[a], seed, r, gegen, gegen_int);

  gsl_rng_free (r);
  
  f.close();
  h.close();
  
  return 0;
}
