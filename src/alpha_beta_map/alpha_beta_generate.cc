#include "setting.hh"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "../../util/enumtree_wf.hh"
using namespace std;

extern void variant_fraction_partition(int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, VLog& vf, VLog& vf_numerator, VLog& vf_denominator, Log& partition);
extern void set_gegen(VVLog &gegen);
extern void set_gegen_integral(VLog &gegen_int, VLog &gegen_int_err);

void write_params(std::ofstream& f, params& pa, hyperparams& hpa, subtypes& tr)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    f << pa.pa[i]->u.eval() << "\t";
  f << endl;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int j=0; j<tr[i].children.size(); ++j)
        f << pa.pa[i]->beta[j].eval() << "\t";
      f << endl;
    }

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    f << pa.pa[i]->xi.eval() << "\t";
  f << endl;
    
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double sum = 0;
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        f << pa.pa[i]->pi[l].eval() << "\t";
      f << endl;

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (int r=1; r<=l; ++r)
            f << pa.pa[i]->kappa[l][r].eval() << "\t";
          f << endl;
        }
      f << endl;
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

double calc_mu(subtypes& st, hyperparams& hpa)
{
  Log denom;
  Log num;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      denom += st[i].n * Log(st[i].total_cn);
      num += st[i].n * st[i].x * Log(st[i].variant_cn);
    }
    
  return (num / denom).eval();
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

void generate_params(params& pa, hyperparams& hpa, subtypes& tr, gsl_rng* rng)
{
  for (int i=0; i<1024; ++i) // for appropriate random number generation
    gsl_rng_uniform(rng);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    pa.pa[i]->u = Log(gsl_rng_uniform(rng));

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int j=0; j<tr[i].children.size(); ++j)
        pa.pa[i]->beta[j] = Log(gsl_rng_uniform(rng));
    }

  calc_t(pa, hpa, tr);
  calc_n(pa, hpa, tr);

  Vdouble v (hpa.MAX_SUBTYPE, 0);
  Vdouble w (hpa.TOTAL_CN, 0);
  
  gsl_ran_dirichlet(rng, hpa.MAX_SUBTYPE, &hpa.gamma[1], &v[0]);
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->xi = Log(v[i-1]);
    }
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      gsl_ran_dirichlet(rng, hpa.TOTAL_CN, &hpa.alpha[1], &w[0]);
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        pa.pa[i]->pi[l] = Log(w[l-1]);

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          gsl_ran_dirichlet(rng, l, &hpa.beta[l][1], &w[0]);
          for (int r=1; r<=l; ++r)
            pa.pa[i]->kappa[l][r] = Log(w[r-1]);
        }
    }
}

void generate_binom(ofstream& f, int M, int n, params& pa, hyperparams& hpa, subtypes& tr, int seed, gsl_rng* rng, VVLog& gegen, VLog& gegen_int)
{
  params pa_cum (hpa);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa_cum.pa[i]->xi = pa.pa[i]->xi;
    }

  for (int i=2; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa_cum.pa[i]->xi = pa_cum.pa[i-1]->xi + pa_cum.pa[i]->xi;
    }
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          pa_cum.pa[i]->pi[l] = pa.pa[i]->pi[l];
        }

      for (int l=2; l<=hpa.TOTAL_CN; ++l)
        {
          pa_cum.pa[i]->pi[l] = pa_cum.pa[i]->pi[l-1] + pa_cum.pa[i]->pi[l];
        }

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (int r=1; r<=l; ++r)
            {
              pa_cum.pa[i]->kappa[l][r] = pa.pa[i]->kappa[l][r];
            }

          for (int r=2; r<=l; ++r)
            {
              pa_cum.pa[i]->kappa[l][r] = pa_cum.pa[i]->kappa[l][r-1] + pa_cum.pa[i]->kappa[l][r];
            }
        }
    }

  int q;
  
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
      
      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        {
          Log x = Log(gsl_rng_uniform(rng));
          for (int l=1; l<=hpa.TOTAL_CN; ++l)
            {
              if (x < pa_cum.pa[i]->pi[l])
                {
                  tr[i].total_cn = l;
                  break;
                }
            }
          
          Log y = Log(gsl_rng_uniform(rng));
          for (int r=1; r<=tr[i].total_cn; ++r)
            {
              if (y < pa_cum.pa[i]->kappa[tr[i].total_cn][r])
                {
                  tr[i].variant_cn = r;
                  break;
                }
            }
        }

      VLog vf(FRACTIONS + 1, Log(0));
      VLog vf_numerator(FRACTIONS + 1, Log(0));
      VLog vf_denominator(FRACTIONS + 1, Log(0));
      Log partition = Log(0);
      
      variant_fraction_partition(0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf, vf_numerator, vf_denominator, partition);

      VLog vf_cum(FRACTIONS + 1, Log(0));
      for (int s=1; s<=FRACTIONS; ++s)
        {
          vf_cum[s] = vf[s];
        }
      
      for (int s=2; s<=FRACTIONS; ++s)
        {
          vf_cum[s] = vf_cum[s-1] + vf_cum[s];
        }

      Log z = Log(gsl_rng_uniform(rng));
      for (int s=1; s<=FRACTIONS; ++s)
        {
          if (z < vf_cum[s])
            {
              tr[q].x = ((double) s) / FRACTIONS;
              break;
            }
        }
      
      double mu = calc_mu(tr, hpa);
      // cerr << "mu: " << mu << endl;
      
      unsigned int m = gsl_ran_binomial(rng, mu, M);
      f << m << "\t" << M << "\t" << q << "\t";
      
      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        {
          f << tr[i].total_cn << "\t" << tr[i].variant_cn << "\t";
        }
      f << endl;
    }
}

int main(int argc, char** argv)
{
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  // _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  cerr << scientific;
  
  if (argc != 10)
    {
      cerr << "usage: ./alpha_beta_generate max_subtype total_cn M n seed (u_pi_kappa outfile) (t_n outfile) (reads outfile) topology" << endl;
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

  ofstream f (argv[6]);
  ofstream g (argv[7]);
  ofstream h (argv[8]);

  int a = atoi(argv[9]);
  
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
  
  generate_params(pa, hpa, tr[a], r);
  
  write_params(f, pa, hpa, tr[a]);
  write_t_n(g, tr[a], hpa);
  
  generate_binom(h, M, n, pa, hpa, tr[a], seed, r, gegen, gegen_int);

  gsl_rng_free (r);
  
  f.close();
  g.close();
  h.close();
  
  return 0;
}
