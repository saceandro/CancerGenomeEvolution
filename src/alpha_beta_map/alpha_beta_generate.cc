#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "../../util/enumtree_prior.hh"
using namespace std;

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
    pa.pa[i].u = Log(gsl_rng_uniform(rng));

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int j=0; j<tr[i].children.size(); ++j)
        pa.pa[i]->beta[j] = Log(gsl_rng_uniform(rng));
    }

  calc_t(pa, hpa, tr);
  calc_n(pa, hpa, tr);
}

void generate_binom(ofstream& f, ofstream& g, int M, int n, params& pa, hyperparams& hpa, subtypes& tr, int seed, gsl_rng* rng)
{
  params pa_cum (hpa);

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

  for (int a=0; a<hpa.MAX_TREE; ++a)
    {
      calc_n(tr[a], hpa); // deleted calc_t() because t is the variables generated. Indeed, n = t for all i because I set labmda to be 1.
    }
  write_t_n(g, tr[2], hpa);
  
  for (int k=0; k<n; ++k)
    {
      int topology = 2;
      Log w = Log(gsl_rng_uniform(rng));
      
      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        {
          Log x = Log(gsl_rng_uniform(rng));
          for (int l=1; l<=hpa.TOTAL_CN; ++l)
            {
              if (x < pa_cum.pa[i]->pi[l])
                {
                  tr[topology][i].total_cn = l;
                  break;
                }
            }
          
          Log y = Log(gsl_rng_uniform(rng));
          for (int r=1; r<=tr[topology][i].total_cn; ++r)
            {
              if (y < pa_cum.pa[i]->kappa[tr[topology][i].total_cn][r])
                {
                  tr[topology][i].variant_cn = r;
                  break;
                }
            }
        }

      double mu = calc_mu(tr[topology], hpa);
      // cerr << "mu: " << mu << endl;
      
      unsigned int m = gsl_ran_binomial(rng, mu, M);
      f << m << "\t" << M << endl;
    }
}

int main(int argc, char** argv)
{
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  cerr << scientific;
  
  if (argc != 9)
    {
      cerr << "usage: ./lda_tree_generate_ms max_subtype total_cn M n seed (u_pi_kappa outfile) (t_n outfile) (reads outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  int M, n, seed, MAX_SUBTYPE, TOTAL_CN, MAX_TREE;
  
  MAX_SUBTYPE = atoi(argv[1]);
  TOTAL_CN = atoi(argv[2]);
  M = atoi(argv[3]);
  n = atoi(argv[4]);
  seed = atoi(argv[5]);

  if (MAX_SUBTYPE != 4)
    {
      cerr << "This program assumes max_subtype to be 4." << endl;
      exit(EXIT_FAILURE);
    }
  
  trees tr;
  trees_cons(tr, MAX_SUBTYPE);
  MAX_TREE = tr.size();

  ofstream f (argv[6]);
  ofstream g (argv[7]);
  ofstream h (argv[8]);
  
  f << scientific;
  g << scientific;
  
  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  // T = gsl_rng_default; // default is mt19937
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed);
  hyperparams hpa (MAX_SUBTYPE, TOTAL_CN, MAX_TREE);
   
  params pa (hpa);
  
  generate_params(pa, hpa, tr, r);
  
  calc_u_from_t(pa, hpa, tr[2]);
  write_params(f, pa, hpa);
  
  generate_binom(h, g, M, n, pa, hpa, tr, seed, r);

  gsl_rng_free (r);
  
  f.close();
  g.close();
  h.close();
  
  return 0;
}
