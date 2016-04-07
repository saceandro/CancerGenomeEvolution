#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "../../util/enumtree.hh"
using namespace std;

void write_params(std::ofstream& f, params& pa, hyperparams& hpa)
{
  for (int a=0; a<hpa.MAX_TREE; ++a)
    f << pa.rho[a].eval() << "\t";
  f << endl << endl;
  
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    f << pa.pa[i]->u.eval() << "\t";
  f << endl << endl;

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

void calc_n(subtypes& sts, hyperparams& hpa)
{
  Log sum;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      sts[i].n = sts[i].t;
      sum += sts[i].t;
    }
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    sts[i].n /= sum;
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
      num += st[i].n * Log(st[i].variant_cn);
    }
    
  return (num / denom).eval();
}

void generate_params(params& pa, hyperparams& hpa, gsl_rng* rng)
{
  Vdouble v (hpa.MAX_TREE, 0);
  
  for (int i=0; i<1024; ++i) // for appropriet random number generation
    gsl_rng_uniform(rng);

  gsl_ran_dirichlet(rng, hpa.MAX_TREE, &hpa.gamma[0], &v[0]);
  
  for (int a=0; a<=hpa.MAX_TREE; ++a)
    pa.rho[a] = Log(v[a]);

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->u = Log(gsl_ran_beta(rng, hpa.be_hpa.first, hpa.be_hpa.second));
    }

  Vdouble w (hpa.TOTAL_CN, 0);
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

void generate_binom(ofstream& f, ofstream& g, int M, int n, params& pa, hyperparams& hpa, trees& tr, int seed, gsl_rng* rng)
{
  params pa_cum (hpa);

  for (int a=0; a<hpa.MAX_TREE; ++a)
    {
      pa_cum.rho[a] = pa.rho[a];
    }
  for (int a=1; a<hpa.MAX_TREE; ++a)
    {
      pa_cum.rho[a] += pa_cum.rho[a-1];
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

  for (int a=0; a<hpa.MAX_TREE; ++a)
    {
      calc_t(pa, hpa, tr[a]);
      calc_n(tr[a], hpa);
      write_t_n(g, tr[a], hpa);
    }
  
  for (int k=0; k<n; ++k)
    {
      int topology;
      Log w = Log(gsl_rng_uniform(rng));
      for (int a=0; a<hpa.MAX_TREE; ++a)
        {
          if (w < pa_cum.rho[a])
            {
              topology = a;
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
      cerr << "usage: ./lda_tree_generate max_subtype total_cn M n seed (rho_u_pi_kappa outfile) (t_n outfile) (reads outfile)" << endl;
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
  
  generate_params(pa, hpa, r);
  write_params(f, pa, hpa);

  generate_binom(h, g, M, n, pa, hpa, tr, seed, r);

  gsl_rng_free (r);
  
  f.close();
  g.close();
  h.close();
  
  return 0;
}
