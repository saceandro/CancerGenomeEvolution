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
#include "../../util/enumtree_wf_n_r_inherited.hh"
using namespace std;

extern void variant_fraction_partition(int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, VLog& vf, VLog& vf_numerator, VLog& vf_denominator, Log& partition);
extern void set_gegen(VVLog &gegen);
extern void set_gegen_integral(VLog &gegen_int, VLog &gegen_int_err);

typedef std::vector<VVLog> VVVLog;

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

double calc_mu(subtypes& tr, hyperparams& hpa)
{
  Log sum = 0;
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      sum += tr[i].n * tr[i].x;
    }

  return sum.eval() / 2.0;
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
  VVVLog vf (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog vf_numerator (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog vf_denominator (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVLog partition (hpa.MAX_SUBTYPE + 1, VLog (hpa.MAX_SUBTYPE + 1, Log(0)));

  for (int q=1; q<=hpa.MAX_SUBTYPE; ++q)
    {
      variant_fraction_partition(0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q][q], vf_numerator[q][q], vf_denominator[q][q], partition[q][q]);
      for (std::vector<subtype*>::iterator ch=tr[q].children.begin(); ch!=tr[q].children.end(); ++ch)
        {
          int ch_index = (*ch)->index;
          variant_fraction_partition(1, q, tr[q].n, tr[q].t, tr[ch_index].t, Log(BETA_TILDA), gegen, gegen_int, vf[q][ch_index], vf_numerator[q][ch_index], vf_denominator[q][ch_index], partition[q][ch_index]);
        }
    }

  VVVLog vf_cum (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));

  for (int q=1; q<=hpa.MAX_SUBTYPE; ++q)
    {
      for (int s=0; s<=FRACTIONS; ++s)
        {
          vf_cum[q][q][s] = vf[q][q][s];
          // cerr << vf_cum[q][q][s].eval() << "\t";
        }
      // cerr << endl;
      
      for (int s=1; s<=FRACTIONS; ++s)
        {
          vf_cum[q][q][s] = vf_cum[q][q][s-1] + vf_cum[q][q][s];
        }
      
      for (std::vector<subtype*>::iterator ch=tr[q].children.begin(); ch!=tr[q].children.end(); ++ch)
        {
          int ch_index = (*ch)->index;
          for (int s=0; s<=FRACTIONS; ++s)
            {
              vf_cum[q][ch_index][s] = vf[q][ch_index][s];
            }
          
          for (int s=1; s<=FRACTIONS; ++s)
            {
              vf_cum[q][ch_index][s] = vf_cum[q][ch_index][s-1] + vf_cum[q][ch_index][s];
            }
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

  int count = 0;
  while (count < n)
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

      bool observed = false;
      if (tr[q].children.size() > 0)
        {
          for (std::vector<subtype*>::iterator ch=tr[q].children.begin(); ch!=tr[q].children.end(); ++ch)
            {
              double inherit_prob = gsl_rng_uniform(rng);
              int y = gsl_ran_bernoulli(rng, inherit_prob);
              if (y == 1)
                {
                  observed = true;
                  mark_inherited(*ch);
                }
            }
        }
      
      int eldest_ch_index = q;
      for (std::vector<subtype*>::iterator ch=tr[q].children.begin(); ch!=tr[q].children.end(); ++ch)
        {
          if ((*ch)->inherited == 1)
            {
              eldest_ch_index = (*ch)->index;
              break;
            }
        }

      Log z = Log(gsl_rng_uniform(rng));

      int s = 1;
      if (eldest_ch_index != q)
        {
          s = 0;
        }
      for (; s<=FRACTIONS; ++s)
        {
          if (z < vf_cum[q][eldest_ch_index][s])
            {
              if (s > 0)
                {
                  observed = true;
                }
              tr[q].x = Log(((double) s) / FRACTIONS);
              break;
            }
        }

      if (observed)
        {
          double mu = calc_mu(tr, hpa);
          // cerr << "mu: " << mu << endl;
          
          unsigned int m = gsl_ran_binomial(rng, mu, M);

          if (m > 0)
            {
              count++;
              f << m << "\t" << M << "\t" << q << "\t";
              for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
                {
                  f << tr[i].inherited << "\t";
                  cout << tr[i].x.eval() << "\t";
                }
              f << endl;
              cout << endl;
            }
        }
      clear_inherited(tr, hpa);
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
