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
#include "../../util/enumtree_wf_n_r.hh"
using namespace std;

extern void variant_fraction_t_partition(Log init_frac, Log n_q, Log t_q, Log t, Log beta_tilda_q, VVLog& gegen, VLog& vf);
extern void set_gegen(VVLog &gegen);

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

void read_params(std::ifstream& f, params& pa, hyperparams& hpa)
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
}

void generate_binom(ofstream& f, int M, int t_fractions, params& pa, hyperparams& hpa, subtypes& tr, int seed, gsl_rng* rng, VVLog& gegen, long long int bp)
{
  double dt = 1.0 / t_fractions;

  int until = 0;
  if (bp >= 7500000000LL) // to avoid time consume
    {
      until = 9;
    }
  
  Log t;
  for (int h = t_fractions - 1; h >= until; --h)
    {
      t = h*dt;

      cout << "t: " << t.eval() << endl;
      
      for (int i=1; tr[i].t > t; ++i)
        {
          // cout << "subtype: " << i << endl;
          Log Nni = Log(CELL_MAX) * tr[i].n;
          Log lNni = Nni.take_log_Log();
          Log Nni_t = ((Log(1) - t / tr[i].t) * lNni).take_exp();
          Log lambda = Nni_t * Log(bp) * pa.pa[i]->r * Log(dt);
          // cout << "lambda: " << lambda.eval() << endl;
          unsigned int s = gsl_ran_poisson(rng, lambda.eval());
          
          VLog vf (FRACTIONS + 1, Log(0));
          variant_fraction_t_partition(Nni_t.inverse(), tr[i].n, tr[i].t, t, BETA_TILDA, gegen, vf);
          // for (int a=0; a<=FRACTIONS; ++a)
          //   cerr << vf[a].eval() << "\t";
          // cerr << endl;

          VLog vf_cum(FRACTIONS + 1, Log(0));
          for (int x=0; x<=FRACTIONS; ++x)
            {
              vf_cum[x] = vf[x];
            }
          for (int x=1; x<=FRACTIONS; ++x)
            {
              vf_cum[x] = vf_cum[x-1] + vf_cum[x];
            }
          
          // cout << "s: " << s << endl;
          for (int k=0; k<s; ++k)
            {
              bool observed = true;
              Log z = Log(gsl_rng_uniform(rng));
              for (int x=0; x<=FRACTIONS; ++x)
                {
                  if (z < vf_cum[x])
                    {
                      if (x == 0)
                        {
                          observed = false;
                        }
                      tr[i].x = Log(((double) x) / FRACTIONS);
                      break;
                    }
                }

              for (std::vector<subtype*>::iterator ch=tr[i].children.begin(); ch!=tr[i].children.end(); ++ch)
                {
                  if ((*ch)->t < t) // bugfix
                    {
                      int y = gsl_ran_bernoulli(rng, Nni_t.inverse().eval());
                      // (*ch)->x = Log(y);
                      if (y == 1)
                        {
                          observed = true;
                          mark_inherited(*ch);
                        }
                    }
                }

              // for (std::vector<subtype*>::iterator ch=tr[i].children.begin(); ch!=tr[i].children.end(); ++ch)
              //   {
              //     cerr << (*ch)->x.eval() << "\t";
              //   }
              // cerr << endl;

              double mu = calc_mu(tr, hpa);
              // cout << mu << "\t";
              // for (int l=1; l<=hpa.MAX_SUBTYPE; ++l)
              //   {
              //     cout << tr[l].x.eval() << "\t";
              //   }
              
              if (observed)
                {
                  // for (int j=1; j<=hpa.MAX_SUBTYPE; ++j)
                  //   {
                  //     cerr << tr[j].inherited << "\t";
                  //   }
                  // cerr << endl;
                  
                  // for (int j=1; j<=hpa.MAX_SUBTYPE; ++j)
                  //   {
                  //     cerr << tr[j].x.eval() << "\t";
                  //   }
                  // cerr << endl;
                  
                  unsigned int m = gsl_ran_binomial(rng, mu, M);
                  if (m == 0)
                    {
                      // cerr << "m == 0!" << endl;
                    }
                  else
                    {
                      // cout << m << "\t" << M;
                      f << m << "\t" << M << "\t" << i << "\t";
                      for (int l=1; l<=hpa.MAX_SUBTYPE; ++l)
                        {
                          f << tr[l].inherited << "\t";
                        }
                      f << endl;
                    }
                }
              // cout << endl;
              clear_inherited(tr, hpa);
            }
          // cout << endl;
        }
      // cerr << endl << endl;
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
      cerr << "usage: ./generate_reads_poisson max_subtype M t_fractions topology bp seed (u_n infile) (reads outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  int M, t_fractions, topology, seed, MAX_SUBTYPE, TOTAL_CN, MAX_TREE;
  
  MAX_SUBTYPE = atoi(argv[1]);
  TOTAL_CN = 2;
  M = atoi(argv[2]);
  t_fractions = atoi(argv[3]);
  topology = atoi(argv[4]);
  long long int bp = atoll(argv[5]);
  cout << "bp: " << bp << endl;
  
  seed = atoi(argv[6]);

  trees tr;
  trees_cons(tr, MAX_SUBTYPE);
  MAX_TREE = tr.size();

  ifstream f (argv[7]);
  ofstream h (argv[8]);
  
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

  hyperparams hpa (MAX_SUBTYPE, TOTAL_CN, MAX_TREE);
   
  params pa (hpa);

  read_params(f, pa, hpa);
  calc_t(pa, hpa, tr[topology]);
  calc_n(pa, hpa, tr[topology]);

  generate_binom(h, M, t_fractions, pa, hpa, tr[topology], seed, r, gegen, bp);
  
  gsl_rng_free (r);
  
  f.close();
  h.close();
  
  return 0;
}
