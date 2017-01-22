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
#include "../../util/enumtree_mutation_em.hh"
using namespace std;

extern void variant_fraction_t_partition(Log init_frac, Log n_q, Log t_q, Log t, Log beta_tilda_q, VVLog& gegen, VLog& vf);
extern void variant_fraction_partition(int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, VLog& vf, VLog& vf_numerator, VLog& vf_denominator, Log& partition);
extern void set_gegen(VVLog &gegen);
extern void set_gegen_integral(VLog &gegen_int, VLog &gegen_int_err);

typedef vector<VVLog> VVVLog;

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

  return sum.eval();
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

void generate_binom(ofstream& f1, ofstream& f2, int M, params& pa, hyperparams& hpa, subtypes& tr, int seed, gsl_rng* rng, VVLog& gegen, VLog& gegen_int, int all_snvs)
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

  vector<unsigned int> snvs (hpa.MAX_SUBTYPE + 1, 0);
  Vdouble xi_double (hpa.MAX_SUBTYPE);
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      xi_double[i-1] = tr[i].xi.eval();
    }
  
  gsl_ran_multinomial(rng, hpa.MAX_SUBTYPE, all_snvs, &xi_double[0], &snvs[1]);
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      cerr << "i = " << i << endl;
      
      Log Nni = Log(CELL_MAX) * tr[i].n;
      Log lNni = Nni.take_log_Log();
      cerr << "snvs: " << snvs[i] << endl;
      
      vector<unsigned int> snv_h (tr[i].omega.size());
      Vdouble omega_double (tr[i].omega.size());
      for (int ch=0; ch<tr[i].omega.size(); ++ch)
        {
          omega_double[ch] = tr[i].omega[ch].eval();
        }
      
      gsl_ran_multinomial(rng, tr[i].omega.size(), snvs[i], &omega_double[0], &snv_h[0]);

      for (int h=0; h<pow(2, tr[i].children.size()); ++h)
        {
          calc_child_x(tr[i], hpa, h);

          int eldest_ch_index = calc_eldest_child_index(tr[i], h);

          vector<unsigned int> snv_h_x (FRACTIONS + 1);
          Vdouble vf_double (FRACTIONS + 1);
          for (int fr=0; fr<=FRACTIONS; ++fr)
            {
              vf_double[fr] = vf[i][eldest_ch_index][fr].eval();
            }
          
          gsl_ran_multinomial(rng, FRACTIONS+1, snv_h[h], &vf_double[0], &snv_h_x[0]);

          for (int x=0; x<=FRACTIONS; ++x)
            {
              tr[i].x = Log(((double) x) / FRACTIONS);

              double mu = calc_mu(tr, hpa);
              
              for (int l=0; l<snv_h_x[x]; ++l)
                {
                  unsigned int m = gsl_ran_binomial(rng, mu, M);

                  if (i==1)
                    {
                      f1 << m << "\t" << M << "\t" << i << "\t";
                      for (int j=1; j<=hpa.MAX_SUBTYPE; ++j)
                        {
                          f1 << tr[j].inherited << "\t";
                        }
                      f1 << tr[i].x.eval() << "\t" << mu << endl;
                    }
                  else
                    {
                      f2 << m << "\t" << M << "\t" << i << "\t";
                      for (int j=1; j<=hpa.MAX_SUBTYPE; ++j)
                        {
                          f2 << tr[j].inherited << "\t";
                        }
                      f2 << tr[i].x.eval() << "\t" << mu << endl;
                    }
                }
            }
          clear_inherited(tr, hpa);
        }
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
      cerr << "usage: ./generate_reads_em_estimate_model_mu_separate_subtype_snv max_subtype M topology (#SNV) seed (u_n infile) (reads subtype1 outfile) (reads subtype2 outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  int M, topology, seed, MAX_SUBTYPE, TOTAL_CN, MAX_TREE;
  
  MAX_SUBTYPE = atoi(argv[1]);
  TOTAL_CN = 2;
  M = atoi(argv[2]);
  topology = atoi(argv[3]);
  int all_snvs = atoi(argv[4]);
  cerr << "all_snvs: " << all_snvs << endl;
  
  seed = atoi(argv[5]);

  trees tr;
  trees_cons(tr, MAX_SUBTYPE);
  MAX_TREE = tr.size();

  ifstream f (argv[6]);
  ofstream read_subtype1 (argv[7]);
  ofstream read_subtype2 (argv[8]);
  
  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  // T = gsl_rng_default; // default is mt19937
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed);
  for (int i=0; i<1024; ++i) // for appropriate random number generation
    gsl_rng_uniform(r);

  VVLog gegen;
  gegen.assign(FRACTIONS+1, VLog(GEGEN_MAX+1, Log(0)));

  set_gegen(gegen);

  VLog gegen_int (FRACTIONS+1, Log(0));
  VLog gegen_int_err (FRACTIONS+1, Log(0));

  set_gegen_integral(gegen_int, gegen_int_err);

  hyperparams hpa (MAX_SUBTYPE, TOTAL_CN, MAX_TREE);
   
  params pa (hpa);

  read_params(f, pa, hpa);
  calc_subtypes(pa, hpa, tr[topology]);

  generate_binom(read_subtype1, read_subtype2, M, pa, hpa, tr[topology], seed, r, gegen, gegen_int, all_snvs);
  gsl_rng_free (r);
  
  f.close();
  read_subtype1.close();
  read_subtype2.close();
  
  return 0;
}
