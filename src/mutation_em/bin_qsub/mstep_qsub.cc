#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf.h>
#include <iomanip>
#include "../../../util/enumtree_mutation_em_multinomial.hh"
#include "../setting.hh"
#include "lbfgsb.hpp"

using namespace std;

#define calc_sigmoid(x) ((tanh((x)/2.0) + 1.0) / 2.0)

typedef void (*myfunc) (int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& partition);

extern void variant_fraction_partition(int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, VLog& vf, VLog& vf_numerator, VLog& vf_denominator, Log& partition);
extern void d_variant_fraction_all(myfunc d_x_variant_fraction, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, VLog& vf, VLog& dvf);
extern void d_t_variant_fraction(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& partition);
extern void d_th_variant_fraction(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& denominator);
extern void d_n_variant_fraction(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& denominator);
extern void set_gegen(VVLog &gegen);
extern void set_gegen_integral(VLog &gegen_int, VLog &gegen_int_err);

typedef vector<VVLog> VVVLog;

double calc_dx_sigmoid(double x)
{
  double y = tanh(x/2.0);
  return (1.0 - y*y) / 4.0;
}

void copy_x_y(const Vdouble& x, params& pa, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->x = x[i-1];
    }
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->y = x[hpa.MAX_SUBTYPE + i-1];
    }
}

void read_x(std::ifstream& f, Vdouble& x, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f >> x[i-1];
    }

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f >> x[hpa.MAX_SUBTYPE + i-1];
    }
}

void calc_params_from_x_y(params& pa, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->u = Log(calc_sigmoid(pa.pa[i]->x));
    }

  Log sum = Log(0);
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->n = Log(pa.pa[i]->y).take_exp();
      sum += pa.pa[i]->n;
    }

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->n = (Log(1) - pa.pa[0]->n) * pa.pa[i]->n / sum;
    }
}

void calc_params_from_x(const Vdouble& x, params& pa, hyperparams& hpa)
{
  copy_x_y(x, pa, hpa);
  calc_params_from_x_y(pa, hpa);
}

void write_params(std::ofstream& f, params& pa, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f << pa.pa[i]->x << "\t";
    }
  f << endl;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f << pa.pa[i]->y << "\t";
    }
  f << endl;

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

void calc_vf_dvf(VVVLog& vf, VVVLog& dtvf, VVVLog& dthvf, VVVLog& dnvf, params& pa, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int)
{
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

  for (int q=1; q<=hpa.MAX_SUBTYPE; ++q)
    {
      d_variant_fraction_all(d_t_variant_fraction, 0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q][q], dtvf[q][q]);
      d_variant_fraction_all(d_n_variant_fraction, 0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q][q], dnvf[q][q]);
      
      for (std::vector<subtype*>::iterator ch=tr[q].children.begin(); ch!=tr[q].children.end(); ++ch)
        {
          int ch_index = (*ch)->index;
          d_variant_fraction_all(d_t_variant_fraction, 1, q, tr[q].n, tr[q].t, tr[ch_index].t, Log(BETA_TILDA), gegen, gegen_int, vf[q][ch_index], dtvf[q][ch_index]);
          d_variant_fraction_all(d_th_variant_fraction, 1, q, tr[q].n, tr[q].t, tr[ch_index].t, Log(BETA_TILDA), gegen, gegen_int, vf[q][ch_index], dthvf[q][ch_index]);
          d_variant_fraction_all(d_n_variant_fraction, 1, q, tr[q].n, tr[q].t, tr[ch_index].t, Log(BETA_TILDA), gegen, gegen_int, vf[q][ch_index], dnvf[q][ch_index]);
        }
    }
}

void write_vf(std::ofstream& f, VVVLog& vf, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int j=1; j<=hpa.MAX_SUBTYPE; ++j)
        {
          for (int s=0; s<=FRACTIONS; ++s)
            {
              f << vf[i][j][s].eval() << endl;
            }
          f << endl;
        }
      f << endl << endl;
    }
  f << endl << endl << endl;
}

void calc_fdf(VLog& du, VLog& dn, Log& llik, hyperparams& hpa, int num_of_split)
{
  char str[1024];
  double a;
  
  for (int sp=0; sp<num_of_split; ++sp)
    {
      int n = sprintf(str, "../du_dn_llik_qsub/%d", sp);

      ifstream f (str);

      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        {
          f >> a;
          du[i] += Log(a);
        }

      for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
        {
          f >> a;
          dn[i] += Log(a);
        }
      
      f >> a;
      llik += Log(a);

      f.close();
    }
  // for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     cerr << du[i].eval() << "\t";
  //   }
  // cerr << endl;
  
  // for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     cerr << dn[i].eval() << "\t";
  //   }
  // cerr << endl;
  // cerr << llik.eval() << endl;
}

Log d_n_s(params& pa, int i, int j) // corrected
{
  if (i == j)
    return pa.pa[i]->n * (Log(1) - pa.pa[i]->n/(Log(1) - pa.pa[0]->n));
  
  return -pa.pa[i]->n * pa.pa[j]->n/(Log(1) - pa.pa[0]->n);
}

struct ComputeFdf {
  ComputeFdf(subtypes& _tr, hyperparams& _hpa, ofstream& _pa_test_f, ofstream& _vf_test_f, ofstream& _llik_f, VVLog& _gegen, VLog& _gegen_int, int _num_of_split) : iter(0), tr(_tr), hpa(_hpa), pa_test_f(_pa_test_f), vf_test_f(_vf_test_f), llik_f(_llik_f), gegen(_gegen), gegen_int(_gegen_int), num_of_split(_num_of_split) {}
  int operator()(const Vdouble& x, double& fn, vector<double>& gr) {
    fn = 0;

    params pa_test (hpa);
    pa_test.pa[0]->n = Log(0);
    calc_params_from_x(x, pa_test, hpa);
    write_params(pa_test_f, pa_test, hpa);
    
    calc_subtypes(pa_test, hpa, tr);

    VVVLog vf_test (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
    VVVLog dtvf_test (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
    VVVLog dthvf_test (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
    VVVLog dnvf_test (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  
    calc_vf_dvf(vf_test, dtvf_test, dthvf_test, dnvf_test, pa_test, hpa, tr, gegen, gegen_int);
    write_vf(vf_test_f, vf_test, hpa);
    write_vf(vf_test_f, dtvf_test, hpa);
    write_vf(vf_test_f, dthvf_test, hpa);
    write_vf(vf_test_f, dnvf_test, hpa);

    std::ostringstream filenum_str;
    filenum_str << num_of_split;
    system(("qsub -N estep -e ../log_qsub/estep.err -o ../log_qsub/estep.out -sync y -tc 200 -t 1-" + filenum_str.str() + ":1 estep.sh").c_str());
    
    VLog du (hpa.MAX_SUBTYPE + 1, Log(0));
    VLog dn (hpa.MAX_SUBTYPE + 1, Log(0));
    Log llik (0);
    calc_fdf(du, dn, llik, hpa, num_of_split);
    fn = llik.eval();
    llik_f << iter << "\t" << fn << endl;
    
    gr.assign(x.size(), 0);

    for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
      {
        gr[i-1] = calc_dx_sigmoid(pa_test.pa[i]->x) * du[i].eval();
      }

    for (int index=1; index<=hpa.MAX_SUBTYPE; ++index)
      {
        Log grad = Log(0);
        
        for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
          {
            grad += d_n_s(pa_test, index, i) * dn[i];
          }
        gr[hpa.MAX_SUBTYPE + index-1] = grad.eval();
      }

    system("mv ../params_qsub/new ../params_qsub/old");
    
    ++iter;
    
    return 0;
  }
  
  int    iter;
  subtypes& tr;
  hyperparams& hpa;
  ofstream& pa_test_f;
  ofstream& vf_test_f;
  ofstream& llik_f;
  VVLog& gegen;
  VLog& gegen_int;
  int num_of_split;
};

int main(int argc, char* argv[]) {
  cout << scientific << setprecision(10);
  cerr << scientific << setprecision(10);
  
  // feenableexcept(FE_INVALID);
  gsl_set_error_handler_off ();
  
  if (argc != 9)
    {
      cerr << "usage: ./mstep_qsub subtype topology num_of_split (x y init infile) (pa_test outfile) (vf_test outfile) (llik outfile) (pa_best outfile)" << endl;
      exit(EXIT_FAILURE);
    }
  
  int n, MAX_SUBTYPE, TOTAL_CN, MAX_TREE;

  MAX_SUBTYPE = atoi(argv[1]);
  TOTAL_CN = 2;
  
  trees trs;
  trees_cons(trs, MAX_SUBTYPE);
  MAX_TREE = trs.size();

  hyperparams hpa (MAX_SUBTYPE, TOTAL_CN, MAX_TREE);

  VVLog gegen;
  gegen.assign(FRACTIONS+1, VLog(GEGEN_MAX+1, Log(0)));
  set_gegen(gegen);

  VLog gegen_int (FRACTIONS+1, Log(0));
  VLog gegen_int_err (FRACTIONS+1, Log(0));
  set_gegen_integral(gegen_int, gegen_int_err);

  int topology = atoi(argv[2]);
  int num_of_split = atoi(argv[3]);

  ifstream init_param_f(argv[4]);
  ofstream pa_test_f(argv[5]);
  ofstream vf_test_f(argv[6]);
  ofstream llik_f(argv[7], ofstream::out | ofstream::app);
  ofstream pa_best_f(argv[8]);

  pa_test_f << scientific << setprecision(10);
  vf_test_f << scientific << setprecision(10);
  llik_f << scientific << setprecision(10);
  pa_best_f << scientific << setprecision(10);
  
  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  // T = gsl_rng_default; // default is mt19937
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, 1);

  Lbfgsb minimizer;
  
  minimizer.set_eps(1.0e-7);
  minimizer.set_maxit(1);
  
  Vdouble x0 (2 * hpa.MAX_SUBTYPE);
  read_x(init_param_f, x0, hpa);
  
  minimizer.minimize(x0, ComputeFdf(trs[topology], hpa, pa_test_f, vf_test_f, llik_f, gegen, gegen_int, num_of_split));
  
  cout << "iter:" << minimizer.iter() << ",best_fn:" << minimizer.best_fn() << ",best_x:";
  
  const vector<double>& x = minimizer.best_x();

  params pa_best(hpa);
  calc_params_from_x(x, pa_best, hpa);

  write_params(pa_best_f, pa_best, hpa);
  
  return 0;
}
