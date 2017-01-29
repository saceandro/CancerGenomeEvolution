#include <iostream>
#include <stdio.h>
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
#define calc_asigmoid(x) (2.0 * atanh(2.0 * (x) - 1.0))

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

void params_rmsd(double& rmsd_u, double& rmsd_n, params& pa1, params& pa2, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double diff = (pa1.pa[i]->u - pa2.pa[i]->u).eval();
      rmsd_u += diff * diff;
    }
  rmsd_u = sqrt(rmsd_u / hpa.MAX_SUBTYPE);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double diff = (pa1.pa[i]->n - pa2.pa[i]->n).eval();
      rmsd_n += diff * diff;
    }
  rmsd_n = sqrt(rmsd_n / hpa.MAX_SUBTYPE);
}

void calc_rmsd(double& rmsd_u, double& rmsd_t, double& rmsd_n, params& pa1, params& pa2, subtypes& tr1, subtypes& tr2, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double diff = (pa1.pa[i]->u - pa2.pa[i]->u).eval();
      rmsd_u += diff * diff;
    }
  rmsd_u = sqrt(rmsd_u / hpa.MAX_SUBTYPE);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double diff = (tr1[i].t - tr2[i].t).eval();
      rmsd_t += diff * diff;
    }
  rmsd_t = sqrt(rmsd_t / hpa.MAX_SUBTYPE);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double diff = (pa1.pa[i]->n - pa2.pa[i]->n).eval();
      rmsd_n += diff * diff;
    }
  rmsd_n = sqrt(rmsd_n / hpa.MAX_SUBTYPE);
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

  pa.pa[0]->n = Log(0);
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

void read_params(std::ifstream& f, params& pa, hyperparams& hpa)
{
  double a;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f >> a;
      pa.pa[i]->x = a;
    }

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f >> a;
      pa.pa[i]->y = a;
    }

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

int calc_fdf(VLog& du, VLog& dn, Log& qfunc, Log& llik, hyperparams& hpa, int num_of_split, string& u1, string& u2, string& n, string& snvs, string& seed, string& est_iter)
{
  double a;
  
  for (int sp=0; sp<num_of_split; ++sp)
    {
      std::ostringstream filenum_str;
      filenum_str << sp;

      ifstream f (("../du_dn_llik_qsub_snv_seed_rand_param/" + u1 + "/" + u2 + "/" + n + "/" + snvs + "/" + seed + "/" + est_iter + "/" + filenum_str.str()).c_str());

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
      if (a > 0) return 1;
      
      qfunc += Log(a);

      f >> a;
      llik += Log(a);

      f.close();
    }
  return 0;
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
  // cerr << qfunc.eval() << endl;
}

Log d_n_s(params& pa, int i, int j) // corrected
{
  if (i == j)
    return pa.pa[i]->n * (Log(1) - pa.pa[i]->n/(Log(1) - pa.pa[0]->n));
  
  return -pa.pa[i]->n * pa.pa[j]->n/(Log(1) - pa.pa[0]->n);
}

struct ComputeFdf {
  ComputeFdf(subtypes& _tr, hyperparams& _hpa, char* _pa_test_file, char* _vf_test_file, VVLog& _gegen, VLog& _gegen_int, int _num_of_split, Log& _llik_final, string& _u1, string& _u2, string& _n, string& _snvs, string& _seed, string& _est_iter, int& _flag) : iter(0), tr(_tr), hpa(_hpa), pa_test_file(_pa_test_file), vf_test_file(_vf_test_file), gegen(_gegen), gegen_int(_gegen_int), num_of_split(_num_of_split), llik_final(_llik_final), u1(_u1), u2(_u2), n(_n), snvs(_snvs), seed(_seed), est_iter(_est_iter), flag(_flag) {}
  
  int operator()(const Vdouble& x, double& fn, vector<double>& gr)
  {
    // pa_test_f.seekp(0);
    // vf_test_f.seekp(0);

    ofstream pa_test_f(pa_test_file);
    ofstream vf_test_f(vf_test_file);
    pa_test_f << scientific << setprecision(10);
    vf_test_f << scientific << setprecision(10);

    params pa_test (hpa);
    calc_params_from_x(x, pa_test, hpa);
    write_params(pa_test_f, pa_test, hpa);
    write_params((ofstream&) cerr, pa_test, hpa);
    pa_test_f.close();
    
    calc_subtypes(pa_test, hpa, tr);

    for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
      {
        if (!((Log(CELL_MAX) * tr[i].n) > Log(1))) // if #CELL_i <= 1
          {
            flag = 1;
            return 1; // end gradient descent
          }
      }
    
    VVVLog vf_test (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
    VVVLog dtvf_test (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
    VVVLog dthvf_test (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
    VVVLog dnvf_test (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  
    calc_vf_dvf(vf_test, dtvf_test, dthvf_test, dnvf_test, pa_test, hpa, tr, gegen, gegen_int);
    write_vf(vf_test_f, vf_test, hpa);
    write_vf(vf_test_f, dtvf_test, hpa);
    write_vf(vf_test_f, dthvf_test, hpa);
    write_vf(vf_test_f, dnvf_test, hpa);
    vf_test_f.close();

    std::ostringstream filenum_str;
    filenum_str << num_of_split;
    system(("qsub -N estep_snv" + snvs + "_seed" + seed + "_rand" + est_iter + "_param" + u1 + "_" + u2 + "_" + n + " -e ../log_qsub_snv_seed_rand_param/estep.err -o ../log_qsub_snv_seed_rand_param/estep.log -sync y -tc 100 -t 1-" + filenum_str.str() + ":1 estep_snv_seed_rand_param.sh " + u1 + " " + u2 + " " + n + " " + snvs + " " + seed + " " + est_iter).c_str());
    
    VLog du (hpa.MAX_SUBTYPE + 1, Log(0));
    VLog dn (hpa.MAX_SUBTYPE + 1, Log(0));
    Log qfunc (0);
    llik_final = Log(0);

    int qfunc_flag = 0;
    qfunc_flag = calc_fdf(du, dn, qfunc, llik_final, hpa, num_of_split, u1, u2, n, snvs, seed, est_iter);
    if (qfunc_flag==1) // vf_new became zero even if vf neq zero
      {
        flag = 2;
        return 1; // end gradient descent
      }
    
    fn = -qfunc.eval(); // minimize!
    cerr << iter << "\t" << -fn << endl;
    cerr << "du_dn_total:" <<endl;
    for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
      {
        cerr << du[i].eval() << "\t";
      }
    cerr << endl;
    for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
      {
        cerr << dn[i].eval() << "\t";
      }
    cerr << endl;
    
    gr.assign(x.size(), 0);

    cerr << "dx_dy_total:" << endl;
    for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
      {
        gr[i-1] = -calc_dx_sigmoid(pa_test.pa[i]->x) * du[i].eval(); // minimize!
        cerr << gr[i-1] << "\t";
      }

    for (int index=1; index<=hpa.MAX_SUBTYPE; ++index)
      {
        Log grad = Log(0);
        
        for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
          {
            grad += d_n_s(pa_test, index, i) * dn[i];
          }
        gr[hpa.MAX_SUBTYPE + index-1] = -grad.eval(); // minimize!
        cerr << gr[hpa.MAX_SUBTYPE + index-1] << "\t";
      }
    cerr << endl << "-------------------------------------------------------------------------------" << endl << endl;
    
    ++iter;
    
    return 0;
  }
  
  int    iter;
  subtypes& tr;
  hyperparams& hpa;
  char* pa_test_file;
  char* vf_test_file;
  VVLog& gegen;
  VLog& gegen_int;
  int num_of_split;
  Log& llik_final;
  string& u1;
  string& u2;
  string& n;
  string& snvs;
  string& seed;
  string& est_iter;
  int& flag;
};

int main(int argc, char* argv[]) {
  cout << scientific << setprecision(10);
  cerr << scientific << setprecision(10);
  
  // feenableexcept(FE_INVALID);
  gsl_set_error_handler_off ();
  
  if (argc != 23)
    {
      cerr << "usage: ./mstep_qsub_snv_seed_rand_param subtype topology num_of_split (x y init infile) (pa_test outfile) (vf_test outfile) (llik outfile) (pa_best outfile) (vf_old outfile) (params log fie) (true params) (rmsd file) (param difference) (u1) (u2) (n) (#SNV) (seed) (iter) (u lower bound) (em max iter) (grad desc max iter)" << endl;
      exit(EXIT_FAILURE);
    }
  
  int MAX_SUBTYPE, TOTAL_CN, MAX_TREE;

  MAX_SUBTYPE = atoi(argv[1]);
  TOTAL_CN = 2;
  
  trees trs;
  trees_cons(trs, MAX_SUBTYPE);
  MAX_TREE = trs.size();

  trees trs_true;
  trees_cons(trs_true, MAX_SUBTYPE);
  
  hyperparams hpa (MAX_SUBTYPE, TOTAL_CN, MAX_TREE);

  VVLog gegen;
  gegen.assign(FRACTIONS+1, VLog(GEGEN_MAX+1, Log(0)));
  set_gegen(gegen);

  VLog gegen_int (FRACTIONS+1, Log(0));
  VLog gegen_int_err (FRACTIONS+1, Log(0));
  set_gegen_integral(gegen_int, gegen_int_err);

  int topology = atoi(argv[2]);
  int num_of_split = atoi(argv[3]);

  char* pa_old_file = argv[4];
  char* pa_test_file = argv[5];
  char* vf_test_file = argv[6];
  ofstream llik_f(argv[7], ofstream::out | ofstream::app);
  ofstream pa_best_f(argv[8]);
  char* vf_old_file = argv[9];
  ofstream pa_log_f(argv[10]);
  ifstream pa_true_f(argv[11]);
  ofstream rmsd_f(argv[12]);
  ofstream params_diff_f(argv[13]);
  string u1 (argv[14]);
  string u2 (argv[15]);
  string n (argv[16]);
  string snvs (argv[17]);
  string seed (argv[18]);
  string est_iter (argv[19]);
  double u_lower = atof(argv[20]);
  int em_max_iter = atoi(argv[21]);
  int grad_desc_max_iter = atoi(argv[22]);

  cerr << "snvs: " << snvs << "\t" << "seed: " << seed << "\t" << "iter: " << est_iter << endl;;
  
  llik_f << scientific << setprecision(10);
  pa_best_f << scientific << setprecision(10);
  pa_log_f << scientific << setprecision(10);
  rmsd_f << scientific << setprecision(10);
  params_diff_f << scientific << setprecision(10);
  
  params pa_true(hpa);
  read_params(pa_true_f, pa_true, hpa);
  calc_subtypes(pa_true, hpa, trs_true[topology]);
  pa_true_f.close();
  
  Lbfgsb minimizer;

  Vdouble lower(2*hpa.MAX_SUBTYPE, 0);
  Vdouble upper(2*hpa.MAX_SUBTYPE, 0);
  Vint bdtype(2*hpa.MAX_SUBTYPE, 0);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      lower[i-1] = calc_asigmoid(u_lower);
      bdtype[i-1] = 1;
    }
  
  minimizer.set_eps(1.0e-5);
  minimizer.set_maxit(grad_desc_max_iter);
  minimizer.set_bounds(lower, upper, bdtype);
  
  Vdouble x0 (2 * hpa.MAX_SUBTYPE);

  for (int i=1; i<=em_max_iter; ++i) // do not use index i to loop subtype in this block!
    {
      ifstream pa_old_f_in(pa_old_file);
      read_x(pa_old_f_in, x0, hpa);
      pa_old_f_in.close();

      params pa_prev(hpa);
      calc_params_from_x(x0, pa_prev, hpa);

      int flag = 0;
      Log llik_final (0);
      minimizer.minimize(x0, ComputeFdf(trs[topology], hpa, pa_test_file, vf_test_file, gegen, gegen_int, num_of_split, llik_final, u1, u2, n, snvs, seed, est_iter, flag));
      if (flag == 1) cerr << "minimizer Nn <= 1 !" << endl;
      else if (flag == 2) cerr << "minimizer numerically cannot calculate error" << endl;
      
      const vector<double>& x = minimizer.best_x();
      
      params pa_old(hpa);
      calc_params_from_x(x, pa_old, hpa);
      ofstream pa_old_f(pa_old_file);
      pa_old_f << scientific << setprecision(10);
      write_params(pa_old_f, pa_old, hpa);
      pa_old_f.close();
      write_params(pa_log_f, pa_old, hpa);
      
      llik_f << i << "\t" << llik_final.eval() << "\t" << minimizer.best_fn() << endl;

      calc_subtypes(pa_old, hpa, trs[topology]);
      
      double rmsd_u = 0;
      double rmsd_t = 0;
      double rmsd_n = 0;
      calc_rmsd(rmsd_u, rmsd_t, rmsd_n, pa_old, pa_true, trs[topology], trs_true[topology], hpa);
      rmsd_f << i << "\t" << rmsd_u << "\t" << rmsd_t << "\t" << rmsd_n << endl;

      double u_diff = 0;
      double n_diff = 0;
      params_rmsd(u_diff, n_diff, pa_old, pa_prev, hpa);
      params_diff_f << i << "\t" << u_diff << "\t" << n_diff << endl;
      cerr << "u_diff: " << u_diff << "\tn_diff" << n_diff << endl;

      cerr << "=========================================================================================" << endl;
      
      for (int j=1; j<=hpa.MAX_SUBTYPE; ++j)
      {
        if (!((Log(CELL_MAX) * trs[topology][j].n) > Log(1))) // if #CELL_j <= 1
          {
            cerr << "Nn <= 1 !" << endl;
            goto finalize_em;
          }

        if ((pa_old.pa[j]->u < Log(1e-5)))
          {
            cerr << "u became too small\tu[" << j << "] = "<< pa_old.pa[j]->u.eval() << endl;
            goto finalize_em;
          }
      }

      if ((u_diff < 1e-6) && (n_diff < 1e-6))
        {
          cerr << "converged. padiff < 1e-6" << endl;
          goto finalize_em;
        }
        
      VVVLog vf_old (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
      VVVLog dtvf_old (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
      VVVLog dthvf_old (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
      VVVLog dnvf_old (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));

      calc_vf_dvf(vf_old, dtvf_old, dthvf_old, dnvf_old, pa_old, hpa, trs[topology], gegen, gegen_int);
      ofstream vf_old_f(vf_old_file);
      vf_old_f << scientific << setprecision(10);
      write_vf(vf_old_f, vf_old, hpa);
      write_vf(vf_old_f, dtvf_old, hpa);
      write_vf(vf_old_f, dthvf_old, hpa);
      write_vf(vf_old_f, dnvf_old, hpa);
      vf_old_f.close();
    }

 finalize_em:
  const vector<double>& y = minimizer.best_x();
  params pa_best(hpa);
  calc_params_from_x(y, pa_best, hpa);
  write_params(pa_best_f, pa_best, hpa);

  llik_f.close();
  pa_best_f.close();
  pa_log_f.close();
  rmsd_f.close();
  params_diff_f.close();
  
  return 0;
}
