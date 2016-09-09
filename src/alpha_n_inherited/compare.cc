#include "setting.hh"
#include <fstream>
#include <algorithm>
#include <numeric>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf.h>
#include "../../util/enumtree_wf_n_r_inherited.hh"
#include <xmmintrin.h>

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

typedef vector<int> Vint;
typedef vector<double> Vdouble;
typedef vector<Vdouble> VVdouble;
typedef vector<bool> Vbool;
typedef vector<Vbool> VVbool;
typedef pair<int, int> READ;
typedef vector<READ*> READS;
typedef vector<VVLog> VVVLog;
typedef vector<int> INHERITEDS;
typedef pair<int*, INHERITEDS*> Q;
typedef vector< Q* > QS;
// typedef vector<int*> QS;

typedef struct _du
{
  int i;
  READS &res;
  QS &qs;
  gsl_vector *x;
  hyperparams &hpa;
  subtypes& tr;
  VVLog gegen;
  VLog gegen_int;
  Log purity;
  
  _du (int _i, READS& _res, QS& _qs, gsl_vector* _x, hyperparams& _hpa, subtypes& _tr, VVLog& _gegen, VLog& _gegen_int, Log& _purity) : i(_i), res(_res), qs(_qs), x(_x), hpa(_hpa), tr(_tr), gegen(_gegen), gegen_int(_gegen_int), purity(_purity) {}
}
  du;

typedef struct _dn
{
  int index;
  READS &res;
  QS &qs;
  gsl_vector *x;
  hyperparams &hpa;
  subtypes& tr;
  VVLog gegen;
  VLog gegen_int;
  Log purity;
  
  _dn (int _index, READS& _res, QS& _qs, gsl_vector* _x, hyperparams& _hpa, subtypes& _tr, VVLog& _gegen, VLog& _gegen_int, Log& _purity) : index(_index), res(_res), qs(_qs), x(_x), hpa(_hpa), tr(_tr), gegen(_gegen), gegen_int(_gegen_int), purity(_purity) {}
}
  dn;

double calc_dx_sigmoid(double x)
{
  double y = tanh(x/2.0);
  return (1.0 - y*y) / 4.0;
}

double sum_vector(Vdouble& v, int s, int e)
{
  double sum = 0;
  for (int i=s; i<=e; ++i)
    sum += v[i];
  return sum;
}

void write_params(std::ofstream& f, params& pa, hyperparams& hpa, subtypes& tr)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f << pa.pa[i]->u.eval() << "\t";
    }
  f << endl;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f << pa.pa[i]->n.eval() << "\t";
    }
  f << endl;
}

void write_t_n(std::ostream& f, subtypes& st, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    f << st[i].t.eval() << "\t";
  f << endl;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    f << st[i].n.eval() << "\t";
  f << endl << endl;
}

void calc_child_x(subtypes& tr, hyperparams&hpa, INHERITEDS& ihs)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      if (ihs[i] == 1)
        {
          tr[i].x = Log(1);
        }
    }
}

void clear_x(subtypes& tr, hyperparams&hpa, INHERITEDS& ihs)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      tr[i].x = Log(0);
    }
}

double calc_mu(subtypes& tr, hyperparams& hpa)
{
  Log sum = Log(0);
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      sum += tr[i].n * tr[i].x;
    }

  return sum.eval() / 2.0;
}

#define log_binomial_pdf(m, mu, M) ((m) * log((mu)) + ((M) - (m)) * log1p(-(mu)) + gsl_sf_lnchoose((M), (m)))

Log d_bin_mu(READ& re, double mu)
{
  Log a (log_binomial_pdf(re.first-1, mu, re.second-1), 1);
  if (re.first < re.second)
    {
      Log b (log_binomial_pdf(re.first, mu, re.second-1), 1);
      return Log(re.second) * (a - b);
    }
  else
    {
      return Log(re.second) * a;
    }
}

Log d_mu_n(subtypes& st, hyperparams& hpa, int i)
{
  return st[i].x / Log(2);
}

void calc_params(const gsl_vector* x, params& pa, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->u = Log(calc_sigmoid(gsl_vector_get(x, i-1)));
    }

  Log sum = Log(0);
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->n = Log(gsl_vector_get(x, hpa.MAX_SUBTYPE + i-1)).take_exp();
      sum += pa.pa[i]->n;
    }

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->n = (Log(1) - pa.pa[0]->n) * pa.pa[i]->n / sum;
    }
}

double calc_llik(READS& res, QS& qs, params& pa, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int)
{
  int K;
  K = res.size();

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
          // variant_fraction_partition(0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q], vf_numerator[q], vf_denominator[q], partition[q]);
        }
    }
  
  Log lik = Log(1);
  
  for (int k=0; k<K; ++k)
    {
      int index = *(qs[k]->first);
      INHERITEDS inh = *(qs[k]->second);

      int eldest_ch_index = index;
      for (std::vector<subtype*>::iterator ch=tr[index].children.begin(); ch!=tr[index].children.end(); ++ch)
        {
          int ch_index = (*ch)->index;
          if (inh[ch_index] == 1)
            {
              eldest_ch_index = ch_index;
              break;
            }
        }
      
      Log lik_k = Log(0);
      for (int s=1; s<=FRACTIONS; ++s)
        {
          tr[index].x = Log(((double) s) / FRACTIONS);
          calc_child_x(tr, hpa, inh);
          double mu = calc_mu(tr, hpa);
          lik_k += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * vf[index][eldest_ch_index][s];
        }
      clear_x(tr, hpa, inh);
      
      lik *= lik_k;
    }

  return lik.take_log();
}

double calc_llik_for_du(double x_i, void* _du)
{
  du* p = (du*) _du;

  gsl_vector_set(p->x, p->i-1, x_i);

  params pa (p->hpa);
  pa.pa[0]->n = Log(1) - p->purity;
  calc_params(p->x, pa, p->hpa);

  calc_t(pa, p->hpa, p->tr);
  calc_n(pa, p->hpa, p->tr);

  // cerr << "calc_llik_for_du" << endl;
  // cerr << "params" << endl;
  // write_params((ofstream&)cerr, pa, p->hpa, p->tr);
  // cerr << "t_n" << endl;
  // write_t_n(cerr, p->tr, p->hpa);

  return calc_llik(p->res, p->qs, pa, p->hpa, p->tr, p->gegen, p->gegen_int);
}

double calc_llik_for_dn(double x_i_j, void* _dn)
{
  dn* p = (dn*) _dn;

  gsl_vector_set(p->x, p->hpa.MAX_SUBTYPE + p->index-1, x_i_j);

  params pa (p->hpa);
  pa.pa[0]->n = Log(1) - p->purity;
  calc_params(p->x, pa, p->hpa);

  calc_t(pa, p->hpa, p->tr);
  calc_n(pa, p->hpa, p->tr);

  // cerr << "calc_llik_for_du" << endl;
  // cerr << "params" << endl;
  // write_params((ofstream&)cerr, pa, p->hpa, p->tr);
  // cerr << "t_n" << endl;
  // write_t_n(cerr, p->tr, p->hpa);

  return calc_llik(p->res, p->qs, pa, p->hpa, p->tr, p->gegen, p->gegen_int);
}

double d_llik(READS& res, QS& qs, params& pa, params& grad, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int)
{
  int K;
  K = res.size();

  VVVLog vf (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dtvf (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dthvf (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dnvf (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));

  for (int q=1; q<=hpa.MAX_SUBTYPE; ++q)
    {
      d_variant_fraction_all(d_t_variant_fraction, 0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q][q], dtvf[q][q]);
      d_variant_fraction_all(d_n_variant_fraction, 0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q][q], dnvf[q][q]);
      
      for (std::vector<subtype*>::iterator ch=tr[q].children.begin(); ch!=tr[q].children.end(); ++ch)
        {
          int ch_index = (*ch)->index;
          d_variant_fraction_all(d_t_variant_fraction, 1, q, tr[q].n, tr[q].t, tr[ch_index].t, Log(BETA_TILDA), gegen, gegen_int, vf[q][ch_index], dtvf[q][ch_index]);
          d_variant_fraction_all(d_th_variant_fraction, 1, q, tr[q].n, tr[q].t, tr[ch_index].t, Log(BETA_TILDA), gegen, gegen_int, vf[q][ch_index], dthvf[q][ch_index]);

          // variant_fraction_partition(1, q, tr[q].n, tr[q].t, tr[ch_index].t, Log(BETA_TILDA), gegen, gegen_int, vf[q][ch_index], vf_numerator[q][ch_index], vf_denominator[q][ch_index], partition[q][ch_index]);
        }
    }

  // for (int q=1; q<=hpa.MAX_SUBTYPE; ++q)
  //   {
  //     d_variant_fraction_all(d_t_variant_fraction, 0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q], dtvf[q]);
  //     d_variant_fraction_all(d_n_variant_fraction, 0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q], dnvf[q]);
  //   }

  Log lik = Log(1);

  VLog d_t (hpa.MAX_SUBTYPE + 1, Log(0));
  
  for (int k=0; k<K; ++k)
    {
      int index = *(qs[k]->first);
      INHERITEDS inh = *(qs[k]->second);

      int eldest_ch_index = index;
      for (std::vector<subtype*>::iterator ch=tr[index].children.begin(); ch!=tr[index].children.end(); ++ch)
        {
          int ch_index = (*ch)->index;
          if (inh[ch_index] == 1)
            {
              eldest_ch_index = ch_index;
              break;
            }
        }

      Log lik_k = Log(0);
      Log d_t_k = Log(0);
      Log d_th_k = Log(0);
      VLog d_n_k (hpa.MAX_SUBTYPE + 1, Log(0));

      for (int s=1; s<=FRACTIONS; ++s)
        {
          tr[index].x = Log(((double) s) / FRACTIONS);
          calc_child_x(tr, hpa, inh);
          double mu = calc_mu(tr, hpa);
          lik_k += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * vf[index][eldest_ch_index][s];
          d_t_k += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * dtvf[index][eldest_ch_index][s];
          d_th_k += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * dthvf[index][eldest_ch_index][s];

          for (int i=1; i<=hpa.MAX_SUBTYPE; ++i) // start from subtype 1
            d_n_k[i] += d_mu_n(tr, hpa, i) * d_bin_mu(*res[k], mu) * vf[index][eldest_ch_index][s];

          d_n_k[index] += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * dnvf[index][eldest_ch_index][s];
        }
      clear_x(tr, hpa, inh);

      d_t[index] += d_t_k / lik_k;
      d_t[eldest_ch_index] += d_th_k / lik_k;
      
      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i) // start from subtype 1
        grad.pa[i]->n += d_n_k[i] / lik_k;
        
      lik *= lik_k;
    }

  // for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
  //   cerr << d_n[i].eval() << "\t";
  // cerr << endl;
  
  for (int j=1; j<=hpa.MAX_SUBTYPE; ++j)
    {
      Log acc = Log(0);
      for (int i=j; i<=hpa.MAX_SUBTYPE; ++i)
        {
          acc += d_t[i] * d_t_u(pa, tr, i, j);
        }
      grad.pa[j]->u += acc;
    }

  return lik.take_log();
}

double calc_rel_err(double x, double y)
{
  return fabs(x - y) / fabs(x);
}

double calc_dx_u_llik_numeric(int i, READS& res, QS& qs, gsl_vector* x, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int, Log purity)
{
  gsl_function F;
  du _du (i, res, qs, x, hpa, tr, gegen, gegen_int, purity);

  F.function = &calc_llik_for_du;
  F.params = &_du;

  double result, abserr;
  gsl_deriv_central(&F, gsl_vector_get(x, i-1), 1e-5, &result, &abserr);
  
  return result;
}

double calc_dx_u_llik_analytic(int i, READS& res, QS& qs, gsl_vector* x, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int, Log purity)
{
  int K = res.size();

  params pa (hpa);
  pa.pa[0]->n = Log(1) - purity;
  calc_params(x, pa, hpa);

  calc_t(pa, hpa, tr);
  calc_n(pa, hpa, tr);

  // cerr << "calc_dx_u_llik_analytic" << endl;
  // cerr << "params" << endl;
  // write_params((ofstream&)cerr, pa, hpa, tr);
  // cerr << "t_n" << endl;
  // write_t_n(cerr, tr, hpa);
  
  params grad (hpa);

  double llik = d_llik(res, qs, pa, grad, hpa, tr, gegen, gegen_int);

  return (Log(calc_dx_sigmoid(gsl_vector_get(x, i-1))) * grad.pa[i]->u).eval();
}

double calc_dx_n_llik_numeric(int index, READS& res, QS& qs, gsl_vector* x, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int, Log purity)
{
  gsl_function F;
  dn _dn (index, res, qs, x, hpa, tr, gegen, gegen_int, purity);

  F.function = &calc_llik_for_dn;
  F.params = &_dn;

  double result, abserr;
  gsl_deriv_central(&F, gsl_vector_get(x, hpa.MAX_SUBTYPE + index-1), 1e-5, &result, &abserr);
  
  return result;
}

Log d_n_s(params& pa, int i, int j) // corrected
{
  if (i == j)
    return pa.pa[i]->n * (Log(1) - pa.pa[i]->n/(Log(1) - pa.pa[0]->n));
  
  return -pa.pa[i]->n * pa.pa[j]->n/(Log(1) - pa.pa[0]->n);
}

double calc_dx_n_llik_analytic(int index, READS& res, QS& qs, gsl_vector* x, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int, Log purity)
{
  int K = res.size();

  params pa (hpa);
  pa.pa[0]->n = Log(1) - purity;
  calc_params(x, pa, hpa);

  calc_t(pa, hpa, tr);
  calc_n(pa, hpa, tr);

  params grad (hpa);

  double llik = d_llik(res, qs, pa, grad, hpa, tr, gegen, gegen_int);

  Log gr = Log(0);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      gr += d_n_s(pa, index, i) * grad.pa[i]->n;
    }
      
  return gr.eval();
}

void gsl_set_random(gsl_vector* x, hyperparams& hpa, gsl_rng* rng)
{
  for (int i=0; i<1024; ++i) // for appropriet random number generation
    gsl_rng_uniform(rng);

  for (int i=0; i<2*hpa.MAX_SUBTYPE; ++i)
    {
      double a = gsl_rng_uniform(rng);
      gsl_vector_set(x, i, a);
    }
}

int main(int argc, char** argv)
{
  cout << scientific;
  cerr << scientific;
  
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  
  gsl_set_error_handler_off ();
  
  if (argc != 11)
    {
      cerr << "usage: ./compare max_subtype n step tumor_purity (reads) (u diff outfile) (n diff outfile) topology u_index beta_index" << endl;
      exit(EXIT_FAILURE);
    }

  READS res;
  QS qs;
  int n, step, MAX_SUBTYPE, TOTAL_CN, MAX_TREE;

  MAX_SUBTYPE = atoi(argv[1]);
  TOTAL_CN = 2;

  trees trs;
  trees_cons(trs, MAX_SUBTYPE);
  MAX_TREE = trs.size();

  hyperparams hpa (MAX_SUBTYPE, TOTAL_CN, MAX_TREE);
  
  n = atoi(argv[2]);
  step = atoi(argv[3]);
  Log purity = Log(atof(argv[4]));
  
  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  // T = gsl_rng_default; // default is mt19937
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, 1);

  ifstream f (argv[5]);
  ofstream ff (argv[6]);
  ofstream g (argv[7]);

  int a = atoi(argv[8]);
  int u_index = atoi(argv[9]);
  int beta_index = atoi(argv[10]);

  ff << scientific;
  g << scientific;

  VVLog gegen;
  gegen.assign(FRACTIONS+1, VLog(GEGEN_MAX+1, Log(0)));

  set_gegen(gegen);

  VLog gegen_int (FRACTIONS+1, Log(0));
  VLog gegen_int_err (FRACTIONS+1, Log(0));

  set_gegen_integral(gegen_int, gegen_int_err);

  gsl_vector* x = gsl_vector_alloc(2*hpa.MAX_SUBTYPE);

  for (int k=0; k<n; ++k)
    {
      READ *re = new READ;
      int *a = new int;
      INHERITEDS *ih = new INHERITEDS (MAX_SUBTYPE + 1, 0);

      f >> re->first >> re->second >> *a;
      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        {
          f >> (*ih)[i];
        }
      
      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        {
          cout << (*ih)[i] << "\t";
        }
      cout << endl;
      
      Q *q = new Q (a, ih);
      res.push_back(re);
      qs.push_back(q);
    }

  for (int i=-step; i<=step; ++i)
    {
      gsl_set_random(x, hpa, r);
      gsl_vector_set(x, u_index-1, 1.0 * ((double)i) / step);
      double num = calc_dx_u_llik_numeric(u_index, res, qs, x, hpa, trs[a], gegen, gegen_int, purity);
      double analytic = calc_dx_u_llik_analytic(u_index, res, qs, x, hpa, trs[a], gegen, gegen_int, purity);

      ff << gsl_vector_get(x, u_index-1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
      // cerr << "-----------------------------------------------------------------------------------------------------------------" << endl;
    }

  for (int i=-step; i<=step; ++i)
    {
      gsl_set_random(x, hpa, r);
      gsl_vector_set(x, hpa.MAX_SUBTYPE + beta_index-1, 1.0 * ((double)i) / step);
      double num = calc_dx_n_llik_numeric(beta_index, res, qs, x, hpa, trs[a], gegen, gegen_int, purity);
      double analytic = calc_dx_n_llik_analytic(beta_index, res, qs, x, hpa, trs[a], gegen, gegen_int, purity);

      if (fabs(num) > 0)
        g << gsl_vector_get(x, hpa.MAX_SUBTYPE + beta_index-1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
      else
        g << gsl_vector_get(x, hpa.MAX_SUBTYPE + beta_index-1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
    }

  for (int i=0; i<n; ++i)
    {
      delete res[i];
      delete qs[i];
    }

  f.close();
  g.close();
  ff.close();

  return 0;
}
