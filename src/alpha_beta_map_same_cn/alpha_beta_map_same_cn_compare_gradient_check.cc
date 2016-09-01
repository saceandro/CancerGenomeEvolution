#include "setting.hh"
#include <fstream>
#include <algorithm>
#include <numeric>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf.h>
#include "../../util/enumtree_wf.hh"
#include <xmmintrin.h>

using namespace std;

#define calc_sigmoid(x) ((tanh((x)/2.0) + 1.0) / 2.0)

typedef void (*myfunc) (int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& partition);

extern void variant_fraction_partition(int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, VLog& vf, VLog& vf_numerator, VLog& vf_denominator, Log& partition);
extern void d_variant_fraction_all(myfunc d_x_variant_fraction, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, VLog& vf, VLog& dvf);
extern void d_t_variant_fraction(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& partition);
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
typedef vector<int*> QS;

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
  
  _du (int _i, READS& _res, QS& _qs, gsl_vector* _x, hyperparams& _hpa, subtypes& _tr, VVLog& _gegen, VLog& _gegen_int) : i(_i), res(_res), qs(_qs), x(_x), hpa(_hpa), tr(_tr), gegen(_gegen), gegen_int(_gegen_int)  {}
}
  du;

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
  f << endl << endl;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int j=0; j<(int)tr[i].children.size(); ++j)
        {
          f << pa.pa[i]->beta[j].eval() << "\t";
        }
      f << endl;
    }
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

double calc_mu(subtypes& tr, hyperparams& hpa, int q)
{
  return (tr[q].x / Log(2) / (Log(1) + tr[0].n/tr[q].n)).eval();
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

Log d_mu_n(subtypes& st, hyperparams& hpa, int i, int q)
{
  Log a = st[0].n + st[q].n;
  
  if (i == 0)
    return -st[q].n * st[q].x / Log(2) / a / a;

  if (i == q)
    return st[0].n * st[q].x / Log(2) / a / a;

  return Log(0);
}

void calc_params(const gsl_vector* x, params& pa, hyperparams& hpa, subtypes& tr)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->u = Log(calc_sigmoid(gsl_vector_get(x, i)));
    }

  int count = 0;
  for (int i=0; i<hpa.MAX_SUBTYPE; ++i)
    {
      for (int j=0; j<(int)tr[i].children.size(); ++j)
        {
          pa.pa[i]->beta[j] = Log(calc_sigmoid(gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + count)));
          count++;
        }
    }
}

double calc_llik(READS& res, QS& qs, params& pa, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int)
{
  int K;
  K = res.size();

  VVLog vf (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0)));
  VVLog vf_numerator (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0)));
  VVLog vf_denominator (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0)));
  VLog partition (hpa.MAX_SUBTYPE + 1, Log(0));

  for (int q=1; q<=hpa.MAX_SUBTYPE; ++q)
    variant_fraction_partition(0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q], vf_numerator[q], vf_denominator[q], partition[q]);

  Log lik = Log(1);
  
  for (int k=0; k<K; ++k)
    {
      Log lik_k = Log(0);
      for (int s=1; s<=FRACTIONS; ++s)
        {
          tr[*qs[k]].x = Log(((double) s) / FRACTIONS);
          double mu = calc_mu(tr, hpa, *qs[k]);
          lik_k += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * vf[*qs[k]][s];
        }
      tr[*qs[k]].x = Log(0);
      
      lik *= lik_k;
    }

  return lik.take_log();
}

double calc_llik_for_du(double x_i, void* _du)
{
  du* p = (du*) _du;

  gsl_vector_set(p->x, p->i, x_i);

  params pa (p->hpa);
  calc_params(p->x, pa, p->hpa, p->tr);

  for (int i=1; i<=p->hpa.MAX_SUBTYPE; ++i)
    p->tr[i].t = pa.pa[i]->u;
  
  // calc_t(pa, p->hpa, p->tr);
  calc_n(pa, p->hpa, p->tr);

  cerr << "calc_llik_for_du" << endl;
  cerr << "params" << endl;
  write_params((ofstream&)cerr, pa, p->hpa, p->tr);
  cerr << "t_n" << endl;
  write_t_n(cerr, p->tr, p->hpa);

  return calc_llik(p->res, p->qs, pa, p->hpa, p->tr, p->gegen, p->gegen_int);
}

double d_llik(READS& res, QS& qs, params& pa, params& grad, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int, myfunc d_x_variant_fraction)
{
  int K;
  K = res.size();

  VVLog vf (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0)));
  VVLog dtvf (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0)));
  
  for (int q=1; q<=hpa.MAX_SUBTYPE; ++q)
    d_variant_fraction_all(d_t_variant_fraction, 0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q], dtvf[q]);

  Log lik = Log(1);

  VLog d_t (hpa.MAX_SUBTYPE + 1, Log(0));
  
  for (int k=0; k<K; ++k)
    {
      Log lik_k = Log(0);
      Log d_t_k = Log(0);

      for (int s=1; s<=FRACTIONS; ++s)
        {
          tr[*qs[k]].x = Log(((double) s) / FRACTIONS);
          double mu = calc_mu(tr, hpa, *qs[k]);
          lik_k += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * vf[*qs[k]][s];
          d_t_k += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * dtvf[*qs[k]][s];
        }
      tr[*qs[k]].x = Log(0);
      d_t_k /= lik_k;

      d_t[*qs[k]] += d_t_k;
      lik *= lik_k;

      // for (int j=1; j<=hpa.MAX_SUBTYPE; ++j)
      //   {
      //     Log acc = Log(0);
      //     for (int i=j; i<=hpa.MAX_SUBTYPE; ++i)
      //       {
      //         acc += d_t[i] * d_t_u(pa, tr, i, j);
      //       }
      //     grad.pa[j]->u += acc;
      //   }
    }

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      grad.pa[i]->u += d_t[i];
    }

  return lik.take_log();
}

double calc_rel_err(double x, double y)
{
  return fabs(x - y) / fabs(x);
}

double calc_dx_u_llik_numeric(int i, READS& res, QS& qs, gsl_vector* x, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int)
{
  gsl_function F;
  du _du (i, res, qs, x, hpa, tr, gegen, gegen_int);

  F.function = &calc_llik_for_du;
  F.params = &_du;

  double result, abserr;
  gsl_deriv_central(&F, gsl_vector_get(x, i), 1e-5, &result, &abserr);
  
  return result;
}

double calc_dx_u_llik_analytic(int i, READS& res, QS& qs, gsl_vector* x, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int)
{
  int K = res.size();

  params pa (hpa);
  calc_params(x, pa, hpa, tr);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    tr[i].t = pa.pa[i]->u;

  // calc_t(pa, hpa, tr);
  calc_n(pa, hpa, tr);

  cerr << "calc_dx_u_llik_analytic" << endl;
  cerr << "params" << endl;
  write_params((ofstream&)cerr, pa, hpa, tr);
  cerr << "t_n" << endl;
  write_t_n(cerr, tr, hpa);
  
  params grad (hpa);

  double llik = d_llik(res, qs, pa, grad, hpa, tr, gegen, gegen_int, d_t_variant_fraction);

  return (Log(calc_dx_sigmoid(gsl_vector_get(x, i))) * grad.pa[i]->u).eval();
}

void gsl_set_random(gsl_vector* x, hyperparams& hpa, gsl_rng* rng)
{
  for (int i=0; i<1024; ++i) // for appropriet random number generation
    gsl_rng_uniform(rng);

  for (int i=0; i<=2*hpa.MAX_SUBTYPE; ++i)
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
  
  if (argc != 8)
    {
      cerr << "usage: ./alpha_beta_map_same_cn_compare_gradient max_subtype n step (reads) (u diff outfile) topology u_index" << endl;
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
  
  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  // T = gsl_rng_default; // default is mt19937
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, 1);

  ifstream f (argv[4]);
  ofstream ff (argv[5]);
  // ofstream g (argv[7]);
  // ofstream h (argv[8]);

  int a = atoi(argv[6]);
  int index = atoi(argv[7]);

  ff << scientific;
  // g << scientific;
  // h << scientific;

  VVLog gegen;
  gegen.assign(FRACTIONS+1, VLog(GEGEN_MAX+1, Log(0)));

  set_gegen(gegen);

  VLog gegen_int (FRACTIONS+1, Log(0));
  VLog gegen_int_err (FRACTIONS+1, Log(0));

  set_gegen_integral(gegen_int, gegen_int_err);

  gsl_vector* x = gsl_vector_alloc(2*hpa.MAX_SUBTYPE + 1);

  for (int k=0; k<n; ++k)
    {
      READ *re = new READ;
      int *q = new int;
      f >> re->first >> re->second >> *q;
      res.push_back(re);
      qs.push_back(q);
    }
  
  for (int i=-step; i<=step; ++i)
    {
      gsl_set_random(x, hpa, r);
      gsl_vector_set(x, index, 1.0 * ((double)i) / step);
      double num = calc_dx_u_llik_numeric(index, res, qs, x, hpa, trs[a], gegen, gegen_int);
      double analytic = calc_dx_u_llik_analytic(index, res, qs, x, hpa, trs[a], gegen, gegen_int);

      ff << gsl_vector_get(x, index) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
      // cerr << "-----------------------------------------------------------------------------------------------------------------" << endl;
    }

  for (int i=0; i<n; ++i)
    {
      delete res[i];
      delete qs[i];
    }

  f.close();
  // g.close();
  // h.close();
  ff.close();

  return 0;
}
