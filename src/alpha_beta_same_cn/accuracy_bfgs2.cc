#include <fstream>
#include <algorithm>
#include <numeric>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf.h>
#include <iomanip>
#include "../../util/enumtree_wf_n.hh"
#include "setting.hh"
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

class diff
{
public:
  READS& res;
  QS& qs;
  hyperparams& hpa;
  subtypes& tr;
  VVLog& gegen;
  VLog& gegen_int;

  diff (READS& _res, QS& _qs, hyperparams& _hpa, subtypes& _tr, VVLog& _gegen, VLog& _gegen_int) : res(_res), qs(_qs), hpa(_hpa), tr(_tr), gegen(_gegen), gegen_int(_gegen_int) {}
};

double calc_accuracy(params& pa, params& pa_ans, hyperparams& hpa)
{
  Log accuracy = Log(0);
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      Log diff = pa.pa[i]->u - pa_ans.pa[i]->u;
      accuracy += diff * diff;
    }

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      Log diff = pa.pa[i]->n - pa_ans.pa[i]->n;
      accuracy += diff * diff;
    }

  return accuracy.eval();
}

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

void write_params(std::ofstream& f, params& pa, hyperparams& hpa)
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

void copy_params(params& pa, params& target, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      target.pa[i]->u = pa.pa[i]->u;
    }
  
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      target.pa[i]->n = pa.pa[i]->n;
    }
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

void calc_params(const gsl_vector* x, params& pa, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->u = Log(calc_sigmoid(gsl_vector_get(x, i-1)));
    }

  Log sum = Log(0);
  
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->n = Log(gsl_vector_get(x, hpa.MAX_SUBTYPE + i)).take_exp();
      sum += pa.pa[i]->n;
    }

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->n /= sum;
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

double d_llik(READS& res, QS& qs, params& pa, params& grad, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int)
{
  int K;
  K = res.size();

  VVLog vf (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0)));
  VVLog dtvf (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0)));
  VVLog dnvf (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0)));
  
  for (int q=1; q<=hpa.MAX_SUBTYPE; ++q)
    {
      d_variant_fraction_all(d_t_variant_fraction, 0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q], dtvf[q]);
      d_variant_fraction_all(d_n_variant_fraction, 0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q], dnvf[q]);
    }

  Log lik = Log(1);

  VLog d_t (hpa.MAX_SUBTYPE + 1, Log(0));
  
  for (int k=0; k<K; ++k)
    {
      Log lik_k = Log(0);
      Log d_t_k = Log(0);
      VLog d_n_k (hpa.MAX_SUBTYPE + 1, Log(0));

      for (int s=1; s<=FRACTIONS; ++s)
        {
          tr[*qs[k]].x = Log(((double) s) / FRACTIONS);
          double mu = calc_mu(tr, hpa, *qs[k]);
          lik_k += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * vf[*qs[k]][s];
          d_t_k += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * dtvf[*qs[k]][s];

          for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
            d_n_k[i] += d_mu_n(tr, hpa, i, *qs[k]) * d_bin_mu(*res[k], mu) * vf[*qs[k]][s];

          d_n_k[*qs[k]] += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * dnvf[*qs[k]][s];
        }
      tr[*qs[k]].x = Log(0);

      d_t[*qs[k]] += d_t_k / lik_k;
      for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
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

double my_f (const gsl_vector *v, void *par)
{
  diff *p = (diff *) par;

  params pa (p->hpa);
  calc_params(v, pa, p->hpa);

  calc_t(pa, p->hpa, p->tr);
  calc_n(pa, p->hpa, p->tr);
  
  double llik = calc_llik(p->res, p->qs, pa, p->hpa, p->tr, p->gegen, p->gegen_int);

  return -llik;
}

Log d_n_s(params& pa, int i, int j)
{
  if (i == j)
    return pa.pa[i]->n * (Log(1.0) - pa.pa[i]->n);

  return -pa.pa[i]->n * pa.pa[j]->n;
}

/* The gradient of f, df = (df/dx). */
void my_df (const gsl_vector *v, void *par, gsl_vector *df)
{
  diff *p = (diff *) par;

  int K = p->res.size();
  
  params pa (p->hpa);
  calc_params(v, pa, p->hpa);

  calc_t(pa, p->hpa, p->tr);
  calc_n(pa, p->hpa, p->tr);

  params grad (p->hpa);

  double llik = d_llik(p->res, p->qs, pa, grad, p->hpa, p->tr, p->gegen, p->gegen_int);

  for (int i=1; i<=p->hpa.MAX_SUBTYPE; ++i)
    {
      double gr = (Log(calc_dx_sigmoid(gsl_vector_get(v, i-1))) * grad.pa[i]->u).eval();
      gsl_vector_set(df, i-1, -gr);
    }

  for (int index=0; index<=p->hpa.MAX_SUBTYPE; ++index)
    {
      Log gr = Log(0);
      
      for (int i=0; i<=p->hpa.MAX_SUBTYPE; ++i)
        {
          gr += d_n_s(pa, index, i) * grad.pa[i]->n;
        }
      gsl_vector_set(df, p->hpa.MAX_SUBTYPE + index, (-gr).eval());
    }
}

/* Compute both f and df together. */
void my_fdf (const gsl_vector *x, void *par, double *f, gsl_vector *df)
{
  *f = my_f(x, par);

  my_df(x, par, df);
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

double minimize(diff& di, params& pa, ofstream& h, gsl_rng* rng)
{
  size_t iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  my_func.n = 2*di.hpa.MAX_SUBTYPE;
  my_func.f = my_f;
  my_func.df = my_df;
  my_func.fdf = my_fdf;
  my_func.params = &di;

  x = gsl_vector_alloc (my_func.n);
  gsl_set_random(x, di.hpa, rng);

  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  s = gsl_multimin_fdfminimizer_alloc (T, my_func.n);

  gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.001, 1e-8);

  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

      if (status)
        {
          if (status == GSL_ENOPROG)
            {
              cout << "No more improvement can be made for current estimate" << endl << endl;
              calc_params(s->x, pa, di.hpa);
            }
          break;
        }
      
      calc_params(s->x, pa, di.hpa);
      write_params((ofstream&)cout, pa, di.hpa);
      cout << -s->f << endl << endl;
      h << iter << "\t" << -s->f << endl;
      
      status = gsl_multimin_test_gradient (s->gradient, 1e-7);
      if (status == GSL_SUCCESS)
        {
          //cerr << "Minimum found at:" << endl;
        }

      //cerr << "iter: " << iter << endl;
      //cerr << "kappa\tgrad: " << endl;
      //for (int r=1; r<=st.TOTAL_CN; ++r)
        //cerr << kappa[r] << "\t" << -gsl_vector_get(s->gradient, r-1) << endl;

      //cerr << "llik: " << -s->f << endl << endl;
    }
  while (status == GSL_CONTINUE && iter < 1000);

  double llik = -s->f;
  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);

  return llik;
}

int main(int argc, char** argv)
{
  cout << scientific;
  cerr << scientific;
  
  // feenableexcept(FE_INVALID);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  gsl_set_error_handler_off ();
  
  if (argc != 10)
    {
      cerr << "usage: ./accuracy subtype topology n (reads infile) (u n outfile) (llik outfile) (llik final outfile) (u n infile) (accuracy outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  READS res;
  QS qs;
  int n, MAX_SUBTYPE, TOTAL_CN, MAX_TREE;

  MAX_SUBTYPE = atoi(argv[1]);
  TOTAL_CN = 2;
  
  trees trs;
  trees_cons(trs, MAX_SUBTYPE);
  MAX_TREE = trs.size();

  hyperparams hpa (MAX_SUBTYPE, TOTAL_CN, MAX_TREE);

  int topology = atoi(argv[2]);
  
  n = atoi(argv[3]);

  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  // T = gsl_rng_default; // default is mt19937
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, 1);

  VVLog gegen;
  gegen.assign(FRACTIONS+1, VLog(GEGEN_MAX+1, Log(0)));

  set_gegen(gegen);

  VLog gegen_int (FRACTIONS+1, Log(0));
  VLog gegen_int_err (FRACTIONS+1, Log(0));

  set_gegen_integral(gegen_int, gegen_int_err);

  ifstream f (argv[4]);
  ofstream g (argv[5]);
  ofstream h (argv[6]);
  ofstream hh (argv[7]);
  ifstream u_n_f (argv[8]);
  ofstream accuracy_f (argv[9]);
           
  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      int *q = new int;
      f >> re->first >> re->second >> *q;
      res.push_back(re);
      qs.push_back(q);
    }
  
  params pa_new (hpa);
  
  g << scientific << setprecision(10);
  h << scientific << setprecision(10);
  hh << scientific << setprecision(10);
  accuracy_f << scientific << setprecision(10);

  params pa_ans(hpa);
  read_params(u_n_f, pa_ans, hpa);
  
  diff di (res, qs, hpa, trs[topology], gegen, gegen_int);
  params pa_best (hpa);
  double llik_best = -DBL_MAX;
  for (int i=0; i<10; ++i)
    {
      cout << "minimize_iter: " << i << endl << endl;
      h << "minimize_iter: " << i << endl << endl;
      hh << "minimize_iter: " << i << endl << endl;
      double llik = minimize(di, pa_new, h, r);
      if (llik > llik_best)
        {
          copy_params(pa_new, pa_best, hpa);
          llik_best = llik;
        }
      hh << "llik: " << llik << endl;
      write_params(hh, pa_new, hpa);
    }
  write_params(g, pa_best, hpa);
  accuracy_f << res.size() << "\t" << calc_accuracy(pa_best, pa_ans, hpa) << endl;
  
  for (int i=0; i<n; ++i)
    {
      delete res[i];
      delete qs[i];
    }
  
  f.close();
  g.close();
  h.close();
  hh.close();
  u_n_f.close();
  accuracy_f.close();

  return 0;
}
