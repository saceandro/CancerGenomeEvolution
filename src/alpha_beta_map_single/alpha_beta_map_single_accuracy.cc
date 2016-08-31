#include <fstream>
#include <algorithm>
#include <numeric>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf.h>
#include <iomanip>
#include "../../util/enumtree_wf.hh"
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

class diff
{
public:
  READS& res;
  hyperparams& hpa;
  subtype& st;
  VVLog& gegen;
  VLog& gegen_int;
  int read_num;

  diff (READS& _res, hyperparams& _hpa, subtype& _st, VVLog& _gegen, VLog& _gegen_int, int _read_num) : res(_res), hpa(_hpa), st(_st), gegen(_gegen), gegen_int(_gegen_int), read_num(_read_num) {}
};

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

void write_params(std::ofstream& f, params& pa, hyperparams& hpa)
{
  f << pa.pa[1]->u.eval() << endl;
}

void copy_params(params& pa, params& target, hyperparams& hpa)
{
  target.pa[1]->u = pa.pa[1]->u;
}

double calc_mu(subtype& st, hyperparams& hpa)
{
  return (st.x / Log(2)).eval();
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

void calc_params(const gsl_vector* x, params& pa, hyperparams& hpa)
{
  pa.pa[1]->u = Log(calc_sigmoid(gsl_vector_get(x, 0)));
}

Log deriv_k(READ& re, subtype& st, hyperparams& hpa, VVLog& gegen, VLog& gegen_int, myfunc d_x_variant_fraction, Log& d_t)
{
  Log denom = Log(0);
  // VLog d_t_num = VLog(hpa.MAX_SUBTYPE + 1, Log(0));
  // VLog d_n_num = VLog(hpa.MAX_SUBTYPE + 1, Log(0));
  
  VLog vf (FRACTIONS + 1, Log(0));
  VLog dtvf (FRACTIONS + 1, Log(0));
      
  d_variant_fraction_all(d_t_variant_fraction, 0, 1, st.n, st.t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf, dtvf);

  for (int s=1; s<=FRACTIONS; ++s)
    {
      st.x = Log(((double) s) / FRACTIONS);
      double mu = calc_mu(st, hpa);
      denom += Log(log_binomial_pdf(re.first, mu, re.second), 1) * vf[s];

      d_t += Log(log_binomial_pdf(re.first, mu, re.second), 1) * dtvf[s];
    }

  d_t /= denom;
  
  return denom;
}

double calc_llik(READS& res, params& pa, hyperparams& hpa, subtype& st, VVLog& gegen, VLog& gegen_int, int read_num)
{
  VLog vf (FRACTIONS + 1, Log(0));
  VLog vf_numerator (FRACTIONS + 1, Log(0));
  VLog vf_denominator (FRACTIONS + 1, Log(0));
  Log partition = Log(0);
  
  variant_fraction_partition(0, 1, st.n, st.t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf, vf_numerator, vf_denominator, partition);

  Log lik = Log(1);
  
  for (int k=0; k<read_num; ++k)
    {
      Log lik_k = Log(0);
      for (int s=1; s<=FRACTIONS; ++s)
        {
          st.x = Log(((double) s) / FRACTIONS);
          double mu = calc_mu(st, hpa);
          lik_k += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * vf[s];
        }
      
      lik *= lik_k;
    }

  // for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     lik *= pa.pa[i]->u.take_pow(hpa.be_hpa_u.first - 1.0) * (Log(1) - pa.pa[i]->u).take_pow(hpa.be_hpa_u.second - 1.0);
  //   }

  return lik.take_log();
}

double d_llik(READS& res, params& pa, params& grad, hyperparams& hpa, subtype& st, VVLog& gegen, VLog& gegen_int, myfunc d_x_variant_fraction, int read_num)
{
  VLog vf (FRACTIONS + 1, Log(0));
  VLog dtvf (FRACTIONS + 1, Log(0));
      
  d_variant_fraction_all(d_t_variant_fraction, 0, 1, st.n, st.t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf, dtvf);

  Log lik = Log(1);

  for (int k=0; k<read_num; ++k)
    {
      Log lik_k = Log(0);
      Log d_t = Log(0);

      for (int s=1; s<=FRACTIONS; ++s)
        {
          st.x = Log(((double) s) / FRACTIONS);
          double mu = calc_mu(st, hpa);
          lik_k += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * vf[s];
          d_t += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * dtvf[s];
        }
      d_t /= lik_k;

      lik *= lik_k;

      grad.pa[1]->u += d_t;
    }

  // use MAP
  // for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     grad.pa[i]->u += Log(hpa.be_hpa_u.first - 1.0) / pa.pa[i]->u - Log(hpa.be_hpa_u.second - 1.0) / (Log(1.0) - pa.pa[i]->u);
  //   }

  // for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     for (int j=0; j<(int)tr[i].children.size(); ++j)
  //       {
  //         grad.pa[i]->beta[j] += Log(hpa.be_hpa_beta.first - 1.0) / pa.pa[i]->beta[j] - Log(hpa.be_hpa_beta.second - 1.0) / (Log(1.0) - pa.pa[i]->beta[j]);
  //       }
  //   }
  
  // for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     for (int l=1; l<=hpa.TOTAL_CN; ++l)
  //       {
  //         grad.pa[i]->pi[l] += Log(hpa.alpha[l] - 1.0);

  //         for (int r=1; r<=l; ++r)
  //           {
  //             grad.pa[i]->kappa[l][r] += Log(hpa.beta[l][r] - 1.0);
  //           }
  //       }
  //   }

  // use MAP
  // for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     lik *= pa.pa[i]->u.take_pow(hpa.be_hpa_u.first - 1.0) * (Log(1) - pa.pa[i]->u).take_pow(hpa.be_hpa_u.second - 1.0);

  //     for (int j=0; j<(int)tr[i].children.size(); ++j)
  //       {
  //         lik *= pa.pa[i]->beta[j].take_pow(hpa.be_hpa_beta.first - 1.0) * (Log(1) - pa.pa[i]->beta[j]).take_pow(hpa.be_hpa_beta.second - 1.0);
  //       }
  //   }

  // for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     for (int l=1; l<=hpa.TOTAL_CN; ++l)
  //       {
  //         lik *= pa.pa[i]->pi[l].take_pow(hpa.alpha[l] - 1.0);
          
  //         for (int r=1; r<=l; ++r)
  //           {
  //             lik *= pa.pa[i]->kappa[l][r].take_pow(hpa.beta[l][r] - 1.0);
  //           }
  //       }
  //   }

  return lik.take_log();
}

double my_f (const gsl_vector *v, void *par)
{
  diff *p = (diff *) par;

  params pa (p->hpa);
  calc_params(v, pa, p->hpa);

  p->st.t = pa.pa[1]->u;
  p->st.n = Log(1);
  
  double llik = calc_llik(p->res, pa, p->hpa, p->st, p->gegen, p->gegen_int, p->read_num);

  return -llik;
}

/* The gradient of f, df = (df/dx). */
void my_df (const gsl_vector *v, void *par, gsl_vector *df)
{
  diff *p = (diff *) par;

  params pa (p->hpa);
  calc_params(v, pa, p->hpa);

  p->st.t = pa.pa[1]->u;
  p->st.n = Log(1);

  params grad (p->hpa);

  double llik = d_llik(p->res, pa, grad, p->hpa, p->st, p->gegen, p->gegen_int, d_t_variant_fraction, p->read_num);

  double gr = (Log(calc_dx_sigmoid(gsl_vector_get(v, 0))) * grad.pa[1]->u).eval();
  gsl_vector_set(df, 0, -gr);
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

  double a = gsl_rng_uniform(rng);
  gsl_vector_set(x, 0, a);
}

double minimize(diff& di, params& pa, ofstream& h, gsl_rng* rng)
{
  size_t iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  my_func.n = 1;
  my_func.f = my_f;
  my_func.df = my_df;
  my_func.fdf = my_fdf;
  my_func.params = &di;

  x = gsl_vector_alloc (my_func.n);
  gsl_set_random(x, di.hpa, rng);

  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc (T, my_func.n);

  gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.001, 1e-7);

  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

      if (status)
        {
          if (status == GSL_ENOPROG)
            {
              // cout << "No more improvement can be made for current estimate" << endl << endl;
              calc_params(s->x, pa, di.hpa);
            }
          break;
        }
      
      calc_params(s->x, pa, di.hpa);
      // write_params((ofstream&)cout, pa, di.hpa);
      // cout << -s->f << endl << endl;
      h << iter << "\t" << -s->f << endl;
      
      status = gsl_multimin_test_gradient (s->gradient, 1e-6);
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
  
  if (argc != 8)
    {
      cerr << "usage: ./alpha_beta_map_single_accuracy n (reads) (u outfile) (llik outfile) (llik final outfile) (ans file) (accuracy file)" << endl;
      exit(EXIT_FAILURE);
    }

  READS res;
  int n, MAX_SUBTYPE, TOTAL_CN, MAX_TREE;

  MAX_SUBTYPE = 1;
  TOTAL_CN = 2;
  
  MAX_TREE = 1;
  
  subtype st (1, 2, 1, 0, Log(0), Log(0), Log(1), Log(0), Log(0), NULL, NULL, std::vector<subtype*>(0));;

  hyperparams hpa (MAX_SUBTYPE, TOTAL_CN, MAX_TREE);
  
  n = atoi(argv[1]);

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

  ifstream f (argv[2]);
  ofstream g (argv[3]);
  ofstream h (argv[4]);
  ofstream hh (argv[5]);
  ifstream ans_f (argv[6]);
  ofstream accuracy_f (argv[7]);

  accuracy_f << scientific << setprecision(10);
  
  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      f >> re->first >> re->second;
      res.push_back(re);
    }

  params pa_ans (hpa);
  double pa_ans_double;
  ans_f >> pa_ans_double;
  pa_ans.pa[1]->u = Log(pa_ans_double);
  
  params pa_new (hpa);
  
  g << scientific << setprecision(10);
  h << scientific << setprecision(10);
  hh << scientific << setprecision(10);

  for (int read_num=16; read_num<n; read_num*=2)
    {
      h << "read_num: " << read_num << endl;
      hh << "read_num: " << read_num << endl;
      
      diff di (res, hpa, st, gegen, gegen_int, read_num);
      params pa_best (hpa);
      double llik_best = -DBL_MAX;
      for (int i=0; i<10; ++i)
        {
          // cout << "minimize_iter: " << i << endl << endl;
          h << "minimize_iter: " << i << endl << endl;
          hh << "minimize_iter: " << i << endl << endl;
          double llik = minimize(di, pa_new, h, r);
          if (llik > llik_best)
            {
              copy_params(pa_new, pa_best, hpa);
              llik_best = llik;
            }
          hh << "llik: " << llik << endl;
          hh << "u: " << pa_new.pa[1]->u.eval() << endl;
        }
      g  << read_num << "\t" << pa_new.pa[1]->u.eval() << "\t" << llik_best << endl;
      accuracy_f << read_num << "\t" << fabs((pa_best.pa[1]->u - pa_ans.pa[1]->u).eval()) << endl;
      h << "----------------------------------------------------------------------" << endl;
      hh << "----------------------------------------------------------------------" << endl;
      g << "----------------------------------------------------------------------" << endl;
    }
  
  for (int i=0; i<n; ++i)
    delete res[i];

  f.close();
  g.close();
  h.close();
  hh.close();

  return 0;
}
