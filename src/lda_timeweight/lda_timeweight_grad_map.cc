#include <fstream>
#include <algorithm>
#include <numeric>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf.h>
#include "../loglib.hh"
using namespace std;

#define calc_sigmoid(x) ((tanh((x)/2.0) + 1.0) / 2.0)

typedef vector<int> Vint;
typedef vector<double> Vdouble;
typedef vector<Vdouble> VVdouble;
typedef pair<int, int> READ;
typedef vector<READ*> READS;

class state 
{
public:
  int k;
  Vint total_cn;
  Vint variant_cn;
  VLog resp_du;
  Log resp_num;
  Log resp;

  // use default constructor
  // state (int, std::vector<int>, std::vector<int>, double);
};

typedef std::vector<state*> states;

class param
{
public:
  Log u;
  Log t;
  Log n;
  VLog pi;
  VVLog kappa;

  param (Log _u, Log _t, Log _n, VLog _pi, VVLog _kappa) : u(_u), t(_t), n(_n), pi(_pi), kappa(_kappa) {}
};

typedef vector<param*> params;

typedef struct _hyperparams
{
  pair<double,double> be_hpa;
  Vdouble alpha;
  VVdouble beta;
  int MAX_SUBTYPE;
  int TOTAL_CN;
}
  hyperparams;

class diff
{
public:
  READS res;
  hyperparams hpa;
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
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    f << pa[i]->u.eval() << "\t";
  f << endl << endl;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    f << pa[i]->t.eval() << "\t";
  f << endl << endl;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    f << pa[i]->n.eval() << "\t";
  f << endl << endl;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double sum = 0;
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        f << pa[i]->pi[l].eval() << "\t";
      f << endl;

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (int r=1; r<=l; ++r)
            f << pa[i]->kappa[l][r].eval() << "\t";
          f << endl;
        }
      f << endl;
    }
}

void copy_params(params& pa, params& target, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    target[i]->u = pa[i]->u;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    target[i]->t = pa[i]->t;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    target[i]->n = pa[i]->n;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        target[i]->pi[l] = pa[i]->pi[l];

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (int r=1; r<=l; ++r)
            target[i]->kappa[l][r] = pa[i]->kappa[l][r];
        }
    }
}

void init_state(state& st, hyperparams& hpa)
{
  st.total_cn.assign(hpa.MAX_SUBTYPE + 1, 0);
  st.variant_cn.assign(hpa.MAX_SUBTYPE + 1, 0);
  st.resp_du.assign(hpa.MAX_SUBTYPE + 1, Log(0));
  st.total_cn[0] = 2;
  st.variant_cn[0] = 0;
}

void init_params(params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    pa.push_back(new param (Log(0), Log(0), Log(0), VLog (hpa.TOTAL_CN + 1, Log(0)), VVLog (hpa.TOTAL_CN + 1, VLog (hpa.TOTAL_CN + 1, Log(0)))));
}

void init_hyperparams(hyperparams& hpa)
{
  hpa.be_hpa.first = 0.1;
  hpa.be_hpa.second = 0.1;
  hpa.alpha.assign(hpa.TOTAL_CN + 1, 0.1);
  hpa.beta.assign(hpa.TOTAL_CN + 1, Vdouble (hpa.TOTAL_CN + 1, 0.1));
}

void delete_params(params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    delete pa[i];
}

void calc_t(params& pa, hyperparams& hpa)
{
  pa[0]->t = pa[0]->u;
  if (hpa.MAX_SUBTYPE >= 1)
    pa[1]->t = pa[1]->u;
  
  for (int i=2; i<=hpa.MAX_SUBTYPE; ++i)
    pa[i]->t = pa[i-1]->t * pa[i]->u;
}

void calc_n(params& pa, hyperparams& hpa)
{
  Log sum;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa[i]->n = pa[i]->t;
      sum += pa[i]->t;
    }
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    pa[i]->n /= sum;
}

double calc_mu(state& st, params& pa, hyperparams& hpa)
{
  Log denom;
  Log num;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      denom += pa[i]->n * Log(st.total_cn[i]);
      num += pa[i]->n * Log(st.variant_cn[i]);
    }
    
  return (num / denom).eval();
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

Log d_n_t(params& pa, hyperparams& hpa, int i, int j)
{
  if (i == j)
    return pa[i]->n * (Log(1) - pa[i]->n)/ pa[i]->t;
  else
    return -pa[i]->n * pa[i]->n / pa[i]->t;
}

Log d_t_u(params& pa, hyperparams& hpa, int j, int y)
{
  if (y == 0)
    {
      if (j == y)
        return Log(1.0);
      else
        return Log(0);
    }
  else if (j < y)
    return Log(0);
  else
    return pa[j]->t / pa[y]->u;
}

Log d_mu_n(state& st, params& pa, hyperparams& hpa, int j)
{
  Log normal;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    normal += pa[i]->n * Log(st.total_cn[i]);

  Log variant;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    variant += pa[i]->n * Log(st.variant_cn[i]);

  Log a = Log(st.variant_cn[j]) / normal;
  Log b = Log(st.total_cn[j]) / normal * variant / normal;

  return a - b;
}

void responsibility_numerator(READ& re, states& sts, state& st, params& pa, hyperparams& hpa)
{
  Log product (1);
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    product *= pa[i]->pi[st.total_cn[i]] * pa[i]->kappa[st.total_cn[i]][st.variant_cn[i]];

  state* new_st = new state;
  init_state(*new_st, hpa);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      new_st->total_cn[i] = st.total_cn[i];
      new_st->variant_cn[i] = st.variant_cn[i];
    }

  double mu = calc_mu(st, pa, hpa);
  new_st->resp_num = product * Log(log_binomial_pdf(re.first, mu, re.second), 1);

  product *= d_bin_mu(re, mu);
  for (int y=0; y<=hpa.MAX_SUBTYPE; ++y)
    {
      Log sum;
      // for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
      //   sum += d_mu_n(st, pa, hpa, i) * d_n_u(pa, hpa, i, y);

      for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
        {
          Log sum2;
          for (int j=y; j<=hpa.MAX_SUBTYPE; ++j)
            sum2 += d_n_t(pa, hpa, i, j) * d_t_u(pa, hpa, j, y);
          sum += d_mu_n(st, pa, hpa, i) * sum2;
        }
      
      new_st->resp_du[y] = product * sum;
    }
  
  sts.push_back(new_st);
}

void responsibility_numerator_all(READ& re, states& sts, state& st, params& pa, hyperparams& hpa, int subtype)
{
  if (subtype < hpa.MAX_SUBTYPE)
    {
      for (st.total_cn[subtype] = 1; st.total_cn[subtype] <= hpa.TOTAL_CN; ++st.total_cn[subtype])
        {
          for (st.variant_cn[subtype] = 1; st.variant_cn[subtype] <= st.total_cn[subtype]; ++st.variant_cn[subtype])
            {
              responsibility_numerator_all(re, sts, st, pa, hpa, subtype + 1);
            }
        }
    }

  else
    {
      for (st.total_cn[subtype] = 1; st.total_cn[subtype] <= hpa.TOTAL_CN; ++st.total_cn[subtype])
        {
          for (st.variant_cn[subtype] = 1; st.variant_cn[subtype] <= st.total_cn[subtype]; ++st.variant_cn[subtype])
            {
              responsibility_numerator(re, sts, st, pa, hpa);
            }
        }
    }
}

Log responsibility_partition(READ& re, states& sts, params& pa, hyperparams& hpa)
{
  Log partition;
  
  for (states::iterator it = sts.begin(); it != sts.end(); ++it)
    {
      // cout << (*it)->resp_num.take_log() << "\t" << (*it)->resp_num.get_sign() << endl;
      partition += (*it)->resp_num;
    }

  for (states::iterator it = sts.begin(); it != sts.end(); ++it)
    {
      (*it)->resp = (*it)->resp_num / partition;
      for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
        {
          (*it)->resp_du[i] /= partition;
        }
    }

  return partition;
}

void delete_states(states& sts)
{
  for (int i=0; i<sts.size(); ++i)
    {
      delete sts[i];
    }
}

void calc_params(const gsl_vector* x, params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa[i]->u = Log(calc_sigmoid(gsl_vector_get(x, i)));
      // pa[i]->u = Log(gsl_vector_get(x, i));
    }
  
  calc_t(pa, hpa);
  calc_n(pa, hpa);
  
  int params_per_subtype = hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double m = gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1));
      
      for (int l=2; l<=hpa.TOTAL_CN; ++l)
        {
          double s = gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + l-1);
          if (m < s) m = s;
        }
  
      Log sum;
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          pa[i]->pi[l] = Log(gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + l-1) - m, 1);
          sum += pa[i]->pi[l];
        }

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        pa[i]->pi[l] /= sum;

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          m = gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + hpa.TOTAL_CN + (l-1) * l / 2);
          for (int r=2; r<=l; ++r)
            {
              double s = gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + hpa.TOTAL_CN + (l-1) * l / 2 + r-1);
              if (m < s) m = s;
            }
      
          sum = Log(0);
      
          for (int r=1; r<=l; ++r)
            {
              pa[i]->kappa[l][r] = Log(gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + hpa.TOTAL_CN + (l-1) * l / 2 + r-1) - m, 1);
              sum += pa[i]->kappa[l][r];
            }

          for (int r=1; r<=l; ++r)
            pa[i]->kappa[l][r] /= sum;
        }
    }
}

#define log_beta_pdf(s, a, b) (-gsl_sf_lnbeta((a), (b)) + ((a) - 1.0) * log((s)) + ((b) - 1.0) * log1p(-(s)) )

double calc_llik(READS& res, params& pa, hyperparams& hpa)
{
  int K;
  K = res.size();
  
  Log lik (1);
  states sts;
  state st;
  init_state(st, hpa);
  
  for (int k=0; k<res.size(); ++k)
    {
      // cerr << "k: " << k << endl;
      
      responsibility_numerator_all(*res[k], sts, st, pa, hpa, 1);
      Log lik_k = responsibility_partition(*res[k], sts, pa, hpa);
      // cout << lik_k.take_log() << endl;
      lik *= lik_k;
      
      delete_states(sts);
      sts.clear();
    }

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    lik *= Log(log_beta_pdf(pa[i]->u.eval(), hpa.be_hpa.first, hpa.be_hpa.second), 1);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          lik *= pa[i]->pi[l].take_pow(hpa.alpha[l] - 1.0);
          
          for (int r=1; r<=l; ++r)
            {
              lik *= pa[i]->kappa[l][r].take_pow(hpa.beta[l][r] - 1.0);
            }
        }
    }

  return lik.take_log();
}

double d_llik(READS& res, params& pa, params& grad_by_param, hyperparams& hpa)
{
  int K;
  K = res.size();
  
  Log lik (1);
  states sts;
  state st;
  init_state(st, hpa);
  
  for (int k=0; k<res.size(); ++k)
    {
      // cerr << "k: " << k << endl;
      
      responsibility_numerator_all(*res[k], sts, st, pa, hpa, 1);
      lik *= responsibility_partition(*res[k], sts, pa, hpa);

      for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
        {
          Log sum;
          for (states::iterator it = sts.begin(); it != sts.end(); ++it)
            {
              sum += (*it)->resp_du[i];
            }
          grad_by_param[i]->u += sum;
        }
      
      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        {
          for (int l=1; l<=hpa.TOTAL_CN; ++l)
            {
              for (states::iterator it = sts.begin(); it != sts.end(); ++it)
                {
                  if ((*it)->total_cn[i] == l)
                    {
                      grad_by_param[i]->pi[l] += (*it)->resp;
                    }
                }
              
              for (int r=1; r<=l; ++r)
                {
                  for (states::iterator it = sts.begin(); it != sts.end(); ++it)
                    {
                      if ((*it)->total_cn[i] == l && (*it)->variant_cn[i] == r)
                        {
                          grad_by_param[i]->kappa[l][r] += (*it)->resp;
                        }
                    }
                }
            }
        }
      delete_states(sts);
      sts.clear();
    }

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    grad_by_param[i]->u += Log(hpa.be_hpa.first - 1.0) / pa[i]->u - Log(hpa.be_hpa.second - 1.0) / (Log(1.0) - pa[i]->u);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          grad_by_param[i]->pi[l] = Log(hpa.alpha[l] - 1.0) + grad_by_param[i]->pi[l];

          for (int r=1; r<=l; ++r)
            {
              grad_by_param[i]->kappa[l][r] = Log(hpa.beta[l][r] - 1.0) + grad_by_param[i]->kappa[l][r];
            }
        }
    }

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    lik *= Log(log_beta_pdf(pa[i]->u.eval(), hpa.be_hpa.first, hpa.be_hpa.second), 1);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          lik *= pa[i]->pi[l].take_pow(hpa.alpha[l] - 1.0);
          
          for (int r=1; r<=l; ++r)
            {
              lik *= pa[i]->kappa[l][r].take_pow(hpa.beta[l][r] - 1.0);
            }
        }
    }

  return lik.take_log();
}

double dx_vec(Vdouble& vec, int l, int s)
{
  if (l == s)
    return vec[l] * (1.0 - vec[l]);

  return -vec[l] * vec[s];
}

// double dx_kappa(VVdouble& kappa, int l, int r, int t)
// {
//   return dx_vec(kappa[l], r, t);
// }

double my_f (const gsl_vector *v, void *par)
{
  diff *p = (diff *) par;

  params pa;
  init_params(pa, p->hpa);
  
  calc_params(v, pa, p->hpa);

  double llik = calc_llik(p->res, pa, p->hpa);

  return -llik;
}

/* The gradient of f, df = (df/dx). */
void my_df (const gsl_vector *v, void *par, gsl_vector *df)
{
  diff *p = (diff *) par;

  int params_per_subtype = p->hpa.TOTAL_CN * (p->hpa.TOTAL_CN + 3) / 2;
  int K = p->res.size();
  
  params pa;
  init_params(pa, p->hpa);
  calc_params(v, pa, p->hpa);
  
  params grad_by_param;
  init_params(grad_by_param, p->hpa);

  double llik = d_llik(p->res, pa, grad_by_param, p->hpa);

  for (int y=0; y<=p->hpa.MAX_SUBTYPE; ++y)
    {
      double grad = (Log(calc_dx_sigmoid(gsl_vector_get(v,y))) * grad_by_param[y]->u).eval();;
      gsl_vector_set(df, y, -grad);
    }
  
  for (int i=1; i<=p->hpa.MAX_SUBTYPE; ++i)
    {
      for (int s=1; s<=p->hpa.TOTAL_CN; ++s)
        {
          double sum_alpha = sum_vector(p->hpa.alpha, 1, p->hpa.TOTAL_CN);
          double grad = (grad_by_param[i]->pi[s] - pa[i]->pi[s] * Log(sum_alpha - p->hpa.TOTAL_CN + K)).eval();

          // double grad = 0;
          // for (int l=1; l<=p->hpa.TOTAL_CN; ++l)
          //   {
          //     grad += dx_vec(pa[i]->pi, l, s) * grad_by_param[i]->pi[l] / pa[i]->pi[l];
          //   }

          gsl_vector_set(df, p->hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + s-1, -grad);
        }

      for (int l=1; l<=p->hpa.TOTAL_CN; ++l)
        {
          for (int t=1; t<=l; ++t)
            {
              double sum_beta = sum_vector(p->hpa.beta[l], 1, l);
              double grad = (grad_by_param[i]->kappa[l][t] - pa[i]->kappa[l][t] * ( Log(sum_beta - l - p->hpa.alpha[l] + 1.0) + grad_by_param[i]->pi[l] )).eval();

              // double grad = 0;
              // for (int r=1; r<=l; ++r)
              //   {
              //     grad += dx_vec(pa[i]->kappa[l], r, t) * grad_by_param[i]->kappa[l][r] / pa[i]->kappa[l][r];
              //   }

              gsl_vector_set(df, p->hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + p->hpa.TOTAL_CN + (l-1) * l / 2 + t-1, -grad);
            }
        }
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
    // gsl_ran_beta(r, hpa.be_hpa.first, hpa.be_hpa.second);
  
  // for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
  //   pa[i]->u = gsl_ran_beta(r, hpa.be_hpa.first, hpa.be_hpa.second);
  // pa[0]->u *= 0.1;

  for (int i = 0; i < hpa.MAX_SUBTYPE + 1 + hpa.MAX_SUBTYPE * hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2; ++i)
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

  my_func.n = di.hpa.MAX_SUBTYPE + 1 + di.hpa.TOTAL_CN * (di.hpa.TOTAL_CN + 3) / 2 * di.hpa.MAX_SUBTYPE;
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
              cout << "No more improvement can be made for current estimate" << endl << endl;
              calc_params(s->x, pa, di.hpa);
            }
          break;
        }
      
      calc_params(s->x, pa, di.hpa);
      write_params((ofstream&)cout, pa, di.hpa);
      cout << -s->f << endl << endl;
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
  
  feenableexcept(FE_INVALID);
  
  if (argc != 8)
    {
      cerr << "usage: ./lda_timeweight_grad max_subtype total_cn n (reads) (n pi kappa outfile) (llik outfile) (llik final outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  diff di;
  int n;

  di.hpa.MAX_SUBTYPE = atoi(argv[1]);
  di.hpa.TOTAL_CN = atoi(argv[2]);
  init_hyperparams(di.hpa);
  n = atoi(argv[3]);

  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  // T = gsl_rng_default; // default is mt19937
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, 1);

  ifstream f (argv[4]);
  ofstream g (argv[5]);
  ofstream h (argv[6]);
  ofstream hh (argv[7]);
  
  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      f >> re->first >> re->second;
      di.res.push_back(re);
    }

  params pa_new;
  init_params(pa_new, di.hpa);

  params pa_best;
  init_params(pa_best, di.hpa);
  
  g << scientific;
  h << scientific;
  hh << scientific;

  double llik_best = -DBL_MAX;
  for (int i=0; i<100; ++i)
    {
      cout << "minimize_iter: " << i << endl << endl;
      h << "minimize_iter: " << i << endl << endl;
      double llik = minimize(di, pa_new, h, r);
      hh << i << "\t" << llik << endl;
      if (llik > llik_best)
        copy_params(pa_new, pa_best, di.hpa);
    }
  write_params(g, pa_best, di.hpa);
  
  for (int i=0; i<n; ++i)
    delete di.res[i];

  delete_params(pa_new, di.hpa);
  delete_params(pa_best, di.hpa);
  
  f.close();
  g.close();
  h.close();
  hh.close();

  return 0;
}
