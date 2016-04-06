#include <fstream>
#include <algorithm>
#include <numeric>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf.h>
#include "../loglib.hh"
using namespace std;

#define calc_sigmoid(x) ((tanh((x)/2.0) + 1.0) / 2.0)

typedef vector<int> Vint;
typedef vector<double> Vdouble;
typedef vector<Vdouble> VVdouble;
typedef vector<bool> Vbool;
typedef vector<Vbool> VVbool;
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

typedef struct _du
{
  READS &res;
  gsl_vector *x;
  int i;
  hyperparams &hpa;
  
  _du (READS& _res, gsl_vector* _x, int _i, hyperparams& _hpa) : res(_res), x(_x), i(_i), hpa(_hpa) {}
}
  du;

typedef struct _dpi
{
  READS &res;
  gsl_vector *x;
  int i;
  int l;
  hyperparams &hpa;
  
  _dpi (READS& _res, gsl_vector* _x, int _i, int _l, hyperparams& _hpa) : res(_res), x(_x), i(_i), l(_l), hpa(_hpa) {}
}
  dpi;

typedef struct _dkappa
{
  READS &res;
  gsl_vector *x;
  int i;
  int l;
  int r;
  hyperparams &hpa;
  
  _dkappa (READS& _res, gsl_vector* _x, int _i, int _l, int _r, hyperparams& _hpa) : res(_res), x(_x), i(_i), l(_l), r(_r), hpa(_hpa) {}
}
  dkappa;

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
  hpa.be_hpa.first = 2.0;
  hpa.be_hpa.second = 2.0;
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
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    pa[i]->t = pa[i-1]->t * pa[i]->u;
}

void calc_n(params& pa, hyperparams& hpa)
{
  Log sum;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa[i]->n = Log(pa[i]->t.eval(), 1);
      sum += pa[i]->n;
    }

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa[i]->n /= sum;
    }
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
  Log b (log_binomial_pdf(re.first, mu, re.second-1), 1);

  return Log(re.second) * (a - b);
}

Log d_n_t(params& pa, hyperparams& hpa, int i, int j)
{
  if (i == j)
    return pa[i]->n * (Log(1) - pa[j]->n);
  else
    return -pa[i]->n * pa[j]->n;
}

Log d_t_u(params& pa, hyperparams& hpa, int j, int y)
{
  if (j < y)
    return Log(0);
  else
    return pa[j]->t / pa[y]->u;
}
Log d_n_u(params& pa, hyperparams& hpa, int i, int y)
{
  Log b;
  for (int j=y; j<=hpa.MAX_SUBTYPE; ++j)
    b += pa[j]->t * pa[j]->n;

  if (i < y)
    return -pa[i]->n / pa[y]->u * b;
  else
    return pa[i]->n / pa[y]->u * (pa[i]->t - b);
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
      for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
        sum += d_mu_n(st, pa, hpa, i) * d_n_u(pa, hpa, i, y);

      // for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
      //   {
      //     Log sum2;
      //     for (int j=y; j<=hpa.MAX_SUBTYPE; ++j)
      //       sum2 += d_n_t(pa, hpa, i, j) * d_t_u(pa, hpa, j, y);
      //     sum += d_mu_n(st, pa, hpa, i) * sum2;
      //   }
      
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

// int calc_bin_pi(state& st, hyperparams& hpa, int i, int l)
// {
//   if (st.total_cn[i] == l)
//     return 1;
//   else
//     return 0;
// }

// int calc_bin_kappa(state& st, hyperparams& hpa, int i, int l, int r)
// {
//   if (st.total_cn[i] == l && st.variant_cn[i] == r)
//     return 1;
//   else
//     return 0;
// }
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

double calc_llik_for_du(double x_i, void* _du)
{
  du* p = (du*) _du;

  gsl_vector_set(p->x, p->i, x_i);

  params pa;
  init_params(pa, p->hpa);
  calc_params(p->x, pa, p->hpa);

  return calc_llik(p->res, pa, p->hpa);
}

double calc_llik_for_dpi(double x_il, void* dp)
{
  dpi* p = (dpi*) dp;
  
  int params_per_subtype = p->hpa.TOTAL_CN * (p->hpa.TOTAL_CN + 3) / 2;
  
  gsl_vector_set(p->x, p->hpa.MAX_SUBTYPE + 1 + params_per_subtype * (p->i - 1) + p->l - 1, x_il);
  
  params pa;
  init_params(pa, p->hpa);
  calc_params(p->x, pa, p->hpa);

  return calc_llik(p->res, pa, p->hpa);
}

double calc_llik_for_dkappa(double x_ilr, void* dk)
{
  dkappa* p = (dkappa*) dk;

  int params_per_subtype = p->hpa.TOTAL_CN * (p->hpa.TOTAL_CN + 3) / 2;

  gsl_vector_set(p->x, p->hpa.MAX_SUBTYPE + 1 + params_per_subtype * (p->i - 1) + p->hpa.TOTAL_CN + (p->l - 1) * p->l / 2 + p->r - 1, x_ilr);
  
  params pa;
  init_params(pa, p->hpa);
  calc_params(p->x, pa, p->hpa);

  return calc_llik(p->res, pa, p->hpa);
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


double calc_rel_err(double x, double y)
{
  return fabs(x - y) / fabs(x);
}

double calc_dx_u_llik_numeric(READS& res, gsl_vector* x, int i, hyperparams& hpa)
{
  gsl_function F;
  du _du (res, x, i, hpa);
  _du.res = res;
  _du.x = x;
  _du.i = i;
  _du.hpa = hpa;

  F.function = &calc_llik_for_du;
  F.params = &_du;

  double result, abserr;
  gsl_deriv_central(&F, gsl_vector_get(x, i), 1e-5, &result, &abserr);
  
  return result;
}

double calc_dx_u_llik_analytic(READS& res, gsl_vector* x, int j, hyperparams& hpa)
{
  int K = res.size();
  
  params pa;
  init_params(pa, hpa);
  calc_params(x, pa, hpa);

  params grad_by_param;
  init_params(grad_by_param, hpa);

  double llik = d_llik(res, pa, grad_by_param, hpa);
  // return calc_dx_sigmoid(gsl_vector_get(x,j)) * grad_by_param[j]->u.eval();
  return (Log(calc_dx_sigmoid(gsl_vector_get(x,j))) * grad_by_param[j]->u).eval();
  // return (Log(calc_dx_sigmoid(pa[j]->u.eval())) * grad_by_param[j]->u).eval();
  // return grad_by_param[j]->u.eval();

  // double grad = 0;
  // for (int l=1; l<=hpa.TOTAL_CN; ++l)
  //   {
  //     grad += dx_vec(pa[i]->pi, l, s) * grad_by_param[i]->pi[l] / pa[i]->pi[l];
  //   }
  
  // return grad;
}

double calc_dx_pi_llik_numeric(READS& res, gsl_vector* x, int i, int l, hyperparams& hpa)
{
  gsl_function F;
  dpi dp (res, x, i, l, hpa);
  dp.res = res;
  dp.x = x;
  dp.i = i;
  dp.l = l;
  dp.hpa = hpa;

  F.function = &calc_llik_for_dpi;
  F.params = &dp;

  double result, abserr;
  int params_per_subtype = hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2;
  gsl_deriv_central(&F, gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + l-1), 1e-5, &result, &abserr);

  return result;
}

double calc_dx_pi_llik_analytic(READS& res, gsl_vector* x, int i, int s, hyperparams& hpa)
{
  int K = res.size();
  
  params pa;
  init_params(pa, hpa);
  calc_params(x, pa, hpa);

  params grad_by_param;
  init_params(grad_by_param, hpa);

  double llik = d_llik(res, pa, grad_by_param, hpa);
  double sum_alpha = sum_vector(hpa.alpha, 1, hpa.TOTAL_CN);
  return (grad_by_param[i]->pi[s] - pa[i]->pi[s] * Log(sum_alpha - hpa.TOTAL_CN + K)).eval();

  // double grad = 0;
  // for (int l=1; l<=hpa.TOTAL_CN; ++l)
  //   {
  //     grad += dx_vec(pa[i]->pi, l, s) * grad_by_param[i]->pi[l] / pa[i]->pi[l];
  //   }
  
  // return grad;
}

double calc_dx_kappa_llik_numeric(READS& res, gsl_vector* x, int i, int l, int r, hyperparams& hpa)
{
  gsl_function F;
  dkappa dk (res, x, i, l, r, hpa);
  dk.res = res;
  dk.x = x;
  dk.i = i;
  dk.l = l;
  dk.r = r;
  dk.hpa = hpa;

  F.function = &calc_llik_for_dkappa;
  F.params = &dk;

  double result, abserr;
  int params_per_subtype = hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2;
  gsl_deriv_central(&F, gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + hpa.TOTAL_CN + (l-1) * l / 2 + r-1), 1e-5, &result, &abserr);
  
  return result;
}

double calc_dx_kappa_llik_analytic(READS& res, gsl_vector* x, int i, int l, int t, hyperparams& hpa)
{
  int K = res.size();
  
  params pa;
  init_params(pa, hpa);
  calc_params(x, pa, hpa);

  params grad_by_param;
  init_params(grad_by_param, hpa);

  double llik = d_llik(res, pa, grad_by_param, hpa);

  double sum_beta = sum_vector(hpa.beta[l], 1, l);
  return (grad_by_param[i]->kappa[l][t] - pa[i]->kappa[l][t] * ( Log(sum_beta - l - hpa.alpha[l] + 1.0) + grad_by_param[i]->pi[l] )).eval();

  // double grad = 0;
  // for (int r=1; r<=l; ++r)
  //   {
  //     grad += dx_vec(pa[i]->kappa[l], r, t) * grad_by_param[i]->kappa[l][r] / pa[i]->kappa[l][r];
  //   }
  
  // return grad;
}

void gsl_set_random(gsl_vector* x, hyperparams& hpa, gsl_rng* rng)
{
  for (int i=0; i<1024; ++i) // for appropriet random number generation
    gsl_rng_uniform(rng);
    // gsl_ran_beta(r, hpa.be_hpa.first, hpa.be_hpa.second);
  
  // for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
  //   pa[i]->u = gsl_ran_beta(r, hpa.be_hpa.first, hpa.be_hpa.second);
  // pa[0]->u *= 0.1;

  for (int i = 0; i<hpa.MAX_SUBTYPE + 1 + hpa.MAX_SUBTYPE * hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2; ++i)
    {
      double a = gsl_rng_uniform(rng);
      gsl_vector_set(x, i, a);
    }
}

int main(int argc, char** argv)
{
  cout << scientific;
  cerr << scientific;
  
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 9)
    {
      cerr << "usage: ./lda_timeweight_compare_diff2 max_subtype total_cn n step (reads) (u diff outfile) (pi diff outfile) (kappa diff outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  READS res;
  hyperparams hpa;
  int n, step;

  hpa.MAX_SUBTYPE = atoi(argv[1]);
  hpa.TOTAL_CN = atoi(argv[2]);
  init_hyperparams(hpa);
  n = atoi(argv[3]);
  step = atoi(argv[4]);
  
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
  ofstream h (argv[8]);

  ff << scientific;
  g << scientific;
  h << scientific;

  gsl_vector* x = gsl_vector_alloc(hpa.MAX_SUBTYPE + 1 + hpa.MAX_SUBTYPE * hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2);

  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      f >> re->first >> re->second;
      res.push_back(re);
    }

  for (int i=-step; i<=step; ++i)
    {
      gsl_set_random(x, hpa, r);
      gsl_vector_set(x, 1, 1.0 * ((double)i) / step);
      double num = calc_dx_u_llik_numeric(res, x, 1, hpa);
      double analytic = calc_dx_u_llik_analytic(res, x, 1, hpa);
      
      if (fabs(num) > 0)
        ff << gsl_vector_get(x, 1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
      else
        ff << gsl_vector_get(x, 1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
    }

  for (int i=-step; i<=step; ++i)
    {
      gsl_set_random(x, hpa, r);
      gsl_vector_set(x, hpa.MAX_SUBTYPE + 1, 1.0 * ((double)i) / step);
      double num = calc_dx_pi_llik_numeric(res, x, 1, 1, hpa);
      double analytic = calc_dx_pi_llik_analytic(res, x, 1, 1, hpa);
      
      if (fabs(num) > 0)
        g << gsl_vector_get(x, hpa.MAX_SUBTYPE + 1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
      else
        g << gsl_vector_get(x, hpa.MAX_SUBTYPE + 1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
    }

  for (int i=-step; i<=step; ++i)
    {
      gsl_set_random(x, hpa, r);
      gsl_vector_set(x, hpa.MAX_SUBTYPE + 1 + hpa.MAX_SUBTYPE * hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2 - 1, 1.0 * ((double)i) / step);
      double num = calc_dx_kappa_llik_numeric(res, x, hpa.MAX_SUBTYPE, hpa.TOTAL_CN, hpa.TOTAL_CN, hpa);
      double analytic = calc_dx_kappa_llik_analytic(res, x, hpa.MAX_SUBTYPE, hpa.TOTAL_CN, hpa.TOTAL_CN, hpa);
      
      if (fabs(num) > 0)
        h << gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + hpa.MAX_SUBTYPE * hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2 - 1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
      else
        h << gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + hpa.MAX_SUBTYPE * hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2 - 1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
  }

  for (int i=0; i<n; ++i)
    delete res[i];

  f.close();
  g.close();

  return 0;
}
