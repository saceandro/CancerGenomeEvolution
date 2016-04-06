#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf.h>
#include "../mt19937-64/mt64.h"
using namespace std;

#define logsum(a, b) (((a) > (b)) ? ((a) + log1p(exp((b)-(a)))) : ((b) + log1p(exp((a)-(b)))))

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
  Vint sign;
  Vdouble resp_dn;
  double resp_num;
  double resp;

  // use default constructor
  // state (int, std::vector<int>, std::vector<int>, double);
};

typedef std::vector<state*> states;

class param
{
public:
  double n;
  int dn_sign;
  bool n_flag;
  Vdouble pi;
  Vbool pi_flag;
  VVdouble kappa;
  VVbool kappa_flag;

  param (double _n, Vdouble _pi, VVdouble _kappa, int _dn_sign, bool _n_flag, Vbool _pi_flag, VVbool _kappa_flag) : n(_n), pi(_pi), kappa(_kappa) , dn_sign(_dn_sign), n_flag(_n_flag), pi_flag(_pi_flag), kappa_flag(_kappa_flag) {}
};

typedef vector<param*> params;

typedef struct _hyperparams
{
  Vdouble gamma;
  Vdouble alpha;
  VVdouble beta;
  int MAX_SUBTYPE;
  int TOTAL_CN;
}
  hyperparams;

typedef struct _dn
{
  READS &res;
  gsl_vector *x;
  int i;
  hyperparams &hpa;
  
  _dn (READS& _res, gsl_vector* _x, int _i, hyperparams& _hpa) : res(_res), x(_x), i(_i), hpa(_hpa) {}
}
  dn;

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

double logsub(double a, double b, int& sign)
{
  if (a < b)
    {
      sign *= -1;
      return b + log1p(-exp(a-b));
    }
  else
    {
      sign *= 1;
      return a + log1p(-exp(b-a));
    }
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
    f << pa[i]->n << "\t";
  f << endl << endl;
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double sum = 0;
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        f << pa[i]->pi[l] << "\t";
      f << endl;

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (int r=1; r<=l; ++r)
            f << pa[i]->kappa[l][r] << "\t";
          f << endl;
        }
      f << endl;
    }
}

void init_state(state& st, hyperparams& hpa)
{
  st.total_cn.assign(hpa.MAX_SUBTYPE + 1, 0);
  st.variant_cn.assign(hpa.MAX_SUBTYPE + 1, 0);
  st.sign.assign(hpa.MAX_SUBTYPE + 1, 1);
  st.resp_dn.assign(hpa.MAX_SUBTYPE + 1, 0);
  st.total_cn[0] = 2;
  st.variant_cn[0] = 0;
}

void init_params(params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    pa.push_back(new param (0, Vdouble (hpa.TOTAL_CN + 1, 0), VVdouble (hpa.TOTAL_CN + 1, Vdouble (hpa.TOTAL_CN + 1, 0)), 1, false, Vbool (hpa.TOTAL_CN + 1, false), VVbool (hpa.TOTAL_CN + 1, Vbool (hpa.TOTAL_CN + 1, false))));
}

void init_hyperparams(hyperparams& hpa)
{
  hpa.gamma.assign(hpa.MAX_SUBTYPE + 1, 0.05);
  hpa.alpha.assign(hpa.TOTAL_CN + 1, 0.05);
  hpa.beta.assign(hpa.TOTAL_CN + 1, Vdouble (hpa.TOTAL_CN + 1, 0.05));
}

void delete_params(params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    delete pa[i];
}

double calc_mu(state& st, params& pa, hyperparams& hpa)
{
  double denom = 0;
  double num = 0;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      denom += pa[i]->n * st.total_cn[i];
      num += pa[i]->n * st.variant_cn[i];
    }
    
  return num / denom;
}

#define log_binomial_pdf(m, mu, M) ((m) * log((mu)) + ((M) - (m)) * log1p(-(mu)) + gsl_sf_lnchoose((M), (m)))

double log_d_bin_mu(READ& re, double mu, int& sign)
{
  double a = log_binomial_pdf(re.first-1, mu, re.second-1);
  double b = log_binomial_pdf(re.first, mu, re.second-1);

  return log(re.second) + logsub(a, b, sign);
}

double log_d_mu_n(state& st, params& pa, hyperparams& hpa, int j)
{
  double normal = 0;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    normal += pa[i]->n * st.total_cn[i];

  double variant = 0;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    variant += pa[i]->n * st.variant_cn[i];

  double a = log(st.variant_cn[j]/normal);
  double b = log(st.total_cn[j]/normal * variant/normal);

  return logsub(a, b, st.sign[j]);
}

void responsibility_numerator(READ& re, states& sts, state& st, params& pa, hyperparams& hpa)
{
  double log_product = 0;
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    log_product += log(pa[i]->pi[st.total_cn[i]] * pa[i]->kappa[st.total_cn[i]][st.variant_cn[i]]);

  state* new_st = new state;
  init_state(*new_st, hpa);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      new_st->total_cn[i] = st.total_cn[i];
      new_st->variant_cn[i] = st.variant_cn[i];
    }

  double mu = calc_mu(st, pa, hpa);
  new_st->resp_num = log_product + log_binomial_pdf(re.first, mu, re.second);

  int sign = 1;
  log_product += log_d_bin_mu(re, mu, sign);
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      new_st->resp_dn[i] = log_product + log_d_mu_n(*new_st, pa, hpa, i);
      new_st->sign[i] *= sign;
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

double responsibility_partition(READ& re, states& sts, params& pa, hyperparams& hpa)
{
  states::iterator it = sts.begin();
  double log_partition = (*it)->resp_num;

  for (++it; it != sts.end(); ++it)
    log_partition = logsum(log_partition, (*it)->resp_num);

  for (it = sts.begin(); it != sts.end(); ++it)
    {
      (*it)->resp = (*it)->resp_num - log_partition;
      for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
        {
          (*it)->resp_dn[i] -= log_partition;
        }
    }

  return log_partition;
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
  double m = gsl_vector_get(x, 0);
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double s = gsl_vector_get(x, i);
      if (m < s) m = s;
    }

  double sum = 0;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa[i]->n = exp(gsl_vector_get(x, i) - m);
      sum += pa[i]->n;
    }
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    pa[i]->n /= sum;
  
  int params_per_subtype = hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2;
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      m = gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1));
      
      for (int l=2; l<=hpa.TOTAL_CN; ++l)
        {
          double s = gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + l-1);
          if (m < s) m = s;
        }
  
      sum = 0;
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          pa[i]->pi[l] = exp(gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + l-1) - m);
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
      
          sum = 0;
      
          for (int r=1; r<=l; ++r)
            {
              pa[i]->kappa[l][r] = exp(gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + hpa.TOTAL_CN + (l-1) * l / 2 + r-1) - m);
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
  
  double llik = 0;
  states sts;
  state st;
  init_state(st, hpa);
  
  for (int k=0; k<res.size(); ++k)
    {
      // cerr << "k: " << k << endl;
      
      responsibility_numerator_all(*res[k], sts, st, pa, hpa, 1);
      double llik_k = responsibility_partition(*res[k], sts, pa, hpa);
      // cout << llik_k << endl;
      llik += llik_k;

      delete_states(sts);
      sts.clear();
    }

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      llik += (hpa.gamma[i] - 1.0) * log(pa[i]->n);
    }
      
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          llik += (hpa.alpha[l] - 1.0) * log(pa[i]->pi[l]);
          
          for (int r=1; r<=l; ++r)
            {
              llik += (hpa.beta[l][r] - 1.0) * log(pa[i]->kappa[l][r]);
            }
        }
    }

  return llik;
}

double calc_llik_for_dn(double x_i, void* _dn)
{
  dn* p = (dn*) _dn;

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
  
  double llik = 0;
  states sts;
  state st;
  init_state(st, hpa);
  
  for (int k=0; k<res.size(); ++k)
    {
      // cerr << "k: " << k << endl;
      
      responsibility_numerator_all(*res[k], sts, st, pa, hpa, 1);
      llik += responsibility_partition(*res[k], sts, pa, hpa);

      for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
        {
          for (states::iterator it = sts.begin(); it != sts.end(); ++it)
            {
              if (!grad_by_param[i]->n_flag)
                {
                  grad_by_param[i]->n = (*it)->resp_dn[i];
                  grad_by_param[i]->dn_sign = (*it)->sign[i];
                  grad_by_param[i]->n_flag = true;
                }
              
              else
                {
                  if (grad_by_param[i]->dn_sign == 1)
                    {
                      if ((*it)->sign[i] == 1)
                        {
                          grad_by_param[i]->n = logsum(grad_by_param[i]->n, (*it)->resp_dn[i]);
                        }
                      else
                        {
                          grad_by_param[i]->n = logsub(grad_by_param[i]->n, (*it)->resp_dn[i], grad_by_param[i]->dn_sign);
                        }
                    }
              
                  else
                    {
                      if ((*it)->sign[i] == -1)
                        {
                          grad_by_param[i]->n = logsum(grad_by_param[i]->n, (*it)->resp_dn[i]);
                        }
                      else
                        {
                          grad_by_param[i]->n = logsub(grad_by_param[i]->n, (*it)->resp_dn[i], grad_by_param[i]->dn_sign);
                        }
                    }
                }
            }
        }
      
      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        {
          for (int l=1; l<=hpa.TOTAL_CN; ++l)
            {
              for (states::iterator it = sts.begin(); it != sts.end(); ++it)
                {
                  if ((*it)->total_cn[i] == l)
                    {
                      if (!grad_by_param[i]->pi_flag[l])
                        {
                          grad_by_param[i]->pi[l] = (*it)->resp;
                          grad_by_param[i]->pi_flag[l] = true;
                        }
                      else
                        {
                          grad_by_param[i]->pi[l] = logsum(grad_by_param[i]->pi[l], (*it)->resp);
                        }
                      
                    }
                }
              
              for (int r=1; r<=l; ++r)
                {
                  for (states::iterator it = sts.begin(); it != sts.end(); ++it)
                    {
                      if ((*it)->total_cn[i] == l && (*it)->variant_cn[i] == r)
                        {
                          if (!grad_by_param[i]->kappa_flag[l][r])
                            {
                              grad_by_param[i]->kappa[l][r] = (*it)->resp;
                              grad_by_param[i]->kappa_flag[l][r] = true;
                            }
                          else
                            {
                              grad_by_param[i]->kappa[l][r] = logsum(grad_by_param[i]->kappa[l][r], (*it)->resp);
                            }
                        }
                    }
                }
            }
        }
      delete_states(sts);
      sts.clear();
    }
  
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      grad_by_param[i]->n += log(pa[i]->n);

      grad_by_param[i]->n = hpa.gamma[i] - 1.0 + grad_by_param[i]->dn_sign * exp(grad_by_param[i]->n);
      
      // if (hpa.gamma[i] < 1)
      //   {
      //     int sign;
      //     grad_by_param[i]->n = logsub(grad_by_param[i]->n, log(1.0 - hpa.gamma[i]), sign);
      //     st->flag[i] *= sign;
      //   }
      // else
      //   {
      //     grad_by_param[i]->n = logsum(grad_by_param[i]->n, log(hpa.gamma[i] - 1.0));
      //   }
    }
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          grad_by_param[i]->pi[l] = hpa.alpha[l] - 1.0 + exp(grad_by_param[i]->pi[l]);

          for (int r=1; r<=l; ++r)
            {
              grad_by_param[i]->kappa[l][r] = hpa.beta[l][r] - 1.0 + exp(grad_by_param[i]->kappa[l][r]);
            }
        }
    }

  return llik;
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

double calc_dx_n_llik_numeric(READS& res, gsl_vector* x, int i, hyperparams& hpa)
{
  gsl_function F;
  dn _dn (res, x, i, hpa);
  _dn.res = res;
  _dn.x = x;
  _dn.i = i;
  _dn.hpa = hpa;

  F.function = &calc_llik_for_dn;
  F.params = &_dn;

  double result, abserr;
  gsl_deriv_central(&F, gsl_vector_get(x, i), 1e-5, &result, &abserr);
  
  return result;
}

double calc_dx_n_llik_analytic(READS& res, gsl_vector* x, int j, hyperparams& hpa)
{
  int K = res.size();
  
  params pa;
  init_params(pa, hpa);
  calc_params(x, pa, hpa);

  params grad_by_param;
  init_params(grad_by_param, hpa);

  double llik = d_llik(res, pa, grad_by_param, hpa);
  double sum_gamma = sum_vector(hpa.gamma, 0, hpa.MAX_SUBTYPE);
  return grad_by_param[j]->n - pa[j]->n * (sum_gamma - hpa.MAX_SUBTYPE - 1);

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
  return grad_by_param[i]->pi[s] - pa[i]->pi[s] * (sum_alpha - hpa.TOTAL_CN + K);

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
  return grad_by_param[i]->kappa[l][t] - pa[i]->kappa[l][t] * (sum_beta - l + grad_by_param[i]->pi[l] - (hpa.alpha[l] - 1.0));

  // double grad = 0;
  // for (int r=1; r<=l; ++r)
  //   {
  //     grad += dx_vec(pa[i]->kappa[l], r, t) * grad_by_param[i]->kappa[l][r] / pa[i]->kappa[l][r];
  //   }
  
  // return grad;
}

void gsl_set_random(gsl_vector* x, hyperparams& hpa)
{
  for (int i=0; i<hpa.MAX_SUBTYPE + 1 + hpa.MAX_SUBTYPE * hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2; ++i)
    {
      double a = genrand64_real2();
      // cout << a << endl;
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
      cerr << "usage: ./lda_mix_compare_gradient max_subtype total_cn n step (reads) (n diff outfile) (pi diff outfile) (kappa diff outfile)" << endl;
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
  
  ifstream f (argv[5]);
  ofstream ff (argv[6]);
  ofstream g (argv[7]);
  ofstream h (argv[8]);

  ff << scientific;
  g << scientific;
  h << scientific;

  gsl_vector* x = gsl_vector_alloc(hpa.MAX_SUBTYPE + 1 + hpa.MAX_SUBTYPE * hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2);

  init_genrand64(1);

  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      f >> re->first >> re->second;
      res.push_back(re);
    }

  for (int i=-step; i<=step; ++i)
    {
      gsl_set_random(x, hpa);
      gsl_vector_set(x, 0, 1.0 * ((double)i) / step);
      double num = calc_dx_n_llik_numeric(res, x, 1, hpa);
      double analytic = calc_dx_n_llik_analytic(res, x, 1, hpa);
      
      if (fabs(num) > 0)
        ff << gsl_vector_get(x, 0) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
      else
        ff << gsl_vector_get(x, 0) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
    }

  for (int i=-step; i<=step; ++i)
    {
      gsl_set_random(x, hpa);
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
      gsl_set_random(x, hpa);
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
