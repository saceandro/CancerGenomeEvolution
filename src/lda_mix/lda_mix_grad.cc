#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <fenv.h>
#include <float.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
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

class diff
{
public:
  READS res;
  hyperparams hpa;
};

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

void copy_params(params& pa, params& target, hyperparams& hpa)
{
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
  hpa.gamma.assign(hpa.MAX_SUBTYPE + 1, 1);
  hpa.alpha.assign(hpa.TOTAL_CN + 1, 1);
  hpa.beta.assign(hpa.TOTAL_CN + 1, Vdouble (hpa.TOTAL_CN + 1, 1));
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
  if (re.first < re.second)
    {
      double b = log_binomial_pdf(re.first, mu, re.second-1);
      return log(re.second) + logsub(a, b, sign);
    }
  else
    {
      return log(re.second) + a;
    }
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

  double sum_gamma = sum_vector(p->hpa.gamma, 0, p->hpa.MAX_SUBTYPE);
  for (int j=0; j<=p->hpa.MAX_SUBTYPE; ++j)
    {
      double grad = grad_by_param[j]->n - pa[j]->n * (sum_gamma - p->hpa.MAX_SUBTYPE - 1);
      gsl_vector_set(df, j, -grad);
    }
  
  for (int i=1; i<=p->hpa.MAX_SUBTYPE; ++i)
    {
      for (int s=1; s<=p->hpa.TOTAL_CN; ++s)
        {
          double sum_alpha = sum_vector(p->hpa.alpha, 1, p->hpa.TOTAL_CN);
          double grad = grad_by_param[i]->pi[s] - pa[i]->pi[s] * (sum_alpha - p->hpa.TOTAL_CN + K);
          
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
              double grad = grad_by_param[i]->kappa[l][t] - pa[i]->kappa[l][t] * (sum_beta - l + grad_by_param[i]->pi[l] - (p->hpa.alpha[l] - 1.0));

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

void gsl_set_random(gsl_vector* x, hyperparams& hpa)
{
  for (int i=0; i<hpa.MAX_SUBTYPE + 1 + hpa.MAX_SUBTYPE * hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2; ++i)
    {
      double a = genrand64_real2();
      // cout << a << endl;
      gsl_vector_set(x, i, a);
    }
}

void write_all_gradient(const gsl_vector* x, hyperparams& hpa)
{
  params pa;
  init_params(pa, hpa);

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
    cout << pa[i]->n / sum << "\t";
  cout << endl << endl;
  
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
        cout << pa[i]->pi[l] / sum << "\t";
      cout << endl;

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
            cout << pa[i]->kappa[l][r] / sum << "\t";
          cout << endl;
        }
      cout << endl;
    }
}

double minimize(diff& di, params& pa, ofstream& h)
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
  gsl_set_random(x, di.hpa);

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
              cout << "gradient:" << endl;
              write_all_gradient(s->gradient, di.hpa);
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
  
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 8)
    {
      cerr << "usage: ./lda_mix_grad max_subtype total_cn n (reads) (n pi kappa outfile) (llik outfile) (llik final outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  init_genrand64(1000);
  for (int i=0; i<2048; ++i)
    genrand64_real2();
  
  diff di;
  int n;

  di.hpa.MAX_SUBTYPE = atoi(argv[1]);
  di.hpa.TOTAL_CN = atoi(argv[2]);
  init_hyperparams(di.hpa);
  n = atoi(argv[3]);
  
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
      double llik = minimize(di, pa_new, h);
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
