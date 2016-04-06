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

int main(int argc, char** argv)
{
  cout << scientific;
  cerr << scientific;
  
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 7)
    {
      cerr << "usage: ./lda_mix_grad max_subtype total_cn n (reads) (n pi kappa infile) (llik outfile)" << endl;
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
  ifstream g (argv[5]);
  ofstream h (argv[6]);
  
  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      f >> re->first >> re->second;
      di.res.push_back(re);
    }

  params pa;
  init_params(pa, di.hpa);

  for (int i=0; i<=di.hpa.MAX_SUBTYPE; ++i)
    g >> pa[i]->n;
  
  for (int i=1; i<=di.hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=di.hpa.TOTAL_CN; ++l)
        {
          g >> pa[i]->pi[l];
        }

      for (int l=1; l<=di.hpa.TOTAL_CN; ++l)
        {
          for (int r=1; r<=l; ++r)
            {
              g >> pa[i]->kappa[l][r];
            }
        }
    }

  h << scientific;
  h << calc_llik(di.res, pa, di.hpa) << endl;
  
  for (int i=0; i<n; ++i)
    delete di.res[i];

  delete_params(pa, di.hpa);
  
  f.close();
  g.close();
  h.close();

  return 0;
}
