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
using namespace std;

typedef vector<double> Vdouble;
typedef vector<Vdouble> VVdouble;
typedef vector<unsigned int> Vuint;
typedef pair<unsigned int, unsigned int> READ;
typedef vector<READ*> READS;

class state 
{
public:
  int k;
  std::vector<int> total_cn;
  std::vector<int> variant_cn;
  double resp_num;
  double resp;

  // use default constructor
  // state (int, std::vector<int>, std::vector<int>, double);
};

typedef std::vector<state*> states;

class params
{
public:
  std::vector<double> pi;
  std::vector<std::vector<double> > kappa;

  // use default constructor
  // params (std::vector<double>, std::vector<std::vector<double> >);
};

typedef struct _hyperparams
{
  unsigned int MAX_SUBTYPE;
  unsigned int TOTAL_CN;
}
  hyperparams;

typedef struct _dpi
{
  READS &res;
  gsl_vector *x;
  unsigned int l;
  hyperparams &hpa;
  
  _dpi (READS& _res, gsl_vector* _x, unsigned int _l, hyperparams& _hpa) : res(_res), x(_x), l(_l), hpa(_hpa) {}
}
  dpi;

typedef struct _dkappa
{
  READS &res;
  gsl_vector *x;
  unsigned int l;
  unsigned int r;
  hyperparams &hpa;
  
  _dkappa (READS& _res, gsl_vector* _x, unsigned int _l, unsigned int _r, hyperparams& _hpa) : res(_res), x(_x), l(_l), r(_r), hpa(_hpa) {}
}
  dkappa;

void write_params(std::ofstream& f, params& pa)
{
  for (int c=1; c<pa.pi.size(); ++c)
    {
      f << pa.pi[c] << "\t";
    }
  f << endl;

  for (int c=1; c<pa.pi.size(); ++c)
    {
      for (int d=1; d<=c; ++d)
        {
          f << pa.kappa[c][d] << "\t";
        }
      f << endl;
    }
}

void init_state(state& st, hyperparams& hpa)
{
  st.total_cn.assign(hpa.MAX_SUBTYPE + 1, 0);
  st.variant_cn.assign(hpa.MAX_SUBTYPE + 1, 0);
}

void init_params(params& pa, hyperparams& hpa)
{
  pa.pi.assign(hpa.TOTAL_CN + 1, 0);
  pa.kappa.assign(hpa.TOTAL_CN + 1, Vdouble (hpa.TOTAL_CN + 1, 0));
}

double calc_mu(state& st, hyperparams& hpa)
{
  unsigned int denom = 0;
  unsigned int num = 0;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      denom += st.total_cn[i];
      num += st.variant_cn[i];
    }
    
  return ((double) num) / denom;
}

void responsibility_numerator(READ& re, states& sts, state& st, params& pa, hyperparams& hpa)
{
  double product = 1;
  
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    product *= pa.pi[st.total_cn[i]] * pa.kappa[st.total_cn[i]][st.variant_cn[i]];

  double mu = calc_mu(st, hpa);
  product *= gsl_ran_binomial_pdf(re.first, mu, re.second);

  state* new_st = new state;
  init_state(*new_st, hpa);

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      new_st->total_cn[i] = st.total_cn[i];
      new_st->variant_cn[i] = st.variant_cn[i];
    }

  // cerr << "resp_num: " << product << endl;
  
  new_st->resp_num = product;
  
  sts.push_back(new_st);
}

void responsibility_numerator_all(READ& re, states& sts, state& st, params& pa, hyperparams& hpa, unsigned int subtype)
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
  double partition = 0;

  for (states::iterator it = sts.begin(); it != sts.end(); ++it)
    partition += (*it)->resp_num;

  for (states::iterator it = sts.begin(); it != sts.end(); ++it)
    (*it)->resp = (*it)->resp_num / partition;

  return partition;
}

void delete_states(states& sts)
{
  for (int i=0; i<sts.size(); ++i)
    {
      delete sts[i];
    }
}

unsigned int calc_bin_pi(state& st, hyperparams& hpa, unsigned int l)
{
  unsigned int bin = 0;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      if (st.total_cn[i] == l)
        bin++;
    }
  return bin;
}

unsigned int calc_bin_kappa(state& st, hyperparams& hpa, unsigned int l, unsigned int r)
{
  unsigned int bin = 0;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      if (st.total_cn[i] == l && st.variant_cn[i] == r)
        bin++;
    }
  return bin;
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
      
      responsibility_numerator_all(*res[k], sts, st, pa, hpa, 0);
      llik += log(responsibility_partition(*res[k], sts, pa, hpa));

      delete_states(sts);
      sts.clear();
    }

  return llik;
}

void calc_params(const gsl_vector* x, params& pa, hyperparams& hpa)
{
  double m = gsl_vector_get(x, 0);
  for (unsigned int l=2; l<=hpa.TOTAL_CN; ++l)
    {
      double s = gsl_vector_get(x, l-1);
      if (m < s) m = s;
    }
  
  double sum = 0;
  for (unsigned int l=1; l<=hpa.TOTAL_CN; ++l)
    {
      pa.pi[l] = exp(gsl_vector_get(x, l-1) - m);
      sum += pa.pi[l];
    }

  for (unsigned int l=1; l<=hpa.TOTAL_CN; ++l)
    pa.pi[l] /= sum;

  for (unsigned int l=1; l<=hpa.TOTAL_CN; ++l)
    {
      m = gsl_vector_get(x, hpa.TOTAL_CN + (l-1) * l / 2);
      for (unsigned int r=2; r<=l; ++r)
        {
          double s = gsl_vector_get(x, hpa.TOTAL_CN + (l-1) * l / 2 + r-1);
          if (m < s) m = s;
        }
      
      sum = 0;
      
      for (unsigned int r=1; r<=l; ++r)
        {
          pa.kappa[l][r] = exp(gsl_vector_get(x, hpa.TOTAL_CN + (l-1) * l / 2 + r-1) - m);
          sum += pa.kappa[l][r];
        }

      for (unsigned int r=1; r<=l; ++r)
        pa.kappa[l][r] /= sum;
    }
}

void init_gsl_vector(gsl_vector *x, hyperparams& hpa)
{
  for (int s=1; s<=hpa.TOTAL_CN; ++s)
    {
      gsl_vector_set(x, s-1, 0);
    }

  for (int l=1; l<=hpa.TOTAL_CN; ++l)
    {
      for (int t=1; t<=l; ++t)
        {
          gsl_vector_set(x, hpa.TOTAL_CN + (l-1) * l / 2 + t-1, 0);
        }
    }
}

double calc_llik_for_dpi(double x_l, void* dp)
{
  dpi* p = (dpi*) dp;
  gsl_vector_set(p->x, p->l -1, x_l);
  
  params pa;
  init_params(pa, p->hpa);
  calc_params(p->x, pa, p->hpa);

  return calc_llik(p->res, pa, p->hpa);
}

double calc_llik_for_dkappa(double x_lr, void* dk)
{
  dkappa* p = (dkappa*) dk;
  gsl_vector_set(p->x, p->hpa.TOTAL_CN + (p->l - 1) * p->l / 2 + p->r - 1, x_lr);
  
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
      
      responsibility_numerator_all(*res[k], sts, st, pa, hpa, 0);
      llik += log(responsibility_partition(*res[k], sts, pa, hpa));

      for (unsigned int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (states::iterator it = sts.begin(); it != sts.end(); ++it)
            {
              int bin_l = calc_bin_pi(**it, hpa, l);
              grad_by_param.pi[l] += (*it)->resp * bin_l;
            }

          for (unsigned int r=1; r<=l; ++r)
            {
              for (states::iterator it = sts.begin(); it != sts.end(); ++it)
                {
                  int bin_lr = calc_bin_kappa(**it, hpa, l, r);
                  grad_by_param.kappa[l][r] += (*it)->resp * bin_lr;
                }
            }
        }
      
      delete_states(sts);
      sts.clear();
    }

  return llik;
}

double dx_vec(Vdouble& vec, unsigned int l, unsigned int s)
{
  if (l == s)
    return vec[l] * (1.0 - vec[l]);

  return -vec[l] * vec[s];
}

// double dx_kappa(VVdouble& kappa, unsigned int l, unsigned int r, unsigned int t)
// {
//   return dx_vec(kappa[l], r, t);
// }


double calc_dx_pi_llik_numeric(READS& res, gsl_vector* x, unsigned int l, hyperparams& hpa)
{
  gsl_function F;
  dpi dp (res, x, l, hpa);
  dp.res = res;
  dp.x = x;
  dp.l = l;
  dp.hpa = hpa;

  F.function = &calc_llik_for_dpi;
  F.params = &dp;

  double result, abserr;
  gsl_deriv_central(&F, gsl_vector_get(x, l-1), 1e-7, &result, &abserr);
  
  return result;
}

double calc_dx_pi_llik_analytic(READS& res, gsl_vector* x, unsigned int s, hyperparams& hpa)
{
  int K = res.size();
  
  params pa;
  init_params(pa, hpa);
  calc_params(x, pa, hpa);

  params grad_by_param;
  init_params(grad_by_param, hpa);

  double llik = d_llik(res, pa, grad_by_param, hpa);

  return grad_by_param.pi[s] - pa.pi[s] * K * (hpa.MAX_SUBTYPE+1);

  // double grad = 0;
  // for (int l=1; l<=hpa.TOTAL_CN; ++l)
  //   {
  //     grad += dx_vec(pa.pi, l, s) * grad_by_param.pi[l] / pa.pi[l];
  //   }
  
  // return grad;
}

double calc_rel_err(double x, double y)
{
  return fabs(x - y) / fabs(x);
}

double calc_dx_kappa_llik_numeric(READS& res, gsl_vector* x, unsigned int l, unsigned int r, hyperparams& hpa)
{
  gsl_function F;
  dkappa dk (res, x, l, r, hpa);
  dk.res = res;
  dk.x = x;
  dk.l = l;
  dk.r = r;
  dk.hpa = hpa;

  F.function = &calc_llik_for_dkappa;
  F.params = &dk;

  double result, abserr;
  gsl_deriv_central(&F, gsl_vector_get(x, hpa.TOTAL_CN + (l-1) * l / 2 + r-1), 1e-7, &result, &abserr);
  
  return result;
}

double calc_dx_kappa_llik_analytic(READS& res, gsl_vector* x, unsigned int l, unsigned int t, hyperparams& hpa)
{
  int K = res.size();
  
  params pa;
  init_params(pa, hpa);
  calc_params(x, pa, hpa);

  params grad_by_param;
  init_params(grad_by_param, hpa);

  double llik = d_llik(res, pa, grad_by_param, hpa);

  return grad_by_param.kappa[l][t] - pa.kappa[l][t] * grad_by_param.pi[l];

  // double grad = 0;
  // for (int r=1; r<=l; ++r)
  //   {
  //     grad += dx_vec(pa.kappa[l], r, t) * grad_by_param.kappa[l][r] / pa.kappa[l][r];
  //   }
  
  // return grad;
}

int main(int argc, char** argv)
{
  cout << scientific;
  // cerr << scientific;
  
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 8)
    {
      cerr << "usage: ./mixture_total_variant_compare_gradient max_subtype total_cn n step (reads) (pi diff outfile) (kappa diff outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  hyperparams hpa;
  READS res;
  unsigned int n;
  int step;
  
  hpa.MAX_SUBTYPE = atoi(argv[1]);
  hpa.TOTAL_CN = atoi(argv[2]);
  n = atoi(argv[3]);
  step = atoi(argv[4]);
  
  ifstream f (argv[5]);
  ofstream g (argv[6]);
  ofstream h (argv[7]);
  
  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      f >> re->first >> re->second;
      res.push_back(re);
    }

  gsl_vector* x = gsl_vector_alloc(hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2);

  g << scientific;
  h << scientific;
  
  for (int i=-step; i<=step; ++i)
    {
      gsl_vector_set_zero(x);
      gsl_vector_set(x, 0, 10.0 * ((double)i) / step);
      double num = calc_dx_pi_llik_numeric(res, x, 1, hpa);
      double analytic = calc_dx_pi_llik_analytic(res, x, 1, hpa);
      
      if (fabs(num) > 0)
        g << gsl_vector_get(x, 0) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
      else
        g << gsl_vector_get(x, 0) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
    }

  for (int i=-step; i<=step; ++i)
    {
      gsl_vector_set_zero(x);
      gsl_vector_set(x, hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2 - 1, 10.0 * ((double)i) / step);
      double num = calc_dx_kappa_llik_numeric(res, x, hpa.TOTAL_CN, hpa.TOTAL_CN, hpa);
      double analytic = calc_dx_kappa_llik_analytic(res, x, hpa.TOTAL_CN, hpa.TOTAL_CN, hpa);
      
      if (fabs(num) > 0)
        h << gsl_vector_get(x, hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2 - 1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
      else
        h << gsl_vector_get(x, hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2 - 1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
    }

  for (int i=0; i<n; ++i)
    delete res[i];

  gsl_vector_free(x);
  
  f.close();
  g.close();
  h.close();

  return 0;
}
