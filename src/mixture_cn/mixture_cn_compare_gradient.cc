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
typedef vector<unsigned int> Vuint;
typedef pair<unsigned int, unsigned int> READ;
typedef vector<READ*> READS;

typedef struct _state
{
  unsigned int MAX_SUBTYPE;
  unsigned int TOTAL_CN;
  READS res;
}
  state;

double calc_sum(Vdouble& x)
{
  double sum = 0;
  for (Vdouble::iterator it = x.begin(); it != x.end(); ++it)
    sum += *it;
  return sum;
}

void write_diff(ofstream& f, Vdouble& x, Vdouble& ans, unsigned int n)
{
  double sum = 0;
  for (int i=1; i<x.size(); ++i)
    {
      sum += pow(x[i] - ans[i], 2.0);
    }
  f << n << "\t" << sum << endl;
}

void calc_kappa(Vdouble& x, Vdouble& kappa, unsigned int TOTAL_CN)
{
  double m = x[1];
  for (unsigned int r=2; r<=TOTAL_CN; ++r)
    {
      if (m < x[r])
        m = x[r];
    }
  
  double sum = 0;
  for (unsigned int r=1; r<=TOTAL_CN; ++r)
    {
      kappa[r] = exp(x[r] - m);
      sum += kappa[r];
    }

  for (unsigned int r=1; r<=TOTAL_CN; ++r)
    kappa[r] /= sum;
}

double dx_kappa(Vdouble& kappa, unsigned int r, unsigned int s)
{
  if (r == s)
    return kappa[r] * (1.0 - kappa[r]);

  return -kappa[r] * kappa[s];
}

double calc_mu(Vuint& a, unsigned int MAX_SUBTYPE, unsigned int TOTAL_CN)
{
  unsigned int sum = 0;
  for (Vuint::iterator it = a.begin(); it != a.end(); ++it)
    sum += *it;
  return ((double) sum) / (MAX_SUBTYPE + 1) / TOTAL_CN;
}

double lik_sub(READ& re, Vuint& cns, Vdouble& kappa, int curr_subtype, unsigned int MAX_SUBTYPE, unsigned int TOTAL_CN)
{
  double sum = 0;
  if (curr_subtype < MAX_SUBTYPE)
    {
      for (int variant_cn = 1; variant_cn <= TOTAL_CN; ++variant_cn)
        {
          cns[curr_subtype] = variant_cn;
          sum += kappa[variant_cn] * lik_sub(re, cns, kappa, curr_subtype + 1, MAX_SUBTYPE, TOTAL_CN);
        }
    }
  
  else
    {
      for (int variant_cn = 1; variant_cn <= TOTAL_CN; ++variant_cn)
        {
          cns[curr_subtype] = variant_cn;
          double mu = calc_mu(cns, MAX_SUBTYPE, TOTAL_CN);
          sum += kappa[variant_cn] * gsl_ran_binomial_pdf(re.first, mu, re.second);
        }
    }
  return sum;
}

double calc_lik(READ& re, Vdouble& kappa, unsigned int MAX_SUBTYPE, unsigned int TOTAL_CN)
{
  Vuint cns (MAX_SUBTYPE+1, 0);
  return lik_sub(re, cns, kappa, 0, MAX_SUBTYPE, TOTAL_CN);
}

double calc_llik(READS& res, Vdouble& kappa, unsigned int MAX_SUBTYPE, unsigned int TOTAL_CN)
{
  Vuint cns (MAX_SUBTYPE+1, 0);
  double lik = 0;

  for (int k=0; k<res.size(); ++k)
    lik += log(lik_sub(*res[k], cns, kappa, 0, MAX_SUBTYPE, TOTAL_CN));
  
  return lik;
}

typedef struct _dkappa
{
  READS &res;
  Vdouble &x;
  unsigned int r;
  unsigned int MAX_SUBTYPE;
  unsigned int TOTAL_CN;

  _dkappa (READS& _res, Vdouble& _x, unsigned int _r, unsigned int _MAX_SUBTYPE, unsigned int _TOTAL_CN) : res(_res), x(_x), r(_r), MAX_SUBTYPE(_MAX_SUBTYPE), TOTAL_CN(_TOTAL_CN) {}
}
  dkappa;

double calc_llik_for_diff(double x_r, void* params)
{
  dkappa* p = (dkappa*) params;
  p->x[p->r] = x_r;
  
  Vdouble kappa (p->TOTAL_CN + 1, 0);
  calc_kappa(p->x, kappa, p->TOTAL_CN);
  
  Vuint cns (p->MAX_SUBTYPE+1, 0);
  double lik = 0;

  for (int k=0; k < p->res.size(); ++k)
    lik += log(lik_sub(*(p->res[k]), cns, kappa, 0, p->MAX_SUBTYPE, p->TOTAL_CN));
  
  return lik;
}

double d_kappa_lik_sub(READ& re, Vuint& cns, Vdouble& kappa, unsigned int r, unsigned int curr_subtype, unsigned int MAX_SUBTYPE, unsigned int TOTAL_CN)
{
  double sum = 0;
  if (curr_subtype < MAX_SUBTYPE)
    {
      for (int variant_cn = 1; variant_cn <= TOTAL_CN; ++variant_cn)
        {
          cns[curr_subtype] = variant_cn;
          sum += d_kappa_lik_sub(re, cns, kappa, r, curr_subtype + 1, MAX_SUBTYPE, TOTAL_CN);
        }
    }
  
  else
    {
      for (int variant_cn = 1; variant_cn <= TOTAL_CN; ++variant_cn)
        {
          cns[curr_subtype] = variant_cn;
          double mu = calc_mu(cns, MAX_SUBTYPE, TOTAL_CN);

          Vuint bin;
          bin.assign(MAX_SUBTYPE+1, 0);
          
          for (int c=0; c<=TOTAL_CN; ++c)
            {
              for (int i=0; i<=MAX_SUBTYPE; ++i)
                {
                  if (cns[i] == c)
                    bin[c]++;
                }
            }

          double prod = 1;
          for (int i=1; i<=TOTAL_CN; ++i)
            {
              if (i == r)
                prod *= bin[i] * pow(kappa[i], bin[i] - 1);
              else
                prod *= pow(kappa[i], bin[i]);
            }

          sum += prod * gsl_ran_binomial_pdf(re.first, mu, re.second);
        }
    }
  return sum;
}

double calc_d_kappa_llik(READS& res, Vdouble& kappa, unsigned int r, unsigned int MAX_SUBTYPE, unsigned int TOTAL_CN)
{
  Vuint cns (MAX_SUBTYPE+1, 0);

  double d_lik = 0;

  for (int k=0; k<res.size(); ++k)
    {
      double lik = calc_lik(*res[k], kappa, MAX_SUBTYPE, TOTAL_CN);
      double d_kappa_lik = d_kappa_lik_sub(*res[k], cns, kappa, r, 0, MAX_SUBTYPE, TOTAL_CN);
      
      d_lik += d_kappa_lik / lik;
    }
  
  return d_lik;
}

double calc_dx_kappa_llik_numeric(READS& res, Vdouble& x, unsigned int r, unsigned int MAX_SUBTYPE, unsigned int TOTAL_CN)
{
  gsl_function F;
  dkappa dk (res, x, r, MAX_SUBTYPE, TOTAL_CN);
  dk.res = res;
  for (int i=1; i<=TOTAL_CN; ++i)
    dk.x[i] = x[i];
  dk.r = r;
  dk.MAX_SUBTYPE = MAX_SUBTYPE;
  dk.TOTAL_CN = TOTAL_CN;

  F.function = &calc_llik_for_diff;
  F.params = &dk;

  double result, abserr;
  gsl_deriv_central(&F, x[r], 1e-7, &result, &abserr);
  
  return result;
}

double calc_dx_kappa_llik_analytic(READS& res, Vdouble& x, unsigned int r, unsigned int MAX_SUBTYPE, unsigned int TOTAL_CN)
{
  Vdouble kappa (TOTAL_CN+1, 0);
  calc_kappa(x, kappa, TOTAL_CN);
  
  Vdouble d_kappa (TOTAL_CN+1, 0);
  
  for (int s=1; s<=TOTAL_CN; ++s)
    d_kappa[s] = calc_d_kappa_llik(res, kappa, s, MAX_SUBTYPE, TOTAL_CN);

  double grad = 0;
  for (int s=1; s<=TOTAL_CN; ++s)
    grad += dx_kappa(kappa, s, r) * d_kappa[s];

  return grad;
}

double calc_rel_err(double x, double y)
{
  return fabs(x - y) / fabs(x);
}

int main(int argc, char** argv)
{
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 6)
    {
      cerr << "usage: ./mixture_cn_compare_gradient max_subtype total_cn n step (infile)" << endl;
      exit(EXIT_FAILURE);
    }

  state st;
  int n;
  int step;

  st.MAX_SUBTYPE = atoi(argv[1]);
  st.TOTAL_CN = atoi(argv[2]);
  n = atoi(argv[3]);
  step = atoi(argv[4]);
  
  ifstream f (argv[5]);

  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      f >> re->first >> re->second;
      st.res.push_back(re);
    }

  Vdouble x (st.TOTAL_CN, 0);

  for (int i=-step; i<=step; ++i)
    {
      x[1] = 10.0 * ((double)i) / step;
      double num = calc_dx_kappa_llik_numeric(st.res, x, 1, st.MAX_SUBTYPE, st.TOTAL_CN);
      double analytic = calc_dx_kappa_llik_analytic(st.res, x, 1, st.MAX_SUBTYPE, st.TOTAL_CN);
      
      if (fabs(num) > 0)
        cout << x[1] << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
      else
        cout << x[1] << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
    }

  for (int i=0; i<n; ++i)
    delete st.res[i];
  
  f.close();
  
  return 0;
}
