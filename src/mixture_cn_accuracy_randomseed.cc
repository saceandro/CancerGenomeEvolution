#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
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

void calc_kappa(const gsl_vector* x, Vdouble& kappa, unsigned int MAX_SUBTYPE)
{
  double m = gsl_vector_get(x, 0);
  for (unsigned int i=1; i<MAX_SUBTYPE; ++i)
    {
      double s = gsl_vector_get(x, i);
      if (m < s)
        m = s;
    }
  
  double sum = 0;
  for (unsigned int k=1; k<=MAX_SUBTYPE; ++k)
    {
      kappa[k] = exp(gsl_vector_get(x, k-1) - m);
      sum += kappa[k];
    }

  for (unsigned int k=1; k<=MAX_SUBTYPE; ++k)
    kappa[k] /= sum;
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


double my_f (const gsl_vector *v, void *params)
{
  state *p = (state *) params;

  Vdouble kappa (p->TOTAL_CN+1, 0);
  calc_kappa(v, kappa, p->MAX_SUBTYPE);

  return -calc_llik(p->res, kappa, p->MAX_SUBTYPE, p->TOTAL_CN);
}

/* The gradient of f, df = (df/dx). */
void my_df (const gsl_vector *v, void *params, gsl_vector *df)
{
  state *p = (state *) params;

  Vdouble kappa (p->TOTAL_CN+1, 0);
  calc_kappa(v, kappa, p->MAX_SUBTYPE);

  Vdouble d_kappa (p->TOTAL_CN+1, 0);
  for (int r=1; r<=p->TOTAL_CN; ++r)
    d_kappa[r] = calc_d_kappa_llik(p->res, kappa, r, p->MAX_SUBTYPE, p->TOTAL_CN);

  for (int s=1; s<=p->TOTAL_CN; ++s)
    {
      double grad = 0;
      for (int r=1; r<=p->TOTAL_CN; ++r)
        grad += dx_kappa(kappa, r, s) * d_kappa[r];
      gsl_vector_set(df, s-1, -grad);
    }
}

/* Compute both f and df together. */
void my_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df)     
{
  *f = my_f(x, params);

  my_df(x, params, df);
}

void minimize(state& st, Vdouble& kappa)
{
  size_t iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  my_func.n = st.TOTAL_CN;
  my_func.f = my_f;
  my_func.df = my_df;
  my_func.fdf = my_fdf;
  my_func.params = &st;

  x = gsl_vector_alloc (st.TOTAL_CN);

  for (int r=1; r<=st.TOTAL_CN; ++r)
    gsl_vector_set (x, r-1, 0); // initiate
  
  //cerr << scientific;
  
  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc (T, st.TOTAL_CN);

  gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.001, 1e-4);

  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

      if (status)
        {
          if (status == GSL_ENOPROG)
            {
              //cerr << "No more improvement can be made for current estimate" << endl << endl;
              calc_kappa(s->x, kappa, st.MAX_SUBTYPE);
            }
          break;
        }

      calc_kappa(s->x, kappa, st.MAX_SUBTYPE);
      
      status = gsl_multimin_test_gradient (s->gradient, 1e-3);
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
  while (status == GSL_CONTINUE && iter < 100);

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);
}

int main(int argc, char** argv)
{
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 8)
    {
      cerr << "usage: ./mixture_cn_accuracy_randomseed max_subtype total_cn n (infile) (kappa_answer_file) (accuracy_file)" << endl;
      exit(EXIT_FAILURE);
    }

  state st;
  unsigned int n;
  
  st.MAX_SUBTYPE = atoi(argv[1]);
  st.TOTAL_CN = atoi(argv[2]);
  n = atoi(argv[3]);
  
  ifstream f (argv[4]);

  ifstream h (argv[5]);
  unsigned int total_cn;
  h >> total_cn;
  Vdouble kappa_ans (total_cn, 0);
  for (int r=1; r<=total_cn; ++r)
    h >> kappa_ans[r];
  
  Vdouble kappa (st.TOTAL_CN + 1, 0);
  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      f >> re->first >> re->second;
      st.res.push_back(re);
    }

  ofstream g (argv[6], ofstream::app);
  
  minimize(st, kappa);
  write_diff(g, kappa, kappa_ans, n);

  for (int i=0; i<n; ++i)
    delete st.res[i];
  
  f.close();
  h.close();
  g.close();
  
  return 0;
}
