#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
using namespace std;

#define calc_sigmoid(x) (((x) < 0) ? (exp(x) / (exp(x) + 1.0)) : (1.0 / (exp(-(x)) + 1.0)))

typedef vector<double> Vdouble;
typedef vector<unsigned int> Vuint;
typedef pair<unsigned int, unsigned int> READ;
typedef vector<READ*> READS;

double calc_dx_sigmoid(double x)
{
  if (x < -10.0)
    return exp(x) / (2.0*exp(x) + 1.0);
  else if (x < 10.0)
    return 1.0 / (exp(-x) + 2.0 + exp(x));
  else
    return exp(-x) / (2.0*exp(-x) + 1.0);
}

double my_f (const gsl_vector *v, void *params)
{
  READS *p = (READS *) params;

  double x = gsl_vector_get(v, 0);

  double mu = calc_sigmoid(x);

  double llik = 0;

  for (READS::iterator it = p->begin(); it != p->end(); ++it)
    {
      llik += log(gsl_ran_binomial_pdf((*it)->first, mu, (*it)->second));
    }
  
  return -llik;
}

/* The gradient of f, df = (df/dx). */
void my_df (const gsl_vector *v, void *params, gsl_vector *df)
{
  READS *p = (READS *)params;

  double x = gsl_vector_get(v, 0);

  double mu = calc_sigmoid(x);

  double grad = 0;

  for (READS::iterator it = p->begin(); it != p->end(); ++it)
    {
      grad += (*it)->first/mu - ((*it)->second - (*it)->first)/(1-mu);
    }
  
  gsl_vector_set(df, 0, -calc_dx_sigmoid(x) * grad);
}


/* Compute both f and df together. */
void my_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df)     
{
  *f = my_f(x, params);

  my_df(x, params, df);
}

void minimize(READS& res, double& mu)
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
  my_func.params = &res;

  x = gsl_vector_alloc (1);
  gsl_vector_set (x, 0, 0); // initiate with x
  
  cerr << scientific;
  
  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc (T, 1);

  gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.001, 1e-4);

  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

      if (status)
        {
          if (status == GSL_ENOPROG)
            {
              cerr << "No more improvement can be made for current estimate" << endl << endl;
              mu = calc_sigmoid(gsl_vector_get(s->x, 0));
            }
          break;
        }

      double x_star = gsl_vector_get(s->x, 0);
      double mu_star = calc_sigmoid(x_star);
      
      status = gsl_multimin_test_gradient (s->gradient, 1e-3);
      if (status == GSL_SUCCESS)
        {
          cerr << "Minimum found at:" << endl;
          mu = mu_star;
        }
      
      cerr << iter << " " << x_star << " " << mu_star << " " << -gsl_vector_get(s->gradient, 0) << " " << -s->f << endl << endl;
    }
  while (status == GSL_CONTINUE && iter < 100);

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);
}

int main(int argc, char** argv)
{
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 3)
    {
      cerr << "usage: ./binomial_inference_multisites n (infile)" << endl;
      exit(EXIT_FAILURE);
    }
  
  ifstream f (argv[2]);
  unsigned int n;
  double mu = -1;
  
  n = atoi(argv[1]);
  
  READS res;
  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      f >> re->first >> re->second;
      res.push_back(re);
    }
  f.close();
  
  minimize(res, mu);
  cout << mu << endl;
  
  for (int i=0; i<n; ++i)
    delete res[i];
  
  return 0;
}
