#include <iostream>
#include <fstream>
#include <math.h>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
using namespace std;

#define calc_sigmoid(x) (((x) < 0) ? (exp(x) / (exp(x) + 1.0)) : (1.0 / (exp(-(x)) + 1.0)))

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
  unsigned int *p = (unsigned int *) params;

  double x = gsl_vector_get(v, 0);

  double mu = calc_sigmoid(x);

  return -log(gsl_ran_binomial_pdf(p[0], mu, p[1]));
}

/* The gradient of f, df = (df/dx). */
void my_df (const gsl_vector *v, void *params, gsl_vector *df)
{
  unsigned int *p = (unsigned int *)params;

  double x = gsl_vector_get(v, 0);

  double mu = calc_sigmoid(x);

  double grad = -calc_dx_sigmoid(x) * (p[0]/mu - (p[1]-p[0])/(1-mu));
  
  gsl_vector_set(df, 0, grad);
}


/* Compute both f and df together. */
void my_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df)     
{
  *f = my_f(x, params);

  my_df(x, params, df);
}

void minimize(int m, int M, double& mu)
{
  size_t iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  unsigned int par[2];
  par[0] = m;
  par[1] = M;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  my_func.n = 1;
  my_func.f = my_f;
  my_func.df = my_df;
  my_func.fdf = my_fdf;
  my_func.params = par;

  x = gsl_vector_alloc (1);
  gsl_vector_set (x, 0, 0); // initiate with x
  
  cerr << scientific;
  
  double x0 = 0;
  double mu0 = calc_sigmoid(x0);
  double grad0 = calc_dx_sigmoid(x0) * (par[0]/mu0 - (par[1]-par[0])/(1-mu0));
  double llk = log(gsl_ran_binomial_pdf(par[0], mu0, par[1]));
  
  cerr << iter << " " << x0 << " " << mu0 << " " << grad0 <<  " " << llk << endl << endl;

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
      double grad_star = calc_dx_sigmoid(x_star) * (par[0]/mu_star - (par[1]-par[0])/(1-mu_star));
      double llk = log(gsl_ran_binomial_pdf(par[0], mu_star, par[1]));
      
      status = gsl_multimin_test_gradient (s->gradient, 1e-3);
      if (status == GSL_SUCCESS)
        {
          cerr << "Minimum found at:" << endl;
          mu = mu_star;
        }
      
      cerr << iter << " " << x_star << " " << mu_star << " " << grad_star << " " << -s->f << endl << endl;
    }
  while (status == GSL_CONTINUE && iter < 100);

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);
}

int main(int argc, char** argv)
{
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 4)
    {
      cerr << "usage: ./binomial_inference (infile) n (outfile)" << endl;
      exit(EXIT_FAILURE);
    }
  
  ifstream f (argv[1]);
  unsigned int n;
  n = atoi(argv[2]);
  ofstream g (argv[3]);
  g << scientific;
  
  unsigned int m, M;
  for (int i=0; i<n; ++i)
    {
      double mu = -1;
      f >> m >> M;
      cerr << "m: " << m << endl;
      cerr << "M: " << M << endl;
      
      minimize(m, M, mu);
      g << ((double) m) / M << "\t" << mu << endl;
    }

  f.close();
  g.close();
  
  return 0;
}
