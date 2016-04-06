#include <iostream>
#include <fstream>
#include <math.h>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
using namespace std;

#define calc_sigmoid(x) ((tanh((x)/2.0) + 1.0) / 2.0)

typedef pair<unsigned int, unsigned int> READ;

double calc_sigmoid_func(double x, void* params)
{
  return (tanh(x/2.0) + 1.0) / 2.0;
}

double calc_dx_sigmoid_numeric(double x)
{
  gsl_function F;

  F.function = &calc_sigmoid_func;
  F.params = NULL;

  double result, abserr;
  gsl_deriv_central(&F, x, 1e-5, &result, &abserr);

  return result;
}

double calc_dx_sigmoid(double x)
{
  double y = tanh(x/2.0);
  return (1.0 - y*y) / 4.0;
}

double calc_llik(double x, void *params)
{
  READ *p = (READ*) params;

  double mu = calc_sigmoid(x);
  
  return log(gsl_ran_binomial_pdf(p->first, mu, p->second));
}

double dx_llik_analytic(double x, READ& re)
{
  double mu = calc_sigmoid(x);
  double dllik_dmu = re.first / mu - (re.second - re.first) / (1 - mu);
  
  return calc_dx_sigmoid(x) * dllik_dmu;
}

double dx_llik_numeric(double x, READ& re)
{
  gsl_function F;
  
  F.function = &calc_llik;
  F.params = &re;

  double result, abserr;
  gsl_deriv_central(&F, x, 1e-8, &result, &abserr);

  return result;
}

double calc_mu_llik(double mu, void *params)
{
  READ *p = (READ*) params;

  return log(gsl_ran_binomial_pdf(p->first, mu, p->second));
}

double dmu_llik_numeric(double mu, READ& re)
{
  gsl_function F;
  
  F.function = &calc_mu_llik;
  F.params = &re;

  double result, abserr;
  gsl_deriv_central(&F, mu, 1e-8, &result, &abserr);

  return result;
}

double dmu_llik_analytic(double mu, READ& re)
{
  return re.first/mu - (re.second - re.first)/(1-mu);
}

double calc_rel_err(double x, double y)
{
  return fabs(x - y) / fabs(x);
}

int main(int argc, char** argv)
{
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 2)
    {
      cerr << "usage: ./binomial_compare_gradient (split)" << endl;
      exit(EXIT_FAILURE);
    }

  int n = atoi(argv[1]);
  
  READ re;
  cin >> re.first >> re.second;
  
  for (int i=-n; i<=n; ++i)
    {
      double x = 15.0 * ((double)i) / n;
      double num = dx_llik_numeric(x, re);
      double analytic = dx_llik_analytic(x, re);

      if (fabs(num) > 0)
        cout << x << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
      else
        cout << x << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
    }

  // for (int i=-n; i<=n; ++i)
  //   {
  //     double x = 100.0 * ((double)i) / n;
  //     double num = calc_dx_sigmoid_numeric(x);
  //     double analytic = calc_dx_sigmoid(x);

  //     if (fabs(num) > 0)
  //       cout << x << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
  //     else
  //       cout << x << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
  //   }

  // for (int i=1; i<n; ++i)
  //   {
  //     double mu = ((double)i) / n;
  //     double num = dmu_llik_numeric(mu, re);
  //     double analytic = dmu_llik_analytic(mu, re);

  //     if (fabs(num) > 0)
  //       cout << mu << "\t" << num << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
  //     else
  //       cout << mu << "\t" << num << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
  //   }

  return 0;
}
