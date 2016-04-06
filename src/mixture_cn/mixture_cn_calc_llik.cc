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

typedef struct _hyperparams
{
  unsigned int MAX_SUBTYPE;
  unsigned int TOTAL_CN;
}
  hyperparams;

double sum_square_error(Vdouble& kappa, Vdouble& kappa_ans)
{
  double kappa_sum = 0;
  for (int r=1; r<kappa.size(); ++r)
    {
      kappa_sum += pow(kappa[r], 2.0);
    }
  
  double sum = 0;
  for (int r=1; r<kappa.size(); ++r)
    {
      sum += pow(kappa[r] - kappa_ans[r], 2.0);
    }
  return sum / kappa_sum;
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

double calc_mu(Vuint& a, hyperparams& hpa)
{
  unsigned int sum = 0;
  for (Vuint::iterator it = a.begin(); it != a.end(); ++it)
    sum += *it;
  return ((double) sum) / (hpa.MAX_SUBTYPE + 1) / hpa.TOTAL_CN;
}

double responsibility_partition_sub(READ& re, Vuint& cns, Vdouble& kappa, hyperparams& hpa, unsigned int subtype)
{
  double sum = 0;

  if (subtype < hpa.MAX_SUBTYPE)
    {
      for (int variant_cn=1; variant_cn<=hpa.TOTAL_CN; ++variant_cn)
        {
          cns[subtype] = variant_cn;
          sum += kappa[variant_cn] * responsibility_partition_sub(re, cns, kappa, hpa, subtype + 1); // caution: ++subtypeにすると、次のループで+1されて計算されてしまう
        }
    }
  
  else
    {
      for (int variant_cn=1; variant_cn<=hpa.TOTAL_CN; ++variant_cn)
        {
          cns[subtype] = variant_cn;
          double mu = calc_mu(cns, hpa);
          sum += kappa[variant_cn] * gsl_ran_binomial_pdf(re.first, mu, re.second);;
        }
    }
  return sum;
}

double calc_llik(READS& res, Vdouble& kappa, hyperparams& hpa)
{
  Vuint cns (hpa.MAX_SUBTYPE+1, 0);
  double llik = 0;

  for (int k=0; k<res.size(); ++k)
    llik += log(responsibility_partition_sub(*res[k], cns, kappa, hpa, 0));
  
  return llik;
}

int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 6)
    {
      cerr << "usage: ./mixture_cn_calc_llik max_subtype total_cn n (infile) (kappa_file)" << endl;
      exit(EXIT_FAILURE);
    }

  hyperparams hpa;
  unsigned int n;
  
  hpa.MAX_SUBTYPE = atoi(argv[1]);
  hpa.TOTAL_CN = atoi(argv[2]);
  n = atoi(argv[3]);
  
  ifstream f (argv[4]);
  ifstream g (argv[5]);

  READS res;
  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      f >> re->first >> re->second;
      res.push_back(re);
    }

  Vdouble kappa (hpa.TOTAL_CN + 1, 0);
  for (int r=1; r<=hpa.TOTAL_CN; ++r)
    g >> kappa[r];

  double llik = calc_llik(res, kappa, hpa);

  cout << scientific;
  cout << "llik: " << llik << endl;

  for (int i=0; i<n; ++i)
    delete res[i];
  
  f.close();
  g.close();

  return 0;
}
