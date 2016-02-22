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

unsigned int calc_bin_r(Vuint& cns, hyperparams& hpa, unsigned int r)
{
  unsigned int bin = 0;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      if (cns[i] == r)
        bin++;
    }
  return bin;
}

double responsibility_weighted_numerator_sub(READ& re, Vuint& cns, Vdouble& kappa, hyperparams& hpa, unsigned int subtype, unsigned int r)
{
  double sum = 0;

  if (subtype < hpa.MAX_SUBTYPE)
    {
      for (int variant_cn=1; variant_cn<=hpa.TOTAL_CN; ++variant_cn)
        {
          cns[subtype] = variant_cn;
          sum += kappa[variant_cn] * responsibility_weighted_numerator_sub(re, cns, kappa, hpa, subtype + 1, r);
        }
    }
  
  else
    {
      for (int variant_cn=1; variant_cn<=hpa.TOTAL_CN; ++variant_cn)
        {
          cns[subtype] = variant_cn;
          double mu = calc_mu(cns, hpa);
          unsigned int bin = calc_bin_r(cns, hpa, r);

          sum += kappa[variant_cn] * gsl_ran_binomial_pdf(re.first, mu, re.second) * bin;
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

double em(READS& res, Vdouble& kappa_old, Vdouble& kappa_new, hyperparams& hpa)
{
  double llik = 0;
  unsigned int K = res.size();
  Vuint cns (hpa.MAX_SUBTYPE + 1, 0);
  Vdouble denominator (res.size(), 0);

  // fill(kappa_new.begin(), kappa_new.end(), 0);

  for (int k=0; k<res.size(); ++k)
    {
      denominator[k] = responsibility_partition_sub(*res[k], cns, kappa_old, hpa, 0);
      llik += log(denominator[k]);
    }
  
  
  for (unsigned int r=1; r<=hpa.TOTAL_CN; ++r)
    {
      kappa_new[r] = 0;
      for (int k=0; k<res.size(); ++k)
        {
          double numerator = responsibility_weighted_numerator_sub(*res[k], cns, kappa_old, hpa, 0, r) / K / (hpa.MAX_SUBTYPE + 1);
          kappa_new[r] += numerator / denominator[k];
        }
    }
  return llik;
}

int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 9)
    {
      cerr << "usage: ./mixture_cn_em max_subtype n (infile) (kappa_answer_file) (accuracy_file) iter tol (log_file)" << endl;
      exit(EXIT_FAILURE);
    }

  hyperparams hpa;
  unsigned int n;
  
  hpa.MAX_SUBTYPE = atoi(argv[1]);
  n = atoi(argv[2]);
  
  ifstream f (argv[3]);
  ifstream h (argv[4]);
  ofstream g (argv[5]);
  ofstream l (argv[8]);
  
  int iter;
  iter = atoi(argv[6]);

  double tol;
  tol = atof(argv[7]);
  
  h >> hpa.TOTAL_CN;
  Vdouble kappa_ans (hpa.TOTAL_CN + 1, 0);
  for (int r=1; r<=hpa.TOTAL_CN; ++r)
    h >> kappa_ans[r];

  Vdouble kappa (hpa.TOTAL_CN + 1, 1.0 / hpa.TOTAL_CN);

  READS res;
  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      f >> re->first >> re->second;
      res.push_back(re);
    }

  l << scientific;
  Vdouble kappa_new (hpa.TOTAL_CN + 1, 0);
  for (int i=0; i<iter; ++i)
    {
      fill(kappa_new.begin(), kappa_new.end(), 0);
      double llik = em(res, kappa, kappa_new, hpa);
      for (int r=1; r<=hpa.TOTAL_CN; ++r)
        {
          l << kappa[r] << " ";
        }
      l << endl;

      l << "llik: " << llik << endl << endl;

      if (sum_square_error(kappa_new, kappa) < tol*tol)
        break;
      
      copy(kappa_new.begin(), kappa_new.end(), kappa.begin());

    }

  for (int r=1; r<=hpa.TOTAL_CN; ++r)
    {
      l << kappa_new[r] << " ";
    }
  l << endl;
  
  l << "llik: " << calc_llik(res, kappa_new, hpa) << endl;

  g << scientific;
  write_diff(g, kappa, kappa_ans, n);

  for (int i=0; i<n; ++i)
    delete res[i];
  
  f.close();
  h.close();
  g.close();
  l.close();

  return 0;
}
