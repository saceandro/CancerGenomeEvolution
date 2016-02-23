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

typedef std::vector<state> states;

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

double calc_mu(state& st, hyperparams& hpa)
{
  unsigned int denom = 0;
  unsigned int num = 0;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      denom += st.l[i];
      num += st.r[i];
    }
    
  return ((double) num) / denom;
}

double responsibility_partition_sub(READ& re, state& st, params& pa, hyperparams& hpa, unsigned int subtype)
{
  double sum1 = 0;
  
  if (subtype < hpa.MAX_SUBTYPE)
    {
      for (st.total_cn[subtype] = 1; st.total_cn[subtype] <= hpa.TOTAL_CN; ++st.total_cn[subtype])
        {
          double sum2 = 0;
          
          for (st.variant_cn[subtype] = 1; st.variant_cn[subtype] <= total_cn; ++st.variant_cn[subtype])
            {
              sum2 += pa.kappa[st.total_cn[subtype]][st.variant_cn[subtype]] * responsibility_partition_sub(re, st, pa, hpa, subtype + 1); // caution: ++subtypeにすると、次のループで+1されて計算されてしまう
            }
          sum1 += pa.pi[st.total_cn[subtype]] * sum2;
        }
    }
  
  else
    {
      for (st.total_cn[subtype] = 1; st.total_cn[subtype] <= hpa.TOTAL_CN; ++st.total_cn[subtype])
        {
          double sum2 = 0;
          
          for (st.variant_cn[subtype] = 1; st.variant_cn[subtype] <= total_cn; ++st.variant_cn[subtype])
            {
              double mu = calc_mu(st, hpa);
              sum2 += pa.kappa[st.total_cn[subtype]][st.variant_cn[subtype]] * gsl_ran_binomial_pdf(re.first, mu, re.second);
            }
          sum1 += pa.pi[st.total_cn[subtype]] * sum2;
        }
    }
  
  return sum1;
}

double responsibility_numerator(READ& re, state& st, params& pa, hyperparams& hpa)
{
  double product = 1;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    product *= pa.pi[st.total_cn[i]] * pa.kappa[st.total_cn[i]][st.variant_cn[i]];

  double mu = calc_mu(st, hpa);
  product *= gsl_ran_binomial_pdf(re.first, mu, re.second);
  
  return product;
}

double calc_responsibility(READ& re, state& st, params& pa, hyperparams& hpa)
{
  
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
      if (st.total_cn[i] == l)
        bin++;
    }
  return bin;
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
