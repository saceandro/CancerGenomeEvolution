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

void write_params(std::ofstream& f, params& pa)
{
  f << endl << "pi" << endl;
  // double sum = 0;
  for (int c=1; c<pa.pi.size(); ++c)
    {
      // sum += pa.pi[c];
      f << pa.pi[c] << "\t";
    }
  f << endl;

  // cout << "pi_sum: " << sum << endl;
  
  f << endl << "kappa" << endl;
  for (int c=1; c<pa.pi.size(); ++c)
    {
      // double sum2 = 0;
      for (int d=1; d<=c; ++d)
        {
          // sum2 += pa.kappa[c][d];
          f << pa.kappa[c][d] << "\t";
        }
      f << endl;
      // cout << "kappa_sum (l = " << c << "): " << sum2 << endl;
    }
  // cout << endl;
  
  return;
}

void write_pi_diff(ofstream& f, Vdouble& x, Vdouble& ans, unsigned int n)
{
  double sum = 0;
  for (int i=1; i<x.size(); ++i)
    {
      sum += pow(x[i] - ans[i], 2.0);
    }
  f << n << "\t" << sum << endl;
}

void write_kappa_diff(ofstream& f, VVdouble& kappa, VVdouble& kappa_ans, unsigned int n)
{
  double sum = 0;
  for (int l=1; l<kappa.size(); ++l)
    {
      for (int r=1; r<=l; ++r)
        {
          sum += pow(kappa[l][r] - kappa_ans[l][r], 2.0);
        }
    }
  sum /= kappa.size() - 1;
  
  f << n << "\t" << sum << endl;
}

void copy_params(params& from, params& target)
{
  for (int l=1; l<from.pi.size(); ++l)
    {
      target.pi[l] = from.pi[l];
      for (int r=1; r<=l; ++r)
        {
          target.kappa[l][r] = from.kappa[l][r];
        }
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

void clear_params(params& pa, hyperparams& hpa)
{
  for (int l=1; l<=hpa.TOTAL_CN; ++l)
    {
      pa.pi[l] = 0;
      for (int r=1; r<=hpa.TOTAL_CN; ++r)
        pa.kappa[l][r] = 0;
    }
}

bool check_tol(params& pa, params& pa_ans, double tol)
{
  bool ret = true;
  
  double pi_sum = 0;
  for (int l=1; l<pa.pi.size(); ++l)
    {
      pi_sum += pow(pa.pi[l], 2.0);
    }
  
  double sum = 0;
  for (int l=1; l<pa.pi.size(); ++l)
    {
      sum += pow(pa.pi[l] - pa_ans.pi[l], 2.0);
    }
  
  ret = ret && (sum / pi_sum < tol*tol);

  for (int l=1; l<pa.pi.size(); ++l)
    {
      double kappa_sum = 0;
      for (int r=1; r<pa.kappa.size(); ++r)
        {
          kappa_sum += pow(pa.kappa[l][r], 2.0);
        }
  
      double diff_sum = 0;
      for (int r=1; r<pa.kappa.size(); ++r)
        {
          diff_sum += pow(pa.kappa[l][r] - pa_ans.kappa[l][r], 2.0);
        }
      
      ret = ret && (diff_sum / kappa_sum < tol*tol);
    }
  return ret;
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

double em(READS& res, params& pa_old, params& pa_new, hyperparams& hpa)
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
      
      responsibility_numerator_all(*res[k], sts, st, pa_old, hpa, 0);
      llik += log(responsibility_partition(*res[k], sts, pa_old, hpa));

      for (unsigned int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (states::iterator it = sts.begin(); it != sts.end(); ++it)
            {
              int bin_l = calc_bin_pi(**it, hpa, l);
              pa_new.pi[l] += (*it)->resp * bin_l;
            }

          for (unsigned int r=1; r<=l; ++r)
            {
              for (states::iterator it = sts.begin(); it != sts.end(); ++it)
                {
                  int bin_lr = calc_bin_kappa(**it, hpa, l, r);
                  pa_new.kappa[l][r] += (*it)->resp * bin_lr;
                }
            }
        }
      
      delete_states(sts);
      sts.clear();
    }

  for (unsigned int l=1; l<=hpa.TOTAL_CN; ++l)
    {
      for (unsigned int r=1; r<=l; ++r)
        {
          pa_new.kappa[l][r] /= pa_new.pi[l];
        }
      pa_new.pi[l] /= K * (hpa.MAX_SUBTYPE + 1);
    }
  
  return llik;
}

int main(int argc, char** argv)
{
  // cerr << scientific;
  
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 11)
    {
      cerr << "usage: ./mixture_total_variant_calc_llik max_subtype n (infile) (kappa_answer_file) (kappa_file) (pi_accuracy_file) (kappa_accuracy_file) (ans_llik_file) (llik_file) total_cn" << endl;
      exit(EXIT_FAILURE);
    }

  hyperparams hpa;
  unsigned int n;
  
  hpa.MAX_SUBTYPE = atoi(argv[1]);
  hpa.TOTAL_CN = atoi(argv[10]);
  n = atoi(argv[2]);
  
  ifstream e (argv[3]);
  ifstream f (argv[4]);
  ifstream g (argv[5]);
  ofstream h (argv[6]);
  ofstream ff (argv[7]);
  ofstream gg (argv[8]);
  ofstream hh (argv[9]);

  READS res;
  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      e >> re->first >> re->second;
      res.push_back(re);
    }

  params pa_ans;
  init_params(pa_ans, hpa);

  params pa;
  init_params(pa, hpa);
  
  for (int l=1; l<=hpa.TOTAL_CN; ++l)
    {
      f >> pa_ans.pi[l];
      g >> pa.pi[l];
    }

  for (int l=1; l<=hpa.TOTAL_CN; ++l)
    {
      for (int r=1; r<=l; ++r)
        {
          f >> pa_ans.kappa[l][r];
          g >> pa.kappa[l][r];
        }
    }

  // write_params((ofstream&) cout, pa);
  
  params pa_new;
  init_params(pa_new, hpa);

  h << scientific;
  ff << scientific;
  gg << scientific;
  hh << scientific;

  write_pi_diff(h, pa.pi, pa_ans.pi, n);
  write_kappa_diff(ff, pa.kappa, pa_ans.kappa, n);

  double ans_llik = em(res, pa_ans, pa_new, hpa);
  gg << ans_llik << endl;
  
  clear_params(pa_new, hpa);
  
  double llik = em(res, pa, pa_new, hpa);
  hh << llik << endl;
  
  for (int i=0; i<n; ++i)
    delete res[i];

  e.close();
  f.close();
  g.close();
  h.close();
  ff.close();
  gg.close();
  hh.close();

  return 0;
}
