#include <fstream>
#include <algorithm>
#include <numeric>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf.h>
#include "../../util/loglib.hh"
using namespace std;

#define calc_sigmoid(x) ((tanh((x)/2.0) + 1.0) / 2.0)

typedef vector<int> Vint;
typedef vector<double> Vdouble;
typedef vector<Vdouble> VVdouble;
typedef pair<int, int> READ;
typedef vector<READ*> READS;

class state 
{
public:
  int k;
  Vint total_cn;
  Vint variant_cn;
  VLog resp_ds;
  Log resp_num;
  Log resp;

  // use default constructor
  // state (int, std::vector<int>, std::vector<int>, double);
};

typedef std::vector<state*> states;

class param
{
public:
  Log s;
  Log n;
  VLog pi;
  VVLog kappa;

  param (Log _s, Log _n, VLog _pi, VVLog _kappa) : s(_s), n(_n), pi(_pi), kappa(_kappa) {}
};

typedef vector<param*> params;

typedef struct _hyperparams
{
  pair<double,double> be_hpa;
  Vdouble alpha;
  VVdouble beta;
  int MAX_SUBTYPE;
  int TOTAL_CN;
}
  hyperparams;

class diff
{
public:
  READS res;
  hyperparams hpa;
};

double calc_dx_sigmoid(double x)
{
  double y = tanh(x/2.0);
  return (1.0 - y*y) / 4.0;
}

double sum_vector(Vdouble& v, int s, int e)
{
  double sum = 0;
  for (int i=s; i<=e; ++i)
    sum += v[i];
  return sum;
}

void write_params(std::ofstream& f, params& pa, hyperparams& hpa)
{
  for (int i=0; i<hpa.MAX_SUBTYPE; ++i)
    f << pa[i]->s.eval() << "\t";
  f << endl << endl;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    f << pa[i]->n.eval() << "\t";
  f << endl << endl;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double sum = 0;
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        f << pa[i]->pi[l].eval() << "\t";
      f << endl;

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (int r=1; r<=l; ++r)
            f << pa[i]->kappa[l][r].eval() << "\t";
          f << endl;
        }
      f << endl;
    }
}

void copy_params(params& pa, params& target, hyperparams& hpa)
{
  for (int i=0; i<hpa.MAX_SUBTYPE; ++i)
    target[i]->s = pa[i]->s;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    target[i]->n = pa[i]->n;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        target[i]->pi[l] = pa[i]->pi[l];

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (int r=1; r<=l; ++r)
            target[i]->kappa[l][r] = pa[i]->kappa[l][r];
        }
    }
}

void init_state(state& st, hyperparams& hpa)
{
  st.total_cn.assign(hpa.MAX_SUBTYPE + 1, 0);
  st.variant_cn.assign(hpa.MAX_SUBTYPE + 1, 0);
  st.resp_ds.assign(hpa.MAX_SUBTYPE, Log(0));
  st.total_cn[0] = 2;
  st.variant_cn[0] = 0;
}

void init_params(params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    pa.push_back(new param (Log(0), Log(0), VLog (hpa.TOTAL_CN + 1, Log(0)), VVLog (hpa.TOTAL_CN + 1, VLog (hpa.TOTAL_CN + 1, Log(0)))));
}

void init_hyperparams(hyperparams& hpa)
{
  hpa.be_hpa.first = 0.1;
  hpa.be_hpa.second = 0.1;
  hpa.alpha.assign(hpa.TOTAL_CN + 1, 0.1);
  hpa.beta.assign(hpa.TOTAL_CN + 1, Vdouble (hpa.TOTAL_CN + 1, 0.1));
}

void delete_params(params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    delete pa[i];
}

void calc_n(params& pa, hyperparams& hpa)
{
  pa[0]->n = Log(1.0);
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    pa[i]->n = pa[i-1]->n * (Log(1.0) - pa[i-1]->s);

  for (int i=0; i<hpa.MAX_SUBTYPE; ++i)
    pa[i]->n *= pa[i]->s;
}

double calc_mu(state& st, params& pa, hyperparams& hpa)
{
  Log denom;
  Log num;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      denom += pa[i]->n * Log(st.total_cn[i]);
      num += pa[i]->n * Log(st.variant_cn[i]);
    }
    
  return (num / denom).eval();
}

#define log_binomial_pdf(m, mu, M) ((m) * log((mu)) + ((M) - (m)) * log1p(-(mu)) + gsl_sf_lnchoose((M), (m)))

Log d_bin_mu(READ& re, double mu)
{
  Log a (log_binomial_pdf(re.first-1, mu, re.second-1), 1);
  if (re.first < re.second)
    {
      Log b (log_binomial_pdf(re.first, mu, re.second-1), 1);
      return Log(re.second) * (a - b);
    }
  else
    {
      return Log(re.second) * a;
    }
}

Log d_n_s(params& pa, hyperparams& hpa, int i, int j)
{
  if (i < j)
    return Log(0);
  else if (i == j)
    return pa[i]->n / pa[j]->s;
  else
    return -pa[i]->n / (Log(1.0) - pa[j]->s);
}

Log d_mu_n(state& st, params& pa, hyperparams& hpa, int j)
{
  Log normal;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    normal += pa[i]->n * Log(st.total_cn[i]);

  Log variant;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    variant += pa[i]->n * Log(st.variant_cn[i]);

  Log a = Log(st.variant_cn[j]) / normal;
  Log b = Log(st.total_cn[j]) / normal * variant / normal;

  return a - b;
}

void responsibility_numerator(READ& re, states& sts, state& st, params& pa, hyperparams& hpa)
{
  Log product (1);
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    product *= pa[i]->pi[st.total_cn[i]] * pa[i]->kappa[st.total_cn[i]][st.variant_cn[i]];

  state* new_st = new state;
  init_state(*new_st, hpa);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      new_st->total_cn[i] = st.total_cn[i];
      new_st->variant_cn[i] = st.variant_cn[i];
    }

  double mu = calc_mu(st, pa, hpa);
  new_st->resp_num = product * Log(log_binomial_pdf(re.first, mu, re.second), 1);

  product *= d_bin_mu(re, mu);
  for (int j=0; j<hpa.MAX_SUBTYPE; ++j)
    {
      Log sum;
      for (int i=j; i<=hpa.MAX_SUBTYPE; ++i)
        {
          sum += d_mu_n(st, pa, hpa, i) * d_n_s(pa, hpa, i, j);
        }
      
      new_st->resp_ds[j] = product * sum;
    }
  
  sts.push_back(new_st);
}

void responsibility_numerator_all(READ& re, states& sts, state& st, params& pa, hyperparams& hpa, int subtype)
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

Log responsibility_partition(READ& re, states& sts, params& pa, hyperparams& hpa)
{
  Log partition;
  
  for (states::iterator it = sts.begin(); it != sts.end(); ++it)
    {
      // cout << (*it)->resp_num.take_log() << "\t" << (*it)->resp_num.get_sign() << endl;
      partition += (*it)->resp_num;
    }

  for (states::iterator it = sts.begin(); it != sts.end(); ++it)
    {
      (*it)->resp = (*it)->resp_num / partition;
      for (int i=0; i<hpa.MAX_SUBTYPE; ++i)
        {
          (*it)->resp_ds[i] /= partition;
        }
    }

  return partition;
}

void delete_states(states& sts)
{
  for (int i=0; i<sts.size(); ++i)
    {
      delete sts[i];
    }
}

#define log_beta_pdf(s, a, b) (-gsl_sf_lnbeta((a), (b)) + ((a) - 1.0) * log((s)) + ((b) - 1.0) * log1p(-(s)) )

double calc_llik(READS& res, params& pa, hyperparams& hpa)
{
  int K;
  K = res.size();
  
  Log lik (1);
  states sts;
  state st;
  init_state(st, hpa);
  
  for (int k=0; k<res.size(); ++k)
    {
      // cerr << "k: " << k << endl;
      
      responsibility_numerator_all(*res[k], sts, st, pa, hpa, 1);
      Log lik_k = responsibility_partition(*res[k], sts, pa, hpa);
      // cout << lik_k.take_log() << endl;
      lik *= lik_k;
      
      delete_states(sts);
      sts.clear();
    }

  for (int i=0; i<hpa.MAX_SUBTYPE; ++i)
    lik *= Log(log_beta_pdf(pa[i]->s.eval(), hpa.be_hpa.first, hpa.be_hpa.second), 1);
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          lik *= pa[i]->pi[l].take_pow(hpa.alpha[l] - 1.0);
          
          for (int r=1; r<=l; ++r)
            {
              lik *= pa[i]->kappa[l][r].take_pow(hpa.beta[l][r] - 1.0);
            }
        }
    }

  return lik.take_log();
}

int main(int argc, char** argv)
{
  cout << scientific;
  cerr << scientific;
  
  feenableexcept(FE_INVALID);
  
  if (argc != 7)
    {
      cerr << "usage: ./lda_stick_map_llik max_subtype total_cn n (reads) (s pi kappa infile) (llik outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  diff di;
  int n;

  di.hpa.MAX_SUBTYPE = atoi(argv[1]);
  di.hpa.TOTAL_CN = atoi(argv[2]);
  init_hyperparams(di.hpa);
  n = atoi(argv[3]);

  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  // T = gsl_rng_default; // default is mt19937
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, 1);

  ifstream f (argv[4]);
  ifstream g (argv[5]);
  ofstream h (argv[6]);
  
  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      f >> re->first >> re->second;
      di.res.push_back(re);
    }

  params pa;
  init_params(pa, di.hpa);

  for (int i=0; i<di.hpa.MAX_SUBTYPE; ++i)
    g >> pa[i]->s;

  for (int i=0; i<=di.hpa.MAX_SUBTYPE; ++i)
    g >> pa[i]->n;

  for (int i=1; i<=di.hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=di.hpa.TOTAL_CN; ++l)
        {
          g >> pa[i]->pi[l];
        }

      for (int l=1; l<=di.hpa.TOTAL_CN; ++l)
        {
          for (int r=1; r<=l; ++r)
            {
              g >> pa[i]->kappa[l][r];
            }
        }
    }

  h << scientific;
  h << calc_llik(di.res, pa, di.hpa) << endl;
  
  for (int i=0; i<n; ++i)
    delete di.res[i];

  delete_params(pa, di.hpa);
  
  f.close();
  g.close();
  h.close();

  return 0;
}
