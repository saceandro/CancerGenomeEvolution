#include <fstream>
#include <algorithm>
#include <numeric>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf.h>
#include <iomanip>
#include "../../util/enumtree_wf.hh"
#include "setting.hh"
#include <xmmintrin.h>
using namespace std;

#define calc_sigmoid(x) ((tanh((x)/2.0) + 1.0) / 2.0)

typedef void (*myfunc) (int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& partition);

extern void variant_fraction_partition(int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, VLog& vf, VLog& vf_numerator, VLog& vf_denominator, Log& partition);
extern void set_gegen(VVLog &gegen);
extern void set_gegen_integral(VLog &gegen_int, VLog &gegen_int_err);

typedef vector<int> Vint;
typedef vector<double> Vdouble;
typedef vector<Vdouble> VVdouble;
typedef vector<bool> Vbool;
typedef vector<Vbool> VVbool;
typedef pair<int, int> READ;
typedef vector<READ*> READS;
typedef vector<VVLog> VVVLog;
typedef vector<int*> QS;

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

void copy_beta(params& pa, params& target, hyperparams& hpa, subtypes&tr)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int j=0; j<(int)tr[i].children.size(); ++j)
        {
          target.pa[i]->beta[j] = pa.pa[i]->beta[j];
        }
    }
}

void read_beta(std::ifstream& f, params& pa, hyperparams& hpa, subtypes& tr)
{
  double a;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int j=0; j<(int)tr[i].children.size(); ++j)
        {
          f >> a;
          pa.pa[i]->beta[j] = Log(a);
        }
    }
}

void read_params(std::ifstream& f, params& pa, hyperparams& hpa, subtypes& tr)
{
  double a;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f >> a;
      pa.pa[i]->u = Log(a);
    }
  
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int j=0; j<(int)tr[i].children.size(); ++j)
        {
          f >> a;
          pa.pa[i]->beta[j] = Log(a);
        }
    }
}

void write_params(std::ofstream& f, params& pa, hyperparams& hpa, subtypes& tr)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f << pa.pa[i]->u.eval() << "\t";
    }
  f << endl << endl;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int j=0; j<(int)tr[i].children.size(); ++j)
        {
          f << pa.pa[i]->beta[j].eval() << "\t";
        }
      f << endl;
    }
}

void copy_params(params& pa, params& target, hyperparams& hpa, subtypes&tr)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      target.pa[i]->u = pa.pa[i]->u;
    }
  
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int j=0; j<(int)tr[i].children.size(); ++j)
        {
          target.pa[i]->beta[j] = pa.pa[i]->beta[j];
        }
    }
}

double calc_mu(subtypes& tr, hyperparams& hpa, int q)
{
  return (tr[q].x / Log(2) / (Log(1) + tr[0].n/tr[q].n)).eval();
}

#define log_binomial_pdf(m, mu, M) ((m) * log((mu)) + ((M) - (m)) * log1p(-(mu)) + gsl_sf_lnchoose((M), (m)))

void calc_params(const gsl_vector* x, params& pa, hyperparams& hpa, subtypes& tr)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->u = Log(calc_sigmoid(gsl_vector_get(x, i-1)));
    }
}

double calc_llik(READS& res, QS& qs, params& pa, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int)
{
  int K;
  K = res.size();

  VVLog vf (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0)));
  VVLog vf_numerator (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0)));
  VVLog vf_denominator (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0)));
  VLog partition (hpa.MAX_SUBTYPE + 1, Log(0));

  for (int q=1; q<=hpa.MAX_SUBTYPE; ++q)
    variant_fraction_partition(0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q], vf_numerator[q], vf_denominator[q], partition[q]);

  Log lik = Log(1);
  
  for (int k=0; k<K; ++k)
    {
      Log lik_k = Log(0);
      for (int s=1; s<=FRACTIONS; ++s)
        {
          tr[*qs[k]].x = Log(((double) s) / FRACTIONS);
          double mu = calc_mu(tr, hpa, *qs[k]);
          lik_k += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * vf[*qs[k]][s];
        }
      tr[*qs[k]].x = Log(0);
      
      lik *= lik_k;
    }

  return lik.take_log();
}

int main(int argc, char** argv)
{
  cout << scientific;
  cerr << scientific;
  
  // feenableexcept(FE_INVALID);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  gsl_set_error_handler_off ();
  
  if (argc != 7)
    {
      cerr << "usage: ./llik #subtype #topology #locus (alpha beta) (reads) (llik outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  READS res;
  QS qs;
  int n, MAX_SUBTYPE, TOTAL_CN, MAX_TREE;

  MAX_SUBTYPE = atoi(argv[1]);
  TOTAL_CN = 2;
  
  trees trs;
  trees_cons(trs, MAX_SUBTYPE);
  MAX_TREE = trs.size();

  hyperparams hpa (MAX_SUBTYPE, TOTAL_CN, MAX_TREE);

  int topology = atoi(argv[2]);
  
  n = atoi(argv[3]);

  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  // T = gsl_rng_default; // default is mt19937
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, 1);

  VVLog gegen;
  gegen.assign(FRACTIONS+1, VLog(GEGEN_MAX+1, Log(0)));

  set_gegen(gegen);

  VLog gegen_int (FRACTIONS+1, Log(0));
  VLog gegen_int_err (FRACTIONS+1, Log(0));

  set_gegen_integral(gegen_int, gegen_int_err);

  ifstream alpha_beta_f (argv[4]);
  ifstream f (argv[5]);
  ofstream g (argv[6]);
  g << scientific << setprecision(10);

  params pa(hpa);
  read_params(alpha_beta_f, pa, hpa, trs[topology]);

  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      int *q = new int;
      f >> re->first >> re->second >> *q;
      res.push_back(re);
      qs.push_back(q);
    }

  calc_t(pa, hpa, trs[topology]);
  calc_n(pa, hpa, trs[topology]);
  g << calc_llik(res, qs, pa, hpa, trs[topology], gegen, gegen_int) << endl;
  
  for (int i=0; i<n; ++i)
    {
      delete res[i];
      delete qs[i];
    }

  alpha_beta_f.close();
  f.close();
  g.close();

  return 0;
}
