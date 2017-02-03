#include <fstream>
#include <algorithm>
#include <numeric>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf.h>
#include <iomanip>
#include "../../../util/enumtree_mutation_em.hh"
#include "../setting.hh"
#include <xmmintrin.h>
using namespace std;

#define calc_sigmoid(x) ((tanh((x)/2.0) + 1.0) / 2.0)
#define calc_asigmoid(x) (2.0 * atanh(2.0 * (x) - 1.0))

typedef void (*myfunc) (int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& partition);

typedef vector<int> Vint;
typedef vector<double> Vdouble;
typedef vector<Vdouble> VVdouble;
typedef vector<bool> Vbool;
typedef vector<Vbool> VVbool;
typedef pair<int, int> READ;
typedef vector<READ*> READS;
typedef vector<VVLog> VVVLog;
typedef vector<int> INHERITEDS;
typedef pair<int*, INHERITEDS*> Q;
typedef vector< Q* > QS;

class diff
{
public:
  READS& res;
  statess& stss;
  hyperparams& hpa;
  subtypes& tr;
  VVLog& gegen;
  VLog& gegen_int;
  Log purity;
  
  diff (READS& _res, statess& _stss, hyperparams& _hpa, subtypes& _tr, VVLog& _gegen, VLog& _gegen_int, Log& _purity) : res(_res), stss(_stss), hpa(_hpa), tr(_tr), gegen(_gegen), gegen_int(_gegen_int), purity(_purity) {}
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

void read_params(std::ifstream& f, params& pa, hyperparams& hpa)
{
  double a;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f >> a;
      pa.pa[i]->u = Log(a);
    }

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f >> a;
      pa.pa[i]->n = Log(a);
    }
}

void read_params_double(std::ifstream& f, gsl_vector* x, hyperparams& hpa)
{
  double a;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f >> a;
      gsl_vector_set(x, i-1, a);
    }

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f >> a;
      gsl_vector_set(x, hpa.MAX_SUBTYPE + i-1, a);
    }
}

void convert_params(params& pa, gsl_vector* x, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      gsl_vector_set(x, i-1, calc_asigmoid(pa.pa[i]->u.eval()));
    }

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      gsl_vector_set(x, hpa.MAX_SUBTYPE + i-1, pa.pa[i]->n.take_log());
    }
}

void write_params(std::ofstream& f, params& pa, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f << pa.pa[i]->u.eval() << "\t";
    }
  f << endl;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f << pa.pa[i]->n.eval() << "\t";
    }
  f << endl;
}

void write_params_double(std::ofstream& f, gsl_vector* x, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f << gsl_vector_get(x, i-1) << "\t";
    }
  f << endl;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f << gsl_vector_get(x, hpa.MAX_SUBTYPE + i-1) << "\t";
    }
  f << endl;
}

void copy_params(params& pa, params& target, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      target.pa[i]->u = pa.pa[i]->u;
    }
  
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      target.pa[i]->n = pa.pa[i]->n;
    }
}

void calc_params(const gsl_vector* x, params& pa, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->u = Log(calc_sigmoid(gsl_vector_get(x, i-1)));
    }

  Log sum = Log(0);
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->n = Log(gsl_vector_get(x, hpa.MAX_SUBTYPE + i-1)).take_exp();
      sum += pa.pa[i]->n;
    }

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->n = (Log(1) - pa.pa[0]->n) * pa.pa[i]->n / sum;
    }
}

void generate_params(params& pa, hyperparams& hpa, subtypes& tr, gsl_rng* rng)
{
  for (int i=0; i<1024; ++i) // for appropriate random number generation
    gsl_rng_uniform(rng);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    pa.pa[i]->u = Log(gsl_rng_uniform(rng));

  Vdouble w (hpa.MAX_SUBTYPE + 1, 0);
  gsl_ran_dirichlet(rng, hpa.MAX_SUBTYPE, &hpa.gamma[1], &w[1]);
  
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->n = Log(w[i]);
    }

  calc_t(pa, hpa, tr);
  calc_n(pa, hpa, tr);
}

int main(int argc, char** argv)
{
  cout << scientific << setprecision(10);
  cerr << scientific;
  
  // feenableexcept(FE_INVALID);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  gsl_set_error_handler_off ();
  
  if (argc != 5)
    {
      cerr << "usage: ./calc_reverse_param subtype topology (x y u n outfile) seed" << endl;
      exit(EXIT_FAILURE);
    }

  int MAX_SUBTYPE, TOTAL_CN, MAX_TREE;

  MAX_SUBTYPE = atoi(argv[1]);
  TOTAL_CN = 2;
  
  trees trs;
  trees_cons(trs, MAX_SUBTYPE);
  MAX_TREE = trs.size();

  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  // T = gsl_rng_default; // default is mt19937
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  int seed = atoi(argv[4]);
  gsl_rng_set(r, seed);
  
  hyperparams hpa (MAX_SUBTYPE, TOTAL_CN, MAX_TREE);

  int topology = atoi(argv[2]);
  
  ofstream pa_gsl_f (argv[3]);

  pa_gsl_f << scientific << setprecision(10);

  gsl_vector* x = gsl_vector_alloc(2*hpa.MAX_SUBTYPE);

  params pa_old(hpa);
  generate_params(pa_old, hpa, trs[topology], r);

  convert_params(pa_old, x, hpa);
  write_params_double(pa_gsl_f, x, hpa);
  write_params(pa_gsl_f, pa_old, hpa);

  pa_gsl_f.close();

  return 0;
}
