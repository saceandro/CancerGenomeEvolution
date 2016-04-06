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

class diff
{
public:
  READS res;
  hyperparams hpa;
};

void write_params(std::ofstream& f, params& pa)
{
  for (int c=1; c<pa.pi.size(); ++c)
    {
      f << pa.pi[c] << "\t";
    }
  f << endl;

  for (int c=1; c<pa.pi.size(); ++c)
    {
      for (int d=1; d<=c; ++d)
        {
          f << pa.kappa[c][d] << "\t";
        }
      f << endl;
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

double calc_llik(READS& res, params& pa, hyperparams& hpa)
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
      
      responsibility_numerator_all(*res[k], sts, st, pa, hpa, 0);
      llik += log(responsibility_partition(*res[k], sts, pa, hpa));

      delete_states(sts);
      sts.clear();
    }

  return llik;
}

double d_llik(READS& res, params& pa, params& grad_by_param, hyperparams& hpa)
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
      
      responsibility_numerator_all(*res[k], sts, st, pa, hpa, 0);
      llik += log(responsibility_partition(*res[k], sts, pa, hpa));

      for (unsigned int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (states::iterator it = sts.begin(); it != sts.end(); ++it)
            {
              int bin_l = calc_bin_pi(**it, hpa, l);
              grad_by_param.pi[l] += (*it)->resp * bin_l;
            }

          for (unsigned int r=1; r<=l; ++r)
            {
              for (states::iterator it = sts.begin(); it != sts.end(); ++it)
                {
                  int bin_lr = calc_bin_kappa(**it, hpa, l, r);
                  grad_by_param.kappa[l][r] += (*it)->resp * bin_lr;
                }
            }
        }
      
      delete_states(sts);
      sts.clear();
    }

  return llik;
}

double dx_vec(Vdouble& vec, unsigned int l, unsigned int s)
{
  if (l == s)
    return vec[l] * (1.0 - vec[l]);

  return -vec[l] * vec[s];
}

// double dx_kappa(VVdouble& kappa, unsigned int l, unsigned int r, unsigned int t)
// {
//   return dx_vec(kappa[l], r, t);
// }

void calc_params(const gsl_vector* x, params& pa, hyperparams& hpa)
{
  double m = gsl_vector_get(x, 0);
  for (unsigned int l=2; l<=hpa.TOTAL_CN; ++l)
    {
      double s = gsl_vector_get(x, l-1);
      if (m < s) m = s;
    }
  
  double sum = 0;
  for (unsigned int l=1; l<=hpa.TOTAL_CN; ++l)
    {
      pa.pi[l] = exp(gsl_vector_get(x, l-1) - m);
      sum += pa.pi[l];
    }

  for (unsigned int l=1; l<=hpa.TOTAL_CN; ++l)
    pa.pi[l] /= sum;

  for (unsigned int l=1; l<=hpa.TOTAL_CN; ++l)
    {
      m = gsl_vector_get(x, hpa.TOTAL_CN + (l-1) * l / 2);
      for (unsigned int r=2; r<=l; ++r)
        {
          double s = gsl_vector_get(x, hpa.TOTAL_CN + (l-1) * l / 2 + r-1);
          if (m < s) m = s;
        }
      
      sum = 0;
      
      for (unsigned int r=1; r<=l; ++r)
        {
          pa.kappa[l][r] = exp(gsl_vector_get(x, hpa.TOTAL_CN + (l-1) * l / 2 + r-1) - m);
          sum += pa.kappa[l][r];
        }

      for (unsigned int r=1; r<=l; ++r)
        pa.kappa[l][r] /= sum;
    }
}

void init_gsl_vector(gsl_vector *x, hyperparams& hpa)
{
  for (int s=1; s<=hpa.TOTAL_CN; ++s)
    {
      gsl_vector_set(x, s-1, 0);
    }

  for (int l=1; l<=hpa.TOTAL_CN; ++l)
    {
      for (int t=1; t<=l; ++t)
        {
          gsl_vector_set(x, hpa.TOTAL_CN + (l-1) * l / 2 + t-1, 0);
        }
    }
}

double my_f (const gsl_vector *v, void *par)
{
  diff *p = (diff *) par;

  params pa;
  init_params(pa, p->hpa);
  
  calc_params(v, pa, p->hpa);

  double llik = calc_llik(p->res, pa, p->hpa);

  return -llik;
}

/* The gradient of f, df = (df/dx). */
void my_df (const gsl_vector *v, void *par, gsl_vector *df)
{
  diff *p = (diff *) par;

  int K = p->res.size();
  
  params pa;
  init_params(pa, p->hpa);
  calc_params(v, pa, p->hpa);
  
  params grad_by_param;
  init_params(grad_by_param, p->hpa);

  double llik = d_llik(p->res, pa, grad_by_param, p->hpa);

  
  for (int s=1; s<=p->hpa.TOTAL_CN; ++s)
    {
      double grad = grad_by_param.pi[s] - pa.pi[s] * K * (p->hpa.MAX_SUBTYPE + 1);
      gsl_vector_set(df, s-1, -grad);
    }

  for (int l=1; l<=p->hpa.TOTAL_CN; ++l)
    {
      for (int t=1; t<=l; ++t)
        {
          double grad = grad_by_param.kappa[l][t] - pa.kappa[l][t] * grad_by_param.pi[l];
          gsl_vector_set(df, p->hpa.TOTAL_CN + (l-1) * l / 2 + t-1, -grad);
        }
    }
}

/* Compute both f and df together. */
void my_fdf (const gsl_vector *x, void *par, double *f, gsl_vector *df)     
{
  *f = my_f(x, par);

  my_df(x, par, df);
}

void minimize(diff& di, params& pa)
{
  size_t iter = 10;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  my_func.n = di.hpa.TOTAL_CN * (di.hpa.TOTAL_CN + 3) / 2;
  my_func.f = my_f;
  my_func.df = my_df;
  my_func.fdf = my_fdf;
  my_func.params = &di;

  x = gsl_vector_alloc (my_func.n);

  init_gsl_vector(x, di.hpa);
  
  //cerr << scientific;
  
  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc (T, my_func.n);

  gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.001, 1e-5);

  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

      if (status)
        {
          if (status == GSL_ENOPROG)
            {
              //cerr << "No more improvement can be made for current estimate" << endl << endl;
              calc_params(s->x, pa, di.hpa);
            }
          break;
        }
      
      calc_params(s->x, pa, di.hpa);
      write_params((ofstream&)cout, pa);
      cout << -s->f << endl << endl;
      
      status = gsl_multimin_test_gradient (s->gradient, 1e-5);
      if (status == GSL_SUCCESS)
        {
          //cerr << "Minimum found at:" << endl;
        }

      //cerr << "iter: " << iter << endl;
      //cerr << "kappa\tgrad: " << endl;
      //for (int r=1; r<=st.TOTAL_CN; ++r)
        //cerr << kappa[r] << "\t" << -gsl_vector_get(s->gradient, r-1) << endl;

      //cerr << "llik: " << -s->f << endl << endl;
    }
  while (status == GSL_CONTINUE && iter < 100);

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);
}

int main(int argc, char** argv)
{
  cout << scientific;
  // cerr << scientific;
  
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 6)
    {
      cerr << "usage: ./mixture_total_variant_grad max_subtype total_cn n (reads) (pi kappa outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  diff di;
  unsigned int n;
  
  di.hpa.MAX_SUBTYPE = atoi(argv[1]);
  di.hpa.TOTAL_CN = atoi(argv[2]);
  n = atoi(argv[3]);
  
  ifstream f (argv[4]);
  ofstream g (argv[5]);
  
  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      f >> re->first >> re->second;
      di.res.push_back(re);
    }

  params pa_new;
  init_params(pa_new, di.hpa);

  g << scientific;
  minimize(di, pa_new);
  write_params(g, pa_new);

  for (int i=0; i<n; ++i)
    delete di.res[i];
  
  f.close();
  g.close();

  return 0;
}
