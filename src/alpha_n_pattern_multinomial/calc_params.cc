#include <fstream>
#include <algorithm>
#include <numeric>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf.h>
#include <iomanip>
#include "../../util/enumtree_wf_n_r_pattern_multinomial.hh"
#include "setting.hh"
#include <xmmintrin.h>
using namespace std;

#define calc_sigmoid(x) ((tanh((x)/2.0) + 1.0) / 2.0)

typedef void (*myfunc) (int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& partition);

extern void variant_fraction_partition(int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, VLog& vf, VLog& vf_numerator, VLog& vf_denominator, Log& partition);
extern void d_variant_fraction_all(myfunc d_x_variant_fraction, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, VLog& vf, VLog& dvf);
extern void d_t_variant_fraction(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& partition);
extern void d_th_variant_fraction(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& denominator);
extern void d_n_variant_fraction(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& denominator);
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
typedef vector<int> INHERITEDS;
typedef pair<int*, INHERITEDS*> Q;
typedef vector< Q* > QS;

class diff
{
public:
  READS& res;
  QS& qs;
  hyperparams& hpa;
  subtypes& tr;
  VVLog& gegen;
  VLog& gegen_int;
  Log purity;
  
  diff (READS& _res, QS& _qs, hyperparams& _hpa, subtypes& _tr, VVLog& _gegen, VLog& _gegen_int, Log& _purity) : res(_res), qs(_qs), hpa(_hpa), tr(_tr), gegen(_gegen), gegen_int(_gegen_int), purity(_purity) {}
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

void write_t_n(std::ostream& f, subtypes& st, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    f << st[i].t.eval() << "\t";
  f << endl;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    f << st[i].n.eval() << "\t";
  f << endl << endl;
}

void calc_child_x(subtype& st, hyperparams&hpa, int pat)
{
  for (int i=0; i<st.children.size(); ++i)
    {
      if ((pat >> (st.children.size()-1-i)) % 2 == 1)
        {
          mark_inherited(st.children[i]);
        }
    }
}

int calc_eldest_child_index(subtype& st, int pat)
{
  for (int i=0; i<st.children.size(); ++i)
    {
      if ((pat >> (st.children.size()-1-i)) % 2 == 1)
        {
          return st.children[i]->index;
        }
    }
  return st.index;
}

void clear_x(subtypes& tr, hyperparams&hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      tr[i].x = Log(0);
    }
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

double calc_mu(subtypes& tr, hyperparams& hpa)
{
  Log sum = Log(0);
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      sum += tr[i].n * tr[i].x;
    }

  return sum.eval() / 2.0;
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

Log d_mu_n(subtypes& st, hyperparams& hpa, int i)
{
  return st[i].x / Log(2);
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

double calc_llik(READS& res, QS& qs, params& pa, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int)
{
  int K;
  K = res.size();

  VVVLog vf (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog vf_numerator (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog vf_denominator (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVLog partition (hpa.MAX_SUBTYPE + 1, VLog (hpa.MAX_SUBTYPE + 1, Log(0)));

  for (int q=1; q<=hpa.MAX_SUBTYPE; ++q)
    {
      variant_fraction_partition(0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q][q], vf_numerator[q][q], vf_denominator[q][q], partition[q][q]);
      for (std::vector<subtype*>::iterator ch=tr[q].children.begin(); ch!=tr[q].children.end(); ++ch)
        {
          int ch_index = (*ch)->index;
          variant_fraction_partition(1, q, tr[q].n, tr[q].t, tr[ch_index].t, Log(BETA_TILDA), gegen, gegen_int, vf[q][ch_index], vf_numerator[q][ch_index], vf_denominator[q][ch_index], partition[q][ch_index]);
          // variant_fraction_partition(0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q], vf_numerator[q], vf_denominator[q], partition[q]);
        }
    }
  
  Log lik = Log(1);
  
  for (int k=0; k<K; ++k)
    {
      int index = *(qs[k]->first);
      Log lik_k = Log(0);
      
      for (int pat=0; pat<pow(2,tr[index].children.size()); ++pat)
        {
          int eldest_ch_index = calc_eldest_child_index(tr[index], pat);

          Log lik_k_pat = Log(0);
          int s = 0;
          if (pat == 0)
            {
              s = 1;
            }
          
          for (; s<=FRACTIONS; ++s)
            {
              tr[index].x = Log(((double) s) / FRACTIONS);
              calc_child_x(tr[index], hpa, pat);
              double mu = calc_mu(tr, hpa);
              lik_k_pat += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * vf[index][eldest_ch_index][s];
            }
          clear_x(tr, hpa);

          lik_k += lik_k_pat / Log(pow(2,tr[index].children.size()));
        }
      lik *= lik_k;
    }
  
  return lik.take_log();
}

double d_llik(READS& res, QS& qs, params& pa, params& grad, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int)
{
  int K;
  K = res.size();

  VVVLog vf (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dtvf (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dthvf (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dnvf (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));

  for (int q=1; q<=hpa.MAX_SUBTYPE; ++q)
    {
      d_variant_fraction_all(d_t_variant_fraction, 0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q][q], dtvf[q][q]);
      d_variant_fraction_all(d_n_variant_fraction, 0, q, tr[q].n, tr[q].t, Log(0), Log(BETA_TILDA), gegen, gegen_int, vf[q][q], dnvf[q][q]);
      
      for (std::vector<subtype*>::iterator ch=tr[q].children.begin(); ch!=tr[q].children.end(); ++ch)
        {
          int ch_index = (*ch)->index;
          d_variant_fraction_all(d_t_variant_fraction, 1, q, tr[q].n, tr[q].t, tr[ch_index].t, Log(BETA_TILDA), gegen, gegen_int, vf[q][ch_index], dtvf[q][ch_index]);
          d_variant_fraction_all(d_th_variant_fraction, 1, q, tr[q].n, tr[q].t, tr[ch_index].t, Log(BETA_TILDA), gegen, gegen_int, vf[q][ch_index], dthvf[q][ch_index]);
          d_variant_fraction_all(d_n_variant_fraction, 1, q, tr[q].n, tr[q].t, tr[ch_index].t, Log(BETA_TILDA), gegen, gegen_int, vf[q][ch_index], dnvf[q][ch_index]);
        }
    }

  Log lik = Log(1);

  VLog d_t (hpa.MAX_SUBTYPE + 1, Log(0));
  
  for (int k=0; k<K; ++k)
    {
      int index = *(qs[k]->first);

      Log lik_k = Log(0);
      Log d_t_k = Log(0);
      VLog d_th_k (hpa.MAX_SUBTYPE + 1, Log(0));
      VLog d_n_k (hpa.MAX_SUBTYPE + 1, Log(0));

      for (int pat=0; pat<pow(2,tr[index].children.size()); ++pat)
        {
          Log lik_k_pat = Log(0);
          Log d_t_k_pat = Log(0);
          VLog d_th_k_pat (hpa.MAX_SUBTYPE + 1, Log(0));
          VLog d_n_k_pat (hpa.MAX_SUBTYPE + 1, Log(0));

          int eldest_ch_index = calc_eldest_child_index(tr[index], pat);

          int s = 0;
          if (pat == 0)
            {
              s = 1;
            }

          for (; s<=FRACTIONS; ++s)
            {
              tr[index].x = Log(((double) s) / FRACTIONS);
              calc_child_x(tr[index], hpa, pat);

              double mu = calc_mu(tr, hpa);
              lik_k_pat += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * vf[index][eldest_ch_index][s];
              d_t_k_pat += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * dtvf[index][eldest_ch_index][s];
              d_th_k_pat[eldest_ch_index] += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * dthvf[index][eldest_ch_index][s];

              for (int i=1; i<=hpa.MAX_SUBTYPE; ++i) // start from subtype 1
                d_n_k_pat[i] += d_mu_n(tr, hpa, i) * d_bin_mu(*res[k], mu) * vf[index][eldest_ch_index][s];

              d_n_k_pat[index] += Log(log_binomial_pdf(res[k]->first, mu, res[k]->second), 1) * dnvf[index][eldest_ch_index][s];
            }
          
          clear_x(tr, hpa); // implemented until here
      
          lik_k += lik_k_pat / Log(pow(2,tr[index].children.size()));
          d_t_k += d_t_k_pat / Log(pow(2,tr[index].children.size()));
          for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
            d_th_k[i] += d_th_k_pat[i] / Log(pow(2,tr[index].children.size()));
          for (int i=1; i<=hpa.MAX_SUBTYPE; ++i) // start from subtype 1
            d_n_k[i] += d_n_k_pat[i] / Log(pow(2,tr[index].children.size()));
        }
      
      d_t[index] += d_t_k / lik_k;
      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        d_t[i] += d_th_k[i] / lik_k;
      
      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i) // start from subtype 1
        {
          grad.pa[i]->n += d_n_k[i] / lik_k;
        }
      
      lik *= lik_k;
    }

  for (int j=1; j<=hpa.MAX_SUBTYPE; ++j)
    {
      Log acc = Log(0);
      for (int i=j; i<=hpa.MAX_SUBTYPE; ++i)
        {
          acc += d_t[i] * d_t_u(pa, tr, i, j);
        }
      grad.pa[j]->u += acc;
    }

  return lik.take_log();
}

double my_f (const gsl_vector *v, void *par)
{
  diff *p = (diff *) par;

  params pa (p->hpa);
  pa.pa[0]->n = Log(1) - p->purity;
  calc_params(v, pa, p->hpa);

  calc_t(pa, p->hpa, p->tr);
  calc_n(pa, p->hpa, p->tr);
  
  double llik = calc_llik(p->res, p->qs, pa, p->hpa, p->tr, p->gegen, p->gegen_int);

  return -llik;
}

Log d_n_s(params& pa, int i, int j) // corrected
{
  if (i == j)
    return pa.pa[i]->n * (Log(1) - pa.pa[i]->n/(Log(1) - pa.pa[0]->n));
  
  return -pa.pa[i]->n * pa.pa[j]->n/(Log(1) - pa.pa[0]->n);
}

/* The gradient of f, df = (df/dx). */
void my_df (const gsl_vector *v, void *par, gsl_vector *df)
{
  diff *p = (diff *) par;

  int K = p->res.size();
  
  params pa (p->hpa);
  pa.pa[0]->n = Log(1) - p->purity;
  calc_params(v, pa, p->hpa);

  calc_t(pa, p->hpa, p->tr);
  calc_n(pa, p->hpa, p->tr);

  params grad (p->hpa);

  double llik = d_llik(p->res, p->qs, pa, grad, p->hpa, p->tr, p->gegen, p->gegen_int);

  for (int i=1; i<=p->hpa.MAX_SUBTYPE; ++i)
    {
      double gr = (Log(calc_dx_sigmoid(gsl_vector_get(v, i-1))) * grad.pa[i]->u).eval();
      gsl_vector_set(df, i-1, -gr);
    }

  for (int index=1; index<=p->hpa.MAX_SUBTYPE; ++index)
    {
      Log gr = Log(0);
      
      for (int i=1; i<=p->hpa.MAX_SUBTYPE; ++i)
        {
          gr += d_n_s(pa, index, i) * grad.pa[i]->n;
        }
      gsl_vector_set(df, p->hpa.MAX_SUBTYPE + index-1, (-gr).eval());
    }
}

/* Compute both f and df together. */
void my_fdf (const gsl_vector *x, void *par, double *f, gsl_vector *df)
{
  *f = my_f(x, par);

  my_df(x, par, df);
}

void gsl_set_random(gsl_vector* x, hyperparams& hpa, gsl_rng* rng)
{
  for (int i=0; i<1024; ++i) // for appropriet random number generation
    gsl_rng_uniform(rng);

  for (int i=0; i<2*hpa.MAX_SUBTYPE; ++i)
    {
      double a = gsl_rng_uniform(rng);
      gsl_vector_set(x, i, a);
    }
}

double minimize(diff& di, params& pa, ofstream& h, gsl_rng* rng)
{
  size_t iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  my_func.n = 2*di.hpa.MAX_SUBTYPE;
  my_func.f = my_f;
  my_func.df = my_df;
  my_func.fdf = my_fdf;
  my_func.params = &di;

  x = gsl_vector_alloc (my_func.n);
  gsl_set_random(x, di.hpa, rng);

  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc (T, my_func.n);

  gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.001, 1e-7);

  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

      if (status)
        {
          if (status == GSL_ENOPROG)
            {
              cout << "No more improvement can be made for current estimate" << endl << endl;
              pa.pa[0]->n = Log(1) - di.purity;
              calc_params(s->x, pa, di.hpa);
            }
          break;
        }

      pa.pa[0]->n = Log(1) - di.purity;
      calc_params(s->x, pa, di.hpa);
      write_params((ofstream&)cout, pa, di.hpa);
      cout << -s->f << endl << endl;
      h << iter << "\t" << -s->f << endl;
      
      status = gsl_multimin_test_gradient (s->gradient, 1e-6);
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
  while (status == GSL_CONTINUE && iter < 1000);

  double llik = -s->f;
  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);

  return llik;
}

int main(int argc, char** argv)
{
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  gsl_set_error_handler_off ();
  cerr << scientific;
  
  if (argc != 5)
    {
      cerr << "usage: ./calc_params max_subtype topology (u_n infile) (t_n outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  int MAX_SUBTYPE, TOTAL_CN, MAX_TREE;
  
  MAX_SUBTYPE = atoi(argv[1]);
  TOTAL_CN = 2;
  int topology = atoi(argv[2]);
  
  trees tr;
  trees_cons(tr, MAX_SUBTYPE);
  MAX_TREE = tr.size();

  ifstream u_n_in (argv[3]);
  ofstream t_n_out (argv[4]);
  t_n_out << scientific << setprecision(10);
  
  hyperparams hpa (MAX_SUBTYPE, TOTAL_CN, MAX_TREE);
   
  params pa (hpa);
  read_params(u_n_in, pa, hpa);
  
  calc_t(pa, hpa, tr[topology]);
  calc_n(pa, hpa, tr[topology]);
  write_t_n(t_n_out, tr[topology], hpa);
  
  u_n_in.close();
  t_n_out.close();
  
  return 0;
}
