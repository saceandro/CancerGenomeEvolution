#include "setting.hh"
#include <iostream>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv.h>
// #include <xmmintrin.h>
#include <fstream>
#include "../../util/loglib_header.hh"

using namespace std;

#define calc_gamma_i(i, n, t, beta_tilda) (Log(i*(i+1)) * beta_tilda * t / Log(2) / n / (Log(log(CELL_MAX)) + n.take_log_Log()))

typedef void (*myfunc) (int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& partition);
  
typedef struct _vfnumdiff_t
{
  int s;
  int h;
  int q;
  Log n_q;
  Log t_q_h;
  Log beta_tilda_q;
  VVLog gegen;
  VLog gegen_int;
}
  vfnumdiff_t;

typedef struct _vfnumdiff_th
{
  int s;
  int h;
  int q;
  Log n_q;
  Log t_q;
  Log beta_tilda_q;
  VVLog gegen;
  VLog gegen_int;
}
  vfnumdiff_th;

typedef struct _vfnumdiff_n
{
  int s;
  int h;
  int q;
  Log t_q;
  Log t_q_h;
  Log beta_tilda_q;
  VVLog gegen;
  VLog gegen_int;
}
  vfnumdiff_n;

void write_vector(ofstream &f, VLog &a, int n)
{
  for (int i=0; i<n; ++i)
    f << a[i].eval() << " ";
  f << endl;
}

void write_matrix(ofstream &f, VVLog &a, int m, int n)
{
  for (int i=0; i<m; ++i)
    {
      for (int j=0; j<n; ++j)
        f << a[i][j].eval() << " ";
      f << endl;
    }
  f << endl;
}

void set_gegen(VVLog &gegen)
{
  for (int s=1; s<=FRACTIONS; ++s)
    {
      double frac = ((double)s) / FRACTIONS;

      for (int i=1; i<=GEGEN_MAX; ++i)
        gegen[s][i] = Log(gsl_sf_gegenpoly_n(i-1, 1.5, 1.0-2.0*frac));
    }
}

double f (double s, void * params) {
  double z = *(double *) params;
  double f = (1.0 - s) * log(s) / pow(1.0 - 2.0*z*s + s*s, 1.5);
  return f;
}

void set_gegen_integral(VLog &gegen_int, VLog &gegen_int_err)
{
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

  double result, error;

  gsl_function F;
  F.function = &f;

  for (int s=1; s<=FRACTIONS; ++s)
    {
      double z = 1.0 - 2.0 * ((double) s) / FRACTIONS;
      
      F.params = &z;

      double gegen_int_double = 0;
      double gegen_int_err_double = 0;
      
      gsl_integration_qags (&F, 0, 1.0, 0, 1e-7, 1000, w, &gegen_int_double, &gegen_int_err_double);
      gegen_int[s] = Log(gegen_int_double);
      gegen_int_err[s] = Log(gegen_int_err_double);
    }
  gsl_integration_workspace_free (w);
}

void variant_fraction(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& partition)
{
  Log x_q = Log(((double) s) / FRACTIONS);
  Log Nnq = Log(CELL_MAX) * n_q;
  Log lNnq = Nnq.take_log_Log();

  if (h == 0)
    {
      Log part_acc = Log(0);
      for (int i=GEGEN_MAX; i>0; --i)
        {
          Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
          
          double e1_diff = gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval());
          if (fabs(e1_diff) <= DBL_MIN)
            continue;

          part_acc += Log((2.0*i + 1.0) / i / (i+1.0)) *
            (Log(g.take_log() + g.eval() + log(e1_diff), 1) - Log(1.0) + (-g * (Nnq - Log(1.0))).take_exp()/Nnq);
        }
      partition = t_q * (Log(1.0) - (Log(1.0) - Log(1.0)/Nnq)/lNnq) + Log(2.0) * n_q / beta_tilda_q * ((lNnq - Log(1.0)) * (Log(1.0) - Log(1.0)/Nnq) + part_acc);
      
      if (s == 0)
        {
          cerr << "err: variant fraction <= 0 even though h=0" << endl;
          exit(EXIT_FAILURE);
        }
      
      else if (s < FRACTIONS)
        {
          Log acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              
              double e1_diff = gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval());
              if (fabs(e1_diff) <= DBL_MIN)
                continue;

              acc += 
                Log((2.0*i + 1.0) / i / (i+1.0)) * gegen[s][i] *
                (Log(g.take_log() + g.eval() + log(e1_diff), 1) - Log(1.0) + Log(1.0)/g + (-g * (Nnq - Log(1.0))).take_exp()/Nnq);
            }
          numerator =  Log(2.0) * n_q / beta_tilda_q * ((Log(1.0) - Log(1.0)/Nnq)/x_q + Log(4.0)*n_q/beta_tilda_q/t_q * lNnq * gegen_int[s] + Log(2.0)*acc);
        }

      else // s == FRACTIONS
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

              double e1_diff = gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval());
              if (fabs(e1_diff) <= DBL_MIN)
                continue;
              
              Log a = Log((2.0*i + 1.0) / i / (i+1.0)) *
                (Log(g.take_log() + g.eval() + log(e1_diff), 1) - Log(1.0) + (-g * (Nnq - Log(1.0))).take_exp()/Nnq);

              if (i % 2 == 1)
                a = -a;

              acc += a;
            }
          
          numerator = t_q * (Log(1.0) - (Log(1.0) - Log(1.0)/Nnq)/lNnq) + Log(2.0) * n_q / beta_tilda_q * (Log(1.0)/Nnq - Log(1.0) + acc);
        }
    }
  
  else // h > 0
    {
      partition = t_q - t_q_h;
      
      Log beki = (t_q_h / t_q * lNnq).take_exp();
      Log beki_1 = ((t_q_h/ t_q  - 1.0) * lNnq).take_exp();
      
      if (0 < s && s < FRACTIONS)
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              acc +=
                Log((2.0*i + 1.0) / i / (i+1.0)) * gegen[s][i] * ( (-g * (beki - Log(1.0))).take_exp() - (-g * (Nnq - Log(1.0))).take_exp() );
            }
          numerator = Log(4.0) * acc / Log(CELL_MAX) / beta_tilda_q;
        }

      else if (s == FRACTIONS)
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              Log a =
                Log((2.0*i + 1.0) / i / (i+1.0)) * ( (-g * (beki - Log(1.0))).take_exp() - (-g * (Nnq - Log(1.0))).take_exp() );
              
              if (i % 2 == 1)
                a = -a;
              
              acc += a;
            }
          numerator = t_q / lNnq * ( Log(1.0) - beki_1 ) + Log(2.0) * acc / Log(CELL_MAX) / beta_tilda_q;
        }

      else // s == 0
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              acc +=
                Log((2.0*i + 1.0) / i / (i+1.0)) * ( (-g * (beki - Log(1.0))).take_exp() - (-g * (Nnq - Log(1.0))).take_exp() );
            }
          numerator = t_q - t_q_h - t_q / lNnq * ( Log(1.0) - beki_1 ) - Log(2.0) * acc / Log(CELL_MAX) / beta_tilda_q;
        }
    }

  if (numerator < Log(0))
    {
      cerr << "err: variant_fraction < 0" << endl;
      exit(EXIT_FAILURE);
    }
}

// double func_t(double t_q, void *params)
// {
//   vfnumdiff_t *p = (vfnumdiff_t*) params;

//   Log numerator = Log(0);
//   Log partition = Log(0);
  
//   variant_fraction(p->s, p->h, p->q, p->n_q, Log(t_q), p->t_q_h, p->beta_tilda_q, p->gegen, p->gegen_int, numerator, partition);
//   return (numerator / partition).eval();
// }

// double func_th(double t_q_h, void *params)
// {
//   vfnumdiff_th *p = (vfnumdiff_th*) params;

//   Log numerator = Log(0);
//   Log partition = Log(0);
  
//   variant_fraction(p->s, p->h, p->q, p->n_q, p->t_q, Log(t_q_h), p->beta_tilda_q, p->gegen, p->gegen_int, numerator, partition);
//   return (numerator / partition).eval();
// }

// void d_t_variant_fraction_numeric(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, double* result, double* abserr)
// {
//   gsl_function F;

//   F.function = &func_t;

//   vfnumdiff_t dt_vf;
//   dt_vf.s = s;
//   dt_vf.h = h;
//   dt_vf.q = q;
//   dt_vf.n_q = n_q;
//   dt_vf.t_q_h = t_q_h;
//   dt_vf.beta_tilda_q = beta_tilda_q;
//   dt_vf.gegen = gegen;
//   dt_vf.gegen_int = gegen_int;

//   F.params = &dt_vf;
              
//   gsl_deriv_backward(&F, t_q.eval(), 1e-7, result, abserr);
// }

// void d_th_variant_fraction_numeric(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, double* result, double* abserr)
// {
//   gsl_function F;

//   F.function = &func_th;

//   vfnumdiff_th dth_vf;
//   dth_vf.s = s;
//   dth_vf.h = h;
//   dth_vf.q = q;
//   dth_vf.n_q = n_q;
//   dth_vf.t_q = t_q;
//   dth_vf.beta_tilda_q = beta_tilda_q;
//   dth_vf.gegen = gegen;
//   dth_vf.gegen_int = gegen_int;

//   F.params = &dth_vf;
              
//   gsl_deriv_backward(&F, t_q_h.eval(), 1e-7, result, abserr);
// }

void variant_fraction_partition(int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, VLog& vf, VLog& vf_numerator, VLog& vf_denominator, Log& partition)
{
  partition = Log(0);
  
  Log interval = Log(1.0) / Log(FRACTIONS - 1);

  if (h == 0)
    {
      for (int s=1; s<=FRACTIONS; ++s)
        {
          Log numerator = Log(0);
          Log denominator = Log(0);
          variant_fraction(s, 0, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, numerator, denominator);
          vf[s] = numerator / denominator;
          vf_numerator[s] = numerator;
          vf_denominator[s] = denominator;
          
          if (s < FRACTIONS)
            vf[s] *= interval;
        }

      for (int s=1; s<=FRACTIONS; ++s)
        {
          partition += vf[s];
        }

      for (int s=1; s<=FRACTIONS; ++s)
        {
          vf[s] /= partition;
        }
    }
  else
    {
      for (int s=0; s<=FRACTIONS; ++s)
        {
          Log numerator = Log(0);
          Log denominator = Log(0);
          variant_fraction(s, 1, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, numerator, denominator);
          vf[s] = numerator / denominator;
          vf_numerator[s] = numerator;
          vf_denominator[s] = denominator;
          
          if (0 < s && s < FRACTIONS)
            vf[s] *= interval;
        }

      for (int s=0; s<=FRACTIONS; ++s)
        {
          partition += vf[s];
        }
  
      for (int s=0; s<=FRACTIONS; ++s)
        {
          vf[s] /= partition;
        }
    }
}

double func_t_normalized(double t_q, void *params)
{
  vfnumdiff_t *p = (vfnumdiff_t*) params;

  VLog vf(FRACTIONS + 1, Log(0));
  VLog vf_numerator(FRACTIONS + 1, Log(0));
  VLog vf_denominator(FRACTIONS + 1, Log(0));
  Log partition = Log(0);
  variant_fraction_partition(p->h, p->q, p->n_q, Log(t_q), p->t_q_h, p->beta_tilda_q, p->gegen, p->gegen_int, vf, vf_numerator, vf_denominator, partition);
  
  return vf[p->s].eval();
}

double func_th_normalized(double t_q_h, void *params)
{
  vfnumdiff_th *p = (vfnumdiff_th*) params;

  VLog vf(FRACTIONS + 1, Log(0));
  VLog vf_numerator(FRACTIONS + 1, Log(0));
  VLog vf_denominator(FRACTIONS + 1, Log(0));
  Log partition = Log(0);
  variant_fraction_partition(p->h, p->q, p->n_q, p->t_q, Log(t_q_h), p->beta_tilda_q, p->gegen, p->gegen_int, vf, vf_numerator, vf_denominator, partition);
  
  return vf[p->s].eval();
}

double func_n_normalized(double n_q, void *params)
{
  vfnumdiff_n *p = (vfnumdiff_n*) params;

  VLog vf(FRACTIONS + 1, Log(0));
  VLog vf_numerator(FRACTIONS + 1, Log(0));
  VLog vf_denominator(FRACTIONS + 1, Log(0));
  Log partition = Log(0);
  variant_fraction_partition(p->h, p->q, Log(n_q), p->t_q, p->t_q_h, p->beta_tilda_q, p->gegen, p->gegen_int, vf, vf_numerator, vf_denominator, partition);
  
  return vf[p->s].eval();
}

void d_t_variant_fraction_numeric_normalized(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, double* result, double* abserr)
{
  gsl_function F;

  F.function = &func_t_normalized;

  vfnumdiff_t dt_vf;
  dt_vf.s = s;
  dt_vf.h = h;
  dt_vf.q = q;
  dt_vf.n_q = n_q;
  dt_vf.t_q_h = t_q_h;
  dt_vf.beta_tilda_q = beta_tilda_q;
  dt_vf.gegen = gegen;
  dt_vf.gegen_int = gegen_int;

  F.params = &dt_vf;
              
  gsl_deriv_backward(&F, t_q.eval(), 1e-7, result, abserr);
}

void d_th_variant_fraction_numeric_normalized(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, double* result, double* abserr)
{
  gsl_function F;

  F.function = &func_th_normalized;

  vfnumdiff_th dth_vf;
  dth_vf.s = s;
  dth_vf.h = h;
  dth_vf.q = q;
  dth_vf.n_q = n_q;
  dth_vf.t_q = t_q;
  dth_vf.beta_tilda_q = beta_tilda_q;
  dth_vf.gegen = gegen;
  dth_vf.gegen_int = gegen_int;

  F.params = &dth_vf;
              
  gsl_deriv_backward(&F, t_q_h.eval(), 1e-7, result, abserr);
}

void d_n_variant_fraction_numeric_normalized(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, double* result, double* abserr)
{
  gsl_function F;

  F.function = &func_n_normalized;

  vfnumdiff_n dn_vf;
  dn_vf.s = s;
  dn_vf.h = h;
  dn_vf.q = q;
  dn_vf.t_q = t_q;
  dn_vf.t_q_h = t_q_h;
  dn_vf.beta_tilda_q = beta_tilda_q;
  dn_vf.gegen = gegen;
  dn_vf.gegen_int = gegen_int;

  F.params = &dn_vf;
              
  gsl_deriv_backward(&F, n_q.eval(), 1e-7, result, abserr);
}

void d_t_variant_fraction(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& partition)
{
  Log x_q = ((double) s) / FRACTIONS;
  Log Nnq = Log(CELL_MAX) * n_q;
  Log lNnq = Nnq.take_log_Log();

  if (h == 0)
    {
      Log part_acc = Log(0);
      for (int i=GEGEN_MAX; i>0; --i)
        {
          Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

          double e1_diff = gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval());
          if (fabs(e1_diff) <= DBL_MIN)
            continue;

          part_acc += Log(2.0*i + 1.0) *
            ((Log(1.0) + g) * g.take_exp() * Log(e1_diff) - Log(1.0) + (-g * (Nnq - Log(1.0))).take_exp()/Nnq);
        }
      partition = Log(1.0) + ( Log(1.0)/Nnq - Log(1.0) + part_acc ) / lNnq;
      
      if (s == 0)
        {
          cerr << "err: variant fraction <= 0 even though h=0" << endl;
          exit(EXIT_FAILURE);
        }
            
      else if (s < FRACTIONS)
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

              double e1_diff = gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval());
              if (fabs(e1_diff) <= DBL_MIN)
                continue;

              acc += 
                Log(2.0*i + 1.0) * gegen[s][i] *
                (Log((Log(1.0) + g).take_log() + g.eval() + log(e1_diff), 1) - Log(1.0) - Log(1.0)/g/g + (-g * (Nnq - Log(1.0))).take_exp()/Nnq); // corrected the sign of 1/g/g
            }
          numerator =  -Log(8.0) * (n_q/beta_tilda_q/t_q).take_pow(2.0) * lNnq * gegen_int[s] + Log(2.0) / lNnq * acc;
        }

      else // s == FRACTIONS
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

              double e1_diff = gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval());
              if (fabs(e1_diff) <= DBL_MIN)
                continue;

              Log a = Log(2.0*i + 1.0) *
                (Log((Log(1.0) + g).take_log() + g.eval() + log(e1_diff), 1) - Log(1.0) + (-g * (Nnq - Log(1.0))).take_exp()/Nnq);

              if (i % 2 == 1)
                a = -a;
              
              acc += a;
            }
          numerator =  Log(1.0) +  (Log(1.0)/Nnq - Log(1.0) + acc) / lNnq;
        }
    }
  
  else // h > 0
    {
      partition = Log(1.0);
      
      Log beki = (t_q_h / t_q * lNnq).take_exp();
      Log beki_1 = ((t_q_h/ t_q  - 1.0) * lNnq).take_exp();
      
      if (0 < s && s < FRACTIONS)
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              acc +=
                Log(2.0*i + 1.0) * gegen[s][i] *
                ( (Log(1.0) / lNnq / Nnq + (t_q_h / t_q - Log(1.0)/lNnq) * beki_1) * (-g * (beki - Log(1.0))).take_exp() + (Log(1.0) - Log(1.0)/Nnq) / lNnq * (-g * (Nnq - Log(1.0))).take_exp() );
            }
          numerator = Log(2.0) * acc;
        }

      else if (s == FRACTIONS)
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              Log a =
                Log(2.0*i + 1.0) *
                ( (Log(1.0) / lNnq / Nnq + (t_q_h / t_q - Log(1.0)/lNnq) * beki_1) * (-g * (beki - Log(1.0))).take_exp() + (Log(1.0) - Log(1.0)/Nnq) / lNnq * (-g * (Nnq - Log(1.0))).take_exp() );
              
              if (i % 2 == 1)
                a = -a;
              
              acc += a;
            }
          numerator = Log(1.0)/lNnq + ( t_q_h / t_q - Log(1.0)/lNnq ) * beki_1 + acc;
        }

      else // s == 0
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              acc +=
                Log(2.0*i + 1.0) *
                ( (Log(1.0) / lNnq / Nnq + (t_q_h / t_q - Log(1.0)/lNnq) * beki_1) * (-g * (beki - Log(1.0))).take_exp() + (Log(1.0) - Log(1.0)/Nnq) / lNnq * (-g * (Nnq - Log(1.0))).take_exp() );
            }
          numerator = Log(1.0) - Log(1.0) / lNnq  - ( t_q_h / t_q - Log(1.0)/lNnq ) * beki_1 - acc;
        }
    }
}

void d_variant_fraction_partition(myfunc d_x_variant_fraction, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, VLog& vf, VLog& vf_numerator, VLog& vf_denominator, Log& vf_partition, VLog& dvf)
{
  Log partition = Log(0);
  
  Log interval = Log(1.0) / Log(FRACTIONS - 1);

  if (h == 0)
    {
      for (int s=1; s<=FRACTIONS; ++s)
        {
          Log numerator = Log(0);
          Log denominator = Log(0);
          d_x_variant_fraction(s, 0, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, numerator, denominator);
          if (s < FRACTIONS)
            {
              dvf[s] = (numerator / vf_denominator[s] - vf_numerator[s] / vf_denominator[s] * denominator / vf_denominator[s]) * interval;
            }
          else
            {
              dvf[s] = numerator / vf_denominator[s] - vf_numerator[s] / vf_denominator[s] * denominator / vf_denominator[s];
            }
        }

      for (int s=1; s<=FRACTIONS; ++s)
        {
          partition += dvf[s];
        }

      for (int s=1; s<=FRACTIONS; ++s)
        {
          dvf[s] = (dvf[s] - vf[s] * partition) / vf_partition;
        }
    }
  else
    {
      for (int s=0; s<=FRACTIONS; ++s)
        {
          Log numerator = Log(0);
          Log denominator = Log(0);
          d_x_variant_fraction(s, 1, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, numerator, denominator);
          if (0 < s && s < FRACTIONS)
            {
              dvf[s] = (numerator / vf_denominator[s] - vf_numerator[s] / vf_denominator[s] * denominator / vf_denominator[s]) * interval;
            }
          else
            {
              dvf[s] = numerator / vf_denominator[s] - vf_numerator[s] / vf_denominator[s] * denominator / vf_denominator[s];
            }
        }

      for (int s=0; s<=FRACTIONS; ++s)
        {
          partition += dvf[s];
        }
  
      for (int s=0; s<=FRACTIONS; ++s)
        {
          dvf[s] = (dvf[s] - vf[s] * partition) / vf_partition;
        }
    }
}

void d_variant_fraction_all(myfunc d_x_variant_fraction, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, VLog& vf, VLog& dvf)
{
  VLog vf_numerator(FRACTIONS + 1, Log(0));
  VLog vf_denominator(FRACTIONS + 1, Log(0));
  Log vf_partition = Log(0);
  
  variant_fraction_partition(h, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, vf, vf_numerator, vf_denominator, vf_partition);
  d_variant_fraction_partition(d_x_variant_fraction, h, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, vf, vf_numerator, vf_denominator, vf_partition, dvf);
}

void d_th_variant_fraction(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& denominator)
{
  denominator = Log(-1);
  
  if (h == 0)
    {
      numerator = Log(0);
      return;
    }
  
  Log x_q = ((double) s) / FRACTIONS;
  Log Nnq = Log(CELL_MAX) * n_q;
  Log lNnq = Nnq.take_log_Log();
  
  // below is case for h > 0
  
  Log beki = (t_q_h / t_q * lNnq).take_exp();
  Log beki_1 = ((t_q_h/ t_q  - 1.0) * lNnq).take_exp();
      
  if (0 < s && s < FRACTIONS)
    {
      Log acc = Log(0);
      for (int i=GEGEN_MAX; i>0; --i)
        {
          Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
          acc +=
            Log(2.0*i + 1.0) * gegen[s][i] *
            (-g * (beki - Log(1.0))).take_exp();
        }
      numerator = -Log(2.0) * beki_1 * acc;
    }

  else if (s == FRACTIONS)
    {
      Log acc = Log(0);
      for (int i=GEGEN_MAX; i>0; --i)
        {
          Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
          Log a =
            Log(2.0*i + 1.0) *
            (-g * (beki - Log(1.0))).take_exp();
              
          if (i % 2 == 1)
            a = -a;
          
          acc += a;
        }
      numerator = -beki_1 * (Log(1.0) + acc);
    }

  else // s == 0
    {
      Log acc = Log(0);
      for (int i=GEGEN_MAX; i>0; --i)
        {
          Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
          acc +=
            Log(2.0*i + 1.0) *
            (-g * (beki - Log(1.0))).take_exp();
        }
      numerator = -Log(1.0) + beki_1 * (Log(1.0) + acc);
    }
}

void d_n_variant_fraction(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& denominator)
{
  Log x_q = ((double) s) / FRACTIONS;
  Log Nnq = Log(CELL_MAX) * n_q;
  Log lNnq = Nnq.take_log_Log();

  if (h == 0)
    {
      Log part_acc = Log(0);
      for (int i=GEGEN_MAX; i>0; --i)
        {
          Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

          double e1_diff = gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval());
          if (fabs(e1_diff) <= DBL_MIN)
            continue;

          part_acc += Log(2.0*i + 1.0) *
            (((Log(1) + lNnq).inverse() + g) * g.take_exp() * Log(e1_diff) - Log(1) + (Log(1) + lNnq.inverse()).inverse()/g + (-g * (Nnq - Log(1))).take_exp()/Nnq);
        }
  
      denominator = Log(2) / beta_tilda_q * (-Nnq.inverse() + lNnq) + t_q / n_q / lNnq / lNnq - t_q/n_q/lNnq * (Log(1) + lNnq.inverse()) * ( Nnq.inverse() + part_acc );
      
      if (s == 0)
        {
          cerr << "err: variant fraction <= 0 even though h=0" << endl;
          exit(EXIT_FAILURE);
        }
            
      else if (s < FRACTIONS)
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

              double e1_diff = gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval());
              if (fabs(e1_diff) <= DBL_MIN)
                continue;
              
              acc += 
                Log(2.0*i + 1.0) * gegen[s][i] *
                (((Log(1) + lNnq).inverse() + g) * g.take_exp() * Log(e1_diff) - Log(1) + (Log(1) + lNnq.inverse()).inverse()/g - (Log(2) + lNnq.inverse())/(Log(1) + lNnq.inverse())/g/g + (-g * (Nnq - Log(1))).take_exp()/Nnq);
            }
          
          numerator = Log(2)/beta_tilda_q * (x_q.inverse() + Log(4)/beta_tilda_q/t_q * n_q * (Log(2)*lNnq + Log(1)) * gegen_int[s])
            -Log(2)*t_q/n_q/lNnq * (Log(1) + lNnq.inverse()) * acc;
        }

      else // s == FRACTIONS
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

              double e1_diff = gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval());
              if (fabs(e1_diff) <= DBL_MIN)
                continue;

              Log a = Log(2.0*i + 1.0) *
                (((Log(1) + lNnq).inverse() + g) * g.take_exp() * Log(e1_diff) - Log(1) + (Log(1) + lNnq.inverse()).inverse()/g + (-g * (Nnq - Log(1))).take_exp()/Nnq);

              if (i % 2 == 1)
                a = -a;
              
              acc += a;
            }
          numerator =  -Log(2)/beta_tilda_q + t_q/n_q/lNnq/lNnq - t_q/n_q/lNnq * (Log(1) + lNnq.inverse()) * (Nnq.inverse() + acc);
        }
    }
  
  else // h > 0
    {
      denominator = Log(0);
      
      Log beki = (t_q_h / t_q * lNnq).take_exp();
      Log beki_1 = ((t_q_h/ t_q  - 1.0) * lNnq).take_exp();
      
      if (0 < s && s < FRACTIONS)
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              acc +=
                Log(2.0*i + 1.0) * gegen[s][i] *
                ( ( (Log(1) + lNnq.inverse()) * (Nnq.inverse() - beki_1) + t_q_h / t_q * beki_1) * (-g * (beki - Log(1))).take_exp()
                  - (Log(1) - (Log(1) + lNnq.inverse()) * (Log(1) - Nnq.inverse())) * (-g * (Nnq - Log(1))).take_exp());
            }
          numerator = -Log(2) / lNnq * t_q / n_q * acc;
        }

      else if (s == FRACTIONS)
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              Log a =
                Log(2.0*i + 1.0) *
                ( ( (Log(1) + lNnq.inverse()) * (Nnq.inverse() - beki_1) + t_q_h / t_q * beki_1) * (-g * (beki - Log(1))).take_exp()
                  - (Log(1) - (Log(1) + lNnq.inverse()) * (Log(1) - Nnq.inverse())) * (-g * (Nnq - Log(1))).take_exp());
              
              if (i % 2 == 1)
                a = -a;
              
              acc += a;
            }
          numerator = -t_q / n_q / lNnq * ( (Log(1) - beki_1) / lNnq + (t_q_h / t_q - Log(1)) * beki_1 + acc );
        }

      else // s == 0
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              acc +=
                Log(2.0*i + 1.0) *
                ( ( (Log(1) + lNnq.inverse()) * (Nnq.inverse() - beki_1) + t_q_h / t_q * beki_1) * (-g * (beki - Log(1))).take_exp()
                  - (Log(1) - (Log(1) + lNnq.inverse()) * (Log(1) - Nnq.inverse())) * (-g * (Nnq - Log(1))).take_exp());
            }
          numerator = t_q / n_q / lNnq * ( (Log(1) - beki_1) / lNnq + (t_q_h / t_q - Log(1)) * beki_1 + acc );
        }
    }
}

// int main(int argc, char **argv)
// {
//   if (argc != 8)
//     {
//       cerr << "usage: ./variant_fraction (outfile) h q n_q t_q t_q_h beta_tilda_q" << endl;
//       return 0;
//     }
  
//   _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  
//   gsl_set_error_handler_off ();

//   VVLog gegen;
//   gegen.assign(FRACTIONS+1, VLog(GEGEN_MAX+1, Log(0)));

//   set_gegen(gegen);

//   VLog gegen_int (FRACTIONS+1, Log(0));
//   VLog gegen_int_err (FRACTIONS+1, Log(0));

//   set_gegen_integral(gegen_int, gegen_int_err);

//   ofstream of(argv[1]);
//   int h = atoi(argv[2]);
//   int q = atoi(argv[3]);
//   Log n_q = Log(atof(argv[4]));
//   Log t_q = Log(atof(argv[5]));
//   Log t_q_h = Log(atof(argv[6]));
//   Log beta_tilda_q = Log(atof(argv[7]));

//   VLog vf(FRACTIONS + 1, Log(0));
//   VLog dtvf(FRACTIONS + 1, Log(0));
//   VLog dthvf(FRACTIONS + 1, Log(0));
//   VLog dnvf(FRACTIONS + 1, Log(0));

//   d_variant_fraction_all(d_t_variant_fraction, h, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, vf, dtvf);
//   d_variant_fraction_all(d_th_variant_fraction, h, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, vf, dthvf);
//   d_variant_fraction_all(d_n_variant_fraction, h, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, vf, dnvf);

//   int s = 0;
//   if (h == 0)
//     {
//       s = 1;
//     }

//   for (; s<=FRACTIONS; ++s)
//     {
//       double dt_result = 0;
//       double dt_abserr = 0;
//       double dth_result = 0;
//       double dth_abserr = 0;
//       double dn_result = 0;
//       double dn_abserr = 0;
//       d_t_variant_fraction_numeric_normalized(s, h, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, &dt_result, &dt_abserr);
//       d_th_variant_fraction_numeric_normalized(s, h, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, &dth_result, &dth_abserr);
//       d_n_variant_fraction_numeric_normalized(s, h, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, &dn_result, &dn_abserr);
//       of << ((double) s) / FRACTIONS << "\t" << vf[s].eval() << "\t" << dtvf[s].eval() << "\t" << dt_result << "\t" << dthvf[s].eval() << "\t" << dth_result << "\t" << dnvf[s].eval() << "\t" << dn_result << endl;
//     }

//   s = 0;
//   if (h == 0)
//     {
//       s = 1;
//     }

//   Log acc = Log(0);
//   Log dt_acc = Log(0);
//   Log dth_acc = Log(0);
//   Log dn_acc = Log(0);
  
//   for (; s<=FRACTIONS; ++s)
//     {
//       acc += vf[s];
//       dt_acc += dtvf[s];
//       dth_acc += dthvf[s];
//       dn_acc += dnvf[s];
//     }
//   of << "sum" << "\t" << acc.eval() << "\t" << dt_acc.eval() << "\t" << "\t" << dth_acc.eval() << "\t" << "\t" << dn_acc.eval() << endl;
  
//   return 0;
// }
