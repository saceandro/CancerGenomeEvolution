#include "setting.hh"
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv.h>
#include <xmmintrin.h>
#include <fstream>
#include "../../util/loglib.hh"

using namespace std;

#define calc_gamma_i(i, n, t, beta_tilda) (Log(i*(i+1)) * beta_tilda * t / Log(2) / n / (Log(log(CELL_MAX)) + n.take_log_Log()))
#define FRACTIONS 10
#define GEGEN_MAX 200
#define BETA_TILDA_MAX 10

typedef double (*myfunc) (int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int);

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
  // cout << "x_q = " << x_q.eval() << endl;
  Log Nnq = Log(CELL_MAX) * n_q;
  Log lNnq = Nnq.take_log_Log();

  if (h == 0)
    {
      Log part_acc = Log(0);
      for (int i=GEGEN_MAX; i>0; --i)
        {
          Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
          
          double e1_diff = gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval());
          if (abs(e1_diff) <= DBL_MIN)
            continue;

          part_acc += Log((2.0*i + 1.0) / i / (i+1.0)) *
            (Log(g.take_log() + g.eval() + log(gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval())), 1) - Log(1.0) + (-g * (Nnq - Log(1.0))).take_exp()/Nnq);
        }
      partition = t_q * (Log(1.0) - (Log(1.0) - Log(1.0)/Nnq)/lNnq) + Log(2.0) * n_q / beta_tilda_q * ((lNnq - Log(1.0)) * (Log(1.0) - Log(1.0)/Nnq) + part_acc);
      // f << "Z = " << partition << endl;
      
      if (s == 0)
        {
          // cout << "err: variant fraction <= 0" << endl;
        }
      
      else if (s < FRACTIONS)
        {
          // cout << "0 < x_q < 1" << endl;
          
          Log acc = 0;
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              
              double e1_diff = gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval());
              if (abs(e1_diff) <= DBL_MIN)
                continue;

              Log a = 
                Log((2.0*i + 1.0) / i / (i+1.0)) * gegen[s][i] *
                (Log(g.take_log() + g.eval() + log(gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval())), 1) - Log(1.0) + Log(1.0)/g + (-g * (Nnq - Log(1.0))).take_exp()/Nnq);
              // cout << a.eval() << endl;

              acc += a;
            }
          numerator =  Log(2.0) * n_q / beta_tilda_q * ((Log(1.0) - Log(1.0)/Nnq)/x_q + Log(4.0)*n_q/beta_tilda_q/t_q * lNnq * gegen_int[s] + Log(2.0)*acc);

          if (numerator.get_sign() < 0)
            {
              // cout << "h0numerator < 0!,\t 2acc = " << (Log(2.0) * acc).eval() << endl;
              // cout << "(1.0 - 1.0/Nnq)/x_q = "  << ((Log(1.0) - Log(1.0)/Nnq)/x_q).eval() << endl;
              // cout << "4.0*n_q/beta_tilda_q/t_q * log(Nnq) * gegen_int[s] = " << (Log(4.0)*n_q/beta_tilda_q/t_q * lNnq * gegen_int[s]).eval() << endl;
              // cout << "4.0*n_q/beta_tilda_q/t_q * log(Nnq)  = " << (Log(4.0)*n_q/beta_tilda_q/t_q * lNnq).eval() << endl;
            }
          
          // cout << "numerator = " << numerator.eval() << endl << endl;

          if (numerator < 0)
            numerator = Log(0);
          
          return;
        }

      else // s == FRACTIONS
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

              double e1_diff = gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval());
              if (abs(e1_diff) <= DBL_MIN)
                continue;
              
              Log a = Log((2.0*i + 1.0) / i / (i+1.0)) *
                (Log(g.take_log() + g.eval() + log(gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval())), 1) - Log(1.0) + (-g * (Nnq - Log(1.0))).take_exp()/Nnq);

              if (i % 2 == 1)
                {
                  // cout << (-a).eval() << endl;
                  acc -= a;
                }
              else
                {
                  // cout << a.eval() << endl;                  
                  acc += a;
                }
            }
          
          numerator = t_q * (Log(1.0) - (Log(1.0) - Log(1.0)/Nnq)/lNnq) + Log(2.0) * n_q / beta_tilda_q * (Log(1.0)/Nnq - Log(1.0) + acc);
          if (numerator.get_sign() < 0)
          //   cout << "h0FRACTIONnumerator < 0!,\t acc = " << acc.eval()  << endl;
          // cout << "numerator = " << numerator.eval() << endl << endl;
          if (numerator.get_sign() < 0)
            numerator = Log(0);
          
          return;
        }
      
    }
  
  else // h > 0
    {
      // f << "Z = " << t_q - t_q_h << endl;
      partition = t_q - t_q_h;
      
      Log beki = (t_q_h / t_q * lNnq).take_exp();
      Log beki_1 = ((t_q_h/ t_q  - 1.0) * lNnq).take_exp();
      
      if (0 < s && s < FRACTIONS)
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              Log a =
                Log((2.0*i + 1.0) / i / (i+1.0)) * gegen[s][i] * ( (-g * (beki - Log(1.0))).take_exp() - (-g * (Nnq - Log(1.0))).take_exp() );
              cout << a.eval() << endl;
              acc += a;
            }
          numerator = Log(4.0) * acc / Log(CELL_MAX) / beta_tilda_q;
          if (numerator.get_sign() < 0)
            cout << "numerator < 0!,\t acc = " << acc.eval()  << endl;
          cout << "numerator = " << numerator.eval() << endl << endl;
          if (numerator.get_sign() < 0)
            numerator = Log(0);
          
          return;
        }

      else if (s == FRACTIONS)
        {
          cout << "x_q = 1" << endl;
          
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              Log a =
                Log((2.0*i + 1.0) / i / (i+1.0)) * ( (-g * (beki - Log(1.0))).take_exp() - (-g * (Nnq - Log(1.0))).take_exp() );
              
              if (i % 2 == 1)
                a = -a;
              
              cout << a.eval() << endl;
              acc += a;
            }
          numerator = t_q / lNnq * ( Log(1.0) - beki_1 ) + Log(2.0) * acc / Log(CELL_MAX) / beta_tilda_q;
          if (numerator.get_sign() < 0)
            cout << "numerator < 0!,\t acc = " << acc.eval() << endl;
          cout << "numerator = " << numerator.eval() << endl << endl;
          if (numerator.get_sign() < 0)
            numerator = Log(0);
          
          return;
        }

      else // s == 0
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              Log a =
                Log((2.0*i + 1.0) / i / (i+1.0)) * ( (-g * (beki - Log(1.0))).take_exp() - (-g * (Nnq - Log(1.0))).take_exp() );
              
              cout << a.eval() << endl;
              acc += a;
            }
          numerator = t_q - t_q_h - t_q / lNnq * ( Log(1.0) - beki_1 ) - Log(2.0) * acc / Log(CELL_MAX) / beta_tilda_q;
          if (numerator.get_sign() < 0)
            cout << "numerator < 0!,\t acc = " << acc.eval() << endl;
          cout << "numerator = " << numerator.eval() << endl << endl;
          if (numerator.get_sign() < 0)
            numerator = Log(0);
          
          return;
        }
    }
  return;
}

double func_t(double t_q, void *params)
{
  vfnumdiff_t *p = (vfnumdiff_t*) params;

  Log numerator = Log(0);
  Log partition = Log(0);
  
  variant_fraction(p->s, p->h, p->q, p->n_q, Log(t_q), p->t_q_h, p->beta_tilda_q, p->gegen, p->gegen_int, numerator, partition);
  return (numerator / partition).eval();
}


void d_t_variant_fraction_numeric(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, double* result, double* abserr)
{
  gsl_function F;

  F.function = &func_t;

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

// void variant_fraction_partition(int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q, VVdouble& gegen, Vdouble& gegen_int, Vdouble& vf, double& partition, ofstream &f)
// {
//   partition = 0;
  
//   double interval = 1.0 / FRACTIONS;

//   if (h == 0)
//     {
//       for (int s=1; s<=FRACTIONS; ++s)
//         {
//           vf[s] = variant_fraction(s, 0, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, f);
//           if (s < FRACTIONS)
//             vf[s] *= interval;
//         }

//       f << "not_partitioned_variant_fraction h0:" << endl;
//       write_vector(f, vf, FRACTIONS+1);
//       cout << endl;      

//       for (int s=1; s<=FRACTIONS; ++s)
//         {
//           partition += vf[s];
//         }

//       f << "vf_partition: " << partition << endl;
  
//       for (int s=1; s<=FRACTIONS; ++s)
//         {
//           vf[s] /= partition;
//         }
//     }
//   else
//     {
//       for (int s=0; s<=FRACTIONS; ++s)
//         {
//           vf[s] = variant_fraction(s, 1, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, f);
//           if (0 < s && s < FRACTIONS)
//             vf[s] *= interval;
//         }

//       f << "not_partitioned_variant_fraction h1:" << endl;
//       write_vector(f, vf, FRACTIONS+1);
//       cout << endl;
      
//       for (int s=0; s<=FRACTIONS; ++s)
  
//       {
//           partition += vf[s];
//         }
//       f << "vf_partition: " << partition << endl;
  
//       for (int s=0; s<=FRACTIONS; ++s)
//         {
//           vf[s] /= partition;
//         }
      
//     }
// }

void d_t_variant_fraction(int s, int h, int q, Log n_q, Log t_q, Log t_q_h, Log beta_tilda_q, VVLog& gegen, VLog& gegen_int, Log& numerator, Log& partition)
{
  Log x_q = ((double) s) / FRACTIONS;
  cout << "x_q = " << x_q.eval() << endl;
  // f << "x_q = " << x_q << endl;
  Log Nnq = Log(CELL_MAX) * n_q;
  Log lNnq = Nnq.take_log_Log();

  if (h == 0)
    {
      Log part_acc = Log(0);
      for (int i=GEGEN_MAX; i>0; --i)
        {
          Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

          double e1_diff = gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval());
          if (abs(e1_diff) <= DBL_MIN)
            continue;

          part_acc += Log(2.0*i + 1.0) *
            ((Log(1.0) + g) * g.take_exp() * Log(gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval())) - Log(1.0) + (-g * (Nnq - Log(1.0))).take_exp()/Nnq);
        }
      partition = Log(1.0) + ( Log(1.0)/Nnq - Log(1.0) + part_acc ) / lNnq;
      // f << "partition = " << partition << endl;
      
      if (s == 0)
        cout << "err: variant fraction <= 0" << endl;
      
      else if (s < FRACTIONS)
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
              // cout << "val(g): " << g.get_val() << endl;
              // cout << "g.eval(): " << g.eval() << endl;
              // cout << scientific << "E1(g) - E1(g*Nnq): " << gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval()) << endl;

              double e1_diff = gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval());
              if (abs(e1_diff) <= DBL_MIN)
                continue;

              cout << "i: " << i << endl;
              
              Log a = 
                Log(2.0*i + 1.0) * gegen[s][i] *
                (Log((Log(1.0) + g).take_log() + g.eval() + log(gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval())), 1) - Log(1.0) - Log(1.0)/g/g + (-g * (Nnq - Log(1.0))).take_exp()/Nnq); // corrected the sign of 1/g/g
              cout << a.eval() << endl << endl;

              acc += a;
            }
          numerator =  -Log(8.0) * (n_q/beta_tilda_q/t_q).take_pow(2.0) * lNnq * gegen_int[s] + Log(2.0) / lNnq * acc;

          cout << "numerator = " << numerator.eval() << endl << endl << endl;
          
          return;
        }

      else // s == FRACTIONS
        {
          Log acc = Log(0);
          for (int i=GEGEN_MAX; i>0; --i)
            {
              Log g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);

              double e1_diff = gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval());
              if (abs(e1_diff) <= DBL_MIN)
                continue;

              cout << "i: " << i << endl;              
              Log a = Log(2.0*i + 1.0) *
                (Log((Log(1.0) + g).take_log() + g.eval() + log(gsl_sf_expint_E1(g.eval()) - gsl_sf_expint_E1((g*Nnq).eval())), 1) - Log(1.0) + (-g * (Nnq - Log(1.0))).take_exp()/Nnq);

              if (i % 2 == 1)
                a = -a;
              
              cout << a.eval() << endl;
              acc += a;
            }
          numerator =  Log(1.0) +  (Log(1.0)/Nnq - Log(1.0) + acc) / lNnq;
          cout << "numerator = " << numerator.eval() << endl << endl;
          
          return;
        }
      
    }
  
  // else // h > 0
  //   {
  //     // f << "partition = 1" << endl;
  //     partition = 1.0;
      
  //     double beki = exp(t_q_h / t_q * lNnq);
  //     double beki_1 = exp((t_q_h/ t_q  - 1.0) * lNnq);
      
  //     if (0 < s && s < FRACTIONS)
  //       {
  //         double acc = 0;
  //         for (int i=GEGEN_MAX; i>0; --i)
  //           {
  //             double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
  //             double a =
  //               (2.0*i + 1.0) * gegen[s][i] *
  //               ( (1.0 / lNnq / Nnq + (t_q_h / t_q - 1.0/lNnq) * beki_1) * exp(-g * (beki - 1.0)) + (1.0 - 1.0/Nnq) / lNnq * exp(-g * (Nnq - 1.0)) );
  //             cout << a << endl;
  //             acc += a;
  //           }
  //         numerator = 2.0 * acc;
          
  //         cout << "numerator = " << numerator << endl << endl;
          
  //         return numerator;
  //       }

  //     else if (s == FRACTIONS)
  //       {
  //         cout << "x_q = 1" << endl;
          
  //         double acc = 0;
  //         for (int i=GEGEN_MAX; i>0; --i)
  //           {
  //             double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
  //             double a =
  //               (2.0*i + 1.0) *
  //               ( (1.0 / lNnq / Nnq + (t_q_h / t_q - 1.0/lNnq) * beki_1) * exp(-g * (beki - 1.0)) + (1.0 - 1.0/Nnq) / lNnq * exp(-g * (Nnq - 1.0)) );
              
  //             if (i % 2 == 1)
  //               a *= -1;
              
  //             cout << a << endl;
  //             acc += a;
  //           }
  //         numerator = 1.0/lNnq + ( t_q_h / t_q - 1.0/lNnq ) * beki_1 + acc;

  //         cout << "numerator = " << numerator << endl << endl;
          
  //         return numerator;
  //       }

  //     else // s == 0
  //       {
  //         double acc = 0;
  //         for (int i=GEGEN_MAX; i>0; --i)
  //           {
  //             double g = calc_gamma_i(i, n_q, t_q, beta_tilda_q);
  //             double a =
  //               (2.0*i + 1.0) *
  //               ( (1.0 / lNnq / Nnq + (t_q_h / t_q - 1.0/lNnq) * beki_1) * exp(-g * (beki - 1.0)) + (1.0 - 1.0/Nnq) / lNnq * exp(-g * (Nnq - 1.0)) );
              
  //             cout << a << endl;
  //             acc += a;
  //           }
  //         numerator = 1.0 - 1.0 / lNnq  - ( t_q_h / t_q - 1.0/lNnq ) * beki_1 - acc;

  //         cout << "numerator = " << numerator << endl << endl;
          
  //         return numerator;
  //       }
  //   }
  // return -1;
  return;
}

// void d_variant_fraction_partition(myfunc d_x_variant_fraction, int h, int q, double n_q, double t_q, double t_q_h, double beta_tilda_q, VVdouble& gegen, Vdouble& gegen_int, Vdouble& vf, Vdouble& dvf, double vf_partition, ofstream &f) // needs 
// {
//   double interval = 1.0 / FRACTIONS;
//   double dvf_partition = 0;

//   if (h == 0)
//     {
//       for (int s=1; s<=FRACTIONS; ++s)
//         {
//           dvf[s] = d_x_variant_fraction(s, 0, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, f);
//           if (s < FRACTIONS)
//             dvf[s] *= interval;
//         }

//       f << "not_partitioned_d_x_variant_fraction h0:" << endl;
//       write_vector(f, dvf, FRACTIONS+1);
//       cout << endl;      

//       for (int s=1; s<=FRACTIONS; ++s)
//         {
//           dvf_partition += dvf[s];
//         }
//       f << "dvf_partition: " << dvf_partition << endl;
  
//       for (int s=1; s<=FRACTIONS; ++s)
//         {
//           dvf[s] = (dvf[s] - vf[s] * dvf_partition)/ vf_partition;
//         }
//     }
//   else
//     {
//       for (int s=0; s<=FRACTIONS; ++s)
//         {
//           dvf[s] = d_x_variant_fraction(s, 1, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, f);
//           if (0 < s && s < FRACTIONS)
//             dvf[s] *= interval;
//         }

//       f << "not_partitioned_d_x_variant_fraction h0:" << endl;
//       write_vector(f, dvf, FRACTIONS+1);
//       cout << endl;      

//       for (int s=0; s<=FRACTIONS; ++s)
//         {
//           dvf_partition += dvf[s];
//         }
//       f << "dvf_partition: " << dvf_partition << endl;
  
//       for (int s=0; s<=FRACTIONS; ++s)
//         {
//           dvf[s] = (dvf[s] - vf[s] * dvf_partition)/ vf_partition;
//           // dvf[s] = interval * (dvf[s] * vf_partition - vf[s] * dvf_partition)/ vf_partition / vf_partition;
//         }
//     }
// }

// int main(int argc, char **argv)
// {
//   _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  
//   gsl_set_error_handler_off ();
    
//   VVdouble gegen;
//   gegen.assign(FRACTIONS+1, Vdouble(GEGEN_MAX+1, 0));

//   set_gegen(gegen);

//   Vdouble gegen_int (FRACTIONS+1, 0);
//   Vdouble gegen_int_err (FRACTIONS+1, 0);

//   set_gegen_integral(gegen_int, gegen_int_err);

//   Vdouble vf_h0 (FRACTIONS+1, 0);
//   Vdouble vf_h1 (FRACTIONS+1, 0);

//   Vdouble dvf_h0 (FRACTIONS+1, 0);
//   Vdouble dvf_h1 (FRACTIONS+1, 0);

//   ofstream f(argv[1]);

//   for (int beta_tilda_disc=1; beta_tilda_disc<=BETA_TILDA_MAX; ++beta_tilda_disc)
//     {
//       double beta_tilda_q = ((double) beta_tilda_disc) / BETA_TILDA_MAX;

//       f << "beta_tilda_q: " << beta_tilda_q << endl;
//       cout << "beta_tilda_q: " << beta_tilda_q << endl;

//       double vf_partition_h0 = 0;
//       double vf_partition_h1 = 0;

//       cout << "h = 0" << endl;
//       f << "h = 0" << endl;
//       variant_fraction_partition(0, 4, 0.1, 0.2, 0.05, beta_tilda_q, gegen, gegen_int, vf_h0, vf_partition_h0, f);
//       // write_vector(f, vf_h0, FRACTIONS+1);
//       d_variant_fraction_partition(d_t_variant_fraction, 0, 4, 0.1, 0.2, 0.05, beta_tilda_q, gegen, gegen_int, vf_h0, dvf_h0, vf_partition_h0, f);
//       f << "d_t_result: ";
//       write_vector(f, dvf_h0, FRACTIONS+1);
//       f << endl;

//       cout << "h > 0" << endl;
//       f << "h > 0" << endl;
//       variant_fraction_partition(1, 4, 0.1, 0.2, 0.05, beta_tilda_q, gegen, gegen_int, vf_h1, vf_partition_h1, f);
//       // write_vector(f, vf_h1, FRACTIONS+1);
//       d_variant_fraction_partition(d_t_variant_fraction, 1, 4, 0.1, 0.2, 0.05, beta_tilda_q, gegen, gegen_int, vf_h1, dvf_h1, vf_partition_h1, f);
//       f << "d_t_result: ";
//       write_vector(f, dvf_h1, FRACTIONS+1);

//       f << endl << endl;
//     }
  
//   f.close();
  
//   return 0;
// }

int main(int argc, char **argv)
{
  if (argc != 8)
    {
      cerr << "usage: ./variant_fraction s h q n_q t_q t_q_h beta_tilda_q" << endl;
      return 0;
    }
  
  
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  
  gsl_set_error_handler_off ();

  VVLog gegen;
  gegen.assign(FRACTIONS+1, VLog(GEGEN_MAX+1, Log(0)));

  set_gegen(gegen);

  VLog gegen_int (FRACTIONS+1, Log(0));
  VLog gegen_int_err (FRACTIONS+1, Log(0));

  set_gegen_integral(gegen_int, gegen_int_err);

  int s = atoi(argv[1]);
  int h = atoi(argv[2]);
  int q = atoi(argv[3]);
  Log n_q = Log(atof(argv[4]));
  Log t_q = Log(atof(argv[5]));
  Log t_q_h = Log(atof(argv[6]));
  Log beta_tilda_q = Log(atof(argv[7]));

  Log numerator = 0;
  Log partition = 0;
  // Log d_t_numerator = 0;
  // Log d_t_partition = 0;
  
  variant_fraction(s, h, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, numerator, partition);
  // d_t_variant_fraction(s, h, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, d_t_numerator, d_t_partition);

  Log var_frac = numerator / partition;
  // Log d_t_var_frac = (d_t_numerator * partition - numerator * d_t_partition) / partition / partition;
  // double result = 0;
  // double abserr = 0;
  
  // d_t_variant_fraction_numeric(s, h, q, n_q, t_q, t_q_h, beta_tilda_q, gegen, gegen_int, &result, &abserr);
  
  cout << "var_frac: " << var_frac.eval() << endl;
  // cout << "d_t_var_frac: " << d_t_var_frac.eval() << endl;
  // cout << "d_t_var_frac (numeric): " << result << endl;
  
  return 0;
}
