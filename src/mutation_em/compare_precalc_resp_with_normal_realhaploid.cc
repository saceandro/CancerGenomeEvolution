#include "setting.hh"
#include <fstream>
#include <algorithm>
#include <numeric>
#include <fenv.h>
#include <bitset>
#include <iomanip>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf.h>
#include "../../util/enumtree_mutation_em_multinomial.hh"
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

typedef struct _du
{
  int i;
  READS &res;
  gsl_vector *x;
  statess& stss;
  hyperparams &hpa;
  subtypes& tr;
  VVLog gegen;
  VLog gegen_int;
  Log purity;
  
  _du (int _i, READS& _res, gsl_vector* _x, statess& _stss, hyperparams& _hpa, subtypes& _tr, VVLog& _gegen, VLog& _gegen_int, Log& _purity) : i(_i), res(_res), x(_x), stss(_stss), hpa(_hpa), tr(_tr), gegen(_gegen), gegen_int(_gegen_int), purity(_purity) {}
}
  du;

typedef struct _dn
{
  int i;
  READS &res;
  gsl_vector *x;
  statess& stss;
  hyperparams &hpa;
  subtypes& tr;
  VVLog gegen;
  VLog gegen_int;
  Log purity;
  
  _dn (int _i, READS& _res, gsl_vector* _x, statess& _stss, hyperparams& _hpa, subtypes& _tr, VVLog& _gegen, VLog& _gegen_int, Log& _purity) : i(_i), res(_res), x(_x), stss(_stss), hpa(_hpa), tr(_tr), gegen(_gegen), gegen_int(_gegen_int), purity(_purity) {}
}
  dn;

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

void write_params(std::ofstream& f, params& pa, hyperparams& hpa, subtypes& tr)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f << pa.pa[i]->u.eval() << "\t";
    }
  f << endl;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
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

// Log mult(state& _state, subtypes& _subtypes, params& pa_old, hyperparams& hpa, READ& re, VVVLog& vf, int eldest_ch_index, int eldest_ch_number)
// {
//   Log mult_var = Log(gsl_sf_lnfact(re.first), 1);
//   Log mult_nor = Log(gsl_sf_lnfact(re.second - re.first), 1);

//   for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
//     {
//       mult_var /= Log(gsl_sf_lnfact(_state.m[i]), 1);
//       mult_nor /= Log(gsl_sf_lnfact(_state.M[i]), 1);
//     }

//   Log prod1(1);
//   Log prod2(1);
  
//   for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
//     prod1 *= (_subtypes[i].n * _subtypes[i].x).take_pow(_state.m[i]);

//   for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
//     prod2 *= (_subtypes[i].n * (Log(1) - _subtypes[i].x)).take_pow(_state.M[i]);

//   return prod1 * prod2 * mult_var * mult_nor;
// }

Log responsibility_num(state& _state, subtypes& _subtypes, params& pa_old, hyperparams& hpa, READ& re, VVVLog& vf, int eldest_ch_index, int eldest_ch_number)
{
  Log prod = _subtypes[_state.q].xi * _subtypes[_state.q].omega[eldest_ch_number] * vf[_state.q][eldest_ch_index][_state.xq];

  Log mult_var = Log(gsl_sf_lnfact(re.first), 1);
  Log mult_nor = Log(gsl_sf_lnfact(re.second - re.first), 1);
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      mult_var /= Log(gsl_sf_lnfact(_state.m[i]), 1);
      mult_nor /= Log(gsl_sf_lnfact(_state.M[i]), 1);
    }

  prod *= mult_var;
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    prod *= (_subtypes[i].n * _subtypes[i].x).take_pow(_state.m[i]);
  prod *= mult_nor;
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    prod *= (_subtypes[i].n * (Log(1) - _subtypes[i].x)).take_pow(_state.M[i]);
  
  return prod;
}

// Log responsibility_partition(subtypes& _subtypes, params& pa_old, hyperparams& hpa, READ& re, VVVLog &vf)
// {
//   state st (1, 1, 0, Vint(hpa.MAX_SUBTYPE+1, 0), Vint(hpa.MAX_SUBTYPE+1, 0), Log(0));
  
//   Log sum = Log(0);

//   for (st.q=1; st.q<=hpa.MAX_SUBTYPE; ++st.q)
//     {
//       Log h_sum = Log(0);
      
//       for (st.h=0; st.h<pow(2,_subtypes[st.q].children.size()); ++st.h)
//         {
//           int eldest_ch_index = calc_eldest_child_index(_subtypes[st.q], st.h);
//           int eldest_ch_number = calc_eldest_child_number(_subtypes[st.q], st.h);

//           st.xq = 0;
//           if (st.h==0) st.xq = 1;

//           Log xq_sum = Log(0);
          
//           for (; st.xq<=FRACTIONS; ++st.xq)
//             {
//               _subtypes[st.q].x = Log(((double) st.xq) / FRACTIONS);
//               calc_child_x(_subtypes[st.q], hpa, st.h);

//               // cout << "x1: " << _subtypes[1].x.eval() << "\tx2: " << _subtypes[2].x.eval() << endl;
//               // cout << "nx1: " << (_subtypes[1].n * _subtypes[1].x).eval() << "\tnx2: " << (_subtypes[2].n * _subtypes[2].x).eval() << endl;
//               Log mu = Log(0);

//               for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
//                 mu += _subtypes[i].n * _subtypes[i].x;

//               Log read_sum (0);
//               for (st.m[1]=0; st.m[1]<=re.first; ++st.m[1])
//                 {
//                   for (st.M[1]=0; st.M[1]<=re.second - re.first; ++st.M[1])
//                     {
//                       st.m[2] = re.first - st.m[1];
//                       st.M[2] = re.second - re.first - st.M[1];

//                       read_sum += mult(st, _subtypes, pa_old, hpa, re, vf, eldest_ch_index, eldest_ch_number);
//                     }
//                 }
//               cout << "mu_pow: " << (mu.take_pow(re.first) * (Log(1)-mu).take_pow(re.second - re.first)).eval() << "\t" << "read_sum: " << read_sum.eval()
//                    << "\t"  << "all_product: " << (_subtypes[st.q].xi * _subtypes[st.q].omega[eldest_ch_number] * vf[st.q][eldest_ch_index][st.xq] * mu.take_pow(re.first) * (Log(1)-mu).take_pow(re.second - re.first)).eval() << endl;
//               read_sum = mu.take_pow(re.first) * (Log(1)-mu).take_pow(re.second - re.first);
              
//               xq_sum += vf[st.q][eldest_ch_index][st.xq] * read_sum;
//             }
//           clear_x(_subtypes, hpa);
          
//           h_sum += _subtypes[st.q].omega[eldest_ch_number] * xq_sum;
//         }
//       sum += _subtypes[st.q].xi * h_sum;
//     }
//   return sum;
// }

Log responsibility_partition(subtypes& _subtypes, params& pa_old, hyperparams& hpa, READ& re, VVVLog &vf)
{
  Log sum = Log(0);

  for (int q=1; q<=hpa.MAX_SUBTYPE; ++q)
    {
      Log h_sum = Log(0);
      
      for (int h=0; h<pow(2,_subtypes[q].children.size()); ++h)
        {
          int eldest_ch_index = calc_eldest_child_index(_subtypes[q], h);
          int eldest_ch_number = calc_eldest_child_number(_subtypes[q], h);

          int xq = 0;
          if (h==0) xq = 1;

          Log xq_sum = Log(0);
          
          for (; xq<=FRACTIONS; ++xq)
            {
              _subtypes[q].x = Log(((double) xq) / FRACTIONS);
              calc_child_x(_subtypes[q], hpa, h);
                  
              Log mu = Log(0);

              for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
                mu += _subtypes[i].n * _subtypes[i].x;

              xq_sum += vf[q][eldest_ch_index][xq] * mu.take_pow(re.first) * (Log(1)-mu).take_pow(re.second - re.first);
            }
          clear_x(_subtypes, hpa);
          
          h_sum += _subtypes[q].omega[eldest_ch_number] * xq_sum;
        }
      sum += _subtypes[q].xi * h_sum;
    }
  return sum;
}

void responsibility_m(states& _states, state _state, subtypes& _subtypes, params& pa_old, hyperparams& hpa, READ& re, int current_m, int current_M, int i, Log& denominator, VVVLog& vf, int eldest_ch_index, int eldest_ch_number)
{
  if (i >= hpa.MAX_SUBTYPE)
    {
      _state.m[i] = current_m;
      _state.M[i] = current_M;
      state* st = new state(_state);
      st->resp = responsibility_num(_state, _subtypes, pa_old, hpa, re, vf, eldest_ch_index, eldest_ch_number) / denominator;
      
      if (!st->resp.iszero())
        _states.push_back(st);
    }

  else
    {
      for (_state.m[i]=0; _state.m[i]<=current_m; ++_state.m[i])
        {
          for (_state.M[i]=0; _state.M[i]<=current_M; ++_state.M[i])
            {
              responsibility_m(_states, _state, _subtypes, pa_old, hpa, re, current_m - _state.m[i], current_M - _state.M[i], i+1, denominator, vf, eldest_ch_index, eldest_ch_number);
            }
        }
    }
}

Log responsibility_all(states& _states, subtypes& _subtypes, params& pa_old, hyperparams& hpa, READ& re, VVVLog& vf)
{
  write_t_n((ofstream&)cout, _subtypes, hpa);
  
  Log denominator = responsibility_partition(_subtypes, pa_old, hpa, re, vf);

  state st (1, 1, 0, Vint(hpa.MAX_SUBTYPE+1, 0), Vint(hpa.MAX_SUBTYPE+1, 0), Log(0));
  
  for (st.q=1; st.q<=hpa.MAX_SUBTYPE; ++st.q)
    {
      for (st.h=0; st.h<pow(2,_subtypes[st.q].children.size()); ++st.h)
        {
          int eldest_ch_index = calc_eldest_child_index(_subtypes[st.q], st.h);
          int eldest_ch_number = calc_eldest_child_number(_subtypes[st.q], st.h);

          st.xq = 0;
          if (st.h==0) st.xq = 1;

          for (; st.xq<=FRACTIONS; ++st.xq)
            {
              _subtypes[st.q].x = Log(((double) st.xq) / FRACTIONS);
              calc_child_x(_subtypes[st.q], hpa, st.h);
              
              responsibility_m(_states, st, _subtypes, pa_old, hpa, re, re.first, re.second - re.first, 1, denominator, vf, eldest_ch_index, eldest_ch_number);
            }
          clear_x(_subtypes, hpa);
        }
    }
}

void calc_vf(VVVLog& vf, VVVLog& dtvf, VVVLog& dthvf, VVVLog& dnvf, params& pa, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int)
{
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
        }
    }

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
}

double q_func(READS& res, statess& _statess, params& pa_new, hyperparams& hpa, subtypes& tr, VVVLog& vf_new, VVLog& gegen, VLog& gegen_int)
{
  int K;
  K = res.size();

  Log llik = Log(0);
  
  for (int k=0; k<K; ++k)
    {
      Log llik_k = Log(0);
      
      for (states::iterator _st = _statess[k]->begin(); _st != _statess[k]->end(); ++_st)
        {
          state* st = *_st;
          int eldest_ch_index = calc_eldest_child_index(tr[st->q], st->h);
          int eldest_ch_number = calc_eldest_child_number(tr[st->q], st->h);
          tr[st->q].x = Log(((double) st->xq) / FRACTIONS);
          calc_child_x(tr[st->q], hpa, st->h);

          llik_k += st->resp * responsibility_num(*st, tr, pa_new, hpa, *res[k], vf_new, eldest_ch_index, eldest_ch_number).take_log_Log();
          
          clear_x(tr, hpa);
        }
      
      llik += llik_k;
    }
  
  return llik.eval();
}

//implemented until here

double q_func_for_du(double x_i, void* _du)
{
  du* p = (du*) _du;

  int K = p->res.size();
  
  gsl_vector_set(p->x, p->i-1, x_i);
  params pa (p->hpa);
  pa.pa[0]->n = Log(1) - p->purity;
  calc_params(p->x, pa, p->hpa); // caluculate new pa
  calc_subtypes(pa, p->hpa, p->tr); // calculate new tr

  VVVLog vf_new (p->hpa.MAX_SUBTYPE + 1, VVLog (p->hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dtvf_new (p->hpa.MAX_SUBTYPE + 1, VVLog (p->hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dthvf_new (p->hpa.MAX_SUBTYPE + 1, VVLog (p->hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dnvf_new (p->hpa.MAX_SUBTYPE + 1, VVLog (p->hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));

  calc_vf(vf_new, dtvf_new, dthvf_new, dnvf_new, pa, p->hpa, p->tr, p->gegen, p->gegen_int); // calc vf_new using pa_new

  double ans = q_func(p->res, p->stss, pa, p->hpa, p->tr, vf_new, p->gegen, p->gegen_int);
  
  return ans;
}

double q_func_for_dn(double x_i_j, void* _dn)
{
  dn* p = (dn*) _dn;

  gsl_vector_set(p->x, p->hpa.MAX_SUBTYPE + p->i-1, x_i_j);
  params pa (p->hpa);
  pa.pa[0]->n = Log(1) - p->purity;
  calc_params(p->x, pa, p->hpa);
  calc_subtypes(pa, p->hpa, p->tr);

  VVVLog vf_new (p->hpa.MAX_SUBTYPE + 1, VVLog (p->hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dtvf_new (p->hpa.MAX_SUBTYPE + 1, VVLog (p->hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dthvf_new (p->hpa.MAX_SUBTYPE + 1, VVLog (p->hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dnvf_new (p->hpa.MAX_SUBTYPE + 1, VVLog (p->hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));

  calc_vf(vf_new, dtvf_new, dthvf_new, dnvf_new, pa, p->hpa, p->tr, p->gegen, p->gegen_int);

  double ans = q_func(p->res, p->stss, pa, p->hpa, p->tr, vf_new, p->gegen, p->gegen_int);
  
  return ans;
}

#define d_xi_v(tr, v_q, q, i) (((q) == (i)) ? ((tr)[q].xi / v_q * (Log(1) - (tr)[q].xi)) : (-(tr)[q].xi / v_q * (tr)[q].xi))

#define d_v_t(tr, v_i, i) (v_i / (tr)[i].t)

#define d_xi_t(tr, v_q, v_i, q, i) (d_xi_v(tr, v_q, q, i) * d_v_t(tr, v_i, i))

#define d_v_n(tr, v_i, i, N) (v_i / (tr)[i].n * (Log(1) - (N * (tr)[i].n).take_log_Log().inverse()))

#define d_xi_n(tr, v_q, v_i, q, i, N) (d_xi_v(tr, v_q, q, i) * d_v_n(tr, v_i, i, N))

double d_llik(READS& res, statess& _statess, params& pa_new, params& grad, hyperparams& hpa, subtypes& tr, VVVLog& vf_new, VVVLog& dtvf_new, VVVLog& dthvf_new, VVVLog& dnvf_new, VVLog& gegen, VLog& gegen_int)
{
  int K;
  K = res.size();

  Log llik = Log(0);
  VLog dt (hpa.MAX_SUBTYPE + 1, 0);

  VLog v (hpa.MAX_SUBTYPE + 1, Log(0));

  for (int k=1; k<=hpa.MAX_SUBTYPE; ++k)
    {
      v[k] = tr[k].t * tr[k].n / (Log(CELL_MAX) * tr[k].n).take_log_Log();
    }
  
  for (int k=0; k<K; ++k)
    {
      Log llik_k = Log(0);
      VLog dt_k (hpa.MAX_SUBTYPE + 1, 0);
      VLog dn_k (hpa.MAX_SUBTYPE + 1, 0);
      
      for (states::iterator _st = _statess[k]->begin(); _st != _statess[k]->end(); ++_st)
        {
          state* st = *_st;
          int eldest_ch_index = calc_eldest_child_index(tr[st->q], st->h);
          int eldest_ch_number = calc_eldest_child_number(tr[st->q], st->h);
          tr[st->q].x = Log(((double) st->xq) / FRACTIONS);
          calc_child_x(tr[st->q], hpa, st->h);

          
          llik_k += st->resp * responsibility_num(*st, tr, pa_new, hpa, *res[k], vf_new, eldest_ch_index, eldest_ch_number).take_log_Log();

          VLog dt_st (hpa.MAX_SUBTYPE + 1, Log(0));
          VLog dn_st (hpa.MAX_SUBTYPE + 1, Log(0));

          for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
            {
              dt_st[i] += d_xi_t(tr, v[st->q], v[i], st->q, i) / tr[st->q].xi + d_lnh_q_t(hpa, tr, i, st->q, st->h);
              dn_st[i] += d_xi_n(tr, v[st->q], v[i], st->q, i, Log(CELL_MAX)) / tr[st->q].xi + Log(st->m[i] + st->M[i]) / tr[i].n;
            }

          dt_st[st->q] += dtvf_new[st->q][eldest_ch_index][st->xq] / vf_new[st->q][eldest_ch_index][st->xq];
          dt_st[eldest_ch_index] += dthvf_new[st->q][eldest_ch_index][st->xq] / vf_new[st->q][eldest_ch_index][st->xq];
          dn_st[st->q] += dnvf_new[st->q][eldest_ch_index][st->xq] / vf_new[st->q][eldest_ch_index][st->xq] + d_s_n(hpa, tr, st->q) *  d_h_q_n(hpa, tr, st->q, st->h) / tr[st->q].omega[eldest_ch_number];

          for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
            {
              dt_k[i] += st->resp * dt_st[i];
              dn_k[i] += st->resp * dn_st[i];
            }
          clear_x(tr, hpa);
        }

      llik += llik_k;

      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i) // start from subtype 1
        {
          dt[i] += dt_k[i];
          grad.pa[i]->n += dn_k[i];
        }
    }
  
  for (int j=1; j<=hpa.MAX_SUBTYPE; ++j)
    {
      Log acc = Log(0);
      for (int i=j; i<=hpa.MAX_SUBTYPE; ++i)
        {
          acc += dt[i] * d_t_u(pa_new, tr, i, j);
        }
      grad.pa[j]->u += acc;
    }

  return llik.eval();
}

double calc_rel_err(double x, double y)
{
  return fabs(x - y) / fabs(x);
}

double calc_dx_u_llik_numeric(int i, READS& res, gsl_vector* x, statess& _stss, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int, Log purity)
{
  gsl_function F;
  du _du (i, res, x, _stss, hpa, tr, gegen, gegen_int, purity);

  F.function = &q_func_for_du;
  F.params = &_du;

  double result, abserr;
  gsl_deriv_central(&F, gsl_vector_get(x, i-1), 1e-5, &result, &abserr);
  
  return result;
}

double calc_dx_u_llik_analytic(int i, READS& res, gsl_vector* x, statess& _stss, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int, Log purity)
{
  int K = res.size();

  params pa (hpa);
  pa.pa[0]->n = Log(1) - purity;
  calc_params(x, pa, hpa);
  calc_subtypes(pa, hpa, tr);
  
  VVVLog vf_new (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dtvf_new (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dthvf_new (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dnvf_new (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));

  calc_vf(vf_new, dtvf_new, dthvf_new, dnvf_new, pa, hpa, tr, gegen, gegen_int); // calc vf_new using pa_new

  params grad (hpa);

  double llik = d_llik(res, _stss, pa, grad, hpa, tr, vf_new, dtvf_new, dthvf_new, dnvf_new, gegen, gegen_int);

  return (Log(calc_dx_sigmoid(gsl_vector_get(x, i-1))) * grad.pa[i]->u).eval();
}

double calc_dx_n_llik_numeric(int i, READS& res, gsl_vector* x, statess& _stss, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int, Log purity)
{
  gsl_function F;
  dn _dn (i, res, x, _stss, hpa, tr, gegen, gegen_int, purity);

  F.function = &q_func_for_dn;
  F.params = &_dn;

  double result, abserr;
  gsl_deriv_central(&F, gsl_vector_get(x, hpa.MAX_SUBTYPE + i-1), 1e-5, &result, &abserr);
  
  return result;
}

Log d_n_s(params& pa, int i, int j) // corrected
{
  if (i == j)
    return pa.pa[i]->n * (Log(1) - pa.pa[i]->n/(Log(1) - pa.pa[0]->n));
  
  return -pa.pa[i]->n * pa.pa[j]->n/(Log(1) - pa.pa[0]->n);
}

double calc_dx_n_llik_analytic(int index, READS& res, gsl_vector* x, statess& _stss, hyperparams& hpa, subtypes& tr, VVLog& gegen, VLog& gegen_int, Log purity)
{
  int K = res.size();

  params pa (hpa);
  pa.pa[0]->n = Log(1) - purity;
  calc_params(x, pa, hpa);
  calc_subtypes(pa, hpa, tr);

  VVVLog vf_new (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dtvf_new (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dthvf_new (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dnvf_new (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));

  calc_vf(vf_new, dtvf_new, dthvf_new, dnvf_new, pa, hpa, tr, gegen, gegen_int);

  params grad (hpa);
  
  double llik = d_llik(res, _stss, pa, grad, hpa, tr, vf_new, dtvf_new, dthvf_new, dnvf_new, gegen, gegen_int);
  
  Log gr = Log(0);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      gr += d_n_s(pa, index, i) * grad.pa[i]->n;
    }

  return gr.eval();
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

int main(int argc, char** argv)
{
  cout << scientific;
  cerr << scientific;
  
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  
  gsl_set_error_handler_off ();
  
  if (argc != 9)
    {
      cerr << "usage: ./compare_precalc_resp_with_normal_u max_subtype n step tumor_purity (x y infile) (reads) (u diff outfile) topology" << endl;
      exit(EXIT_FAILURE);
    }

  READS res;
  int n, step, MAX_SUBTYPE, TOTAL_CN, MAX_TREE;

  MAX_SUBTYPE = atoi(argv[1]);
  TOTAL_CN = 2;

  trees trs;
  trees_cons(trs, MAX_SUBTYPE);
  MAX_TREE = trs.size();

  hyperparams hpa (MAX_SUBTYPE, TOTAL_CN, MAX_TREE);
  
  n = atoi(argv[2]);
  step = atoi(argv[3]);
  Log purity = Log(atof(argv[4]));
  
  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  // T = gsl_rng_default; // default is mt19937
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, 1);

  ifstream x_y_old_f (argv[5]);
  ifstream f (argv[6]);
  ofstream ff (argv[7]);
  // ofstream g (argv[7]);

  int a = atoi(argv[8]);

  ff << scientific << setprecision(10);
  // g << scientific << setprecision(10);

  VVLog gegen;
  gegen.assign(FRACTIONS+1, VLog(GEGEN_MAX+1, Log(0)));

  set_gegen(gegen);

  VLog gegen_int (FRACTIONS+1, Log(0));
  VLog gegen_int_err (FRACTIONS+1, Log(0));

  set_gegen_integral(gegen_int, gegen_int_err);

  gsl_vector* x = gsl_vector_alloc(2*hpa.MAX_SUBTYPE);

  for (int k=0; k<n; ++k)
    {
      READ *re = new READ;
      int a;
      INHERITEDS ih (MAX_SUBTYPE + 1, 0);
      double fra;

      f >> re->first >> re->second >> a;
      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        {
          f >> ih[i];
        }
      f >> fra;
      
      res.push_back(re);
    }

  gsl_vector* y = gsl_vector_alloc(2*hpa.MAX_SUBTYPE);
  read_params_double(x_y_old_f, y, hpa);
  for (int i=0; i<2*hpa.MAX_SUBTYPE; ++i)
    {
      gsl_vector_set(x, i, gsl_vector_get(y, i));
    }
  // gsl_set_random(y, hpa, r);

  params pa_old(hpa);
  pa_old.pa[0]->n = Log(1) - purity;
  calc_params(y, pa_old, hpa);
  calc_subtypes(pa_old, hpa, trs[a]);
  
  VVVLog vf (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dtvf (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dthvf (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dnvf (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));

  calc_vf(vf, dtvf, dthvf, dnvf, pa_old, hpa, trs[a], gegen, gegen_int); // calculate vf using pa_old

  statess _statess;
  for (int k=0; k<n; ++k)
    {
      states* _states = new states;
      responsibility_all(*_states, trs[a], pa_old, hpa, *(res[k]), vf); // calculate states using pa_old

      Log sum (0);
      cerr << "number of hidden state: " << _states->size() << endl;
      for (int s=0; s<_states->size(); ++s)
        {
          sum += (*_states)[s]->resp.eval();
          cerr << (*_states)[s]->resp.eval() << endl;
        }
      cerr << endl;
      cerr << sum.eval() << endl;

      for (int s=0; s<_states->size(); ++s)
        {
          if (((*_states)[s]->q == 2) && ((*_states)[s]->h == 0) && ((*_states)[s]->xq == 5))
            {
              for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
                {
                  cerr << (*_states)[s]->m[i] << "\t";
                }
              for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
                {
                  cerr << (*_states)[s]->M[i] << "\t";
                }
              cerr << endl;
            }
        }
      _statess.push_back(_states);
    }
  
  for (int u_index=1; u_index<=hpa.MAX_SUBTYPE; ++u_index)
    {
      ff << "u_index = " << u_index << endl;
  
      for (int i=-step; i<=step; ++i)
        {
          // gsl_set_random(x, hpa, r);
          gsl_vector_set(x, u_index-1, 1.0 * ((double)i) / step);
          double analytic = calc_dx_u_llik_analytic(u_index, res, x, _statess, hpa, trs[a], gegen, gegen_int, purity);
          double num = calc_dx_u_llik_numeric(u_index, res, x, _statess, hpa, trs[a], gegen, gegen_int, purity);
          ff << gsl_vector_get(x, u_index-1) << "\t" << analytic << "\t" << num << "\t" << fabs(analytic - num) << "\t" << calc_rel_err(analytic, num) << endl;
        }
      ff << endl;
    }

  for (int n_index=1; n_index<=hpa.MAX_SUBTYPE; ++n_index)
    {
      ff << "n_index = " << n_index << endl;
      
      for (int i=-step; i<=step; ++i)
        {
          // gsl_set_random(x, hpa, r);
          gsl_vector_set(x, hpa.MAX_SUBTYPE + n_index-1, 1.0 * ((double)i) / step);
          double num = calc_dx_n_llik_numeric(n_index, res, x, _statess, hpa, trs[a], gegen, gegen_int, purity);
          double analytic = calc_dx_n_llik_analytic(n_index, res, x, _statess, hpa, trs[a], gegen, gegen_int, purity);

          // if (fabs(num) > 0)
          ff << gsl_vector_get(x, hpa.MAX_SUBTYPE + n_index-1) << "\t" << analytic << "\t" << num << "\t" << fabs(analytic - num) << "\t" << calc_rel_err(analytic, num) << endl;
          // else
          //   g << gsl_vector_get(x, hpa.MAX_SUBTYPE + n_index-1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
        }
      ff << endl;
    }
  
  for (int k=0; k<n; ++k)
    {
      for (int s=0; s<_statess[k]->size(); ++s) // memory leak! l=0; l<(_statess[k])->size() is accurate
        {
          delete (*(_statess[k]))[s];
        }
      delete _statess[k];
    }

  for (int i=0; i<n; ++i)
    {
      delete res[i];
    }

  f.close();
  ff.close();

  return 0;
}
