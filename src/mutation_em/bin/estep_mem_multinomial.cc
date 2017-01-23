#include <fstream>
#include <algorithm>
#include <numeric>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf.h>
#include <iomanip>
#include "../../../util/enumtree_mutation_em_multinomial.hh"
#include "../setting.hh"
#include <xmmintrin.h>
using namespace std;

#define calc_sigmoid(x) ((tanh((x)/2.0) + 1.0) / 2.0)
#define calc_asigmoid(x) (2.0 * atanh(2.0 * (x) - 1.0))

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

void read_vf(std::ifstream& f, VVVLog& vf, hyperparams& hpa)
{
  double a;
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int j=1; j<=hpa.MAX_SUBTYPE; ++j)
        {
          for (int s=0; s<=FRACTIONS; ++s)
            {
              f >> a;
              vf[i][j][s] = Log(a);
            }
        }
    }
}

void write_vf(std::ofstream& f, VVVLog& vf, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int j=1; j<=hpa.MAX_SUBTYPE; ++j)
        {
          for (int s=0; s<=FRACTIONS; ++s)
            {
              f << (int) vf[i][j][s].get_sign() << "\t" << vf[i][j][s].get_val() << endl;
            }
          f << endl;
        }
      f << endl << endl;
    }
  f << endl << endl << endl;
}

void read_params(std::ifstream& f, params& pa, hyperparams& hpa)
{
  double a;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f >> a;
      pa.pa[i]->x = a;
    }

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      f >> a;
      pa.pa[i]->y = a;
    }

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
    {
      Log a = _subtypes[i].n * _subtypes[i].x;
      if (a < Log(0)) a = Log(0);
      else if (a > Log(1)) a = Log(1);
      
      prod *= a.take_pow(_state.m[i]);
    }
  
  prod *= mult_nor;
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      Log a = _subtypes[i].n * (Log(1) - _subtypes[i].x);
      if (a < Log(0)) a = Log(0);
      else if (a > Log(1)) a = Log(1);
      
      prod *= a.take_pow(_state.M[i]);
    }
  
  return prod;
}

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

              if (mu < Log(0)) mu = Log(0);
              else if (mu > Log(1)) mu = Log(1);
              
              xq_sum += vf[q][eldest_ch_index][xq] * mu.take_pow(re.first) * (Log(1)-mu).take_pow(re.second - re.first);
            }
          clear_x(_subtypes, hpa);
          
          h_sum += _subtypes[q].omega[eldest_ch_number] * xq_sum;
        }
      sum += _subtypes[q].xi * h_sum;
    }
  return sum;
}

Log calc_lik(subtypes& _subtypes, params& pa, hyperparams& hpa, READS& res, VVVLog& vf)
{
  Log lik (1);
  
  for (int k=0; k<res.size(); ++k)
    {
      Log lik_k = responsibility_partition(_subtypes, pa, hpa, *(res[k]), vf);
      lik *= lik_k;
    }
  
  return lik;
}

void responsibility_m(states& _states, state _state, subtypes& _subtypes, params& pa_old, hyperparams& hpa, READ& re, int current_m, int current_M, int i, Log& denominator, VVVLog& vf, int eldest_ch_index, int eldest_ch_number)
{
  if (i >= hpa.MAX_SUBTYPE)
    {
      _state.m[i] = current_m;
      _state.M[i] = current_M;
      state* st = new state(_state);
      st->resp = responsibility_num(_state, _subtypes, pa_old, hpa, re, vf, eldest_ch_index, eldest_ch_number) / denominator;
      
      // if (!st->resp.iszero()) // might be negative due to the accumulation of numerical error
      if (Log(0) < st->resp)
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

#define d_xi_v(tr, v_q, q, i) (((q) == (i)) ? ((tr)[q].xi / v_q * (Log(1) - (tr)[q].xi)) : (-(tr)[q].xi / v_q * (tr)[q].xi))

#define d_v_t(tr, v_i, i) (v_i / (tr)[i].t)

#define d_xi_t(tr, v_q, v_i, q, i) (d_xi_v(tr, v_q, q, i) * d_v_t(tr, v_i, i))

#define d_v_n(tr, v_i, i, N) (v_i / (tr)[i].n * (Log(1) - (N * (tr)[i].n).take_log_Log().inverse()))

#define d_xi_n(tr, v_q, v_i, q, i, N) (d_xi_v(tr, v_q, q, i) * d_v_n(tr, v_i, i, N))

double d_llik(READS& res, params& pa_old, VVVLog& vf_old, params& pa_new, params& grad, hyperparams& hpa, subtypes& tr_old, subtypes& tr, VVVLog& vf_new, VVVLog& dtvf_new, VVVLog& dthvf_new, VVVLog& dnvf_new)
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
      states _states;
      responsibility_all(_states, tr_old, pa_old, hpa, *(res[k]), vf_old); // calculate states using pa_old bug! Need to use tr calculated based on pa_old!

      Log llik_k = Log(0);
      VLog dt_k (hpa.MAX_SUBTYPE + 1, 0);
      VLog dn_k (hpa.MAX_SUBTYPE + 1, 0);

      for (states::iterator _st = _states.begin(); _st != _states.end(); ++_st)
        {
          state* st = *_st;
          int eldest_ch_index = calc_eldest_child_index(tr[st->q], st->h);
          int eldest_ch_number = calc_eldest_child_number(tr[st->q], st->h);
          tr[st->q].x = Log(((double) st->xq) / FRACTIONS);
          calc_child_x(tr[st->q], hpa, st->h);

          Log l_resp_n = responsibility_num(*st, tr, pa_new, hpa, *res[k], vf_new, eldest_ch_index, eldest_ch_number);
          if (l_resp_n < Log(0)) l_resp_n = Log(0);
          
          llik_k += st->resp * responsibility_num(*st, tr, pa_new, hpa, *res[k], vf_new, eldest_ch_index, eldest_ch_number).take_log_Log();

          VLog dt_st (hpa.MAX_SUBTYPE + 1, Log(0));
          VLog dn_st (hpa.MAX_SUBTYPE + 1, Log(0));

          for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
            {
              dt_st[i] += d_xi_t(tr, v[st->q], v[i], st->q, i) / tr[st->q].xi + d_lnh_q_t(hpa, tr, i, st->q, st->h);
              dn_st[i] += d_xi_n(tr, v[st->q], v[i], st->q, i, Log(CELL_MAX)) / tr[st->q].xi + Log(st->m[i] + st->M[i]) / tr[i].n;
            }

          if (!(vf_new[st->q][eldest_ch_index][st->xq] > Log(0))) // if vf_new[st->q][eldest_ch_index][st->xq] = 0, even if vf[st->q][eldest_ch_index][st->xq] > 0
            {
              return 1; // 1 is the sign of param estimation fails ()
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

      for (int l=0; l<_states.size(); ++l)
        {
          delete _states[l];
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

Log d_n_s(params& pa, int i, int j) // corrected
{
  if (i == j)
    return pa.pa[i]->n * (Log(1) - pa.pa[i]->n/(Log(1) - pa.pa[0]->n));
  
  return -pa.pa[i]->n * pa.pa[j]->n/(Log(1) - pa.pa[0]->n);
}

int main(int argc, char** argv)
{
  cout << scientific << setprecision(10);
  cerr << scientific;
  
  // feenableexcept(FE_INVALID);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  gsl_set_error_handler_off ();
  
  if (argc != 10)
    {
      cerr << "usage: ./estep_mem_multinomial subtype topology n (u n old infile) (u n test infile) (vf_old infile) (vf_test infile) (reads infile) (du dn llik outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  READS res;
  int n, MAX_SUBTYPE, TOTAL_CN, MAX_TREE;

  MAX_SUBTYPE = atoi(argv[1]);
  TOTAL_CN = 2;
  
  trees trs;
  trees_cons(trs, MAX_SUBTYPE);
  MAX_TREE = trs.size();

  trees trs_old;
  trees_cons(trs_old, MAX_SUBTYPE);
  
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

  ifstream pa_old_f(argv[4]);
  ifstream pa_test_f(argv[5]);
  ifstream vf_old_f(argv[6]);
  ifstream vf_test_f(argv[7]);
  ifstream f (argv[8]);
  ofstream g (argv[9]);

  g << scientific << setprecision(10);

  params pa_old(hpa);
  read_params(pa_old_f, pa_old, hpa);
  
  VVVLog vf_old (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  read_vf(vf_old_f, vf_old, hpa);
  calc_subtypes(pa_old, hpa, trs_old[topology]);

  for (int k=0; k<n; ++k)
    {
      READ *re = new READ;
      int a;
      INHERITEDS ih (MAX_SUBTYPE + 1, 0);
      double d;

      f >> re->first >> re->second >> a;
      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        {
          f >> ih[i];
        }
      f >> d; // read tr[q].x
      f >> d; // read mu
      
      res.push_back(re);
    }

  params pa_test(hpa);
  read_params(pa_test_f, pa_test, hpa);
  calc_subtypes(pa_test, hpa, trs[topology]);

  VVVLog vf_test (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dtvf_test (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dthvf_test (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  VVVLog dnvf_test (hpa.MAX_SUBTYPE + 1, VVLog (hpa.MAX_SUBTYPE + 1, VLog (FRACTIONS + 1, Log(0))));
  read_vf(vf_test_f, vf_test, hpa);
  read_vf(vf_test_f, dtvf_test, hpa);
  read_vf(vf_test_f, dthvf_test, hpa);
  read_vf(vf_test_f, dnvf_test, hpa);

  params grad (hpa);
  double qfunc = d_llik(res, pa_old, vf_old, pa_test, grad, hpa, trs_old[topology], trs[topology], vf_test, dtvf_test, dthvf_test, dnvf_test);
  Log lik = calc_lik(trs[topology], pa_test, hpa, res, vf_test);
  write_params(g, grad, hpa);
  g << qfunc << endl;
  g << lik.take_log() << endl;
  
  for (int i=0; i<n; ++i)
    {
      delete res[i];
    }

  pa_old_f.close();
  pa_test_f.close();
  vf_old_f.close();
  vf_test_f.close();
  f.close();
  g.close();

  return 0;
}
