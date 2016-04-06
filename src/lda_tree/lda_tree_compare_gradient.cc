#include <fstream>
#include <algorithm>
#include <numeric>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf.h>
#include "../enumtree.hh"
using namespace std;

#define napier exp(1)
// #define logsum(a, b) (((a) > (b)) ? ((a) + log1p(exp((b)-(a)))) : ((b) + log1p(exp((a)-(b)))))
#define calc_sigmoid(x) ((tanh((x)/2.0) + 1.0) / 2.0)

typedef vector<int> Vint;
typedef vector<double> Vdouble;
typedef vector<Vdouble> VVdouble;
typedef vector<bool> Vbool;
typedef vector<Vbool> VVbool;
typedef pair<int, int> READ;
typedef vector<READ*> READS;
typedef vector<VVLog> VVVLog;

// class state 
// {
// public:
//   int k;
//   Vint total_cn;
//   Vint variant_cn;
//   VLog resp_du;
//   Log resp;

//   // use default constructor
//   // state (int, std::vector<int>, std::vector<int>, double);
// };

class state 
{
public:
  int locus;
  subtypes st;
  Log resp;

  // use default constructor
};

typedef std::vector<state*> states;

typedef struct _drho
{
  READS &res;
  gsl_vector *x;
  int a;
  hyperparams &hpa;
  
  _drho (READS& _res, gsl_vector* _x, int _a, hyperparams& _hpa) : res(_res), x(_x), a(_a), hpa(_hpa) {}
}
  drho;

typedef struct _du
{
  READS &res;
  gsl_vector *x;
  int i;
  hyperparams &hpa;
  
  _du (READS& _res, gsl_vector* _x, int _i, hyperparams& _hpa) : res(_res), x(_x), i(_i), hpa(_hpa) {}
}
  du;

typedef struct _dpi
{
  READS &res;
  gsl_vector *x;
  int i;
  int l;
  hyperparams &hpa;
  
  _dpi (READS& _res, gsl_vector* _x, int _i, int _l, hyperparams& _hpa) : res(_res), x(_x), i(_i), l(_l), hpa(_hpa) {}
}
  dpi;

typedef struct _dkappa
{
  READS &res;
  gsl_vector *x;
  int i;
  int l;
  int r;
  hyperparams &hpa;
  
  _dkappa (READS& _res, gsl_vector* _x, int _i, int _l, int _r, hyperparams& _hpa) : res(_res), x(_x), i(_i), l(_l), r(_r), hpa(_hpa) {}
}
  dkappa;

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
  for (int a=0; a<hpa.MAX_TREE; ++a)
    f << pa.rho[a].eval() << "\t";
  f << endl << endl;
  
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    f << pa.pa[i]->u.eval() << "\t";
  f << endl << endl;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    f << pa.pa[i]->t.eval() << "\t";
  f << endl << endl;

  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    f << pa.pa[i]->n.eval() << "\t";
  f << endl << endl;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double sum = 0;
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        f << pa.pa[i]->pi[l].eval() << "\t";
      f << endl;

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (int r=1; r<=l; ++r)
            f << pa.pa[i]->kappa[l][r].eval() << "\t";
          f << endl;
        }
      f << endl;
    }
}

void init_state(state& st, hyperparams& hpa)
{
  st.st.assign(hpa.MAX_SUBTYPE + 1, subtype (0, 0, 0, Log(0), NULL, NULL, vector<subtype*> (0)));
  st.st[0].total_cn = 2;
  st.st[0].variant_cn = 0;
}

void delete_params(params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    delete pa.pa[i];
}

void calc_n(params& pa, hyperparams& hpa)
{
  Log sum;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->n = pa.pa[i]->t;
      sum += pa.pa[i]->t;
    }
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    pa.pa[i]->n /= sum;
}

double calc_mu(subtypes& st, params& pa, hyperparams& hpa)
{
  Log denom;
  Log num;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      denom += pa.pa[i]->n * Log(st[i].total_cn);
      num += pa.pa[i]->n * Log(st[i].variant_cn);
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

Log d_n_t(params& pa, hyperparams& hpa, int i, int j)
{
  if (i == j)
    return pa.pa[i]->n * (Log(1) - pa.pa[i]->n)/ pa.pa[i]->t;
  else
    return -pa.pa[i]->n * pa.pa[i]->n / pa.pa[i]->t;
}

Log d_mu_n(subtypes& st, params& pa, hyperparams& hpa, int j)
{
  Log normal;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    normal += pa.pa[i]->n * Log(st[i].total_cn);

  Log variant;
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    variant += pa.pa[i]->n * Log(st[i].variant_cn);

  Log a = Log(st[j].variant_cn) / normal;
  Log b = Log(st[j].total_cn) / normal * variant / normal;

  return a - b;
}

Log d_rho_x(params& pa, hyperparams& hpa, int a, int b)
{
  if (a == b)
    return pa.rho[a] * (Log(1.0) - pa.rho[b]);
  else
    return -pa.rho[a] * pa.rho[b];
}

void responsibility_numerator(READ& re, states& sts, subtypes& st, params& pa, hyperparams& hpa)
{
  Log product (1);
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    product *= pa.pa[i]->pi[st[i].total_cn] * pa.pa[i]->kappa[st[i].total_cn][st[i].variant_cn];

  state* new_state = new state;
  init_state(*new_state, hpa);
  copy(new_state->st, st);

  double mu = calc_mu(st, pa, hpa);
  new_state->resp = product * Log(log_binomial_pdf(re.first, mu, re.second), 1);

  product *= d_bin_mu(re, mu);
  for (int y=0; y<=hpa.MAX_SUBTYPE; ++y)
    {
      Log sum;
      // for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
      //   sum += d_mu_n(st, pa, hpa, i) * d_n_u(pa, hpa, i, y);

      for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
        {
          Log sum2;
          for (int j=0; j<=hpa.MAX_SUBTYPE; ++j) // modified to be j=0
            sum2 += d_n_t(pa, hpa, i, j) * d_t_u(pa, hpa, st, j, y);
          sum += d_mu_n(st, pa, hpa, i) * sum2;
        }
      
      new_state->st[y].resp_du = product * sum;
    }
  
  sts.push_back(new_state);
}

void responsibility_numerator_all(READ& re, states& sts, subtypes& st, params& pa, hyperparams& hpa, int i)
{
  if (i < hpa.MAX_SUBTYPE)
    {
      for (st[i].total_cn = 1; st[i].total_cn <= hpa.TOTAL_CN; ++st[i].total_cn)
        {
          for (st[i].variant_cn = 1; st[i].variant_cn <= st[i].total_cn; ++st[i].variant_cn)
            {
              responsibility_numerator_all(re, sts, st, pa, hpa, i + 1);
            }
        }
    }

  else
    {
      for (st[i].total_cn = 1; st[i].total_cn <= hpa.TOTAL_CN; ++st[i].total_cn)
        {
          for (st[i].variant_cn = 1; st[i].variant_cn <= st[i].total_cn; ++st[i].variant_cn)
            {
              responsibility_numerator(re, sts, st, pa, hpa);
            }
        }
    }
}

void responsibility_numerator_tree(READ& re, vector<states>& sts, trees& tr, params& pa, hyperparams& hpa)
{
  for (int a=0; a<tr.size(); ++a)
    {
      calc_t(pa, hpa, tr[a]);
      calc_n(pa, hpa);

      responsibility_numerator_all(re, sts[a], tr[a], pa, hpa, 1);
    }
}

Log responsibility_partition(vector<vector<states> >& sts, params& pa, hyperparams& hpa)
{
  Log partition (0);
  VLog u_partition (hpa.MAX_SUBTYPE + 1, 0);
  
  for (int a=0; a<sts[0].size(); ++a)
    {
      Log partition_a (1);
      VLog u_partition_a(hpa.MAX_SUBTYPE + 1, 1);
      for (int k=0; k<sts.size(); ++k)
        {
          Log partition_k (0);
          VLog u_partition_k(hpa.MAX_SUBTYPE + 1, 0);
          for (states::iterator it = sts[k][a].begin(); it != sts[k][a].end(); ++it)
            {
              // cout << (*it)->resp.take_log() << "\t" << (*it)->resp.get_sign() << endl;
              partition_k += (*it)->resp;
              
              for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
                {
                  u_partition_k[i] += (*it)->st[i].resp_du;
                }
            }
          partition_a *= partition_k;
          for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
            {
              u_partition_a[i] *= u_partition_k[i];
            }
        }
      partition += pa.rho[a] * partition_a;
      for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
        {
          u_partition[i] += pa.rho[a] * u_partition_a[i];
        }
    }

  for (int a=0; a<sts[0].size(); ++a)
    {
      for (int k=0; k<sts.size(); ++k)
        {
          for (states::iterator it = sts[k][a].begin(); it != sts[k][a].end(); ++it)
            {
              (*it)->resp /= partition;
              for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
                {
                  (*it)->st[i].resp_du /= partition;
                }
            }
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

// int calc_bin_pi(state& st, hyperparams& hpa, int i, int l)
// {
//   if (st.total_cn[i] == l)
//     return 1;
//   else
//     return 0;
// }

// int calc_bin_kappa(state& st, hyperparams& hpa, int i, int l, int r)
// {
//   if (st.total_cn[i] == l && st.variant_cn[i] == r)
//     return 1;
//   else
//     return 0;
// }
void calc_params(const gsl_vector* x, params& pa, hyperparams& hpa, subtypes& st)
{
  Log sum_rho (0);
  for (int a=0; a<hpa.MAX_TREE; ++a)
    {
      pa.rho[a] = Log(calc_sigmoid(gsl_vector_get(x, a)));
      sum_rho += pa.rho[a];
    }

  for (int a=0; a<hpa.MAX_TREE; ++a)
    pa.rho[a] /= sum_rho;

  
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->u = Log(calc_sigmoid(gsl_vector_get(x, hpa.MAX_TREE + i)));
      // pa.pa[i]->u = Log(gsl_vector_get(x, i));
    }
  
  calc_t(pa, hpa, st);
  calc_n(pa, hpa);
  
  int params_per_subtype = hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double m = gsl_vector_get(x, hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1));
      
      for (int l=2; l<=hpa.TOTAL_CN; ++l)
        {
          double s = gsl_vector_get(x, hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + l-1);
          if (m < s) m = s;
        }
  
      Log sum;
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          pa.pa[i]->pi[l] = Log(gsl_vector_get(x, hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + l-1) - m, 1);
          sum += pa.pa[i]->pi[l];
        }

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        pa.pa[i]->pi[l] /= sum;

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          m = gsl_vector_get(x, hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + hpa.TOTAL_CN + (l-1) * l / 2);
          for (int r=2; r<=l; ++r)
            {
              double s = gsl_vector_get(x, hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + hpa.TOTAL_CN + (l-1) * l / 2 + r-1);
              if (m < s) m = s;
            }
      
          sum = Log(0);
      
          for (int r=1; r<=l; ++r)
            {
              pa.pa[i]->kappa[l][r] = Log(gsl_vector_get(x, hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + hpa.TOTAL_CN + (l-1) * l / 2 + r-1) - m, 1);
              sum += pa.pa[i]->kappa[l][r];
            }

          for (int r=1; r<=l; ++r)
            pa.pa[i]->kappa[l][r] /= sum;
        }
    }
}

double calc_llik(READS& res, params& pa, hyperparams& hpa)
{
  int K;
  K = res.size();
  
  // Log lik (1);
  vector<vector<states> > sts (K, vector<states> (hpa.MAX_TREE));

  trees trs;
  trees_cons(trs, hpa);

  for (int k=0; k<K; ++k)
    {
      // cerr << "k: " << k << endl;
      
      responsibility_numerator_tree(*res[k], sts[k], trs, pa, hpa);
      // Log lik_k = responsibility_partition(*res[k], sts[k], pa, hpa);
      // cout << lik_k.take_log() << endl;
      // lik *= lik_k;
      
      // delete_states(sts);
      // sts.clear();
    }

  Log lik = responsibility_partition(sts, pa, hpa);

  for (int k=0; k<K; ++k)
    {
      for (int a=0; a<sts[0].size(); ++a)
        {
          delete_states(sts[k][a]);
        }
    }

  // for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     for (int l=1; l<=hpa.TOTAL_CN; ++l)
  //       {
  //         lik *= pa.pa[i]->pi[l].take_pow(hpa.alpha[l] - 1.0);
          
  //         for (int r=1; r<=l; ++r)
  //           {
  //             lik *= pa.pa[i]->kappa[l][r].take_pow(hpa.beta[l][r] - 1.0);
  //           }
  //       }
  //   }

  return lik.take_log();
}

double calc_llik_for_drho(double x_a, void* _drho)
{
  drho* p = (drho*) _drho;

  gsl_vector_set(p->x, p->a, x_a);

  trees trs;
  trees_cons(trs, p->hpa);

  params pa;
  init_params(pa, p->hpa);
  calc_params(p->x, pa, p->hpa, trs[0]);

  return calc_llik(p->res, pa, p->hpa);
}

double calc_llik_for_du(double x_i, void* _du)
{
  du* p = (du*) _du;

  gsl_vector_set(p->x, p->hpa.MAX_TREE + p->i, x_i);

  trees trs;
  trees_cons(trs, p->hpa);

  params pa;
  init_params(pa, p->hpa);
  calc_params(p->x, pa, p->hpa, trs[0]);

  return calc_llik(p->res, pa, p->hpa);
}

double calc_llik_for_dpi(double x_il, void* dp)
{
  dpi* p = (dpi*) dp;
  
  int params_per_subtype = p->hpa.TOTAL_CN * (p->hpa.TOTAL_CN + 3) / 2;
  
  gsl_vector_set(p->x, p->hpa.MAX_TREE + p->hpa.MAX_SUBTYPE + 1 + params_per_subtype * (p->i - 1) + p->l - 1, x_il);
  
  trees trs;
  trees_cons(trs, p->hpa);

  params pa;
  init_params(pa, p->hpa);
  calc_params(p->x, pa, p->hpa, trs[0]);

  return calc_llik(p->res, pa, p->hpa);
}

double calc_llik_for_dkappa(double x_ilr, void* dk)
{
  dkappa* p = (dkappa*) dk;

  int params_per_subtype = p->hpa.TOTAL_CN * (p->hpa.TOTAL_CN + 3) / 2;

  gsl_vector_set(p->x, p->hpa.MAX_TREE + p->hpa.MAX_SUBTYPE + 1 + params_per_subtype * (p->i - 1) + p->hpa.TOTAL_CN + (p->l - 1) * p->l / 2 + p->r - 1, x_ilr);

  trees trs;
  trees_cons(trs, p->hpa);

  params pa;
  init_params(pa, p->hpa);
  calc_params(p->x, pa, p->hpa, trs[0]);

  return calc_llik(p->res, pa, p->hpa);
}

double d_llik(READS& res, params& pa, params& grad_by_param, hyperparams& hpa)
{
  int K;
  K = res.size();

  // Log lik (1);
  vector<vector<states> > sts (K, vector<states> (hpa.MAX_TREE));

  trees trs;
  trees_cons(trs, hpa);

  for (int k=0; k<K; ++k)
    {
      // cerr << "k: " << k << endl;
      
      responsibility_numerator_tree(*res[k], sts[k], trs, pa, hpa);
    }

  Log lik = responsibility_partition(sts, pa, hpa);

  // for (int a=0; a<hpa.MAX_TREE; ++a)
  //   {
  //     grad_by_param.rho[a] = Log(1.0);
  //     VLog rho_numerator_k (hpa.MAX_TREE, Log(0));
      
  //     for (int k=0; k<K; ++k)
  //       {
  //         // Log rho_numerator_k(0);
  //         for (states::iterator it = sts[k][a].begin(); it != sts[k][a].end(); ++it)
  //           {
  //             // cout << (*it)->resp.take_log() << "\t" << (*it)->resp.get_sign() << endl;
  //             rho_numerator_k[a] += (*it)->resp;
  //           }
  //         grad_by_param.rho[a] *= rho_numerator_k[a];
  //       }
  //   }

  // for (int a=0; a<hpa.MAX_TREE; ++a)
  //   {
  //     VLog u_numerator_a(hpa.MAX_SUBTYPE + 1, Log(1.0));
  //     for (int k=0; k<K; ++k)
  //       {
  //         VLog u_numerator_k(hpa.MAX_SUBTYPE + 1, Log(0));
  //         for (states::iterator it = sts[k][a].begin(); it != sts[k][a].end(); ++it)
  //           {
  //             // cout << (*it)->resp.take_log() << "\t" << (*it)->resp.get_sign() << endl;
  //             for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
  //               {
  //                 u_numerator_k[i] += (*it)->st[i].resp_du;
  //               }
  //           }

  //         for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
  //           {
  //             u_numerator_a[i] *= u_numerator_k[i];
  //           }
  //       }

  //     for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
  //       {
  //         grad_by_param.pa[i]->u += pa.rho[a] * u_numerator_a[i];
  //       }
  //   }


  for (int a=0; a<hpa.MAX_TREE; ++a)
    {
      grad_by_param.rho[a] = Log(1.0);
      VLog u_numerator_a(hpa.MAX_SUBTYPE + 1, Log(1.0));
      VVLog pi_numerator_a(hpa.MAX_SUBTYPE + 1, VLog (hpa.TOTAL_CN + 1, Log(1.0)));
      VVVLog kappa_numerator_a(hpa.MAX_SUBTYPE + 1, VVLog (hpa.TOTAL_CN + 1, VLog (hpa.TOTAL_CN + 1, Log(1.0))));
      for (int k=0; k<K; ++k)
        {
          VLog rho_numerator_k (hpa.MAX_TREE, Log(0));
          VLog u_numerator_k(hpa.MAX_SUBTYPE + 1, Log(0));
          VVLog pi_numerator_k(hpa.MAX_SUBTYPE + 1, VLog (hpa.TOTAL_CN + 1, Log(0)));
          VVVLog kappa_numerator_k(hpa.MAX_SUBTYPE + 1, VVLog (hpa.TOTAL_CN + 1, VLog (hpa.TOTAL_CN + 1, Log(0))));
          for (states::iterator it = sts[k][a].begin(); it != sts[k][a].end(); ++it)
            {
              // cout << (*it)->resp.take_log() << "\t" << (*it)->resp.get_sign() << endl;

              rho_numerator_k[a] += (*it)->resp;
              
              for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
                {
                  u_numerator_k[i] += (*it)->st[i].resp_du;
                }

              for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
                {
                  for (int l=1; l<=hpa.TOTAL_CN; ++l)
                    {
                      if ((*it)->st[i].total_cn == l)
                        {
                          pi_numerator_k[i][l] += (*it)->resp;

                          for (int r=1; r<=l; ++r)
                            {
                              if ((*it)->st[i].variant_cn == r)
                                {
                                  kappa_numerator_k[i][l][r] += (*it)->resp;
                                }
                            }
                        }
                    }
                }
            }

          grad_by_param.rho[a] *= rho_numerator_k[a];
          
          for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
            {
              u_numerator_a[i] *= u_numerator_k[i];
            }
          
          for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
            {
              for (int l=1; l<=hpa.TOTAL_CN; ++l)
                {
                  pi_numerator_a[i][l] *= pi_numerator_k[i][l];

                  for (int r=1; r<=l; ++r)
                    {
                      kappa_numerator_a[i][l][r] *= kappa_numerator_k[i][l][r];
                    }
                }
            }
          delete_states(sts[k][a]);
        }

      for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
        {
          grad_by_param.pa[i]->u += pa.rho[a] * u_numerator_a[i];
        }

      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        {
          for (int l=1; l<=hpa.TOTAL_CN; ++l)
            {
              grad_by_param.pa[i]->pi[l] += pa.rho[a] * pi_numerator_a[i][l];

              for (int r=1; r<=l; ++r)
                {
                  grad_by_param.pa[i]->kappa[l][r] += pa.rho[a] * kappa_numerator_a[i][l][r];
                }
            }
        }
    }

  // for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     for (int l=1; l<=hpa.TOTAL_CN; ++l)
  //       {
  //         grad_by_param.pa[i]->pi[l] = Log(hpa.alpha[l] - 1.0) + grad_by_param.pa[i]->pi[l];

  //         for (int r=1; r<=l; ++r)
  //           {
  //             grad_by_param.pa[i]->kappa[l][r] = Log(hpa.beta[l][r] - 1.0) + grad_by_param.pa[i]->kappa[l][r];
  //           }
  //       }
  //   }

  // for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     for (int l=1; l<=hpa.TOTAL_CN; ++l)
  //       {
  //         lik *= pa.pa[i]->pi[l].take_pow(hpa.alpha[l] - 1.0);
          
  //         for (int r=1; r<=l; ++r)
  //           {
  //             lik *= pa.pa[i]->kappa[l][r].take_pow(hpa.beta[l][r] - 1.0);
  //           }
  //       }
  //   }

  return lik.take_log();
}

double dx_vec(Vdouble& vec, int l, int s)
{
  if (l == s)
    return vec[l] * (1.0 - vec[l]);

  return -vec[l] * vec[s];
}

// double dx_kappa(VVdouble& kappa, int l, int r, int t)
// {
//   return dx_vec(kappa[l], r, t);
// }


double calc_rel_err(double x, double y)
{
  return fabs(x - y) / fabs(x);
}

double calc_dx_rho_llik_numeric(READS& res, gsl_vector* x, int a, hyperparams& hpa)
{
  gsl_function F;
  drho _drho (res, x, a, hpa);
  _drho.res = res;
  _drho.x = x;
  _drho.a = a;
  _drho.hpa = hpa;

  F.function = &calc_llik_for_drho;
  F.params = &_drho;

  double result, abserr;
  gsl_deriv_central(&F, gsl_vector_get(x, a), 1e-5, &result, &abserr);
  
  return result;
}

double calc_dx_rho_llik_analytic(READS& res, gsl_vector* x, int b, hyperparams& hpa)
{
  int K = res.size();

  trees trs;
  trees_cons(trs, hpa);

  params pa;
  init_params(pa, hpa);
  calc_params(x, pa, hpa, trs[0]);

  params grad_by_param;
  init_params(grad_by_param, hpa);

  double llik = d_llik(res, pa, grad_by_param, hpa);

  Log grad (0);
  cerr << "max_tree: " << hpa.MAX_TREE << endl;
  
  for (int a=0; a<hpa.MAX_TREE; ++a)
    {
      cerr << d_rho_x(pa, hpa, a, b).get_val() << "\t" << d_rho_x(pa, hpa, a, b).get_sign() << "\t" << grad_by_param.rho[a].get_val() << "\t" << grad_by_param.rho[a].get_sign() << endl;
      grad += d_rho_x(pa, hpa, a, b) * grad_by_param.rho[a];
    }

  if (grad.iszero())
    cerr << "dx_rho is 0!" << endl;
  return grad.take_log();
}

double calc_dx_u_llik_numeric(READS& res, gsl_vector* x, int i, hyperparams& hpa)
{
  gsl_function F;
  du _du (res, x, i, hpa);
  _du.res = res;
  _du.x = x;
  _du.i = i;
  _du.hpa = hpa;

  F.function = &calc_llik_for_du;
  F.params = &_du;

  double result, abserr;
  gsl_deriv_central(&F, gsl_vector_get(x, hpa.MAX_TREE + i), 1e-5, &result, &abserr);
  
  return result;
}

double calc_dx_u_llik_analytic(READS& res, gsl_vector* x, int j, hyperparams& hpa)
{
  int K = res.size();

  trees trs;
  trees_cons(trs, hpa);

  params pa;
  init_params(pa, hpa);
  calc_params(x, pa, hpa, trs[0]);

  params grad_by_param;
  init_params(grad_by_param, hpa);

  double llik = d_llik(res, pa, grad_by_param, hpa);
  // return calc_dx_sigmoid(gsl_vector_get(x,j)) * grad_by_param[j]->u.eval();
  return (Log(calc_dx_sigmoid(gsl_vector_get(x, hpa.MAX_TREE + j))) * grad_by_param.pa[j]->u).eval();
  // return (Log(calc_dx_sigmoid(pa[j]->u.eval())) * grad_by_param[j]->u).eval();
  // return grad_by_param[j]->u.eval();

  // double grad = 0;
  // for (int l=1; l<=hpa.TOTAL_CN; ++l)
  //   {
  //     grad += dx_vec(pa[i]->pi, l, s) * grad_by_param[i]->pi[l] / pa[i]->pi[l];
  //   }
  
  // return grad;
}

double calc_dx_pi_llik_numeric(READS& res, gsl_vector* x, int i, int l, hyperparams& hpa)
{
  gsl_function F;
  dpi dp (res, x, i, l, hpa);
  dp.res = res;
  dp.x = x;
  dp.i = i;
  dp.l = l;
  dp.hpa = hpa;

  F.function = &calc_llik_for_dpi;
  F.params = &dp;

  double result, abserr;
  int params_per_subtype = hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2;
  gsl_deriv_central(&F, gsl_vector_get(x, hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + l-1), 1e-5, &result, &abserr);

  return result;
}

double calc_dx_pi_llik_analytic(READS& res, gsl_vector* x, int i, int s, hyperparams& hpa)
{
  int K = res.size();

  trees trs;
  trees_cons(trs, hpa);

  params pa;
  init_params(pa, hpa);
  calc_params(x, pa, hpa, trs[0]);

  params grad_by_param;
  init_params(grad_by_param, hpa);

  double llik = d_llik(res, pa, grad_by_param, hpa);
  double sum_alpha = sum_vector(hpa.alpha, 1, hpa.TOTAL_CN);
  return (grad_by_param.pa[i]->pi[s] - pa.pa[i]->pi[s] * Log(sum_alpha - hpa.TOTAL_CN + K)).eval();

  // double grad = 0;
  // for (int l=1; l<=hpa.TOTAL_CN; ++l)
  //   {
  //     grad += dx_vec(pa[i]->pi, l, s) * grad_by_param[i]->pi[l] / pa[i]->pi[l];
  //   }
  
  // return grad;
}

double calc_dx_kappa_llik_numeric(READS& res, gsl_vector* x, int i, int l, int r, hyperparams& hpa)
{
  gsl_function F;
  dkappa dk (res, x, i, l, r, hpa);
  dk.res = res;
  dk.x = x;
  dk.i = i;
  dk.l = l;
  dk.r = r;
  dk.hpa = hpa;

  F.function = &calc_llik_for_dkappa;
  F.params = &dk;

  double result, abserr;
  int params_per_subtype = hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2;
  gsl_deriv_central(&F, gsl_vector_get(x, hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + hpa.TOTAL_CN + (l-1) * l / 2 + r-1), 1e-5, &result, &abserr);
  
  return result;
}

double calc_dx_kappa_llik_analytic(READS& res, gsl_vector* x, int i, int l, int t, hyperparams& hpa)
{
  int K = res.size();

  trees trs;
  trees_cons(trs, hpa);

  params pa;
  init_params(pa, hpa);
  calc_params(x, pa, hpa, trs[0]);

  params grad_by_param;
  init_params(grad_by_param, hpa);

  double llik = d_llik(res, pa, grad_by_param, hpa);

  double sum_beta = sum_vector(hpa.beta[l], 1, l);
  return (grad_by_param.pa[i]->kappa[l][t] - pa.pa[i]->kappa[l][t] * ( Log(sum_beta - l - hpa.alpha[l] + 1.0) + grad_by_param.pa[i]->pi[l] )).eval();

  // double grad = 0;
  // for (int r=1; r<=l; ++r)
  //   {
  //     grad += dx_vec(pa[i]->kappa[l], r, t) * grad_by_param[i]->kappa[l][r] / pa[i]->kappa[l][r];
  //   }
  
  // return grad;
}

void gsl_set_random(gsl_vector* x, hyperparams& hpa, gsl_rng* rng)
{
  for (int i=0; i<1024; ++i) // for appropriet random number generation
    gsl_rng_uniform(rng);
    // gsl_ran_beta(r, hpa.be_hpa.first, hpa.be_hpa.second);
  
  // for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
  //   pa[i]->u = gsl_ran_beta(r, hpa.be_hpa.first, hpa.be_hpa.second);
  // pa[0]->u *= 0.1;

  for (int i = 0; i<hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1 + hpa.MAX_SUBTYPE * hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2; ++i)
    {
      double a = gsl_rng_uniform(rng);
      gsl_vector_set(x, i, a);
    }
}

int main(int argc, char** argv)
{
  cout << scientific;
  cerr << scientific;
  
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 10)
    {
      cerr << "usage: ./lda_tree_compare_gradient max_subtype total_cn n step (reads) (rho diff outfile) (u diff outfile) (pi diff outfile) (kappa diff outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  READS res;
  hyperparams hpa;
  int n, step;

  hpa.MAX_SUBTYPE = atoi(argv[1]);
  hpa.TOTAL_CN = atoi(argv[2]);
  init_hyperparams(hpa);
  n = atoi(argv[3]);
  step = atoi(argv[4]);
  
  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  // T = gsl_rng_default; // default is mt19937
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, 1);

  ifstream f (argv[5]);
  ofstream q (argv[6]);
  ofstream ff (argv[7]);
  ofstream g (argv[8]);
  ofstream h (argv[9]);

  q << scientific;
  ff << scientific;
  g << scientific;
  h << scientific;

  gsl_vector* x = gsl_vector_alloc(hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1 + hpa.MAX_SUBTYPE * hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2);

  for (int i=0; i<n; ++i)
    {
      READ *re = new READ;
      f >> re->first >> re->second;
      res.push_back(re);
    }

  for (int i=-step; i<=step; ++i)
    {
      gsl_set_random(x, hpa, r);
      gsl_vector_set(x, 0, 1.0 * ((double)i) / step);
      double num = calc_dx_rho_llik_numeric(res, x, 0, hpa);
      double analytic = calc_dx_rho_llik_analytic(res, x, 0, hpa);
      
      if (fabs(num) > 0)
        q << gsl_vector_get(x, 0) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
      else
        q << gsl_vector_get(x, 0) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
    }

  for (int i=-step; i<=step; ++i)
    {
      gsl_set_random(x, hpa, r);
      gsl_vector_set(x, hpa.MAX_TREE + 1, 1.0 * ((double)i) / step);
      double num = calc_dx_u_llik_numeric(res, x, 1, hpa);
      double analytic = calc_dx_u_llik_analytic(res, x, 1, hpa);
      
      if (fabs(num) > 0)
        ff << gsl_vector_get(x, hpa.MAX_TREE + 1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
      else
        ff << gsl_vector_get(x, hpa.MAX_TREE + 1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
    }

  for (int i=-step; i<=step; ++i)
    {
      gsl_set_random(x, hpa, r);
      gsl_vector_set(x, hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1, 1.0 * ((double)i) / step);
      double num = calc_dx_pi_llik_numeric(res, x, 1, 1, hpa);
      double analytic = calc_dx_pi_llik_analytic(res, x, 1, 1, hpa);
      
      if (fabs(num) > 0)
        g << gsl_vector_get(x, hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
      else
        g << gsl_vector_get(x, hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
    }

  for (int i=-step; i<=step; ++i)
    {
      gsl_set_random(x, hpa, r);
      gsl_vector_set(x, hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1 + hpa.MAX_SUBTYPE * hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2 - 1, 1.0 * ((double)i) / step);
      double num = calc_dx_kappa_llik_numeric(res, x, hpa.MAX_SUBTYPE, hpa.TOTAL_CN, hpa.TOTAL_CN, hpa);
      double analytic = calc_dx_kappa_llik_analytic(res, x, hpa.MAX_SUBTYPE, hpa.TOTAL_CN, hpa.TOTAL_CN, hpa);
      
      if (fabs(num) > 0)
        h << gsl_vector_get(x, hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1 + hpa.MAX_SUBTYPE * hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2 - 1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << calc_rel_err(num, analytic) << endl;
      else
        h << gsl_vector_get(x, hpa.MAX_TREE + hpa.MAX_SUBTYPE + 1 + hpa.MAX_SUBTYPE * hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2 - 1) << "\t" << num << "\t" << analytic << "\t" << fabs(num - analytic) << "\t" << -1 << endl;
  }

  for (int i=0; i<n; ++i)
    delete res[i];

  f.close();
  g.close();
  h.close();
  ff.close();
  q.close();

  return 0;
}