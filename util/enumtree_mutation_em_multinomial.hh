#ifndef ENUMTREE_H_
#define ENUMTREE_H_

#include <fstream>
#include "loglib_header.hh"
#include <iostream>
#include "../src/alpha_n_subtype_expand_pattern_multinomial/setting.hh"

typedef std::vector<int> Vint;
typedef std::vector<unsigned int> Vuint;
typedef std::vector<Vint> VVint;
typedef std::vector<bool> Vbool;
typedef std::vector <Vbool> VVbool;
typedef std::vector <VVbool> VVVbool;
typedef std::vector<Vint*> Vints;
typedef std::vector<double> Vdouble;
typedef std::vector<Vdouble> VVdouble;

class subtype
{
public:
  int index;
  int total_cn;
  int variant_cn;
  int brother_index;
  int inherited;
  Log t;
  Log nu;
  Log n;
  Log xi;
  VLog omega;
  VLog zeta;
  Log x;
  Log growth;
  subtype* parent;
  subtype* above;
  std::vector<subtype*> children;

  subtype () {}
  subtype (int _index, int _total_cn, int _variant_cn, int _brother_index, int _inherited, Log _t, Log _nu, Log _n, Log _xi, VLog _omega, VLog _zeta, Log _x, Log _growth, subtype* _parent, subtype* _above, std::vector<subtype*> _children) : index(_index), total_cn(_total_cn), variant_cn(_variant_cn), brother_index(_brother_index), inherited(_inherited), t(_t), nu(_nu), n(_n), xi(_xi), zeta(_zeta), omega(_omega), x(_x), growth(_growth), parent(_parent), above(_above), children(_children) {}
};

typedef std::vector<subtype> subtypes;

void copy(subtypes& x, const subtypes& y)
{
  for (int i=0; i<(int)x.size(); ++i)
    {
      x[i].index = y[i].index;
      x[i].total_cn = y[i].total_cn;
      x[i].variant_cn = y[i].variant_cn;
      x[i].brother_index = y[i].brother_index;
      x[i].inherited = y[i].inherited;
      x[i].t = y[i].t;
      x[i].nu = y[i].nu;
      x[i].n = y[i].n;
      x[i].xi = y[i].xi;
      x[i].omega = y[i].omega;
      x[i].zeta = y[i].zeta;
      x[i].x = y[i].x; // bug fix
      x[i].growth = y[i].growth;
      
      if (y[i].parent == NULL)
        x[i].parent = NULL;
      else
        x[i].parent = &x[y[i].parent->index]; // copy index, not the pointer itself !
      
      if (y[i].above == NULL)
        x[i].above = NULL;
      else
        x[i].above = &x[y[i].above->index]; // copy index, not the pointer itself !

      x[i].children.assign(y[i].children.size(), NULL);
      for (int j=0; j<(int)y[i].children.size(); ++j)
        x[i].children[j] = &x[y[i].children[j]->index]; // copy index, not the pointer itself !
    }
}

typedef std::vector<subtypes> trees;

class hyperparams
{
public:
  std::pair<double,double> be_hpa_u;
  std::pair<double,double> be_hpa_beta;
  Vdouble alpha;
  VVdouble beta;
  Vdouble gamma;
  int MAX_SUBTYPE;
  int TOTAL_CN;
  int MAX_TREE;

  hyperparams () {}
  hyperparams (int, int, int);
};

class param
{
public:
  Log u;
  Log r;
  Log n;
  Log xi;
  VLog pi;
  VVLog kappa;

  param (Log _u, Log _r, Log _n, Log _xi, VLog _pi, VVLog _kappa) : u(_u), r(_r), n(_n), xi(_xi), pi(_pi), kappa(_kappa) {}
};

class params
{
public:
  std::vector<param*> pa;

  params (hyperparams&);
  ~params ();
};

params::params(hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    pa.push_back(new param (Log(0), Log((1e-6) / 6.0), Log(0), Log(0), VLog (hpa.TOTAL_CN + 1, Log(0)), VVLog (hpa.TOTAL_CN + 1, VLog (hpa.TOTAL_CN + 1, Log(0))))); // set mutation rate to be (1e-6) / 6.0 for all subtype
}

params::~params()
{
  for (int i=0; i<(int)pa.size(); ++i)
    {
      delete pa[i];
    }
}

hyperparams::hyperparams(int _MAX_SUBTYPE, int _TOTAL_CN, int _MAX_TREE) : MAX_SUBTYPE(_MAX_SUBTYPE), TOTAL_CN(_TOTAL_CN), MAX_TREE(_MAX_TREE)
{
  be_hpa_u.first = 1.0;
  be_hpa_u.second = 1.0;
  be_hpa_beta.first = 1.0;
  be_hpa_beta.second = 1.0;

  alpha.assign(_TOTAL_CN + 1, 0.1);
  beta.assign(_TOTAL_CN + 1, Vdouble (_TOTAL_CN + 1, 0.1));
  gamma.assign(_MAX_SUBTYPE + 1, 1.0); // changed to be uniform
}

class state
{
public:
  int q;
  int h;
  int xq;
  Vint m;
  Vint M;
  Log resp;

  state (int _q, int _h, int _xq, Vint _m, Vint _M, Log _resp) : q(_q), h(_h), xq(_xq), m(_m), M(_M), resp(_resp) {}
};

typedef std::vector<state*> states;
typedef std::vector<states*> statess;

// class state2
// {
// public:
//   int q;
//   int h;
//   int xq;
//   Vint i;
//   Log resp;

//   state2 (int _q, int _h, int _xq, Vint _i, Log _resp) : q(_q), h(_h), xq(_xq), i(_i), resp(_resp) {}
// };

void write_Vint(std::ofstream& f, Vint& v)
{
  for (int i=0; i<(int)v.size(); ++i)
    f << v[i] << " ";
  f << std::endl;
}

void write_Vints(std::ofstream& f, Vints& acc)
{
  for (int i=0; i<(int)acc.size(); ++i)
    write_Vint(f, *acc[i]);
}

void rooted_ordered_tree(Vint* dfs, int n, Vints& acc)
{
  if (n <= (int)dfs->size())
    {
      // write_Vint((std::ofstream&)std::cerr, *dfs);
      acc.push_back(dfs);
      return;
    }

  int k = dfs->size();
  for (int p = k; 0 < p; p = (*dfs)[p-1])
    {
      Vint* v = new Vint;
      v->assign(k + 1, 0);
      copy(dfs->begin(), dfs->end(), v->begin());
      (*v)[k] = p;
      rooted_ordered_tree(v, n, acc);
    }
}

void child_matrix(VVVbool& a, Vints& acc)
{
  for (int t=0; t<(int)acc.size(); ++t)
    {
      for (int i=1; i<(int)acc[0]->size(); ++i)
        {
          a[t][(*acc[t])[i]-1][i] = true;
        }
    }
}

void write_bool_matrix(VVVbool& a)
{
  for (int t=0; t<(int)a.size(); ++t)
    {
      for (int i=0; i<(int)a[t].size(); ++i)
        {
          for (int j=0; j<(int)a[t][i].size(); ++j)
            {
              if (a[t][i][j])
                std::cerr << 1 << " ";
              else
                std::cerr << 0 << " ";
            }
          std::cerr << std::endl;
        }
      std::cerr << std::endl;
    }
  std::cerr << "-----------------------------------" << std::endl;
}

void trees_cons(trees& trs, Vints& acc)
{
  for (int t=0; t<(int)acc.size(); ++t)
    {
      for (int i=1; i<=(int)acc[t]->size(); ++i)
        {
          trs[t][i].omega.push_back(Log(0));
        }
    }
  
  for (int t=0; t<(int)acc.size(); ++t)
    {
      trs[t][0].index = 0;
      trs[t][0].total_cn = 2;
      trs[t][0].variant_cn = 0;
      trs[t][0].parent = NULL;
      trs[t][0].above = NULL;

      for (int i=1; i<=(int)acc[t]->size(); ++i)
        {
          trs[t][i].index = i;
          trs[t][i].parent = &trs[t][(*acc[t])[i-1]];
          
          if (trs[t][i].parent->children.size() == 0)
            trs[t][i].above = trs[t][i].parent;
          else
            trs[t][i].above = trs[t][i].parent->children.back();

          trs[t][i].brother_index = trs[t][i].parent->children.size();

          trs[t][(*acc[t])[i-1]].children.push_back(&trs[t][i]); // elder child is pushed back first
          trs[t][(*acc[t])[i-1]].omega.push_back(Log(0));
        }
    }
}

void trees_cons(trees& trs, int max_subtype)
{
  // assumes max_subtype >= 1
  
  Vint* dfs = new Vint;
  dfs->assign(1, 0);

  Vints acc;
  rooted_ordered_tree(dfs, max_subtype, acc);
  // write_Vints((std::ofstream&)std::cerr, acc);

  // VVVbool a (acc.size(), VVbool(max_subtype, Vbool(max_subtype, false)));
  
  // child_matrix(a, acc);
  // write_bool_matrix(a);
  
  trs.assign(acc.size(), subtypes (max_subtype + 1, subtype (0, 0, 0, 0, 0, Log(0), Log(0), Log(0), Log(0), VLog(0, 0), VLog(FRACTIONS + 1, 0), Log(0), Log(0), NULL, NULL, std::vector< subtype* > (0, NULL) )));
  trees_cons(trs, acc);
}

void traversal(subtype* p)
{
  std::cerr << p->index << std::endl;
  if (p->index > 0)
    {
      std::cerr << p->index << "'s parent " << p->parent->index << std::endl;
      std::cerr << p->index << "'s above " << p->above->index << std::endl;
    }

  for (int j = 0; j < (int)p->children.size(); ++j)
    {
      std::cerr << p->index << "'s " << j+1 << " th child ";
      traversal(p->children[j]);
    }
}

void mark_inherited(subtype* p)
{
  p->inherited = 1;
  p->x = Log(1);
  for (int j = 0; j < (int)p->children.size(); ++j)
    {
      // p->children[j]->x = Log(1);
      mark_inherited(p->children[j]);
    }
}

void clear_inherited(subtypes& tr, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      tr[i].inherited = 0;
      tr[i].x = Log(0);
    }
}

void calc_t(params& pa, hyperparams& hpa, subtypes& sts)
{
  sts[0].t = Log(1);
  
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      sts[i].t = sts[sts[i].above->index].t * pa.pa[i]->u;
    }
}

bool above_time(subtypes& sts, int j, int y)
{
  if (j == 0 || y == 0)
    return false;
  
  for (subtype* p = &sts[j]; p->index != 0; p = p->above)
    {
      if (p->index == y)
        return true;
    }
  return false;
}

bool is_ancestor(subtypes& sts, int i, int j)
{
  if (j == 0)
    return true;

  for (subtype* p = &sts[i]; p->index != 0; p = p->parent)
    {
      if (p->index == j)
        return true;
    }
  return false;
}

int ancestor_brother_index(subtypes& sts, int i, int j)
{
  if (j == 0)
    return (i == j);

  subtype* q = &sts[i];
  for (subtype* p = &sts[i]; p->index != 0; q = p, p = p->parent)
    {
      if (p->index == j)
        return q->brother_index;
    }
  return -1;
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

int calc_eldest_child_number(subtype& st, int pat)
{
  for (int i=0; i<st.children.size(); ++i)
    {
      if ((pat >> (st.children.size()-1-i)) % 2 == 1)
        {
          return i + 1;
        }
    }
  return 0;
}

void clear_x(subtypes& tr, hyperparams&hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      tr[i].x = Log(0);
    }
}

Log d_t_u(params& pa, subtypes& sts, int j, int y)
{
  if (above_time(sts, j, y))
    return sts[j].t / pa.pa[y]->u;
  else
    return Log(0);
}

void calc_n(params& pa, hyperparams& hpa, subtypes& sts)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      sts[i].n = pa.pa[i]->n; //implemented here
    }
}

void calc_xi(hyperparams& hpa, subtypes& sts)
{
  Log sum = Log(0);
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      Log lNni = (Log(CELL_MAX) * sts[i].n).take_log_Log();
      sts[i].xi = sts[i].t * sts[i].n / lNni;
      sum += sts[i].xi;
    }
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      sts[i].xi /= sum;
    }
}

void calc_omega(hyperparams& hpa, subtypes& tr)
{
  VLog t_ratio (hpa.MAX_SUBTYPE+1, 0);

  for (int q=1; q<=hpa.MAX_SUBTYPE; ++q)
    {
      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        t_ratio[i] = tr[i].t / tr[q].t;

      Log Nnq = Log(CELL_MAX) * tr[q].n;
      Log s = (Nnq - Log(1)) / Nnq.take_log_Log();

      Log sum_th (0);
      for (int h=1; h<pow(2,tr[q].children.size()); ++h)
        {
          int ic = calc_eldest_child_index(tr[q], h);
          sum_th += t_ratio[ic];
        }
      
      tr[q].omega[0] = (s - Log(1)) / (s + Log(pow(2,tr[q].children.size()) - 2) - sum_th);
      
      for (int j=1; j<=(int)tr[q].children.size(); ++j) // debuged from < to <=
        {
          int ic = tr[q].children[j-1]->index;
          tr[q].omega[j] = (Log(1) - t_ratio[ic]) / (s + Log(pow(2,tr[q].children.size()) - 2) - sum_th);
        }
    }
}

void calc_subtypes(params& pa, hyperparams& hpa, subtypes& tr)
{
  calc_t(pa, hpa, tr);
  calc_n(pa, hpa, tr);
  calc_xi(hpa, tr);
  calc_omega(hpa, tr);
}

Log d_h_q_t(hyperparams& hpa, subtypes& tr, int i, int q, int h)
{
  Log Nnq = Log(CELL_MAX) * tr[q].n;
  Log s = (Nnq - Log(1)) / Nnq.take_log_Log();

  Log sum_th (0);
  Vint count (hpa.MAX_SUBTYPE + 1, 0);
  for (int _h=1; _h<pow(2,tr[q].children.size()); ++_h)
    {
      int ic = calc_eldest_child_index(tr[q], _h);
      sum_th += tr[ic].t;
      count[ic]++;
    }

  if (i < q)
    {
      return Log(0);
    }

  if (i == q) // and case by h=0 or h > 1
    {
      if (h == 0)
        {
          return -(s - Log(1)) * sum_th / ((s + Log(pow(2,tr[q].children.size()) - 2)) * tr[q].t - sum_th).take_pow(2);
        }
      else
        {
          int ih = calc_eldest_child_index(tr[q], h);
          return ((s + Log(pow(2,tr[q].children.size()) - 2)) * tr[ih].t - sum_th) / ((s + Log(pow(2,tr[q].children.size()) - 2)) * tr[q].t - sum_th).take_pow(2);
        }
    }
  
  if (i > q) // and case by h=0 or h > 1
    {
      VLog t_tild (hpa.MAX_SUBTYPE + 1, Log(0));
      for (int j=1; j<=hpa.MAX_SUBTYPE; ++j)
        {
          t_tild[j] = tr[j].t / tr[q].t;
        }
      
      Log sum_th_tild (0);
      for (int _h=1; _h<pow(2,tr[q].children.size()); ++_h)
        {
          int ic = calc_eldest_child_index(tr[q], _h);
          sum_th_tild += t_tild[ic];
        }
      
      if (h == 0)
        {
          return (s - Log(1)) / (s + Log(pow(2,tr[q].children.size()) - 2) - sum_th_tild).take_pow(2) * count[i] / tr[q].t;
        }
      else
        {
          int ic = calc_eldest_child_index(tr[q], h);
          Log extra (0);
          if (i == ic) extra = Log(1);
          return Log(1) / tr[q].t / (s + Log(pow(2,tr[q].children.size()) - 2) - sum_th_tild) * (-extra - (Log(1) - t_tild[ic]) * sum_th_tild / (s + Log(pow(2,tr[q].children.size()) - 2) - sum_th_tild));
        }
    }
}

Log d_lnh_q_t(hyperparams& hpa, subtypes& tr, int i, int q, int h)
{
  Log Nnq = Log(CELL_MAX) * tr[q].n;
  Log s = (Nnq - Log(1)) / Nnq.take_log_Log();

  Log sum_th (0);
  Vint count (hpa.MAX_SUBTYPE + 1, 0);
  for (int _h=1; _h<pow(2,tr[q].children.size()); ++_h)
    {
      int ic = calc_eldest_child_index(tr[q], _h);
      sum_th += tr[ic].t;
      count[ic]++;
    }

  if (i < q)
    {
      return Log(0);
    }

  if (i == q) // and case by h=0 or h > 1
    {
      if (h == 0)
        {
          return -sum_th / tr[q].t / ((s + Log(pow(2,tr[q].children.size()) - 2)) * tr[q].t - sum_th);
        }
      else
        {
          int ih = calc_eldest_child_index(tr[q], h);
          return (tr[q].t - tr[ih].t).inverse() - (s + Log(pow(2,tr[q].children.size()) - 2)) / ((s + Log(pow(2,tr[q].children.size()) - 2)) * tr[q].t - sum_th);
        }
    }
  
  if (i > q) // and case by h=0 or h > 1
    {
      if (h == 0)
        {
          return Log(count[i]) / ((s + Log(pow(2,tr[q].children.size()) - 2)) * tr[q].t - sum_th);
        }
      else
        {
          int ih = calc_eldest_child_index(tr[q], h);
          Log extra (0);
          if (i == ih) extra = Log(1);
          // std::cerr << "count[" << i << "] = " << count[i] << std::endl;
          return -extra / (tr[q].t - tr[ih].t) + Log(count[i]) / ((s + Log(pow(2,tr[q].children.size()) - 2)) * tr[q].t - sum_th);
        }
    }
}

Log d_s_n(hyperparams& hpa, subtypes& tr, int i)
{
  Log Nni = Log(CELL_MAX) * tr[i].n;
  
  return Log(CELL_MAX) / Nni.take_log_Log() * (Log(1) - (Log(1) - Nni.inverse()) / Nni.take_log_Log());
}

Log d_h_q_n(hyperparams& hpa, subtypes& tr, int q, int h)
{
  Log Nnq = Log(CELL_MAX) * tr[q].n;
  Log s = (Nnq - Log(1)) / Nnq.take_log_Log();
  
  VLog t_ratio (hpa.MAX_SUBTYPE+1, 0);
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    t_ratio[i] = tr[i].t / tr[q].t;

  Log sum_th (0);
  for (int _h=1; _h<pow(2,tr[q].children.size()); ++_h)
    {
      int ic = calc_eldest_child_index(tr[q], _h);
      sum_th += t_ratio[ic];
    }

  if (h == 0)
    {
      return (Log(pow(2,tr[q].children.size()) - 1) - sum_th) / (s + Log(pow(2,tr[q].children.size()) - 2) - sum_th).take_pow(2);
    }

  else
    {
      int ih = calc_eldest_child_index(tr[q], h);
      return -(Log(1) - t_ratio[ih]) / (s + Log(pow(2,tr[q].children.size()) - 2) - sum_th).take_pow(2);
    }
}

void calc_growth(params& pa, hyperparams& hpa, subtypes& st)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      st[i].growth = (Log(CELL_MAX) * st[i].n).take_log_Log() / st[i].t;
    }
}

#endif  // ENUMTREE_H_
