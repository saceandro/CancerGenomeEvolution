#ifndef ENUMTREE_H_
#define ENUMTREE_H_

#include <fstream>
#include "loglib.hh"

typedef std::vector<int> Vint;
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
  Log t;
  Log n;
  Log resp_du;
  subtype* parent;
  subtype* above;
  std::vector<subtype*> children;

  subtype () {}
  subtype (int _index, int _total_cn, int _variant_cn, Log _t, Log _n, Log _resp_du, subtype* _parent, subtype* _above, std::vector<subtype*> _children) : index(_index), total_cn(_total_cn), variant_cn(_variant_cn), t(_t), n(_n), resp_du(_resp_du), parent(_parent), above(_above), children(_children) {}
};

typedef std::vector<subtype> subtypes;

void copy(subtypes& x, const subtypes& y)
{
  for (int i=0; i<(int)x.size(); ++i)
    {
      x[i].total_cn = y[i].total_cn;
      x[i].variant_cn = y[i].variant_cn;
      x[i].t = y[i].t;
      x[i].n = y[i].n;
      x[i].resp_du = y[i].resp_du;
      
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
  std::pair<double,double> be_hpa;
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
  VLog pi;
  VVLog kappa;

  param (Log _u, VLog _pi, VVLog _kappa) : u(_u), pi(_pi), kappa(_kappa) {}
};

class params
{
public:
  std::vector<param*> pa;
  VLog rho;

  params (hyperparams&);
};

params::params(hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    pa.push_back(new param (Log(0), VLog (hpa.TOTAL_CN + 1, Log(0)), VVLog (hpa.TOTAL_CN + 1, VLog (hpa.TOTAL_CN + 1, Log(0)))));
  rho.assign(hpa.MAX_TREE, Log(0));
}

hyperparams::hyperparams(int _MAX_SUBTYPE, int _TOTAL_CN, int _MAX_TREE) : MAX_SUBTYPE(_MAX_SUBTYPE), TOTAL_CN(_TOTAL_CN), MAX_TREE(_MAX_TREE)
{
  be_hpa.first = 0.1;
  be_hpa.second = 0.1;
  alpha.assign(_TOTAL_CN + 1, 0.1);
  beta.assign(_TOTAL_CN + 1, Vdouble (_TOTAL_CN + 1, 0.1));
  gamma.assign(_MAX_TREE, 0.1);
}

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

          trs[t][(*acc[t])[i-1]].children.push_back(&trs[t][i]); // elder child is pushed back first
        }
    }
}

void trees_cons(trees& trs, hyperparams& hpa)
{
  // assumes hpa.MAX_SUBTYPE >= 1
  
  Vint* dfs = new Vint;
  dfs->assign(1, 0);

  Vints acc;
  rooted_ordered_tree(dfs, hpa.MAX_SUBTYPE, acc);
  // write_Vints((std::ofstream&)std::cerr, acc);

  // VVVbool a (acc.size(), VVbool(hpa.MAX_SUBTYPE, Vbool(hpa.MAX_SUBTYPE, false)));
  
  // child_matrix(a, acc);
  // write_bool_matrix(a);
  
  trs.assign(acc.size(), subtypes (hpa.MAX_SUBTYPE + 1, subtype (0, 0, 0, Log(0), Log(0), Log(0), NULL, NULL, std::vector< subtype* > (0, NULL) )));
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

void calc_t(params& pa, hyperparams& hpa, subtypes& sts)
{
  sts[0].t = pa.pa[0]->u;
  if (hpa.MAX_SUBTYPE >= 1)
    sts[1].t = pa.pa[1]->u;
  
  for (int i=2; i<=hpa.MAX_SUBTYPE; ++i)
    {
      sts[i].t = sts[sts[i].above->index].t * pa.pa[i]->u;
    }
}

bool above_time(subtypes& sts, int j, int y)
{
  if (y == 0)
    return (j == y);
  
  for (subtype* p = &sts[j]; p->index != 0; p = p->above)
    {
      if (p->index == y)
        return true;
    }
  return false;
}

Log d_t_u(params& pa, subtypes& sts, int j, int y)
{
  if (above_time(sts, j, y))
    return sts[j].t / pa.pa[y]->u;
  else
    return Log(0);
}

#endif  // ENUMTREE_H_
