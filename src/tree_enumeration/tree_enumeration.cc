#include <fstream>
#include <algorithm>
#include "../loglib.hh"
using namespace std;

typedef vector<int> Vint;
typedef vector<Vint> VVint;
typedef vector<bool> Vbool;
typedef vector <Vbool> VVbool;
typedef vector <VVbool> VVVbool;
typedef vector<Vint*> Vints;
typedef vector<double> Vdouble;
typedef vector<Vdouble> VVdouble;

class subtype
{
public:
  int k;
  int i;
  int total_cn;
  int variant_cn;
  Log resp_ds;
  Log resp;
  vector<subtype*> children;
  
  // use default constructor
};

typedef vector<subtype> subtypes;
typedef vector<subtypes> trees;

class param
{
public:
  Log s;
  Log n;
  VLog pi;
  VVLog kappa;

  param (Log _s, Log _n, VLog _pi, VVLog _kappa) : s(_s), n(_n), pi(_pi), kappa(_kappa) {}
};

typedef vector<param*> params;

typedef struct _hyperparams
{
  pair<double,double> be_hpa;
  Vdouble alpha;
  VVdouble beta;
  int MAX_SUBTYPE;
  int TOTAL_CN;
}
  hyperparams;

void init_hyperparams(hyperparams& hpa)
{
  hpa.be_hpa.first = 0.1;
  hpa.be_hpa.second = 0.1;
  hpa.alpha.assign(hpa.TOTAL_CN + 1, 0.1);
  hpa.beta.assign(hpa.TOTAL_CN + 1, Vdouble (hpa.TOTAL_CN + 1, 0.1));
}

void write_Vint(ofstream& f, Vint& v)
{
  for (int i=0; i<v.size(); ++i)
    f << v[i] << " ";
  f << endl;
}

void write_Vints(ofstream& f, Vints& acc)
{
  for (int i=0; i<acc.size(); ++i)
    write_Vint(f, *acc[i]);
}


void rooted_ordered_tree(Vint* dfs, int n, Vints& acc)
{
  if (n <= dfs->size())
    {
      // write_Vint((ofstream&)cout, *dfs);
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
  for (int t=0; t<acc.size(); ++t)
    {
      for (int i=1; i<acc[0]->size(); ++i)
        {
          a[t][(*acc[t])[i]-1][i] = true;
        }
    }
}

void write_bool_matrix(VVVbool& a)
{
  for (int t=0; t<a.size(); ++t)
    {
      for (int i=0; i<a[t].size(); ++i)
        {
          for (int j=0; j<a[t][i].size(); ++j)
            {
              if (a[t][i][j])
                cerr << 1 << " ";
              else
                cerr << 0 << " ";
            }
          cerr << endl;
        }
      cerr << endl;
    }
  cerr << "-----------------------------------" << endl;
}

void trees_cons(trees& trs, VVVbool& a)
{
  for (int t=0; t<a.size(); ++t)
    {
      for (int i=0; i<=a[t].size(); ++i)
        {
          trs[t][i].i = i;
        }

      if (a[t].size() > 0)
        trs[t][0].children.push_back(&trs[t][1]);

      for (int i=0; i<a[t].size(); ++i)
        {
          for (int j=i+1; j<a[t][i].size(); ++j)
            {
              if (a[t][i][j])
                trs[t][i+1].children.push_back(&trs[t][j+1]);
            }
        }
    }
}

void trees_cons(trees& trs, hyperparams& hpa)
{
  Vint* dfs = new Vint;
  dfs->assign(1, 0);

  Vints acc;
  rooted_ordered_tree(dfs, hpa.MAX_SUBTYPE, acc);
  // write_Vints((ofstream&)cout, acc);

  VVVbool a (acc.size(), VVbool(hpa.MAX_SUBTYPE, Vbool(hpa.MAX_SUBTYPE, false)));
  
  child_matrix(a, acc);
  write_bool_matrix(a);
  
  trs.assign(acc.size(), subtypes (hpa.MAX_SUBTYPE + 1));
  trees_cons(trs, a);
}

void traversal(subtype* p)
{
  cerr << p->i << endl;

  for (int j = 0; j < p->children.size(); ++j)
    {
      cerr << p->i << "'s " << j+1 << " th child ";
      traversal(p->children[j]);
    }
}

int main(int argc, char** argv)
{
  if (argc != 2)
    {
      cerr << "usage: ./tree_enumeration n" << endl;
      exit(EXIT_FAILURE);
    }

  hyperparams hpa;
  hpa.MAX_SUBTYPE = atoi(argv[1]);
  hpa.TOTAL_CN = 3;
  init_hyperparams(hpa);

  trees trs;
  trees_cons(trs, hpa);
  
  for (int t=0; t<trs.size(); ++t)
    {
      traversal(&trs[t][0]);
      cerr << "---------------------------------------" << endl;
    }
  
  return 0;
}
