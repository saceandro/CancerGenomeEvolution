#include <fstream>
#include <algorithm>
#include "../loglib.hh"
using namespace std;

typedef vector<int> Vint;
typedef vector<Vint*> Vints;
typedef vector<double> Vdouble;
typedef vector<Vdouble> VVdouble;

class state 
{
public:
  int k;
  int i;
  int total_cn;
  int variant_cn;
  Log resp_ds;
  Log resp;

  // use default constructor
};

typedef std::vector<state*> states;

class state_node
{
public:
  vector<state*> childs;
}
  ;

typedef vector<state_node*> tree;
typedef vector<tree*> trees;

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

void init_trees(trees& trs, hyperparams& hpa)
{
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

int main(int argc, char** argv)
{
  if (argc != 2)
    {
      cerr << "usage: ./tree_enumeration n" << endl;
      exit(EXIT_FAILURE);
    }

  int n;
  n = atoi(argv[1]);
  
  Vint* dfs = new Vint;
  dfs->assign(1, 0);

  Vints acc;
  rooted_ordered_tree(dfs, n, acc);

  // write_Vints((ofstream&)cout, acc);

  hyperparams hpa;
  hpa.MAX_SUBTYPE = 2;
  hpa.TOTAL_CN = 3;
  init_hyperparams(hpa);
  
  trees trs;

  for (int t=0; t<acc.size(); ++t)
    {
      tree* tr = new tree(acc[0]->size(), new state_node);


      for (int i=0; i<acc[0]->size(); ++i)
        {
          state* st = new state;
          st->i = i+1;
          cout << (*acc[t])[i] << endl;
          cout << (*tr)[(*acc[t])[i]] << endl;
          
          (*tr)[(*acc[t])[i]]->childs.push_back(st);
        }
      trs.push_back(tr);
    }

  for (int t=0; t<acc.size(); ++t)
    {
      for (int i=0; i<acc[0]->size(); ++i)
        {
          for (states::iterator it = (*trs[t])[i]->childs.begin(); it != (*trs[t])[i]->childs.end(); ++it)
            cout << (*it)->i << endl;
          cout << endl;
        }
      cout << endl;
    }
  
  return 0;
}
