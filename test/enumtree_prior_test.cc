#include <algorithm>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include "../util/enumtree_prior.hh"
using namespace std;

#define calc_sigmoid(x) ((tanh((x)/2.0) + 1.0) / 2.0)

void calc_u(const gsl_vector* x, params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->u = Log(calc_sigmoid(gsl_vector_get(x, i)));
    }
}

bool strictly_greater_than(double i, double j)
{
  return (i > j);
}

int main(int argc, char** argv)
{
  if (argc != 3)
    {
      cerr << "usage: ./tree_enumeration max_subtype seed" << endl;
      exit(EXIT_FAILURE);
    }

  int seed = atoi(argv[2]);
  const gsl_rng_type * T;

  gsl_rng * r;

  gsl_rng_env_setup();

  // T = gsl_rng_default; // default is mt19937
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed);

  for (int i=0; i<1024; ++i) // for appropriet random number generation
    gsl_rng_uniform(r);

  int MAX_SUBTYPE, TOTAL_CN, MAX_TREE;
  MAX_SUBTYPE = atoi(argv[1]);
  TOTAL_CN = 3;

  if (MAX_SUBTYPE > 0)
    {
      trees trs;
      trees_cons(trs, MAX_SUBTYPE);

      trees trs2 (2, subtypes (MAX_SUBTYPE + 1, subtype()));
      copy(trs2[0], trs[0]);
      copy(trs2[1], trs[2]);
      MAX_TREE = trs2.size();

      hyperparams hpa (MAX_SUBTYPE, TOTAL_CN, MAX_TREE);
      
      // subtype copy assignment operator is ok
      // trs[0][4] = trs[0][3];
      // cerr << trs[0][4].parent->index << endl;
      // trs[0][3].parent = trs[0][1].parent;
      // cerr << trs[0][4].parent->index << endl;

      // subtypes copy assignment is clone !
      // traversal(&trs[4][0]);
      // cerr << endl;
      
      // trs[4] = trs[0];
      // traversal(&trs[4][0]);
      // cerr << endl;
      
      // for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
      //   {
      //     trs[0][i] = trs[1][i];
      //   }
      // traversal(&trs[4][0]);
      // cerr << endl;

      // subtypes copy assignment is clone !
      // traversal(&trs[4][0]);
      // cerr << endl;
      
      // copy(trs[4], trs[0]);

      // traversal(&trs[4][0]);
      // cerr << endl;

      // traversal(&trs[0][0]);
      // cerr << endl;

      // copy(trs[0], trs[1]);

      // traversal(&trs[4][0]);
      // cerr << endl;

      // traversal(&trs[0][0]);
      // cerr << endl;
      
      for (int t=0; t<(int)trs2.size(); ++t)
        {
          traversal(&trs2[t][0]);
          cerr << "---------------------------------------" << endl;
        }

      const gsl_vector* x = gsl_vector_calloc(hpa.MAX_SUBTYPE + 1);
      
      params pa (hpa);
      calc_u(x, pa, hpa);
      
      for (int t=0; t<(int)trs2.size(); ++t)
        {
          calc_t(pa, hpa, trs2[t]);

          for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
            cerr << i << ":\t" << trs2[t][i].t.eval() << endl;
          cerr << "---------------------------------------" << endl;

          for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
            {
              // if (i==0)
              //   cerr << "1\t";
              // else
              //   cerr << "0\t";
              for (int j=0; j<=hpa.MAX_SUBTYPE; ++j)
                {
                  // cerr << above_time(trs2[t], i, j) << "\t";
                  cerr << d_t_u(pa, trs2[t], i, j).eval() << "\t";
                }
              cerr << endl;
            }
          cerr << "=======================================" << endl;
        }

      // Vdouble n (hpa.MAX_SUBTYPE + 1, 0);
      // Vdouble n_prior (hpa.MAX_SUBTYPE + 1, 0.1);
      // gsl_ran_dirichlet(r, hpa.MAX_SUBTYPE + 1, &n_prior[0], &n[0]);
      // sort(n.begin()+1, n.end(), strictly_greater_than);
      // double tmp = n[4];
      // n[4] = n[3];
      // n[3] = tmp;
      
      // cerr << "t:" << endl;
      // for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
      //   cerr << i << ":\t" << n[i] << endl;

      // cerr << "=======================================" << endl;
      // cerr << "---------------------------------------" << endl;
      
      // for (int t=0; t<hpa.MAX_TREE; ++t)
      //   {
      //     for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
      //       {
      //         trs[t][i].t = Log(n[i]);
      //       }
          
      //     calc_u_from_t(pa, hpa, trs[t]);

      //     cerr << "calculated u:" << endl;
      //     for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
      //       cerr << i << ":\t" << pa.pa[i]->u.eval() << endl;
      //     cerr << "---------------------------------------" << endl;

      //     calc_t(pa, hpa, trs[t]);
      //     cerr << "recalculated t:" << endl;
      //     for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
      //       cerr << i << ":\t" << trs[t][i].t.eval() << endl;
      //     cerr << "=======================================" << endl;
      //   }
    }
  
  return 0;
}
