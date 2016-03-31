#include <gsl/gsl_vector.h>
#include "enumtree.hh"
using namespace std;

#define calc_sigmoid(x) ((tanh((x)/2.0) + 1.0) / 2.0)

void calc_params(const gsl_vector* x, params& pa, hyperparams& hpa, subtypes& sts)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa[i]->u = Log(calc_sigmoid(gsl_vector_get(x, i)));
      // pa[i]->u = Log(gsl_vector_get(x, i));
    }
  
  calc_t(pa, hpa, sts);
  // calc_n(pa, hpa);
  
  int params_per_subtype = hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2;

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double m = gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1));
      
      for (int l=2; l<=hpa.TOTAL_CN; ++l)
        {
          double s = gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + l-1);
          if (m < s) m = s;
        }
  
      Log sum;
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          pa[i]->pi[l] = Log(gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + l-1) - m, 1);
          sum += pa[i]->pi[l];
        }

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        pa[i]->pi[l] /= sum;

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          m = gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + hpa.TOTAL_CN + (l-1) * l / 2);
          for (int r=2; r<=l; ++r)
            {
              double s = gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + hpa.TOTAL_CN + (l-1) * l / 2 + r-1);
              if (m < s) m = s;
            }
      
          sum = Log(0);
      
          for (int r=1; r<=l; ++r)
            {
              pa[i]->kappa[l][r] = Log(gsl_vector_get(x, hpa.MAX_SUBTYPE + 1 + params_per_subtype * (i-1) + hpa.TOTAL_CN + (l-1) * l / 2 + r-1) - m, 1);
              sum += pa[i]->kappa[l][r];
            }

          for (int r=1; r<=l; ++r)
            pa[i]->kappa[l][r] /= sum;
        }
    }
}

int main(int argc, char** argv)
{
  if (argc != 2)
    {
      cerr << "usage: ./tree_enumeration max_subtype" << endl;
      exit(EXIT_FAILURE);
    }

  hyperparams hpa;
  hpa.MAX_SUBTYPE = atoi(argv[1]);
  hpa.TOTAL_CN = 3;
  init_hyperparams(hpa);

  if (hpa.MAX_SUBTYPE > 0)
    {
      trees trs;
      trees_cons(trs, hpa);

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
      traversal(&trs[4][0]);
      cerr << endl;
      
      copy(trs[4], trs[0]);

      traversal(&trs[4][0]);
      cerr << endl;

      traversal(&trs[0][0]);
      cerr << endl;

      copy(trs[0], trs[1]);

      traversal(&trs[4][0]);
      cerr << endl;

      traversal(&trs[0][0]);
      cerr << endl;
      
      // for (int t=0; t<trs.size(); ++t)
      //   {
      //     traversal(&trs[t][0]);
      //     cerr << "---------------------------------------" << endl;
      //   }

      // const gsl_vector* x = gsl_vector_calloc(hpa.MAX_SUBTYPE + 1 + hpa.MAX_SUBTYPE * hpa.TOTAL_CN * (hpa.TOTAL_CN + 3) / 2);

      // params pa;
      // init_params(pa, hpa);

      // for (int t=0; t<trs.size(); ++t)
      //   {
      //     calc_params(x, pa, hpa, trs[t]);
      //     calc_t(pa, hpa, trs[t]);

      //     for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
      //       cerr << i << ":\t" << pa[i]->t.eval() << endl;
      //     cerr << "---------------------------------------" << endl;

      //     for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
      //       {
      //         // if (i==0)
      //         //   cerr << "1\t";
      //         // else
      //         //   cerr << "0\t";
      //         for (int j=0; j<=hpa.MAX_SUBTYPE; ++j)
      //           {
      //             // cerr << above_time(trs[t], i, j) << "\t";
      //             cerr << d_t_u(pa, hpa, trs[t], i, j).eval() << "\t";
      //           }
      //         cerr << endl;
      //       }
      //     cerr << "=======================================" << endl;
      //   }
    }
  
  return 0;
}
