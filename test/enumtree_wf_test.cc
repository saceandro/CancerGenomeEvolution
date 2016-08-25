#include <gsl/gsl_vector.h>
#include "../util/enumtree_wf.hh"
using namespace std;

#define calc_sigmoid(x) ((tanh((x)/2.0) + 1.0) / 2.0)

void calc_params(const gsl_vector* x, params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      pa.pa[i]->u = Log(calc_sigmoid(gsl_vector_get(x, i)));
    }

  // pa.pa[0]->beta[0] = Log(calc_sigmoid(gsl_vector_get(x, hpa.MAX_SUBTYPE + 1)));

  // for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     for (int j=0; i<hpa.MAX_SUBTYPE-1; ++i)
  //       {
  //         pa.pa[i]->beta[j] = Log(calc_sigmoid(gsl_vector_get(x, hpa.MAX_SUBTYPE + 2 + (hpa.MAX_SUBTYPE - 1) * (i-1) + j)));
  //       }
  //   }
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int j=0; j<hpa.MAX_SUBTYPE-1; ++j)
        {
          pa.pa[i]->beta[j] = Log(0.5);
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


  int MAX_SUBTYPE = atoi(argv[1]);
  int TOTAL_CN = 3;

  if (MAX_SUBTYPE > 0)
    {
      trees trs;
      trees_cons(trs, MAX_SUBTYPE);
      int MAX_TREE = trs.size();
      
      hyperparams hpa(MAX_SUBTYPE, TOTAL_CN, MAX_TREE);
      
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
      
      for (int t=0; t<(int)trs.size(); ++t)
        {
          traversal(&trs[t][0]);
          cerr << "---------------------------------------" << endl;
        }

      const gsl_vector* x = gsl_vector_calloc(hpa.MAX_SUBTYPE + 2 + (hpa.MAX_SUBTYPE - 1) * (hpa.MAX_SUBTYPE - 1) + hpa.MAX_SUBTYPE-1);

      params pa (hpa);
      calc_params(x, pa, hpa);
      
      for (int t=0; t<(int)trs.size(); ++t)
        {
          calc_t(pa, hpa, trs[t]);

          for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
            cerr << i << ":\t" << trs[t][i].t.eval() << endl;
          cerr << "---------------------------------------" << endl;

          for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
            {
              // if (i==0)
              //   cerr << "1\t";
              // else
              //   cerr << "0\t";
              for (int j=0; j<=hpa.MAX_SUBTYPE; ++j)
                {
                  // cerr << above_time(trs[t], i, j) << "\t";
                  cerr << d_t_u(pa, trs[t], i, j).eval() << "\t";
                }
              cerr << endl;
            }
          cerr << "=======================================" << endl;
        }
      cerr << endl;
      
      for (int t=0; t<(int)trs.size(); ++t)
        {
          calc_n(pa, hpa, trs[t]);

          for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
            cerr << i << ":\t" << trs[t][i].n.eval() << endl;
          cerr << "---------------------------------------" << endl;

          for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
            {
              // if (i==0)
              //   cerr << "1\t";
              // else
              //   cerr << "0\t";
              for (int j=0; j<=hpa.MAX_SUBTYPE; ++j)
                {
                  // cerr << above_time(trs[t], i, j) << "\t";
                  for (int k=0; k<trs[t][j].children.size(); ++k)
                    {
                      cerr << d_n_beta(pa, trs[t], i, j, k).eval() << "\t";
                    }
                  cerr << endl;
                }
              cerr << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" << endl;
            }
          cerr << "=======================================" << endl;
        }
    }
  
  return 0;
}
