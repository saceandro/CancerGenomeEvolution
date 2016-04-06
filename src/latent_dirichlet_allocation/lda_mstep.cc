#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
using namespace std;

typedef vector<int> Vint;
typedef vector<double> Vdouble;
typedef vector<Vdouble> VVdouble;
typedef pair<int, int> READ;
typedef vector<READ*> READS;

class state 
{
public:
  int k;
  Vint total_cn;
  Vint variant_cn;
  double resp_num;
  double resp;

  // use default constructor
  // state (int, std::vector<int>, std::vector<int>, double);
};

typedef std::vector<state*> states;

class param
{
public:
  Vdouble pi;
  VVdouble kappa;

  param (Vdouble _pi, VVdouble _kappa) : pi(_pi), kappa(_kappa) {}
};

typedef vector<param*> params;

typedef struct _hyperparams
{
  Vdouble alpha;
  VVdouble beta;
  int MAX_SUBTYPE;
  int TOTAL_CN;
}
  hyperparams;

void write_params(std::ofstream& f, params& pa, hyperparams& hpa)
{
  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      double sum = 0;
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        f << pa[i]->pi[l] << "\t";
      f << endl;

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (int r=1; r<=l; ++r)
            f << pa[i]->kappa[l][r] << "\t";
          f << endl;
        }
      f << endl;
    }
}

void init_state(state& st, hyperparams& hpa)
{
  st.total_cn.assign(hpa.MAX_SUBTYPE + 1, 0);
  st.variant_cn.assign(hpa.MAX_SUBTYPE + 1, 0);
  st.total_cn[0] = 2;
  st.variant_cn[0] = 0;
}

void init_params(params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    pa.push_back(new param (Vdouble (hpa.TOTAL_CN + 1, 0), VVdouble (hpa.TOTAL_CN + 1, Vdouble (hpa.TOTAL_CN + 1, 0))));
}

void init_hyperparams(hyperparams& hpa)
{
  // different from generating code, set to be uniform (same as maximum liklihood)
  hpa.alpha.assign(hpa.TOTAL_CN + 1, 1.0);
  hpa.beta.assign(hpa.TOTAL_CN + 1, Vdouble (hpa.TOTAL_CN + 1, 1.0));
}

void delete_params(params& pa, hyperparams& hpa)
{
  for (int i=0; i<=hpa.MAX_SUBTYPE; ++i)
    delete pa[i];
}

double mstep(params& pa_old, params& pa_new, hyperparams& hpa, int n, int num_of_split)
{
  double llik_sum = 0;
  double pi, kappa, llik;
  char str[1024];
  
  for (int k=0; k<num_of_split; ++k)
    {
      int n = sprintf(str, "../params/%02d", k);
      cout << str << endl;
      
      ifstream f (str);

      for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
        {
          for (int l=1; l<=hpa.TOTAL_CN; ++l)
            {
              f >> pi;
              pa_new[i]->pi[l] += pi;
            }

          for (int l=1; l<=hpa.TOTAL_CN; ++l)
            {
              for (int r=1; r<=l; ++r)
                {
                  f >> kappa;
                  pa_new[i]->kappa[l][r] += kappa;
                }
            }
        }
      f >> llik;
      llik_sum += llik;
    }

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          pa_new[i]->pi[l] += hpa.alpha[l] - 1.0;
          for (int r=1; r<=l; ++r)
            {
              pa_new[i]->kappa[l][r] += hpa.beta[l][r] - 1.0;
              pa_new[i]->kappa[l][r] /= l*(hpa.beta[l][r] - 1.0) + pa_new[i]->pi[l] - (hpa.alpha[l] - 1.0);
            }
          pa_new[i]->pi[l] /= hpa.TOTAL_CN * (hpa.alpha[l] - 1.0) + n * num_of_split;
        }
    }

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          llik_sum += (hpa.alpha[l] - 1.0) * log(pa_old[i]->pi[l]);
          for (int r=1; r<=l; ++r)
            {
              llik_sum += (hpa.beta[l][r] - 1.0) * log(pa_old[i]->kappa[l][r]);
            }
        }
    }

  // for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     for (int l=1; l<=hpa.TOTAL_CN; ++l)
  //       {
  //         for (int r=1; r<=l; ++r)
  //           {
  //             pa_new[i]->kappa[l][r] /= pa_new[i]->pi[l];
  //           }
  //         pa_new[i]->pi[l] /= n * num_of_split;
  //       }
  //   }

  return llik_sum;
}

int main(int argc, char** argv)
{
  cout << scientific;
  cerr << scientific;
  
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  if (argc != 8)
    {
      cerr << "usage: ./lda_mstep max_subtype total_cn n num_of_splits (pi kappa infile) (pi kappa outfile) (log_file)" << endl;
      exit(EXIT_FAILURE);
    }

  hyperparams hpa;
  int n, num_of_split;

  hpa.MAX_SUBTYPE = atoi(argv[1]);
  hpa.TOTAL_CN = atoi(argv[2]);
  init_hyperparams(hpa);
  n = atoi(argv[3]);
  num_of_split = atoi(argv[4]);
  
  ifstream f (argv[5]);
  ofstream g (argv[6]);
  ofstream h (argv[7], ofstream::out | ofstream::app);
  g << scientific;
  h << scientific;

  params pa;
  init_params(pa, hpa);

  for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
    {
      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          f >> pa[i]->pi[l];
        }

      for (int l=1; l<=hpa.TOTAL_CN; ++l)
        {
          for (int r=1; r<=l; ++r)
            {
              f >> pa[i]->kappa[l][r];
            }
        }
    }

  params pa_new;
  init_params(pa_new, hpa);

  double llik = mstep(pa, pa_new, hpa, n, num_of_split);
  write_params(g, pa_new, hpa);
  write_params(h, pa_new, hpa);
  h << llik << endl << endl;

  delete_params(pa, hpa);
  delete_params(pa_new, hpa);
  
  f.close();
  g.close();
  h.close();

  return 0;
}
