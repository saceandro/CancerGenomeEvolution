#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
using namespace std;

#define MAX_SUBTYPE 4
#define TOTAL_CN 4
#define calc_sigmoid(x) (((x) < 0) ? (exp(x) / (exp(x) + 1.0)) : (1.0 / (exp(-(x)) + 1.0)))

typedef vector<double> Vdouble;
typedef vector<unsigned int> Vuint;

class state
{
public:
  state& nxt;
  unsigned int variant_cn;
};

typedef pair<unsigned int, unsigned int> read;

double softmax(Vdouble& x, unsigned int r)
{
  double m = *max_element(x.begin(), x.end());
  double sum = 0;
  for (unsigned int i=0; i<x.size(); ++i)
    sum += exp(x[i] - m);
  
  return exp(x[r] - m) / sum;
}

double dx_softmax(Vdouble& x, unsigned int r, unsigned int s)
{
  double m = *max_element(x.begin(), x.end());
  double sum = 0;
  for (int i=0; i<x.size(); ++i)
    sum += exp(x[i] - m);
  
  double ans = 0;
  