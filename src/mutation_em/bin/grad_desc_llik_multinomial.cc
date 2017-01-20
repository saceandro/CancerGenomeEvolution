#include <fstream>
#include <algorithm>
#include <numeric>
#include <fenv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf.h>
#include <iomanip>
#include "../../../util/enumtree_mutation_em_multinomial.hh"
#include "../setting.hh"
#include <xmmintrin.h>
using namespace std;

void calc_f(Log& llik, int num_of_split)
{
  char str[1024];
  double a;
  
  for (int sp=0; sp<num_of_split; ++sp)
    {
      int n = sprintf(str, "../du_dn_llik/%d", sp);

      ifstream f (str);

      f >> a;
      llik += Log(a);

      f.close();
    }
  // for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     cerr << du[i].eval() << "\t";
  //   }
  // cerr << endl;
  
  // for (int i=1; i<=hpa.MAX_SUBTYPE; ++i)
  //   {
  //     cerr << dn[i].eval() << "\t";
  //   }
  // cerr << endl;
  // cerr << llik.eval() << endl;
}

int main(int argc, char** argv)
{
  cout << scientific << setprecision(10);
  cerr << scientific;
  
  // feenableexcept(FE_INVALID);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  gsl_set_error_handler_off ();
  
  if (argc != 4)
    {
      cerr << "usage: ./grad_desc_llik num_of_split (u) (llik outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  int num_of_split = atoi(argv[1]);
  double u = atof(argv[2]);
  ofstream llik_f (argv[3], ofstream::out | ofstream::app);

  llik_f << scientific << setprecision(10);

  Log llik (0);
  calc_f(llik, num_of_split);

  llik_f << u << "\t" << llik.eval() << endl;

  llik_f.close();

  return 0;
}
