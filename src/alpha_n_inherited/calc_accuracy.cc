#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <xmmintrin.h>
using namespace std;

typedef vector<double> Vdouble;

double calc_rmsd(Vdouble& u, Vdouble& u_ans, int start, int end)
{
  double rmsd = 0;
  
  for (int i=start; i<=end; ++i)
    {
      double diff = u[i] - u_ans[i];
      rmsd += diff * diff;
    }

  return sqrt(rmsd / (end - start + 1));
}

void read_vector(std::ifstream& f, Vdouble& x, int start, int end)
{
  for (int i=start; i<=end; ++i)
    {
      f >> x[i];
    }
}

int main(int argc, char** argv)
{
  cout << scientific;
  cerr << scientific;
  
  // feenableexcept(FE_INVALID);
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  
  if (argc != 9)
    {
      cerr << "usage: ./accuracy subtype (u n ans infile) (t ans infile) (u n infile) (t infile) (u rmsd outfile) (t rmsd outfile) (n rmsd outfile)" << endl;
      exit(EXIT_FAILURE);
    }

  int MAX_SUBTYPE;

  MAX_SUBTYPE = atoi(argv[1]);

  ifstream u_n_ans_in(argv[2]);
  ifstream t_ans_in(argv[3]);
  ifstream u_n_in (argv[4]);
  ifstream t_in (argv[5]);
  ofstream u_out (argv[6]);
  ofstream t_out (argv[7]);
  ofstream n_out (argv[8]);
           
  u_out << scientific << setprecision(10);
  t_out << scientific << setprecision(10);
  n_out << scientific << setprecision(10);

  Vdouble u_ans (MAX_SUBTYPE + 1, 0);
  Vdouble t_ans (MAX_SUBTYPE + 1, 0);
  Vdouble n_ans (MAX_SUBTYPE + 1, 0);

  Vdouble u (MAX_SUBTYPE + 1, 0);
  Vdouble t (MAX_SUBTYPE + 1, 0);
  Vdouble n (MAX_SUBTYPE + 1, 0);

  read_vector(u_n_ans_in, u_ans, 1, MAX_SUBTYPE);
  read_vector(t_ans_in, t_ans, 0, MAX_SUBTYPE);
  read_vector(u_n_ans_in, n_ans, 0, MAX_SUBTYPE);

  read_vector(u_n_in, u, 1, MAX_SUBTYPE);
  read_vector(t_in, t, 0, MAX_SUBTYPE);
  read_vector(u_n_in, n, 0, MAX_SUBTYPE);

  u_out << calc_rmsd(u, u_ans, 1, MAX_SUBTYPE) << endl;
  t_out << calc_rmsd(t, t_ans, 1, MAX_SUBTYPE) << endl;
  n_out << calc_rmsd(n, n_ans, 1, MAX_SUBTYPE) << endl;

  u_n_ans_in.close();
  t_ans_in.close();
  u_n_in.close();
  t_in.close();
  u_out.close();
  t_out.close();
  n_out.close();

  return 0;
}
