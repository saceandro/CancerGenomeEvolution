#ifndef LOGLIB_H_
#define LOGLIB_H_

#include <vector>

#define DBL_MAX_1 (0x7feffffffffffffeULL)

typedef union _double_union
{
  double d;
  unsigned long long int c;
}
  double_union;

class Log
{
  double val;
  signed char sign;
  // double_union dm_1;

public:
  Log ();
  Log (double, signed char);
  Log (double);
  
  double eval();
  bool iszero();
  signed char get_sign();
  double get_val();
  double take_log();
  Log take_log_Log();
  Log take_pow(double);
  Log take_exp();
  Log inverse();

  Log operator -();
  Log operator +(const Log&);
  Log operator +=(const Log&);
  Log operator -(const Log&);
  Log operator -=(const Log&);
  Log operator *(const Log&);
  Log operator *=(const Log&);
  Log operator /(const Log&);
  Log operator /=(const Log&);
  bool operator <(const Log&);
  bool operator >(const Log&);
};

typedef std::vector<Log> VLog;
typedef std::vector<VLog> VVLog;

#endif  // LOGLIB_H_
