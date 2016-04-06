#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <float.h>

class Log
{
  double val;
  int sign;
  
public:
  Log ();
  Log (double, int);
  Log (double);
  
  double eval();
  bool iszero();
  int get_sign();
  double take_log();
  Log take_pow(double);

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

Log::Log ()
{
  val = -DBL_MAX;
  sign = 0;
}

Log::Log (double val, int sign)
{
  this->val = val;
  this->sign = sign;
}

Log::Log (double a)
{
  if (a > 0)
    {
      val = log(a);
      sign = 1;
    }

  else if (a < 0)
    {
      val = log(-a);
      sign = -1;
    }
  
  else
    {
      val = -DBL_MAX;
      sign = 0;
    }
}

double Log::eval()
{
  return sign * exp(val);
}

bool Log::iszero()
{
  return (sign == 0);
}

int Log::get_sign()
{
  return sign;
}

double Log::take_log()
{
  if (this->sign != 1)
    {
      std::cout << "take_log(<=0)" << std::endl;
      std::cout << val << "\t" << sign << std::endl << std::endl;
      exit(EXIT_FAILURE);
    }
  
  return val;
}

Log Log::take_pow(double b)
{
  if (this->sign == -1)
    {
      std::cout << "take_pow(<0)" << std::endl;
      exit(EXIT_FAILURE);
    }
    
  Log temp;
  temp.val = this->val * b;
  temp.sign = this->sign;
  
  return temp;
}

Log Log::operator -()
{
  Log temp;
  temp.val = this->val;
  temp.sign = -this->sign;
  return temp;
}

Log Log::operator +(const Log& b)
{
  Log temp;

  if (this->val > b.val)
    {
      temp.val = this->val + log1p(this->sign * b.sign * exp(b.val - this->val));
      temp.sign = this->sign;
    }
  else
    {
      temp.val = b.val + log1p(this->sign * b.sign * exp(this->val - b.val));
      temp.sign = b.sign;
    }

  return temp;
}

Log Log::operator +=(const Log& x)
{
  if (this->val > x.val)
    {
      this->val += log1p(this->sign * x.sign * exp(x.val - this->val));
    }
  else
    {
      this->val = x.val + log1p(this->sign * x.sign * exp(this->val - x.val));
      this->sign = x.sign;
    }
  return *this;
}

Log Log::operator -(const Log& b)
{
  Log temp;

  if (this->val > b.val)
    {
      temp.val = this->val + log1p(this->sign * (-b.sign) * exp(b.val - this->val));
      temp.sign = this->sign;
    }
  else
    {
      temp.val = b.val + log1p(this->sign * (-b.sign) * exp(this->val - b.val));
      temp.sign = -b.sign;
    }

  return temp;
}

Log Log::operator -=(const Log& x)
{
  if (this->val > x.val)
    {
      this->val += log1p(this->sign * (-x.sign) * exp(x.val - this->val));
    }
  else
    {
      this->val = x.val + log1p(this->sign * (-x.sign) * exp(this->val - x.val));
      this->sign = -x.sign;
    }
  return *this;
}

Log Log::operator *(const Log& x)
{
  Log temp;
  temp.val = this->val + x.val;
  temp.sign = this->sign * x.sign;
  
  return temp;
}

Log Log::operator *=(const Log& x)
{
  this->val += x.val;
  this->sign *= x.sign;
  
  return *this;
}

Log Log::operator /(const Log& x)
{
  Log temp;

  temp.val = this->val - x.val;
  temp.sign = this->sign * x.sign;
  
  return temp;
}

Log Log::operator /=(const Log& x)
{
  this->val -= x.val;
  this->sign *= x.sign;
  
  return *this;
}

bool Log::operator <(const Log& x)
{
  // if ((this->sign == 1) && (x.sign == 1))
  //   return (this->val < x.val);
  // else if ((this->sign == 1) && (x.sign == -1))
  //   return false;
  // else if ((this->sign == -1) && (x.sign == 1))
  //   return true;
  // else
  //   return (this->val > x.val);
  if (this->sign == 1)
    {
      if (x.sign == 1)
        return (this->val < x.val);
      else
        return false;
    }
  else if (this->sign == 0)
    {
      if (x.sign == 1)
        return true;
      else
        return false;
    }
  else
    {
      if (x.sign == -1)
        return (this->val > x.val);
      else
        return true;
    }
}

bool Log::operator >(const Log& x)
{
  // if ((this->sign == 1) && (x.sign == 1))
  //   return (this->val > x.val);
  // else if ((this->sign == 1) && (x.sign == -1))
  //   return true;
  // else if ((this->sign == -1) && (x.sign == 1))
  //   return false;
  // else
  //   return (this->val < x.val);
  if (this->sign == 1)
    {
      if (x.sign == 1)
        return (this->val > x.val);
      else
        return true;
    }
  else if (this->sign == 0)
    {
      if (x.sign == -1)
        return true;
      else
        return false;
    }
  else
    {
      if (x.sign == -1)
        return (this->val < x.val);
      else
        return false;
    }
}

typedef std::vector<Log> VLog;
typedef std::vector<VLog> VVLog;

Log operator *(VLog& x, VLog& y)
{
  Log temp;
  for (int i=0; i<x.size(); ++i)
    temp += x[i] * y[i];
  return temp;
}
