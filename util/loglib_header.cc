#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <float.h>
#include "loglib_header.hh"

Log::Log ()
{
  val = -DBL_MAX;
  sign = 0;
  // dm_1.c = DBL_MAX_1;
}

Log::Log (double val, signed char sign)
{
  this->val = val;
  this->sign = sign;
  // dm_1.c = DBL_MAX_1;
}

Log::Log (double a)
{
  // dm_1.c = DBL_MAX_1;
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
  if (sign == 0)
    return 0;
  else
    return sign * exp(val);
}

bool Log::iszero()
{
  return (sign == 0);
}

signed char Log::get_sign()
{
  return sign;
}

double Log::get_val()
{
  return val;
}

double Log::take_log()
{
  if (this->sign != 1)
    {
      std::cout << "take_log(<=0)" << std::endl;
      std::cout << val << "\t" << (int) sign << std::endl << std::endl;
      abort();
    }
  
  return val;
}

Log Log::take_log_Log()
{
  if (this->sign != 1)
    {
      std::cout << "take_log(<=0)" << std::endl;
      std::cout << val << "\t" << (int) sign << std::endl << std::endl;
      abort();
    }
  
  Log temp = Log(val);
  
  return temp;
}

// Log Log::take_pow(double b)
// {
//   if (this->sign == -1)
//     {
//       std::cout << "take_pow(<0)" << std::endl;
//       abort();
//     }
    
//   Log temp;
//   if (this->iszero())
//     {
//       temp = Log(0); // bug! 0^0 = 1, thus temp = Log(1)
//     }
//   else
//     {
//       temp.val = this->val * b;
//       if (temp.val <= -DBL_MAX )
//         temp = Log(0); // bug! e^{-DBL_MAX} = Log(1)
//       else
//         temp.sign = 1;
//     }
  
//   return temp;
// }

Log Log::take_pow(double b)
{
  if (this->sign == -1)
    {
      std::cout << "take_pow(<0)" << std::endl;
      abort();
    }
    
  Log temp;
  if (this->iszero())
    {
      if (b == 0)
        temp = Log(1);
      else
        temp = Log(0);
    }
  else
    {
      temp.val = this->val * b;
      if (temp.val <= -DBL_MAX )
        temp = Log(1);
      else
        temp.sign = 1;
    }
  
  return temp;
}

Log Log::take_exp()
{
  Log temp;
  
  temp.val = this->eval();
  
  if (temp.val <= -DBL_MAX)
    {
      temp = Log(0);
    }
  else if (temp.val >= DBL_MAX)
    {
      std::cout << "err: take_exp overflow!" << std::endl;
      abort();
    }
  else
    {
      temp.sign = 1;
    }
  
  return temp;
}

Log Log::inverse()
{
  if (this->sign == 0)
    {
      std::cout << "inverse(0)" << std::endl;
      abort();
    }
  
  Log temp;
  temp.val = -this->val;
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

  if (this->iszero())
    {
      temp.val = b.val;
      temp.sign = b.sign;
    }
  else if (b.sign == 0)
    {
      temp.val = this->val;
      temp.sign = this->sign;
    }
  else if (this->val > b.val)
    {
      temp.val = this->val + log1p(this->sign * b.sign * exp(b.val - this->val));
      if (temp.val <= -DBL_MAX)
        temp = Log(0);
      else
        temp.sign = this->sign;
    }
  else
    {
      temp.val = b.val + log1p(this->sign * b.sign * exp(this->val - b.val));
      if (temp.val <= -DBL_MAX)
        temp = Log(0);
      else
        temp.sign = b.sign;
    }

  return temp;
}

Log Log::operator +=(const Log& x)
{
  if (this->iszero())
    {
      this->val = x.val;
      this->sign = x.sign;
    }
  else if (x.sign == 0)
    {
    }
  else if (this->val > x.val)
    {
      this->val += log1p(this->sign * x.sign * exp(x.val - this->val));
      if (this->val <= -DBL_MAX)
        {
          this->val = -DBL_MAX;
          this->sign = 0;
        }
    }
  else
    {
      this->val = x.val + log1p(this->sign * x.sign * exp(this->val - x.val));
      if (this->val <= -DBL_MAX)
        {
          this->val = -DBL_MAX;
          this->sign = 0;
        }
      else
        {
          this->sign = x.sign;
        }
    }
  return *this;
}

Log Log::operator -(const Log& b)
{
  Log temp;

  if (this->iszero())
    {
      temp.val = b.val;
      temp.sign = -b.sign;
    }
  else if (b.sign == 0)
    {
      temp.val = this->val;
      temp.sign = this->sign;
    }
  else if (this->val > b.val)
    {
      temp.val = this->val + log1p(this->sign * (-b.sign) * exp(b.val - this->val));
      if (temp.val <= -DBL_MAX)
        temp = Log(0);
      else
        temp.sign = this->sign;
    }
  else
    {
      temp.val = b.val + log1p(this->sign * (-b.sign) * exp(this->val - b.val));
      if (temp.val <= -DBL_MAX)
        temp = Log(0);
      else
        temp.sign = -b.sign;
    }

  return temp;
}

Log Log::operator -=(const Log& x)
{
  if (this->iszero())
    {
      this->val = x.val;
      this->sign = -x.sign;
    }
  else if (x.sign == 0)
    {
    }
  else if (this->val > x.val)
    {
      this->val += log1p(this->sign * (-x.sign) * exp(x.val - this->val));
      if (this->val <= -DBL_MAX)
        {
          this->val = -DBL_MAX;
          this->sign = 0;
        }
    }
  else
    {
      this->val = x.val + log1p(this->sign * (-x.sign) * exp(this->val - x.val));
      if (this->val <= -DBL_MAX)
        {
          this->val = -DBL_MAX;
          this->sign = 0;
        }
      else
        this->sign = -x.sign;
    }
  return *this;
}

Log Log::operator *(const Log& x)
{
  Log temp;

  if (this->iszero() || (x.sign == 0))
    {
      temp = Log(0);
    }
  else
    {
      temp.val = this->val + x.val;
      if (temp.val <= -DBL_MAX)
        temp = Log(0);
      else
        temp.sign = this->sign * x.sign;
    }
  
  return temp;
}

Log Log::operator *=(const Log& x)
{
  if (this->iszero())
    {
    }
  else if (x.sign == 0)
    {
      this->val = -DBL_MAX;
      this->sign = 0;
    }
  else
    {
      this->val += x.val;
      if (this->val <= -DBL_MAX)
        {
          this->val = -DBL_MAX;
          this->sign = 0;
        }
      else
        this->sign *= x.sign;
    }
  
  return *this;
}

Log Log::operator /(const Log& x)
{
  Log temp;

  if (x.sign == 0)
    {
      std::cout << "division by zero!" << std::endl;
      abort();
    }
  else if (this->iszero())
    {
      temp = Log(0);
    }
  else
    {
      temp.val = this->val - x.val;
      if (temp.val <= -DBL_MAX)
        temp = Log(0);
      else
        temp.sign = this->sign * x.sign;
    }
  
  return temp;
}

Log Log::operator /=(const Log& x)
{
  if (x.sign == 0)
    {
      std::cout << "division by zero!" << std::endl;
      abort();
    }
  else if (this->iszero())
    {
    }
  else
    {
      this->val -= x.val;
      if (this->val <= -DBL_MAX)
        {
          this->val = -DBL_MAX;
          this->sign = 0;
        }
      else
        this->sign *= x.sign;
    }
  
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

Log operator *(VLog& x, VLog& y)
{
  Log temp;
  for (int i=0; i<(int)x.size(); ++i)
    temp += x[i] * y[i];
  return temp;
}
