//!#####################################################################
//! \file Interval_Data.h
//!#####################################################################
// Class Interval_Data
//######################################################################
#ifndef __Interval_Data__
#define __Interval_Data__

class interval_data
{
  public:
    double x;         /* original input */
    double lb;        /* lower bound */
    double ub;        /* upper bound */

    interval_data() {}

    interval_data(const double& x_input,const double& lb_input,const double& ub_input)
        :x(x_input),lb(lb_input),ub(ub_input)
    {}

    ~interval_data() {}

    interval_data& operator=(const interval_data& rhs)
    {
        if(this==&rhs) return *this;

        x = rhs.x;
        lb = rhs.lb;
        ub = rhs.ub;
        return *this;
    }
};
#endif
