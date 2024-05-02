//!#####################################################################
//! \file Sample_Data.h
//!#####################################################################
// Class Sample_Data
//######################################################################
#ifndef __Sample_Data__
#define __Sample_Data__

#include "Interval_Data.h"

class sample_data
{
  public:
    double x;           /* original input */
    double lb;          /* lower bound */
    double ub;          /* upper bound */
    double orig_lb;     /* original lower bound */
    double orig_ub;     /* original upper bound */

    sample_data() {}

    sample_data(const double& x_input,const double& lb_input,const double& ub_input)
        :x(x_input),lb(lb_input),ub(ub_input),orig_lb(lb_input),orig_ub(ub_input)
    {}

    ~sample_data() {}

    sample_data& operator=(const interval_data& rhs)
    {
        x = rhs.x;
        lb = rhs.lb;
        ub = rhs.ub;
        orig_lb = rhs.lb;
        orig_ub = rhs.ub;
        return *this;
    }

    sample_data& operator=(const sample_data& rhs)
    {
        if(this==&rhs) return *this;

        x = rhs.x;
        lb = rhs.lb;
        ub = rhs.ub;
        orig_lb = rhs.orig_lb;
        orig_ub = rhs.orig_ub;
        return *this;
    }
};
#endif
