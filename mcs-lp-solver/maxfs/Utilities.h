//!#####################################################################
//! \file Utilities.h
//!#####################################################################
// Class Utilities
//######################################################################
#ifndef __Utilities__
#define __Utilities__

#include "Interval_Data.h"

class Utilities
{
  public:

//######################################################################
    static unsigned long Number_Of_Intervals(char** argv);
    static void Read_Intervals_From_File(char** argv,interval_data* intervals);
//######################################################################
};
#endif
