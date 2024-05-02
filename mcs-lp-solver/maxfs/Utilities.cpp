//!#####################################################################
//! \file Utilities.cpp
//!#####################################################################
#include "Utilities.h"
#include <cassert>
#include <fstream>
//######################################################################
// Number_Of_Intervals
//######################################################################
unsigned long Utilities::
Number_Of_Intervals(char** argv)
{
    FILE* fp = fopen(argv[1], "r");
    assert(fp != nullptr);

    fseek(fp, 0, SEEK_END);
    unsigned long nentries = ftell(fp);
    nentries /= (3*sizeof(double));
    printf("number of intervals in file = %lu\n", nentries);
    fclose(fp);

    return nentries;
}
//######################################################################
// Read_Intervals_From_Files
//######################################################################
void Utilities::
Read_Intervals_From_File(char** argv,interval_data* intervals)
{
    FILE* fp = fopen(argv[1], "r");
    assert(fp != nullptr);

    fseek(fp, 0, SEEK_END);
    unsigned long nentries = ftell(fp);
    nentries /= (3*sizeof(double));
    fseek(fp, 0, SEEK_SET);

    for(unsigned long k = 0; k < nentries; ++k){
        double data_entry[3];
        size_t bytes = fread(data_entry, sizeof(double), 3, fp);
        intervals[k] = interval_data(data_entry[0],data_entry[1],data_entry[2]);
    }
    fclose(fp);
}
//######################################################################
