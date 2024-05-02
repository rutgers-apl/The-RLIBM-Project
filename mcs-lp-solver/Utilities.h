#include "Interval_Data.h"
#include <fstream>

class Utilities
{
  public:

    static unsigned long Number_Of_Entries(FILE* fp)
    {
        fseek(fp, 0, SEEK_END);
        unsigned long nentries = ftell(fp);
        nentries = nentries/(3*sizeof(double));
        fseek(fp, 0, SEEK_SET);

        return nentries;
    }

    static void Read_Intervals(FILE* fp, Interval_Data* intervals, unsigned long nentries)
    {
        for (unsigned long i = 0; i < nentries; i++){
            double data_entry[3];
            size_t bytes = fread(data_entry, sizeof(double), 3, fp);
            intervals[i].x = data_entry[0];
            intervals[i].lb = data_entry[1];
            intervals[i].ub = data_entry[2];
        }
    }
};
