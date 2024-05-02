#include "Utilities.h"
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/property_map.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Convex_hull_traits_adapter_2<K,CGAL::Pointer_property_map<Point_2>::type> Convex_hull_traits_2;

void Compute_Extreme_Points(Interval_Data* intervals, unsigned long nentries, std::set<size_t>& index_set,
                            Point_2* points, int* powers, int powers_size, bool bound_flag)
{
    if(bound_flag) std::cout<<"Processing upper halfspaces..."<<std::endl;
    else std::cout<<"Processing lower halfspaces..."<<std::endl;

    for(int k=0; k<powers_size; ++k) if(powers[k]>0){                 // loops over all powers
        for(size_t i=0; i<nentries; ++i){
            K::FT value(intervals[i].x),product(1.0),bound = bound_flag ? K::FT(intervals[i].ub) : K::FT(intervals[i].lb);
            for(int t=0; t<powers[k]; ++t) product*=value;
            points[i] = Point_2(product,bound);
        }
        std::vector<std::size_t> indices(nentries),out;
        std::iota(indices.begin(), indices.end(),0);
        CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
                            Convex_hull_traits_2(CGAL::make_property_map(points)));

        for(std::size_t i : out)
            if(index_set.count(i)==0) index_set.insert(i);
        std::cout<<"Processing power "<<powers[k]<<" the convex hull now has "<<index_set.size()<<" points."<<std::endl;
        break;
    }
}

int main(int argc, char** argv)
{
    if(argc != 4){
        printf("Usage: %s <interval file> <name of the function's configuration file> <function name>\n", argv[0]);
        return 1;
    }

    FILE* fp = fopen(argv[1], "r");
    assert(fp != nullptr);

    unsigned long nentries = Utilities::Number_Of_Entries(fp);
    printf("number of intervals = %lu\n", nentries);

    Interval_Data* intervals = (Interval_Data*) calloc(nentries, sizeof(Interval_Data));
    Utilities::Read_Intervals(fp, intervals, nentries);
    fclose(fp);

    FILE* powers_file = fopen(argv[2], "r");
    assert(powers_file != nullptr);

    int N_RLIBM_PIECES = 1;
    int powers_size = 1;
    
    int s_ret_value = fscanf(powers_file, "%d\n", &N_RLIBM_PIECES);
    assert(s_ret_value == 1);
    
    s_ret_value = fscanf(powers_file, "%d", &powers_size);
    assert(s_ret_value == 1);
  
    int* powers = (int*) calloc(powers_size, sizeof(int));
    for(int i = 0; i < powers_size; i++){
      s_ret_value = fscanf(powers_file, "%d", &powers[i]);
      assert(s_ret_value == 1);
    }
    fclose(powers_file);
        
    Point_2* points = (Point_2*) calloc(nentries, sizeof(Point_2));
    std::set<size_t> upper_index_set, lower_index_set;

    // compute extreme constraints
    Compute_Extreme_Points(intervals, nentries, upper_index_set, points, powers, powers_size, true);
    Compute_Extreme_Points(intervals, nentries, lower_index_set, points, powers, powers_size, false);

    printf("# of extreme constraints: %ld\n", upper_index_set.size()+lower_index_set.size());

    char filename[100] ="";
    sprintf(filename, "%s_upper_hull_indices.txt", argv[3]);
    
    fp = fopen(filename,"w");
    fprintf(fp, "%ld\n", upper_index_set.size());
    for(std::set<size_t>::iterator it=upper_index_set.begin(); it!=upper_index_set.end(); ++it){
        fprintf(fp, "%ld\n", *it);
    }
    fclose(fp);

    sprintf(filename, "%s_lower_hull_indices.txt", argv[3]);
    fp = fopen(filename,"w");
    fprintf(fp, "%ld\n", lower_index_set.size());
    for(std::set<size_t>::iterator it=lower_index_set.begin(); it!=lower_index_set.end(); ++it){
        fprintf(fp, "%ld\n", *it);
    }
    fclose(fp);

    return 0;
}
