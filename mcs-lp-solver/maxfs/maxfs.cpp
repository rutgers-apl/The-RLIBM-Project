#include "maxfs.h"

//Rational coeffs[6];

//double rlibm_poly_evaluation(double x, polynomial* poly){
//    Rational ret_val = 0.0;
//
//    // simulated Horner's method
//    for(int i = poly->termsize-1; i> 0; i--){
//        ret_val = ret_val + coeffs[i];
//        Rational xmul(1.0);
//        for(int j = 0; j < (poly->power[i] - poly->power[i-1]); j++){
//            xmul = xmul * x;
//        }
//        ret_val = ret_val * xmul;	  
//    }
//    ret_val = ret_val + coeffs[0];
//    
//    for(int j = 0; j < poly->power[0]; j++){
//        ret_val = ret_val * x;
//    }
//
//    return mpq_get_d(*(ret_val.getMpqPtr_w()));
//}


void rlibm_print_polyinfo(polynomial* p){
    if(p->termsize == 0){
        printf("Polynomial has no terms!\n");
        exit(0);
    }

    printf("Polynomial: y=%a x^(%d)", p->coeffs[0], p->power[0]);
    for(int j = 1; j < p->termsize; j++){
        printf(" + %a x^(%d)", p->coeffs[j], p->power[j]);
    }
    printf("\n");
}

size_t rlibm_compute_violated_indices(size_t* violated_indices, interval_data* intervals, size_t nentries, polynomial* poly)
{
    size_t num_violated_indices = 0;

    for(size_t i = 0; i < nentries; i++){
        double y = rlibm_poly_evaluation(intervals[i].x, poly);
        if( y < intervals[i].lb || y > intervals[i].ub){
            violated_indices[num_violated_indices] = i;
            num_violated_indices++;
        }
    }

    return num_violated_indices;
}

bool rlibm_validate_and_fix_intervals(sample_data* sintervals, size_t ssize, polynomial* poly)
{
    bool return_val = true;

    for(size_t i = 0; i < ssize; i++){
        double y = rlibm_poly_evaluation(sintervals[i].x, poly);

        if(y < sintervals[i].orig_lb){
            return_val = false;
            double_x lbx;
            lbx.d = sintervals[i].lb;
            if(lbx.d >= 0){
                lbx.x = lbx.x + 1;
            }
            else{
                lbx.x = lbx.x - 1 ;
            }
            sintervals[i].lb = lbx.d;
        }
        else if(y > sintervals[i].orig_ub){
            return_val = false;
            double_x ubx;
            ubx.d = sintervals[i].ub;
            if(ubx.d >= 0){
                ubx.x = ubx.x - 1;
            }
            else{
                ubx.x = ubx.x + 1;
            }
            sintervals[i].ub = ubx.d;
        }
    }

    return return_val;
}

double Compute_Loss(interval_data* intervals, unsigned long* sp_indices, unsigned long sp_indices_size,
                    polynomial* p, double *u, double *w, double *s, double *r, double alpha)
{
    double loss = 0.;

    for(unsigned long k=0; k < sp_indices_size; ++k){
        unsigned long index = sp_indices[k];

        // compute F
        loss += u[k];
        loss += w[k];

        // compute Q
        double y = rlibm_poly_evaluation(intervals[index].x, p);
        loss += alpha*(s[k] - u[k]*(y - intervals[index].ub));
        loss += alpha*(r[k] - w[k]*(intervals[index].lb - y));
    }

    return loss;
}

bool Compute_Maximum_Feasible_Subset(interval_data* intervals, unsigned long* sp_indices, unsigned long sp_indices_size,
                                     sample_data* basis_intervals, unsigned long basis_indices_size, polynomial* p, double alpha)
{
    double* u = (double*) calloc(sp_indices_size, sizeof(double));      // indicator variable for upper bound
    double* w = (double*) calloc(sp_indices_size, sizeof(double));      // indicator variable for lower bound
    double* s = (double*) calloc(sp_indices_size, sizeof(double));      // slack variable for upper bound
    double* r = (double*) calloc(sp_indices_size, sizeof(double));      // slack variable for lower bound

    // initialize indicator/slack variables
    for(unsigned k = 0; k < sp_indices_size; ++k){
        u[k] = 1.;
        w[k] = 1.;
        s[k] = 0.;
        r[k] = 0.;
    }

    //double loss = Compute_Loss(intervals, sp_indices, sp_indices_size, p, u, w, s, r, alpha);
    //printf("Initial Loss: %.70e\n", loss);

    SoPlex mysoplex;
    mysoplex.setBoolParam(SoPlex::RATFACJUMP,true);
    mysoplex.setIntParam(SoPlex::SOLVEMODE,2);
    mysoplex.setIntParam(SoPlex::CHECKMODE,2);
    mysoplex.setIntParam(SoPlex::SYNCMODE,1);
    mysoplex.setIntParam(SoPlex::READMODE,1);
    mysoplex.setRealParam(SoPlex::FEASTOL,0.0);
    mysoplex.setRealParam(SoPlex::OPTTOL,0.0);
    mysoplex.setRealParam(SoPlex::EPSILON_ZERO,0.0);
    mysoplex.setRealParam(SoPlex::EPSILON_FACTORIZATION,0.0);
    mysoplex.setRealParam(SoPlex::EPSILON_UPDATE,0.0);
    mysoplex.setRealParam(SoPlex::EPSILON_PIVOT,0.0);
    mysoplex.setIntParam(SoPlex::VERBOSITY,0);
    mysoplex.setRealParam(SoPlex::TIMELIMIT,5*60);

    // set objective sense
    mysoplex.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);

    DSVectorRational dummycol(0);

    // objective function for slack variables
    for(unsigned long k = 0; k < sp_indices_size; ++k){
        auto column1 = LPColRational(1.0, dummycol, infinity, 0.);
        mysoplex.addColRational(column1);
        auto column2 = LPColRational(1.0, dummycol, infinity, 0.);
        mysoplex.addColRational(column2);
    }

    // objective function for polynomial coefficients
    for(int i = 0; i < p->termsize; i++){
        Rational coeff(0.0);
        for(unsigned long k = 0; k < sp_indices_size; ++k){
            unsigned long index = sp_indices[k];
            Rational xValR(intervals[index].x);

            Rational toAdd(1.0);
            for(int t = 0; t < p->power[i]; ++t) toAdd*=xValR;

            coeff += Rational(w[k] - u[k])*toAdd;
        }
        auto column = LPColRational(coeff, dummycol, infinity, -infinity);
        mysoplex.addColRational(column);
    }

    // add linear constraints for special cases
    for(unsigned long k = 0; k < sp_indices_size; ++k){
        unsigned long index = sp_indices[k];
        DSVectorRational row1(2*sp_indices_size + p->termsize);
        DSVectorRational row2(2*sp_indices_size + p->termsize);
        Rational xValR(intervals[index].x);

        row1.add(2*k, Rational(1.0));       // s
        row2.add(2*k+1, Rational(1.0));     // r

        for(int j=0; j < p->termsize;j++){
            Rational toAdd(1.0);
            for(int t = 0; t < p->power[j]; t++) toAdd*=xValR;

            row1.add(2*sp_indices_size + j, -toAdd);
            row2.add(2*sp_indices_size + j, toAdd);
        }

        double lbnd= intervals[index].lb;
        double ubnd= intervals[index].ub;
        mysoplex.addRowRational(LPRowRational(-ubnd,row1,infinity));
        mysoplex.addRowRational(LPRowRational(lbnd,row2,infinity));
    }

    // add linear constraints for basis intervals
    for(unsigned long k = 0; k < basis_indices_size; ++k){
        Rational xValR(basis_intervals[k].x);
        DSVectorRational row1(2*sp_indices_size + p->termsize);

        for(int j = 0; j < p->termsize; j++){
            Rational toAdd(1.0);
            for(int t = 0; t < p->power[j]; t++) toAdd*=xValR;

            row1.add(2*sp_indices_size + j,toAdd);
        }

        double lbnd= basis_intervals[k].lb;
        double ubnd= basis_intervals[k].ub;
        mysoplex.addRowRational(LPRowRational(lbnd,row1,ubnd));
    }

    /* solve LP */
    SPxSolver::Status stat;
    stat=mysoplex.optimize();

    bool solution_found = false;
    
    /* get solution */
    if(stat==SPxSolver::OPTIMAL){
        solution_found = true;
        DVectorRational prim(2*sp_indices_size + p->termsize);
        mysoplex.getPrimalRational(prim);

        for(unsigned long k = 0; k < sp_indices_size; ++k){
            s[k] = mpq_get_d(*(prim[2*k].getMpqPtr_w()));
            r[k] = mpq_get_d(*(prim[2*k+1].getMpqPtr_w()));
        }
      
        for(int i = 0; i < p->termsize; i++){
            p->coeffs[i] = mpq_get_d(*(prim[2*sp_indices_size + i].getMpqPtr_w()));
            //coeffs[i] = prim[2*sp_indices_size + i];
        }
    }

    if(solution_found){
        //for(unsigned long k = 0; k < sp_indices_size; ++k){
        //    unsigned long index = sp_indices[k];
        //    double y = rlibm_poly_evaluation(intervals[index].x, p);
        //    double lbnd= intervals[index].lb;
        //    double ubnd= intervals[index].ub;

        //    u[k] = (1. - alpha*(y - ubnd) <= 0.)?1.:0.;
        //    w[k] = (1. - alpha*(lbnd - y) <= 0.)?1.:0.;
        //}

        //double loss_new = Compute_Loss(intervals, sp_indices, sp_indices_size, p, u, w, s, r, alpha);
        //printf("Final Loss: %.70e\n", loss_new);
    }
    else printf("No solution found.\n");

    free(u);
    free(w);
    free(s);
    free(r);

    return solution_found;
}

my_function function_to_process= NONE;

int main(int argc, char** argv)
{
    if(argc != 6){
        // input original interval file and basis file for the reduced set
        // (i.e., without special cases)
        printf("Usage: %s <interval file> <basis file> <violated indices file> <function config file> <function name>\n", argv[0]);
        return 1;
    }

    set_function_process(argv);

    // count the number of entries
    unsigned long nentries = Utilities::Number_Of_Intervals(argv);

    // allocate memory for intervals
    interval_data* intervals = (interval_data*) calloc(nentries, sizeof(interval_data));

    // read intervals from file
    Utilities::Read_Intervals_From_File(argv, intervals);

    size_t sp_indices_size = 0;
    unsigned long* sp_indices = (unsigned long*) calloc(nentries, sizeof(unsigned long));

    FILE* fp = fopen(argv[3], "r");
    assert(fp != nullptr);
    int retval = fscanf(fp, "%ld", &sp_indices_size);
    for(size_t i=0; i<sp_indices_size; ++i){
        size_t index;
        retval = fscanf(fp, "%ld", &index);
        sp_indices[i] = index;
    }
    fclose(fp);

    std::fstream inFile(argv[2], std::ios::in);
    unsigned long nindices = 0;
    inFile >> nindices;
    printf("Number of basis indices: %lu\n", nindices);
    unsigned long org_nindices = nindices;

    // allocate memory for intervals
    sample_data* basis_intervals = (sample_data*) calloc(nentries, sizeof(sample_data));
    unsigned long* basis_indices = (unsigned long*) calloc(nentries, sizeof(unsigned long));

    for(unsigned long k = 0; k < nindices; ++k) inFile >> basis_indices[k];
    inFile.close();

    FILE* powers_file = fopen(argv[4], "r");
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
    
    polynomial* p = (polynomial*) calloc(1, sizeof(polynomial));
    p->termsize = powers_size;
    p->power = powers;
    p->coeffs = (double*) calloc(powers_size, sizeof(double));

#if 0
    // log
    p->coeffs[0] = 9.9999999999999555910790149937383830547332763671875000000000000000000000e-01;
    p->coeffs[1] = -5.0000000184725479357439326122403144836425781250000000000000000000000000e-01;
    p->coeffs[2] = 3.3333510938414240287741563406598288565874099731445312500000000000000000e-01;
    p->coeffs[3] = -2.5068839702791845258644798377645201981067657470703125000000000000000000e-01;
    p->coeffs[4] = 3.1664776292877067076503294629219453781843185424804687500000000000000000e-01;
    p->coeffs[5] = -7.0189102481577148040514657623134553432464599609375000000000000000000000e+00; 
#endif

#if 0
    // log10
    p->coeffs[0] = 4.3429448190459402079355299974849913269281387329101562500000000000000000e-01;
    p->coeffs[1] = -2.1714724281710182096638561688450863584876060485839843750000000000000000e-01;
    p->coeffs[2] = 1.4476567429437961487259656223614001646637916564941406250000000000000000e-01;
    p->coeffs[3] = -1.0872812888267488484750344923668308183550834655761718750000000000000000e-01;
    p->coeffs[4] = 9.6039176808723147882318471602047793567180633544921875000000000000000000e-02;
    p->coeffs[5] = 0.0;
#endif

#if 0
    // exp
    p->coeffs[0] = 1.0000000000000002220446049250313080847263336181640625000000000000000000e+00;
    p->coeffs[1] = 1.0000000000002970956813896918902173638343811035156250000000000000000000e+00;
    p->coeffs[2] = 4.9999999999114697057933653923100791871547698974609375000000000000000000e-01;
    p->coeffs[3] = 1.6666664683520868162069916706968797370791435241699218750000000000000000e-01;
    p->coeffs[4] = 4.1666814544402248909893415884653222747147083282470703125000000000000000e-02;
    p->coeffs[5] = 8.5181762320615077299645889752355287782847881317138671875000000000000000e-03;
#endif

    // allocate memory for violated indices
    unsigned long* violated_indices = (unsigned long*) calloc(nentries, sizeof(unsigned long));
    unsigned long nviolated_indices = 0, min_violated_indices = nentries;
    double alpha = 1.;

    bool iterate = true;

    while(iterate){
        iterate = false;

        for(unsigned long k = 0; k < nindices; ++k) basis_intervals[k] = intervals[basis_indices[k]];

        for(int count = 0; count < MAX_TRIES; ++count){
            bool solution_found = Compute_Maximum_Feasible_Subset(intervals, sp_indices, sp_indices_size, basis_intervals, nindices, p, alpha);
            if(solution_found && rlibm_validate_and_fix_intervals(basis_intervals, nindices, p)) break;
        }

        nviolated_indices = rlibm_compute_violated_indices(violated_indices, intervals, nentries, p);
        printf("Number of violated intervals: %lu\n", nviolated_indices);

	if(nviolated_indices <40){
	  for(int i = 0; i<nviolated_indices; i++){
	    size_t index = violated_indices[i];
	    printf("input x=%a, lb=%a, ub=%a\n", intervals[index].x, intervals[index].lb, intervals[index].ub);
	  }

          if(min_violated_indices > nviolated_indices){
            min_violated_indices = nviolated_indices;
            printf("Writing violated indices...\n");

	    char filename[100];
	    sprintf(filename, "%s_violated_inputs", argv[5]);
	    FILE *fpv = fopen(filename, "w");
	    for(size_t m=0; m < nviolated_indices; ++m){
	      double data[1];
	      data[0] = intervals[violated_indices[m]].x;
	      
	      size_t bytes = fwrite(data, sizeof(double), 1, fpv);
	    }
	    fclose(fpv);	  
          }
	}

        rlibm_print_polyinfo(p);

        std::unordered_set<unsigned long> basis_set, sp_set;
        for(unsigned long k = 0; k < nindices; ++k) basis_set.insert(basis_indices[k]);
        for(unsigned long k = 0; k < sp_indices_size; ++k) sp_set.insert(sp_indices[k]);

        printf("Number of bases: %lu, Number of special cases: %lu\n", nindices, sp_indices_size);
        unsigned long current_sp_indices_size = sp_indices_size;

        for(unsigned long k = 0; k < nviolated_indices; ++k){
            if(basis_set.count(violated_indices[k]) == 0 && sp_set.count(violated_indices[k]) == 0){
                iterate = true;
                if(current_sp_indices_size > 0)
                    basis_indices[nindices++] = violated_indices[k];
                else sp_indices[sp_indices_size++] = violated_indices[k];
            }
            else if(basis_set.count(violated_indices[k]) != 0 && sp_set.count(violated_indices[k]) == 0){
                iterate = true;
                sp_indices[sp_indices_size++] = violated_indices[k];

                // remove violated index from basis
                unsigned long t = 0;
                for(; t < nindices; ++t) if(basis_indices[t] == violated_indices[k]) break;
                for(unsigned long j = t; j < nindices-1; ++j) basis_indices[j] = basis_indices[j+1];
                --nindices;
            }
            else if(nindices > org_nindices){
                iterate = true;
                while(nindices != org_nindices){
                    sp_indices[sp_indices_size++] = basis_indices[--nindices];
                }
            }
        }

        printf("New number of bases: %lu, New number of special cases: %lu\n", nindices, sp_indices_size);
    }

    rlibm_print_polyinfo(p);

    return 0;
}
