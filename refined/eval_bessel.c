#include "logbesselk.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>

int main() {
   
    // generate a sequence of x and nu in which x from 0.001 to 141 with interval 1 and nu from 0.01 to
    // 0.5 with interval 0.025, from 0.5 to 4 with interval 0.5, from 4 to 22 with interval 2
    // int x_count = 100; // x goes from 0.001 to 141 with interval 1
    // int nu_count = 100; // nu has three ranges: 21 values in [0.01, 0.5], 6 values in [0.5, 4], 9 values in [4, 22]
    
    // int total_count = x_count * nu_count;
    // int total = total_count * 2;
    // double *flattened_list = (double *)malloc(total * sizeof(double));
    // double *result = (double*)malloc(total_count * sizeof(double));

    // // Generate x from 0.001 to 141 with interval 1
    // int index = 0;
    // for (int i = 0; i < x_count; ++i) {
    //     double x = 0.001 + i * 0.01;
    //     for (int j = 0; j < nu_count; ++j){
    //         double nu = 0.001 + j * 0.01;
    //         flattened_list[index++] = x;
    //         flattened_list[index++] = nu;
    //     }
    //     // // Generate nu in the first range [0.001, 0.5] with interval 0.025
    //     // for (int j = 0; j < 20; ++j) {
    //     //     double nu = 0.001 + j * 0.025;
    //     //     flattened_list[index++] = x;
    //     //     flattened_list[index++] = nu;
    //     // }

    //     // // Generate nu in the second range [0.5, 4] with interval 0.5
    //     // for (int j = 0; j < 7; ++j) {
    //     //     double nu = 0.5 + j * 0.5;
    //     //     flattened_list[index++] = x;
    //     //     flattened_list[index++] = nu;
    //     // }

    //     // // Generate nu in the third range [4, 22] with interval 2
    //     // for (int j = 0; j < 9; ++j) {
    //     //     double nu = 4 + j * 2;
    //     //     flattened_list[index++] = x;
    //     //     flattened_list[index++] = nu;
    //     // }
    // }
 
    // // Call the CUDA function
    // // the CUDA function should have separate x and nu arrays
    // double *x = (double *)malloc(total_count * sizeof(double));
    // for (int i = 0; i < total_count; ++i) {
    //     x[i] = flattened_list[2 * i];
    // }
    // double *nu = (double *)malloc(total_count * sizeof(double));
    // for (int i = 0; i < total_count; ++i) {
    //     nu[i] = flattened_list[2 * i + 1];
    // }
    // BesselK_CUDA(x, nu, result, total_count);

    // // log the results to a file in the format of x, nu, gpu_val with header
    // FILE *f = fopen("new_refined_newnew.csv", "w");
    // fprintf(f, "x, nu, gpu_val\n");
    // for (int i = 0; i < total_count; ++i) {
    //     fprintf(f, "%lf, %lf, %.17f\n", x[i], nu[i], log(result[i]));
    // }
    // fclose(f);
    // free(x);
    // free(nu);
    // free(result);
    // free(flattened_list);

    // Generate the log value of Beseelk_CUA (x, nu) for x = 0.1 to 1, nu = 10
    // int x_count = 10;
    // double *x = (double *)malloc(x_count * sizeof(double));
    // double *nu = (double *)malloc(x_count * sizeof(double));
    // double *result = (double *)malloc(x_count * sizeof(double));
    // for (int i = 0; i < x_count; ++i) {
    //     x[i] = 0.15;
    //     nu[i] = 10 + i * 1.0;
    // }
    // BesselK_CUDA(x, nu, result, x_count);
    // // logarithm the results
    // for (int i = 0; i < x_count; ++i) {
    //     result[i] = log(result[i]);
    // }

    // for (int i = 0; i < x_count; ++i) {
    //     printf("x = %lf, nu = %lf, log(BesselK_CUDA(x, nu)) = %.17f\n", x[i], nu[i], log(result[i]));
    // }

    // define a double array which saves 1, 2, 4
//     double tt[10] = {38.01072697839095,
// 42.90364172999614,
// 47.89185530010114,
// 52.96707172588464,
// 58.12232430429789,
// 63.351679704360485,
// 68.6500238554363,
// 74.01290317971952,
// 79.43640436867013,
// 84.91706167381598
// }; 
    
    // calcuate the relative error which is log(1 + |x - y| / epsilon)
    // double epsilon = 2.2204460492503131e-16;
    // for (int i = 0; i < x_count; ++i) {
    //     double relative_error = log(1 + fabs(result[i] - tt[i]) / epsilon);
    //     printf("x = %lf, nu = %lf, relative_error = %.17f\n", x[i], nu[i], relative_error);
    // }

    // free(x);
    // free(nu);
    // free(result);    
    
    // Calculate array sizes based on the ranges
    int x_count = (0.1 - 0.001) / 0.005 + 1;  // From 0.001 to 0.1 with step 0.005
    int nu_count = (5.0 - 0.001) / 0.25 + 1;  // From 0.001 to 5 with step 0.25

    int total_count = x_count * nu_count;
    int total = total_count * 2;
    double *flattened_list = (double *)malloc(total * sizeof(double));
    double *result = (double *)malloc(total_count * sizeof(double));

    // Generate x and nu values
    int index = 0;
    for (int i = 0; i < x_count; ++i) {
        double x = 0.001 + i * 0.005;  // Step of 0.005 for x
        for (int j = 0; j < nu_count; ++j) {
            double nu = 0.001 + j * 0.25;  // Step of 0.25 for nu
            flattened_list[index++] = x;
            flattened_list[index++] = nu;
        }
    }

    // Separate x and nu arrays for CUDA
    double *x = (double *)malloc(total_count * sizeof(double));
    double *nu = (double *)malloc(total_count * sizeof(double));
    for (int i = 0; i < total_count; ++i) {
        x[i] = flattened_list[2 * i];
        nu[i] = flattened_list[2 * i + 1];
    }

    // Call the CUDA function
    BesselK_CUDA(x, nu, result, total_count);

    // Log results to CSV file
    FILE *f = fopen("newnew_refined.csv", "w");
    fprintf(f, "x,nu,refined_val\n");  // Header
    for (int i = 0; i < total_count; ++i) {
        fprintf(f, "%.15f,%.15f,%.15f\n", x[i], nu[i], log(result[i]));
    }
    fclose(f);

    // Free memory
    free(x);
    free(nu);
    free(result);
    free(flattened_list);



    return 0;
}
