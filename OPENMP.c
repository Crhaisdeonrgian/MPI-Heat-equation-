#include <math.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <stdlib.h>


int main(int argc, char** argv) {
        int numthreads=4;
        const long long n = 10;
        double maxtime = 1;
        double x_n =1 / (double)n;
        double time_n = 0.3 * x_n * x_n;
        double tic = omp_get_wtime();
        double *Un;
        double *Un1;
        Un = (double*)malloc((n+2) * sizeof(double));
        Un1 = (double*)malloc((n+2) * sizeof(double));
        Un[n+1]=1;
        Un1[n+1]=1;
        int i=0;
        double timenow=0;
        #pragma omp parallel num_threads(numthreads) 
        {
                for (timenow = 0; timenow <= maxtime; timenow += time_n) {

                        #pragma omp for schedule(dynamic, 1)
                                for (i = 1; i < n+ 1; ++i) {
                                        Un1[i] = Un[i] + 0.3 * (Un[i - 1] - 2 * Un[i] + Un[i + 1]);
                        }

                        #pragma omp master
                        {
                                for (i = 1; i < n + 1; i++) {
                        Un[i] = Un1[i];
                }
                        }
                        #pragma omp barrier
                }
        }



        for (i = 1; i < n+1; i++) {
                                printf("№:%lf  |  %lf \n",i*x_n, Un[i]);
        }

        double toc = omp_get_wtime();
        double timelast = toc - tic;
         printf("Время: %lf",timelast);
        free(Un);
        free(Un1);
        return 0;

}
