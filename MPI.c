#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include "mpi.h"
#include <stdlib.h>


int main(int argc, char** argv) {
                int numthreads;
        int rank;
        MPI_Init(&argc, &argv);
                MPI_Comm_size(MPI_COMM_WORLD, &numthreads);
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const long long n = atoll(argv[1]);
        double maxtime = atof(argv[2]);
        double x_n =1 / (double)n;
        double time_n = 0.3 * x_n * x_n;

        long long my_n = n / numthreads;
        if(n % numthreads)
                        my_n++;
        long long last_n = n - my_n * ((numthreads) - 1);
        if (rank == numthreads - 1) {
                my_n = last_n;
        }

        double tic = MPI_Wtime();
        double *Un;
        double *Un1;
        Un = (double*)malloc((my_n+2) * sizeof(double));
        Un1 = (double*)malloc((my_n+2) * sizeof(double));
        if (rank == numthreads - 1) {
                Un[my_n+1] = 1;
        }
        int i=0;
        double timenow=0;
        for (timenow = 0; timenow <= maxtime; timenow += time_n) {
                double left, right;

                if (!(rank % 2)) {
                        left = Un[1];
                        right = Un[my_n];
                        if (rank != 0) {
                                MPI_Send(&left, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
                        }

                        if (rank != numthreads - 1) {
                                MPI_Send(&right, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
                        }
                }
                else {
                        MPI_Recv(&left, 1, MPI_DOUBLE, rank- 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                        if (rank == numthreads - 1) {
                                right = 1;
                        } else {
                                MPI_Recv(&right, 1, MPI_DOUBLE, rank+ 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        }
                        Un[0] = left;
                        Un[my_n+1] = right;
                }

                if (rank % 2) {
                        left = Un[1];
                        right = Un[my_n];

                        MPI_Send(&left, 1, MPI_DOUBLE, rank- 1, 0, MPI_COMM_WORLD);

                        if (rank != numthreads - 1) {
                                MPI_Send(&right, 1, MPI_DOUBLE, rank+ 1, 0, MPI_COMM_WORLD);
                        }
                }
                else {
                        if (rank == 0) {
                                left = 0;
                        }
                        else {
                                MPI_Recv(&left, 1, MPI_DOUBLE, rank- 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                        }

                        if (rank == numthreads- 1) {
                                right = 1;
                        }
                        else {
                                MPI_Recv(&right, 1, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        }
                        Un[0] = left;
                        Un[my_n+1] = right;
                }


                for (i = 1; i < my_n + 1; i++) {
                        Un1[i] = Un[i] + 0.3 * (Un[i - 1] - 2 * Un[i] + Un[i + 1]);
                }

                for (i = 1; i < my_n + 1; i++) {
                        Un[i] = Un1[i];
                }
        }

        for (i = 1; i < my_n+1; i++) {
                                printf("№: %d  | %lf  |  %lf \n",rank,rank* my_n*x_n + i*x_n, Un[i]);
        }


        if (rank == 0) {
                for (i = 0; i != n; i++) {
                }
                double toc = MPI_Wtime();
                double timelast = toc - tic;
                printf("Время: %lf",timelast);
        }
        free(Un);
        free(Un1);
        MPI_Finalize();
        return 0;
}
