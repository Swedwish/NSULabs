#include <iostream>
#include <math.h>
#include <cstdio>
#include <stdlib.h>
#include <cstring>
#include <wchar.h>
#include "omp.h"

bool check(double* A, double* b, double* x, int N) {
    double* ans = vec_x_mat(A, x, N);
    double eps = 0.0001;
    for (int i = 0; i < N; i++) {
        if (!((ans[i] + 100*eps > b[i]) && (ans[i] - 100*eps < b[i])))
            return 0;
    }
    return 1;
}


double* solve(int N, double* x, double* A, double* b) {
    double alpha, beta, sum1 = 0, sum2 = 0, eps = 0.00001, *vxc = new double[N], * vxm = new double[N], *r_prev = new double[N], * z = new double[N], *r = new double[N], const1;
#pragma omp parallel //num_threads(4)
    {
#pragma omp for  
    for (int i = 0; i < N; i++) {
        vxm[i] = 0;
        for (int j = 0; j < N; j++) {
            vxm[i] += x[j] * A[i * N + j];
        }
    }

#pragma omp for  
    for (int i = 0; i < N; i++) {
        r[i] = b[i] - vxm[i];
    }

    //double* r = plusminus(b, vec_x_mat(A, x, N), N, -1);
#pragma omp for  
    for (int i = 0; i < N; i++) {
        sum1 += b[i] * b[i];
    }
#pragma omp single
    const1 = sqrt(sum1);
    //double const1 = module(b, N);

#pragma omp single
    memcpy(z, r, sizeof(double) * N);

        while (module(r, N) / const1 > eps) {
#pragma omp single
            memcpy(r_prev, r, sizeof(double) * N);

#pragma omp single
            sum1 = 0;

#pragma omp for   reduction(+:sum1)
            for (int i = 0; i < N; i++) {
                sum1 += r[i] * r[i];
            }

#pragma omp for  
            for (int i = 0; i < N; i++) {
                vxm[i] = 0;
                for (int j = 0; j < N; j++) {
                    vxm[i] += z[j] * A[i * N + j];
                }
            }

#pragma omp single
            sum2 = 0;
#pragma omp for   reduction(+:sum2)
            for (int i = 0; i < N; i++) {
                sum2 += vxm[i] * z[i];
            }

#pragma omp single
            alpha = sum1 / sum2;
            //alpha = scalar(r, r, N) / (scalar(vec_x_mat(A, z, N), z, N));

#pragma omp for  
            for (int i = 0; i < N; i++) {
                vxc[i] = z[i] * alpha;
            }

#pragma omp for  
            for (int i = 0; i < N; i++) {
                x[i] = x[i] + vxc[i];
            }
            //x = plusminus(x, vxc = vec_x_const(z, alpha, N), N, 1);
            
#pragma omp for  
            for (int i = 0; i < N; i++) {
                vxm[i] = 0;
                for (int j = 0; j < N; j++) {
                    vxm[i] += z[j] * A[i * N + j];
                }
                vxm[i] *= alpha;
            }
            
#pragma omp for  
            for (int i = 0; i < N; i++) {
                r[i] = r[i] - vxm[i];
            }
            //r = plusminus(r, vxm = vec_x_mat(A, z, N, alpha), N, -1);

#pragma omp single
            sum1 = 0;

#pragma omp for   reduction(+:sum1)
            for (int i = 0; i < N; i++) {
                sum1 += r[i] * r[i];
            }

#pragma omp single
            sum2 = 0;

#pragma omp for   reduction(+:sum2)
            for (int i = 0; i < N; i++) {
                sum2 += r_prev[i] * r_prev[i];
            }

#pragma omp single
            beta = sum1 / sum2;
            //beta = scalar(r, r, N) / scalar(r_prev, r_prev, N);

#pragma omp for  
            for (int i = 0; i < N; i++) {
                vxc[i] = z[i] * beta;
            }

#pragma omp for  
            for (int i = 0; i < N; i++) {
                z[i] = r[i] + vxc[i];
            }
            //z = plusminus(r, vxc = vec_x_const(z, beta, N), N, 1);
        }
    }
    return x;
}


int main() {
    int N = 2160;
    double* x = mk_x0(N);
    double* A = random_sym_matrix(N);
    double* b = mk_b(N);
    double start = omp_get_wtime();
    x = solve(N, x, A, b);
    double end = omp_get_wtime();
    printf("%f\n", end - start);
    //print_mat(1, N, x);
    if (check(A, b, x, N) == 0) {
        std::cout << "Incorrect";
    }

    return 0;
}

