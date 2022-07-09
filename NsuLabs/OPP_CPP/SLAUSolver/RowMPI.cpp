#include <iostream>
#include <math.h>
#include <cstdio>
#include <stdlib.h>
#include <cstring>
#include <wchar.h>
#include "mpi.h"
#include <cassert>

void random_sym_matrix(double* res, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if (i != j) {
                res[N * i + j] = rand() % 201 - 100;
                res[N * j + i] = res[N * i + j];
                 
            }
            else {
                res[N * i + j] = rand() % 201 + 100;
            }
        }
    }
}

void mk_b(double* b, int N) {
    for (int i = 0; i < N; i++) {
        b[i] = rand() % 201 - 100;
    }
}


void print_mat(int N, int M, double* A) {
    int rank, numtasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cout << rank << ":";
    for (int i = 0; i < N * M; i++) {
        std::cout << A[i] << " ";
        if ((i + 1) % M == 0) std::cout << "\n";
    }
}

double* vec_x_mat(double* A, double* x, int N, int start, int end,double mult = 1) {
    double* res = new double[N];
    for (int i = start; i < end; i++) res[i] = 0;
    //
    for (int i = start; i < end; i++) {
        for (int j = 0; j < N; j++) {
            res[i] += x[j] * A[i * N + j];
        }
        res[i] *= mult;
    }
    return res;
}

void mk_x0(double* x,int N) {
    for (int i = 0; i < N; i++) {
        x[i] = 0;
    }
}


double* plus(double* vec1, double* vec2, int N, int start, int end, int sign, int needAllGather = 1) {
    double* res = new double[N];
    double* res1 = new double[N];
    for (int i = start; i < end; i++) {
        res[i] = vec1[i] + vec2[i]*sign;
    }
    
    if ((start != 0) || (end != N) && needAllGather == 1) {
        MPI_Allgather(res + start, end - start, MPI_DOUBLE, res1, end - start, MPI_DOUBLE, MPI_COMM_WORLD);
        return res1;
    }
    
    return res;
}

double scalar(double* vec1, double* vec2, int N, int start, int end) {
    double sum = 0;
    for (int i = start; i < end; i++) {
        sum += vec1[i] * vec2[i];
    }
    double sum1 = 0;
    MPI_Allreduce(&sum, &sum1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sum1;
}

double module(double* A, int N, int start, int end) {
    return sqrt(scalar(A, A, N, start, end));
}

double* vec_x_const(double* vec1, double cons, int N, int start,int end) {
    double* res1 = new double[N];
    double* res = new double[N];
    for (int i = start; i < end; i++) {
        res[i] = vec1[i] * cons;
    }
    MPI_Allgather(res + start, end - start, MPI_DOUBLE, res1, end - start, MPI_DOUBLE, MPI_COMM_WORLD);
    return res1;
}

double* solve(int N, double* x, double* A, double* b, int start, int end) {
    double div = module(b, N, start, end);
    int rank, numtasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    double eps = 0.00001;
    double* r = plus(b, vec_x_mat(A, x, N, start, end), N, start, end,-1);
    //print_mat(1, N, r);
    double* z = new double[N];
    double alpha;
    double beta;
    double* r_prev = new double[N];
    memcpy(z, r, sizeof(double) * N);
    int c = 0; 
    while ((module(r, N, start, end) / div) > eps) {
        c++;
        memcpy(r_prev, r, sizeof(double) * N); 
        double* res = vec_x_mat(A, z, N, start, end);   
        double* res1 = new double[N]; 
        MPI_Allgather(res + start, end - start, MPI_DOUBLE, res1, end - start, MPI_DOUBLE, MPI_COMM_WORLD); 
        alpha = scalar(r, r, N, start, end) / (scalar(res1, z, N, start, end)); 
        x = plus(x, vec_x_const(z, alpha, N, start, end), N, 0 , N, 1, 0); 
        r = plus(r, vec_x_mat(A, z, N, start, end, alpha), N, start, end,-1); 
        beta = scalar(r, r, N,start, end) / scalar(r_prev, r_prev, N, start, end); 
        z = plus(r, vec_x_const(z, beta, N, start,end), N, 0, N,1); 
        
        if (c == 10000) break;
    }
    return x;
}

bool check(double* A, double* b, double* x, int N) {
    double* ans = vec_x_mat(A, x, N, 0, N);
    double eps = 0.0001;
    for (int i = 0; i < 0; i++) {
        if (!((ans[i] + eps > b[i]) && (ans[i] - eps < b[i])))
            return 0;
    }
    return 1;
}

int main(int argc, char ** argv) {
    int N = 1080; //2160
    
    int rank, numtasks;
    MPI_Init(&argc, &argv); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    double* A = new double[N * N];
    double* A1 = new double[N * N];
    double* x = new double[N];
    double* x1 = new double[N];
    double* b = new double[N];
    double* b1 = new double[N];
    if (rank == 0) {
        mk_x0(x,N);
        random_sym_matrix(A,N);
        mk_b(b,N);
    }
    int start = (N / numtasks) * rank;
    int end = (N / numtasks) * (rank + 1);
    MPI_Scatter(A, (N * N)/numtasks, MPI_DOUBLE, A1+start*N, (N * N) / numtasks, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD); (x, end - start, MPI_DOUBLE, x1 + start, end - start, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double start1 = MPI_Wtime(); 
    x = solve(N, x, A1, b, start, end); 
    double end1 = MPI_Wtime();
    /*if (rank == 0) {
        print_mat(1, N, x);
    }*/
    if (rank == 0) {
        std::cout << end1 - start1<<"\n";
    }
    if (rank == 0) {
        if (check(A, b, x, N) == 0) {
            std::cout << "Incorrect";
        }
    }
    MPI_Finalize();
    return 0;
}

