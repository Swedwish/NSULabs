#include <iostream>
#include <math.h>
#include <cstdio>
#include <stdlib.h>
#include <cstring>
#include <wchar.h>


double *random_sym_matrix(int N) {
    double *res = new double[N * N];
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if (i != j) {
                res[N * i + j] = rand() % 201 - 100;
                res[N * j + i] = res[N * i + j];

            } else {
                res[N * i + j] = rand() % 201 + 100;
            }
        }
    }
    return res;
}

double *mk_b(int N) {
    double *res = new double[N];
    for (int i = 0; i < N; i++){
        res[i] = rand() % 201 - 100;
    }
    return res;
}


void print_mat(int N,int M, double *A) {
    for (int i = 0; i < N * M; i++) {
        std::cout << A[i] << " ";
        if ((i + 1) % M == 0) std::cout << "\n";
    }
}

double *vec_x_mat(double *A, double *x, int N, double mult = 1) {
    double *res = new double[N];
    for (int i = 0; i < N; i++) res[i] = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++){
            res[i]+= x[j]*A[i*N+j];
        }
        res[i]*=mult;
    }
    return res;
}

double* mk_x0(int N){
    double *res = new double[N];
    for (int i = 0; i < N; i++){
        res[i] = 0;
    }
    return res;
}

double* minus(double* vec1, double* vec2, int N){
    double* res = new double[N];
    for (int i = 0; i < N; i++){
        res[i] = vec1[i] - vec2[i];
    }
    return res;
}

double* plus(double* vec1, double* vec2, int N){
    double* res = new double[N];
    for (int i = 0; i < N; i++){
        res[i] = vec1[i] + vec2[i];
    }
    return res;
}

double scalar(double *vec1, double *vec2, int N) {
    double sum = 0;
    for (int i = 0; i < N; i++){
        sum+= vec1[i]*vec2[i];
    }
    return sum;
}

double module(double* A, int N){
    return sqrt(scalar(A, A, N));
}

double* vec_x_const(double* vec1, double cons, int N){
    double* res = new double[N];
    for (int i = 0; i < N; i++){
        res[i] = vec1[i] * cons;
    }
    return res;
}

double* solve(int N, double* x, double* A, double* b){
    double eps = 0.00001;
    double* r = minus(b, vec_x_mat(A, x, N), N);
    double *z = new double[N];
    double alpha;
    double beta;
    double* r_prev = new double[N];
    memcpy(z, r, sizeof(double) * N);
    while (module(r,N)/module(b,N)>eps){
        memcpy(r_prev,r,sizeof(double)*N);
        alpha = scalar(r, r, N) / (scalar(vec_x_mat(A, z, N), z, N));
        x = plus(x, vec_x_const(z, alpha,N),N);
        r = minus(r, vec_x_mat(A,z,N,alpha),N);
        beta = scalar(r,r,N)/scalar(r_prev,r_prev,N);
        z = plus(r, vec_x_const(z,beta,N),N);
    }
    return x;
}


int main() {
    int N = 2160;
    double* x = mk_x0(N);
    double* A = random_sym_matrix(N);
    double* b = mk_b(N);
    x = solve(N,x,A,b);
    print_mat(1,N,x);
    return 0;
}
