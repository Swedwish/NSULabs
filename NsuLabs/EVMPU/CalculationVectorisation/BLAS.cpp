#include <iostream>
#include <xmmintrin.h>
#include <stdlib.h>
#include <cblas.h>
float* Transp(const float* A, int N) {    //Транспонирование матрицы
    float * T = new float [N*N];
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            T[i*N + j] = A[j*N + i];
    return T;
}

float* mk_I(int N){
    float * res = new float [N*N];
    for (int i = 0; i < N*N; i++){
        if (i%N == i/N) res[i] = 1;
        else res[i] = 0;
    }
    return res;
}

float max_sum(const float* A, int N){
    float temp_sum = 0, max_sum = 0;
    for (int i = 0; i<N*N; i++){
        temp_sum+=A[i];
        if (i%N == N-1){
            if (temp_sum>max_sum){
                max_sum = temp_sum;
            }
            temp_sum = 0;
        }
    }
    return max_sum;
}


float* mk_B(float* A, float* I, int N, float A1, float Ainf) {
    float* B = new float [N*N];
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N, N, N, (float)(1 / A1*Ainf), A, N, I, N, 0, B, N); //CblasRowMajor-матрицы хранятся по строкам, А транспонирована, B не транспонирована, N, N, N размеры матриц, коэффицент альфа C= α A*B + βC, s-float
    return B;
}

float* mk_inverted(float * I,float* R, float* B, int N, int M) {
    float* prevRDeg;
    float* curRDeg;
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, R, N, I, N, 1, I, N);//I = sum(I, R, N);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, R, N, R, N, 0, curRDeg, N);//curRDeg = multi_SSE(R, R, N);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, curRDeg, N, I, N, 1, I, N);//I = sum(I, curRDeg, N);
    prevRDeg = curRDeg;
    for (int i = 2; i < M; ++i) {
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, R, N, prevRDeg, N, 0, curRDeg, N);//curRDeg = multi_SSE(R, prevRDeg, N);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, curRDeg, N, I, N, 1, I, N);//I = sum(I, curRDeg, N);
        prevRDeg = curRDeg;
    }
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, I, N, B, N, 0, I, N);//multi_SSE(I, B, N);
    return I;
}

float* mk_R(float* I, float* B, float* A, int N) {
    float* R = mk_I(N);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, -1, B, N, A, N, 1, R, N); //R+=-BA=>R=-BA+I
    return R;
}

int main() {
    int N = 2048, M = 10;
    float* A = new float [N*N];
    for (int i = 0; i < N*N; ++i)
        A[i] = rand()%1000+1;
    float* AT = Transp(A,N);
    float A1 = max_sum(A,N);
    float Ainf = max_sum(AT,N);
    float* I = mk_I(N);
    float* B = mk_B(A,I,N,A1,Ainf);
    float* R = mk_R(I, B, A, N);
    float* AInv = mk_inverted(I, R, B, N, M);
    return 0;
}

