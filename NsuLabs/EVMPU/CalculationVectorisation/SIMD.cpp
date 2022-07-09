#include <iostream>
#include <xmmintrin.h>
#include <stdlib.h>
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

float* multi_SSE(float* matrix_1, float* matrix_2, int N)
{
    __m128 vector_2;//переменная для хранения значения второй матрицы
    __m128 result_vector;//переменная для записи итоговой марицы
    float* result_matrix = (float*)_mm_malloc(N * N * sizeof(float), 16);//выделение памяти с выравниванием
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            __m128 element_vector = _mm_set1_ps(matrix_1[i * N + j]);//4 позиции в одно значение
            for (int k = 0; k < N; k += 4) {
                vector_2 = _mm_load_ps(&matrix_2[j * N + k]);//4 значение по адресу
                result_vector = _mm_load_ps(&result_matrix[i * N + k]);
                result_vector = _mm_add_ps(result_vector, _mm_mul_ps(element_vector, vector_2));
                _mm_store_ps(&result_matrix[i * N + k], result_vector);
            }
        }
    }
    return result_matrix;
}


float* sum(const float* A, const float* B, int N, int z = 1) {
    float* res = new float [N*N];
    z/=abs(z);
    for (int i = 0; i < N*N; ++i)
        res[i] =A[i] + z*B[i];
    return res;
}

float* mk_B(float A1,float Ainf, float * AT, int N){
    float* B = new float [N*N];
    for (int i = 0; i < N*N; ++i)
        B[i] = AT[i] / A1*Ainf;
    return B;
}

float* mk_R(float* I, float* B, float* A, int N){
    return sum(I, multi_SSE(A,B,N),N,-1);
}

float* mk_inverted(float * I,float* R, float* B, int N, int M) {
    float* prevRDeg;
    float* curRDeg;
    I = sum(I, R, N);
    curRDeg = multi_SSE(R, R, N);
    I = sum(I, curRDeg, N);
    prevRDeg = curRDeg;
    for (int i = 2; i < M; ++i) {
        curRDeg = multi_SSE(R, prevRDeg, N);
        I = sum(I, curRDeg, N);
        prevRDeg = curRDeg;
    }
    return multi_SSE(I, B, N);
}

int main() {
    int N = 2048, M = 10;
    float* A = new float [N*N];
    for (int i = 0; i < N*N; ++i)
        A[i] = rand()%1000+1;
    float* AT = Transp(A,N);
    float A1 = max_sum(A,N);
    float Ainf = max_sum(AT,N);
    float* B = mk_B(A1,Ainf,AT,N);
    float* I = mk_I(N);
    float* R = mk_R(I, B, A, N);
    float* AInv = mk_inverted(I, R, B, N, M);
    return 0;
}
