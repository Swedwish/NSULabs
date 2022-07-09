#include <iostream>
#include <stdlib.h>
#include <ctime>

#define offset 16*1024*1024/4 //16 МБ/sizeof(int)
#define size (3*1024+256+32)*1024/4 //cash size


void* multi(float* A, float* B,float* res, int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            res[i*N + j] = 0;
            for (int k = 0; k < N; ++k)
                res[i*N + j] += A[i*N + k] * B[k*N + j];
        }
    }
}


union ticks{
    unsigned long long t64;
    struct s32 { long th, tl; } t32;
} start, end;


unsigned long long time(int* arr, int N) {
    int k = 0;
    asm("rdtsc\n":"=a"(start.t32.th),"=d"(start.t32.tl));
    for(int i = 0; i<size*N; ++i) {
        k = arr[k];
    }
    asm("rdtsc\n":"=a"(end.t32.th),"=d"(end.t32.tl));
    if (k == 1234) std::cout << "Wow!";
    return (end.t64-start.t64)/(N*size);
}

void fill(int* arr, int N) {
    for (int i = 0; i < N - 1; ++i) {
        for (int j = 0; j < size/N; ++j) {
            arr[offset * i + j] = (i + 1) * offset + j;
        }
    }
    for (int i = 0; i < size/N; ++i) {
        arr[offset * (N - 1) + i] = i + 1;
    }
    arr[offset * (N - 1) + size/N - 1] = 0;
}


void calculate (int* arr) {
    for(int i = 1; i <= 32; ++i) {
        fill(arr, i);
        unsigned long long minTime = time(arr, i);
        for (int j = 0; j < 5; ++j) {
            minTime = std::min(minTime, time(arr, i));
        }
        std::cout << i << " " << minTime << std::endl;
    }
}

int main() {
    int N = 500;
    float* A = new float [N*N];
    float* B = new float [N*N];
    float* res = new float [N*N];
    multi(A, B, res, N);
    if(res[0] == 100)  std::cout << "Wow";
    delete[] A;
    delete[] B;
    delete[] res;
    int* arr = new int [offset*32];
    calculate(arr);
    return 0;
}
