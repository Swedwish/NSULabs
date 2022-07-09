#include <iostream>
#include <stdlib.h>
#include <ctime>

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

void heat(int* arr, int N, int K) {
    int k = 0;
    for (int i = 0; i < N * K; ++i) {
        k = arr[k];
    }
    if (k == 1234) std::cout << "Wow!";
}


unsigned long long time(int* arr, int N, int K) {
    heat(arr, N, K);
    int k = 0;
    asm("rdtsc\n":"=a"(start.t32.th),"=d"(start.t32.tl));
    for(int i = 0; i<N*K; ++i) {
        k = arr[k];
    }
    asm("rdtsc\n":"=a"(end.t32.th),"=d"(end.t32.tl));
    if (k == 1234) std::cout << "Wow!";
    return (end.t64-start.t64)/(N*K);
}

void directFill(int* arr, int N) {
    for (int i = 0; i < N-1; ++i) {
        arr[i] = i+1;
    }
    arr[N-1] = 0;
}

void reverseFill(int* arr, int N) {
    for (int i = N-1; i > 0; --i) {
        arr[i] = i-1;
    }
    arr[0] = N-1;
}

void randomFill(int* arr, int N) {
    srand(time(NULL));
    int* fill = new int[N];
    for (int i = 0; i < N-1; ++i) {
        fill[i] = i + 1;
    }
    int unused = 0;
    int unusedTotal = N-1;
    for(int i = 0; i < N-1; ++i)  {
        int index = rand() % (unusedTotal-1);
        arr[unused] = fill[index];
        unused = arr[unused];
        std::swap(fill[index], fill[unusedTotal-1]);
        unusedTotal--;
    }
    delete[] fill;
    arr[unused] = 0;
}



void calculate (int Nmin, int Nmax, int K) {
    for (int i = Nmin; i <= Nmax; i += (i/10)) {
        unsigned long long minTime;
        int* direct = new int [i];
        int* reverse = new int [i];
        int* random = new int [i];
        double mem = (double)i*4/1024;
        std::cout << mem << "; ";
        directFill(direct, i);
        minTime = time(direct, i, K);
        for (int j = 0; j < 10; ++j) {
            minTime = std::min(minTime, time(direct, i, K));
        }
        std::cout << minTime<< "; ";
        reverseFill(reverse, i);
        minTime = time(reverse, i, K);
        for (int j = 0; j < 10; ++j) {
            minTime = std::min(minTime, time(reverse, i, K));
        }
        std::cout << minTime<< "; ";
        randomFill(random, i);
        minTime = time(random, i, K);
        for (int j = 0; j < 10; ++j) {
            minTime = std::min(minTime, time(random, i, K));
        }
        std::cout << minTime<< std::endl;
        delete[] direct;
        delete[] reverse;
        delete[] random;
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
    int Nmin = 256;
    int Nmax = 16 * 1024 * 256;
    int K = 10;
    calculate(Nmin, Nmax, K);
    return 0;
}
