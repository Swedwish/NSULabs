#include <iostream>
#include <math.h>
#include <cstdio>
#include <stdlib.h>
#include <cstring>
#include <wchar.h>
#include "mpi.h"
#include <cassert>

char isCellAlive(char* arr, int x, int y, int N, int M) {
    int sum = 0;
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            if (!(i == 0 && j == 0))
                sum += arr[(x + i + M) % M + (y + j) * M];
        }
    }

    if (arr[y * M + x] == 0) {
        if (sum == 3)
            return 1;
        else
            return 0;
    }
    else if (arr[y * M + x] == 1) {
        if (sum > 3 || sum < 2) {
            return 0;
        }
        else
            return 1;
    }
    else {
        printf("Error in isCellAlive.");
        MPI_Finalize();
        exit(6);
    }
}

void print_mat(int N, int M, char* A) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cout << rank << ":\n";
    for (int i = 0; i < N * M; i++) {
        printf("%d ", A[i]);
        if ((i + 1) % M == 0) std::cout << "\n";
    }
}

void liveCycle(char* arr, char* result, int N, int M) {
    for (int i = 1; i < N - 1; i++) {
        for (int j = 0; j < M; j++) {
            result[(i - 1) * M + j] = isCellAlive(arr, j, i,N, M);
            //print_mat(5, 12, result); MPI_Finalize(); exit(0);
        }
    }
}




int check(int count, char** history, int N, int M, int numthreads, int rank) {
    char flag = 0, buf = 0, exit = -1;
    for (int i = 0; i < count - 1; i++) {
        flag = 1;
        for (int j = 0; j < N * M; j++) {
            if (history[count][j] != history[i][j]) {
                flag = -1;
                break;
            }
        }
        if (flag != -1) {
            break;
        }
    }
    MPI_Allreduce(&flag, &buf, 1, MPI_CHAR, MPI_SUM, MPI_COMM_WORLD);
    if (buf == numthreads) {
        exit = 0;
    }
    return exit;
}

void setGlider(char* field, int N, int M) {
    field[1] = field[M + 2] = field[2 * M] = field[2 * M + 1] = field[2 * M + 2] = 1;
}


int main(int argc, char ** argv) {
    int N = 128; 
    int M = 112; 

    

    int rank, numtasks;
    char* field = new char[N * M];
    for (int i = 0; i < N * M; i++) {
        field[i] = 0;
    }
    setGlider(field, N, M);
    MPI_Init(&argc, &argv); 
    double start = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    char* oldFieldPart = new char [((N / numtasks + 2) * M)];
    char* newFieldPart = new char [((N / numtasks + 2) * M)];
    int counter = 0;
    char** history = new char* [100000];
    memcpy(&oldFieldPart[M], &field[(N / numtasks) * rank * M], (N / numtasks) * M * sizeof(char));
    history[0] = new char[N / numtasks * M];
   
    memcpy(history[0], &oldFieldPart[M], N / numtasks * M * sizeof(char));

    //print_mat(2, 280, oldFieldPart); MPI_Finalize(); exit(0);


    MPI_Request req1, req2, req3, req4;
    N = N / numtasks;
    
    while (check(counter, history, N, M, numtasks, rank) == -1) {
        if (numtasks != 1) {
            MPI_Irecv(oldFieldPart, M, MPI_CHAR, (rank + numtasks - 1) % numtasks, MPI_ANY_TAG, MPI_COMM_WORLD, &req1);
            MPI_Irecv(&oldFieldPart[(N + 1) * M], M, MPI_CHAR, (rank + 1) % numtasks, MPI_ANY_TAG, MPI_COMM_WORLD, &req2);
            MPI_Isend(&oldFieldPart[M], M, MPI_CHAR, (rank + numtasks - 1) % numtasks, 1, MPI_COMM_WORLD, &req3);
            MPI_Isend(&oldFieldPart[N * M], M, MPI_CHAR, (rank + 1) % numtasks, 1, MPI_COMM_WORLD, &req4);
        }
        else {
            memcpy(oldFieldPart, &oldFieldPart[N*M], M * sizeof(char));
            memcpy(&oldFieldPart[(N + 1) * M], &oldFieldPart[M], M * sizeof(char));
        }

        liveCycle(&oldFieldPart[M], &newFieldPart[2 * M], N, M);
         
        
        if (numtasks != 1) {
            MPI_Wait(&req1, MPI_STATUSES_IGNORE);
            MPI_Wait(&req2, MPI_STATUSES_IGNORE);



            liveCycle(oldFieldPart, &newFieldPart[M], 3, M);

            liveCycle(&oldFieldPart[(N - 1) * M], &newFieldPart[(N)*M], 3, M);

            //print_mat(3, 8, newFieldPart); MPI_Finalize(); exit(0);

            MPI_Wait(&req3, MPI_STATUSES_IGNORE);
            MPI_Wait(&req4, MPI_STATUSES_IGNORE);
        }
        counter++;
        history[counter] = new char[N * M];
        memcpy(history[counter], &newFieldPart[M], N * M * sizeof(char));
        memcpy(&oldFieldPart[M],&newFieldPart[M], N * M * sizeof(char));
        
        //if (counter == 1200) { print_mat(3, 8, history[32]); MPI_Finalize(); exit(0); }
        //if (counter == 2) { print_mat(1, 8, history[1]); MPI_Finalize(); exit(0); }
        //delete[](oldFieldPart);
        //oldFieldPart = newFieldPart;
    }

    double end = MPI_Wtime();

    if (rank == 0) {
        printf("time:%f, iterations:%d", end - start, counter);
    }
    
    MPI_Finalize();
    return 0;
}
