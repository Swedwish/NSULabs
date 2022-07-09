#include <iostream>
#include <math.h>
#include <cstdio>
#include <stdlib.h>
#include <cstring>
#include <wchar.h>
#include "mpi.h"
#include <cassert>
#include <iostream>

double* createMatrix(int d1, int d2) {
    double* res = new double[d1 * d2];
    for (int i = 0; i < d1 * d2; i++) {
        res[i] = rand() % 100 - 50;
    }
    return res;
}

void print_mat(int N, int M, double* A) {
    for (int i = 0; i < N * M; i++) {
        std::cout << A[i] << " ";
        if ((i + 1) % M == 0) std::cout << "\n";
    }
}

int main(int argc, char** argv){
    int n1 = 3600, n2 = 3600, n3 = 3600, p1 = 2, p2 = 2, dims[2] = { p2,p1 }, periods[2] = { 0,0 }, reorder = 1;
    MPI_Datatype col_type, col;
    double* A = new double[n1 * n2], * B = new double[n3 * n2], * C = new double[n1 * n3];
    MPI_Comm comm2d;
    double start, end;
    int rank, numtasks, coords[2];
    const char* a = nullptr;
    //print_mat(n1, n2, B);
    MPI_Init(&argc,&argv);

    start = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    //MPI_Dims_create(numtasks, 2, dims);
    if (numtasks != 8) {              //проверка на кол-во процессов
        printf("numtasks!=18");
        MPI_Finalize();
        exit(0);
    }

    p1 = dims[0] = 2; p2 = dims[1] = 4; //изменение размеров решетки
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm2d);
    
    MPI_Comm_rank(comm2d, &rank);
    
    //MPI_Cart_get(comm2d, 2, dims, periods, coords);
    MPI_Type_vector(n2, n3 / p2, n3, MPI_DOUBLE, &col);
    MPI_Type_create_resized(col, 0, n3/p2 * sizeof(double), &col_type);
    MPI_Type_commit(&col_type);
    int prevy, prevx, nexty, nextx;
    double* rows = new double[n2 * (n1 / p1)];
    double* cols = new double[n2 * (n3 / p2)]; 
    if (rank == 0) {
        delete A, B;
        A = createMatrix(n1, n2);
        B = createMatrix(n2, n3);
        //print_mat(n2, n1, A);
    }


    MPI_Comm ROW, COL;
    int remain_dims[2] = { 0, 1 };
    MPI_Cart_sub(comm2d, remain_dims, &ROW);
    remain_dims[0] = 1;
    remain_dims[1] = 0;
    MPI_Cart_sub(comm2d, remain_dims, &COL);
    MPI_Cart_shift(comm2d, 0, 1, &prevx, &nextx);
    MPI_Cart_shift(comm2d, 1, 1, &prevy, &nexty);

    if (coords[1] == 0) {
        MPI_Scatter(B, 1, col_type, cols, n2 * (n3 / p2), MPI_DOUBLE, 0, COL); 
    }
    if (coords[0] == 0) {
        MPI_Scatter(A, n2 * (n1 / p1), MPI_DOUBLE, rows, n2 * (n1 / p1), MPI_DOUBLE, 0, ROW);
    } 



    MPI_Bcast(rows, n2 * (n1 / p1), MPI_DOUBLE, 0, COL);//Раздаем строки A остальным процессам
    MPI_Bcast(cols, n2 * (n3 / p2), MPI_DOUBLE, 0, ROW);//Раздаем столбцы B остальным процессам



    double* res = new double[(n1 * n3)/(p1*p2)];
    double* result = new double[(n1 * n3)];
    double sum;
    for (int i = 0; i < n1 / p1; i++) {
        for (int j = 0; j < n3 / p2; j++) {
            sum = 0;
            for (int k = 0; k < n2; k++) {
                sum += rows[n2 * i + k] * cols[(n3 / p2) * k + j];
            }
            res[i * (n3 / p2) + j] = sum;
        }
    }

    MPI_Gather(res, (n1 * n3) / (p1 * p2), MPI_DOUBLE, result, (n1 * n3) / (p1 * p2), MPI_DOUBLE, 0, comm2d);

    MPI_Type_free(&col_type);
    if (rank == 0) {
        double* realResult = new double[(n1 * n3)];
        for (int i = 0; i < p2; i++) {
            for (int j = 0; j < p1; j++) {
                for (int k = 0; k < (n1 * n3) / (p1 * p2); k++) {
                    realResult[k / (n3 / p2) * n3 + i * (n3 / p2) + j * n3 * (n1 / p1) + k % (n3 / p2)] = result[((n1 * n3) / (p1 * p2)) * (j + i * p1) + k];
                    //realResult[Строчка в блоке + оступ слева к началу+ отступ сверху к
                    //началу+внутри строки] =
                    //result[Сдвиг к началу блока + k]
                }
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    if (rank == 0) {
        printf("%f", end - start);
    }
    MPI_Finalize();
    return 0;
}