#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>
#include <iostream>

using namespace std;

#define LIST_COUNT 10
#define L 10000
#define TASKS_PER_PROCESS 1000
#define MIN_TASKS_TO_SHARE 10

double globalRes = 0;
double globalRes_sum = 0;
int procRank, procSize;
int *tasks;
int tasks_remaining;
double totalDisbalance = 0;

pthread_mutex_t mutex_tasks;
pthread_mutex_t mutex_tasks_remaining;
pthread_t thread_receiver;

void setWeight(int *tasks, int taskCount, int iterCounter) {
    pthread_mutex_lock(&mutex_tasks);
    for (int i = 0; i < taskCount; ++i) {
        tasks[i] = abs(50 - i % 100) * abs(procRank - (iterCounter % procSize)) * L;
    }
    pthread_mutex_unlock(&mutex_tasks);
}

void doTasks(int *tasks, int *tasks_executed) {
    pthread_mutex_lock(&mutex_tasks_remaining);

    for (int i = 0; tasks_remaining; ++i, --tasks_remaining) {
        pthread_mutex_unlock(&mutex_tasks_remaining);
        pthread_mutex_lock(&mutex_tasks);
        int weight = tasks[i];
        pthread_mutex_unlock(&mutex_tasks);

        for (int j = 0; j < weight; ++j) {
            globalRes += sin(j);
        }

        ++(*tasks_executed);
        pthread_mutex_lock(&mutex_tasks_remaining);
    }
    pthread_mutex_unlock(&mutex_tasks_remaining);
}

void *receiverThreadStart(void *args) {
    int tasks_to_share;
    int rankRequestedTasks;

    while (true) {
        MPI_Recv(&rankRequestedTasks, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //process rank that requests tasks

        if (rankRequestedTasks == procRank)
            break;

        pthread_mutex_lock(&mutex_tasks_remaining);
        if (tasks_remaining >= MIN_TASKS_TO_SHARE) {
            tasks_to_share = tasks_remaining / 2;
            tasks_remaining -= tasks_to_share;

            MPI_Send(&tasks_to_share, 1, MPI_INT, rankRequestedTasks, 1, MPI_COMM_WORLD);

            pthread_mutex_lock(&mutex_tasks);
            MPI_Send(tasks + tasks_remaining - 1, tasks_to_share, MPI_INT, rankRequestedTasks, 1, MPI_COMM_WORLD);
            pthread_mutex_unlock(&mutex_tasks);
        } else {
            tasks_to_share = 0;
            MPI_Send(&tasks_to_share, 1, MPI_INT, rankRequestedTasks, 1, MPI_COMM_WORLD);
        }
        pthread_mutex_unlock(&mutex_tasks_remaining);
    }
    return NULL;
}

void *workerThreadStart(void *args) {
    tasks = new int[TASKS_PER_PROCESS];

    double startTime, totalTime;
    double minTime, maxTime;

    for (int iterCounter = 0; iterCounter < LIST_COUNT; ++iterCounter) {
        setWeight(tasks, TASKS_PER_PROCESS, iterCounter);

        pthread_mutex_lock(&mutex_tasks_remaining);
        tasks_remaining = TASKS_PER_PROCESS;
        pthread_mutex_unlock(&mutex_tasks_remaining);
        int tasks_executed = 0;
        int tasks_received;

        MPI_Barrier(MPI_COMM_WORLD);
        startTime = MPI_Wtime();

        doTasks(tasks, &tasks_executed);

        for (int i = 0; i < procSize; ++i) { //requesting tasks
            if (i == procRank)
                continue;

            MPI_Send(&procRank, 1, MPI_INT, i, 0, MPI_COMM_WORLD); //waits for extra tasks
            MPI_Recv(&tasks_received, 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //number of extra tasks

            if (tasks_received > 0) {
                MPI_Recv(tasks, tasks_received, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //extra tasks

                pthread_mutex_lock(&mutex_tasks_remaining);
                tasks_remaining = tasks_received;
                pthread_mutex_unlock(&mutex_tasks_remaining);

                doTasks(tasks, &tasks_executed);
            }
        }

        totalTime = MPI_Wtime() - startTime;

        MPI_Allreduce(&totalTime, &minTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&totalTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        totalDisbalance += (maxTime - minTime) / maxTime;

        if (procRank == 0) {
            cout << endl << "Iteration " << iterCounter + 1 << endl;
            cout << "Disbalance time: " << maxTime - minTime << endl;
            cout << "Disbalance percentage: " << (maxTime - minTime) / maxTime*100 << "%" << endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        for (int i = 0; i < procSize; i++) {
            if (procRank == i) {
                cout << "Process " << procRank+1 << endl;
                cout << "Executed tasks: " << tasks_executed << endl;
                cout << "Result of iteration: " << globalRes << endl;
                cout << "Time of iteration: " << totalTime << endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    MPI_Send(&procRank, 1, MPI_INT, procRank, 0, MPI_COMM_WORLD);
    MPI_Allreduce(&globalRes, &globalRes_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    delete[] tasks;
    return NULL;
}

void createThreads() {
    pthread_attr_t attr; //object sets the attributes of thread
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE); //joinable

    pthread_create(&thread_receiver, &attr, receiverThreadStart, NULL); //create joinable thread fo receiving requests
    pthread_attr_destroy(&attr);

    workerThreadStart(NULL);

    pthread_join(thread_receiver, NULL); //wait for thread_receiver finish
}

int main(int argc, char **argv) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);

    pthread_mutex_init(&mutex_tasks, NULL);
    pthread_mutex_init(&mutex_tasks_remaining, NULL);

    double startTime = MPI_Wtime();
    createThreads();
    double totalTime = MPI_Wtime() - startTime;

    pthread_mutex_destroy(&mutex_tasks);
    pthread_mutex_destroy(&mutex_tasks_remaining);

    if (procRank == 0) {
        cout << endl << "Time taken: " << totalTime << "sec" << endl;
        cout << "Result of all iterations: " << globalRes_sum << endl;
        cout << "Total disbalance: " << totalDisbalance/LIST_COUNT*100 << "%" << endl;
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
