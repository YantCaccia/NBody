#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mycollective.h"

int compute_sendreceivecount_withoffset(int sendbufsize, int send_recv_count[], int offsets[], int size)
{

    int elementsPerProcess = sendbufsize / (size);
    int remainder = sendbufsize % (size);

    // Inizializzo array
    send_recv_count[0] = elementsPerProcess;
    offsets[0] = 0;
    for (int i = 1; i < size; i++)
    {
        send_recv_count[i] = elementsPerProcess;
        if (remainder > 0)
        {
            send_recv_count[i]++;
            remainder--;
        }
        offsets[i] = offsets[i - 1] + send_recv_count[i - 1];
    }
}

void debugPrint(const char *stringToPrint, int j, Body p[], int rank)
{
    if (rank == 0 && DEBUG)
    {
        printf("%s:\n", stringToPrint);
        for (int i = 0; i < j; i++)
        {
            printf("p[%d].x: %f\tp[%d].vx: %f\n", i, p[i].x, i, p[i].vx);
        }
        printf("\n");
    }
}

void writeResults(float totalTime, Body p[], int j, int rank, const char *filename)
{
    if (rank == 0)
    {
        FILE *fp;
        fp = fopen(filename, "w");
        fprintf(fp, "Total time: %fs\n\n", totalTime);
        for (int i = 0; i < j; i++)
        {
            fprintf(fp, "p[%d].x: %5.3f\tp[%d].y: %5.3f\tp[%d].z: %5.3f\np[%d].vx: %5.3f\tp[%d].vy: %5.3f\tp[%d].vz: %5.3f\n\n", i, p[i].x, i, p[i].y, i, p[i].z, i, p[i].vx, i, p[i].vy, i, p[i].vz);
        }
        fclose(fp);
    }
}

void randomizeBodies(Body p[], int n)
{
    for (int i = 0; i < n; i++)
    {
        p[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p[i].vx = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p[i].vy = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p[i].vz = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    }
}

void bodiesPosition(Body p[], int nBodies, float dt)
{
    for (int i = 0; i < nBodies; i++)
    { // integrate position
        p[i].x += p[i].vx * dt;
        p[i].y += p[i].vy * dt;
        p[i].z += p[i].vz * dt;
    }
}