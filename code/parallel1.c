#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include "./mycollective/mycollective.h"

void bodiesForce(Body p[], int n, int myPartStart, int myPartNumberOfElements, float dt)
{
    for (int i = myPartStart; i < myPartStart + myPartNumberOfElements; i++)
    {
        float Fx = 0.0f;
        float Fy = 0.0f;
        float Fz = 0.0f;

        for (int j = 0; j < n; j++)
        {
            float dx = p[j].x - p[i].x;
            float dy = p[j].y - p[i].y;
            float dz = p[j].z - p[i].z;
            float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
            float invDist = 1.0f / sqrtf(distSqr);
            float invDist3 = invDist * invDist * invDist;

            Fx += dx * invDist3;
            Fy += dy * invDist3;
            Fz += dz * invDist3;
        }

        // Update velocity
        p[i].vx += dt * Fx;
        p[i].vy += dt * Fy;
        p[i].vz += dt * Fz;
    }
}

void computePositions(Body p[], int nBodies, int nIters, MPI_Datatype datatype, float dt, int offsets[], int sendrecvcount[], int rank)
{
    for (int iter = 0; iter < nIters; iter++)
    {

        // Bcast complete array to slaves
        MPI_Bcast(p, nBodies, datatype, 0, MPI_COMM_WORLD);

        // Debug
        debugPrint("After Bcast", sendrecvcount[rank], p, rank);

        // Compute forces, each slave just his part
        bodiesForce(p, nBodies, offsets[rank], sendrecvcount[rank], dt);

        // Gather results from slaves to master
        MPI_Gatherv(&p[offsets[rank]], sendrecvcount[rank], datatype, p, sendrecvcount, offsets, datatype, 0, MPI_COMM_WORLD);

        // Debug
        debugPrint("After bodiesForce", nBodies, p, rank);

        // If master
        if (rank == 0)
        {
            bodiesPosition(p, nBodies, dt);
        }
    }
}

int main(int argc, char **argv)
{

    // Init MPI
    MPI_Init(&argc, &argv);

    // Compute size and rank
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Take parameters from input
    int nBodies = argc > 1 ? atoi(argv[1]) : DEFAULT_NBODIES;
    int nIters = argc > 2 ? atoi(argv[2]) : DEFAULT_NITERS;
    const char *resultsFileName = argc > 3 ? argv[3] : DEFAULT_FILENAME;

    // Time step
    const float dt = 0.01f;

    // Compute sendcount and offsets
    int sendrecvcount[size];
    int offsets[size];
    compute_sendreceivecount_withoffset(nBodies, sendrecvcount, offsets, size);

    // Create new MPI_Datatype for Body
    MPI_Datatype MPI_Body;
    MPI_Type_contiguous(6, MPI_FLOAT, &MPI_Body);
    MPI_Type_commit(&MPI_Body);

    // Bodies buffer
    Body p[nBodies];

    // Start counting time
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    // Init bodies pos & vel
    if (rank == 0)
    {
        randomizeBodies(p, nBodies);
    }

    // Fine 1a parte in comune

    // Debug
    debugPrint("After init", nBodies, p, rank);

    // Main computational part
    computePositions(p, nBodies, nIters, MPI_Body, dt, offsets, sendrecvcount, rank);

    // Inizio 2a parte in comune

    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();

    double totalTime = (end - start);

    debugPrint("After all", nBodies, p, rank);

    // Write results in a file
    writeResults(totalTime, p, nBodies, rank, resultsFileName);

    // Finalize
    MPI_Type_free(&MPI_Body);
    MPI_Finalize();
}