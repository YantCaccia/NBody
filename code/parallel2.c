#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include "./mycollective/mycollective.h"

void bodiesForce(Body p[], int nBodies, Body myBodies[], int nLocalBodies, float dt)
{
    for (int i = 0; i < nLocalBodies; i++)
    {
        float Fx = 0.0f;
        float Fy = 0.0f;
        float Fz = 0.0f;

        for (int j = 0; j < nBodies; j++)
        {
            float dx = p[j].x - myBodies[i].x;
            float dy = p[j].y - myBodies[i].y;
            float dz = p[j].z - myBodies[i].z;
            float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
            float invDist = 1.0f / sqrtf(distSqr);
            float invDist3 = invDist * invDist * invDist;

            Fx += dx * invDist3;
            Fy += dy * invDist3;
            Fz += dz * invDist3;
        }

        // Update velocity
        myBodies[i].vx += dt * Fx;
        myBodies[i].vy += dt * Fy;
        myBodies[i].vz += dt * Fz;
    }
}

void computePositions(Body p[], Body myBodies[], int nLocalBodies, int nIters, MPI_Datatype datatype, float dt, int offsets[], int sendrecvcount[], int rank, int size)
{
    for (int iter = 0; iter < nIters; iter++)
    {

        MPI_Request requests[size];

        // Compute forces, each slave just his part
        bodiesForce(myBodies, nLocalBodies, myBodies, nLocalBodies, dt);
        debugPrint("After local bodiesForce", nLocalBodies, myBodies, rank);

        // Non-blocking bcast to others slaves +
        // Receive from other slaves
        for (int i = 0; i < size; i++)
        {
            void *buffer;
            int count;

            // I love one-liners
            (i == rank) ? (buffer = myBodies, count = nLocalBodies) : (buffer = &p[offsets[i]], count = sendrecvcount[i]);

            // 1 will be to send myBodies, size - 1 will be to receive other processes' myBodies
            MPI_Ibcast(buffer, count, datatype, i, MPI_COMM_WORLD, &requests[i]);
        }

        // Compute forces with local bodies from other processes as they arrive
        for (int i = 0; i < size; i++)
        {

            // It will be equal to the rank of the process for which the req is completed
            int index;
            MPI_Status status;
            MPI_Waitany(size, requests, &index, &status);

            if (index != rank)
            {
                // Compute forces on myBodies with forces by other processes' myBodies
                bodiesForce(&p[offsets[index]], sendrecvcount[index], myBodies, nLocalBodies, dt);
            }
        }

        // Compute myBodies position
        bodiesPosition(myBodies, nLocalBodies, dt);
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

    // Bodies buffers
    Body p[nBodies];                    // All bodies
    Body myBodies[sendrecvcount[rank]]; // Local bodies

    // Start counting time
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    // Init bodies pos & vel
    if (rank == 0)
    {
        randomizeBodies(p, nBodies);
    }

    // From master to slaves, each slave just his part
    MPI_Scatterv(p, sendrecvcount, offsets, MPI_Body, myBodies, sendrecvcount[rank], MPI_Body, 0, MPI_COMM_WORLD);

    // Debug
    debugPrint("After init", nBodies, p, rank);

    // Main computational part
    computePositions(p, myBodies, sendrecvcount[rank], nIters, MPI_Body, dt, offsets, sendrecvcount, rank, size);

    // Results to the master
    MPI_Gatherv(myBodies, sendrecvcount[rank], MPI_Body, p, sendrecvcount, offsets, MPI_Body, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();

    double totalTime = (end - start);

    // Debug
    debugPrint("After all", nBodies, p, rank);

    // Write results in a file
    writeResults(totalTime, p, nBodies, rank, resultsFileName);

    // Finalize
    MPI_Type_free(&MPI_Body);
    MPI_Finalize();

}