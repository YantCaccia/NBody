#define SOFTENING 1e-9f
#define DEBUG 0
#define DEFAULT_NBODIES 10000
#define DEFAULT_NITERS 5
#define DEFAULT_FILENAME "../results.txt"

typedef struct
{
    float x, y, z, vx, vy, vz;
} Body;

int compute_sendreceivecount_withoffset(int sendbufsize, int send_recv_count[], int offsets[], int size);
void debugPrint(const char *stringToPrint, int j, Body p[], int rank);
void randomizeBodies(Body p[], int n);
void bodiesPosition(Body p[], int nBodies, float dt);
void writeResults(float totalTime, Body p[], int j, int rank, const char *filename);