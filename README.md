# Introduzione
Nel problema N-Body, ci troviamo di fronte alla sfida di determinare le posizioni e le velocità di una serie di particelle interagenti nel corso del tempo. Immaginiamo, ad esempio, un astrofisico che desideri conoscere le posizioni e le velocità di un gruppo di stelle, oppure un chimico interessato alle posizioni e alle velocità di una collezione di molecole o atomi.

Per affrontare questo compito ci si serve di un risolutore, un programma in grado di trovare la soluzione al problema simulando il comportamento delle particelle. L'input richiesto comprende la massa, la posizione e la velocità di ciascuna particella all'inizio della simulazione, mentre l'output tipicamente fornisce le posizioni e le velocità di ogni particella in una sequenza di tempi o al termine di un determinato intervallo di tempo specificati dall'utente.

Ai fini del progetto si è scelto di implementare il risolutore servendosi della soluzione banale quadratica rispetto al numero di particelle. Tuttavia si noti che esistono anche altre soluzioni più efficienti, ad esempio quella che utilizza l'algoritmo di Barnes-Hut.

Il risolutore è stato testato sulla piattaforma Google Cloud su un cluster di 4 istanze e2-standard-8. Ciascuna istanza offre 32GB di RAM e 8vCPU; di queste, solo 4 sono realmente core fisici, per un totale di 16 core effettivi per l'intero cluster.

# Le soluzioni proposte
Il risolutore è stato implementato 2 volte.
Entrambe le versioni parallelizzano la versione sequenziale del risolutore fornitoci contestualmente all'assegnazione del progetto.
La differenza tra le 2 soluzioni risiede nelle modalità di comunicazione e sincronizzazione tra processi.
 
## Soluzione 1
In questa soluzione si è seguito l'approccio più semplice e lineare possibile.

Ad ogni iterazione, l'intero array contente le informazioni su tutti i bodies del problema viene mandato in broadcast a tutti i processi che partecipano alla computazione della soluzione.

Ciascuno di essi è in grado di conoscere quale è la parte di questo array che deve computare. Per ogni body appartenente alla propria parte di array, il processo aggiorna le velocities.

Le informazioni così calcolate vengono raccolte dal processo master, che infine aggiorna la posizione di ciascun body.

## Soluzione 2
La Soluzione 2 è più creativa.

Inizialmente il processo master suddivide l'array di body in tante parti quanti sono i processi che partecipano alla computazione. Ad ogni processo viene inviata una parte di array.

Ad ogni iterazione, ciascun processo:
* invia i propri bodies in broadcast a tutti gli altri processi;
* calcola le velocities dei propri bodies considerando solo le forze che intercorrono tra di essi;
* riceve da ogni processo n i bodies del processo n;
* ad ogni ricezione dal processo n, aggiorna le velocities dei propri bodies considerando solo le forze che intercorrono tra essi e i bodies del processo n;
* completate tutte le ricezioni e le relative computazioni delle velocities, aggiorna le posizioni dei propri bodies.

Al termine delle n iterazioni (numero deciso dall'utente), il processo master raccoglie tutti i bodies per la stampa dei risultati.

# Dettagli implementativi
L'implementazione delle due soluzioni è stata scritta in C, servendosi della libreria OpenMPI per parallelizzare le operazioni. Essa è un'implementazione open source di MPI sviluppata e mantenuta da un consorzio di partner di ricerca e industriali. MPI è uno standard di message-passing standardizzato e portabile, progettato per funzionare su architetture di calcolo parallelo. Esso definisce la sintassi e la semantica di librerie utilizzate per scrivere programmi di message-passing portabili in C e altri linguaggi di programmazione.

Alcuni snippets / funzioni vengono utilizzati da entrambe le soluzioni. Vediamole:

## Body structure
```
typedef struct
{
    float x, y, z, vx, vy, vz;
} Body;
```
È la struttura che rappresenta il singolo Body.
Ciascun Body è formato dalle *posizioni* (.x, .y, .z) e dalle *velocities* (.vx, .vy, .vz).

## int compute_sendreceivecount_withoffset
```
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
```
Questa funzione calcola il numero di oggetti da computare per ciascun processo.
Il calcolo viene eseguito distribuendo in maniera equa lo stesso numero di oggetti su ciascun processo.
Se ciò non è possibile (cioè se la divisione nOggetti/nProcessi non ha resto 0) i rimanenti oggetti vengono distribuiti in maniera equa sui primi m processi, con m pari al resto della divisione di cui sopra.

Nel codice, *sendbufsize* è il numero totale di oggetti da processare e *size* è il numero di processi che partecipano al calcolo.
I due array *send_recv_count* e *offsets* vengono utilizzati per memorizzare l'output della funzione:
* send_recv_count[i] conterrà il numero di oggetti da processare dal processo con *rank* i;
* offsets[i] conterrà l'indice del primo oggetto da processare dal processo con *rank* i nell'array contenente tutti gli oggetti; viene utilizzato come parametro *displs* nelle chiamate a funzioni di comunicazione collettiva messe a disposizione da MPI.


## randomizeBodies
```
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
```
Questa funzione inizializza in maniera randomica n bodies all'interno dell'array p.

## bodiesPosition
```
void bodiesPosition(Body p[], int nBodies, float dt)
{
    for (int i = 0; i < nBodies; i++)
    { // integrate position
        p[i].x += p[i].vx * dt;
        p[i].y += p[i].vy * dt;
        p[i].z += p[i].vz * dt;
    }
}
```
Questa funzione aggiorna le posizioni di *nBodies* bodies all'interno dell'array p.
Viene chiamata al termine di ciascuna iterazione, ma:
* Nella Soluzione 1, viene chiamata dal processo master che aggiorna le posizioni di tutti i bodies del problema;
* Nella Soluzione 2, viene chiamata da ciascun processo che aggiorna le posizioni solo dei propri bodies.

## main
Gran parte della funzione main è uguale nelle due soluzioni.
```
// Init MPI
MPI_Init(&argc, &argv);

// Compute size and rank
int rank, size;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
```
Si inizializza OpenMPI e e si salvano il numero di processi che partecipano alla computazione (*size*) e il rank del processo locale (*rank*).

```
// Take parameters from input
int nBodies = argc > 1 ? atoi(argv[1]) : DEFAULT_NBODIES;
int nIters = argc > 2 ? atoi(argv[2]) : DEFAULT_NITERS;
const char *resultsFileName = argc > 3 ? argv[3] : DEFAULT_FILENAME;
```
Si legge dagli argomenti da riga di comando alcuni parametri utili per l'esecuzione:
* *nBodies*: numero di bodies in input al problema
* *nIters*: numero di iterazioni da computare; ciascuna iterazione calcola le velocities e le posizioni dei bodies ad un intervallo di 0.01s rispetto all'iterazione precedente;
* *resultsFileName*: il nome del file all'interno del quale scrivere i risultati finali dell'esecuzione.

```
// Compute sendcount and offsets
int sendrecvcount[size];
int offsets[size];
compute_sendreceivecount_withoffset(nBodies, sendrecvcount, offsets, size);
```
Si calcola per ciascun processo il numero di bodies da computare e l'indice di partenza nell'array di bodies completo.

```
// Create new MPI_Datatype for Body
MPI_Datatype MPI_Body;
MPI_Type_contiguous(6, MPI_FLOAT, &MPI_Body);
MPI_Type_commit(&MPI_Body);
```
Si crea una nuovo MPI_Datatype che rappresenta la struttura Body. La scelta di creare un MPI_Datatype ha permesso di utilizzare le funzioni di comunicazione collettiva di MPI con maggior semplicità logica e sintattica.

```
// Init bodies pos & vel
if (rank == 0)
{
    randomizeBodies(p, nBodies);
}
```
Si inizializza l'array di Body *p* dal processo master.

## computePositions 
La grande differenza implementativa tra le due soluzioni risiede nella funzione *computePositions*.

### Soluzione 1
```
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
```
In questa versione di *computePositions* ad ogni iterazione viene eseguito il Broadcast dal master a tutti gli altri processi di tutto l'array completo di bodies.
Ogni processo chiama poi la *bodiesForce*, che aggiorna le velocities *solo dei suoi bodies* ma avendo a disposizione *tutti* i bodies.
I bodies così aggiornati vengono ripresi attraverso la Gather dal processo master.
Esso poi aggiorna le *posizioni* di tutti i bodies. 

### Soluzione 2
```
void computePositions(Body p[], Body myBodies[], int nLocalBodies, int nIters, MPI_Datatype datatype, float dt, int offsets[], int sendrecvcount[], int rank, int size)
{
    for (int iter = 0; iter < nIters; iter++)
    {

        MPI_Request requests[size];

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

        // Compute forces, each slave just his part
        bodiesForce(myBodies, nLocalBodies, myBodies, nLocalBodies, dt);
        debugPrint("After local bodiesForce", nLocalBodies, myBodies, rank);

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
```
In questa versione di *computePositions* ciascun processo effettua la Broadcast in modalità *non bloccante* verso gli altri processi, inviando il proprio array di Bodies. Questo perchè ciascun processo ha a disposizione solo la propria parte di bodies e non è a conoscenza dell'array completo di bodies (riempito e valido solo nel processo master). Inoltre, ciascun processo aspetta di ricevere, sempre in modalità non bloccante, la corrsipondente Broadcast dagli altri processi.
Nel frattempo, calcola le velocities dei suoi bodies tenendo in considerazione solo le forze che intercorrono tra di essi. In seguito, ad ogni ricezione dei bodies dagli altri processi, aggiorna le velocities dei propri bodies considerando le forze che intercorrono tra essi e i bodies appena ricevurti.
Infine, completate tutte le ricezioni e le relative computazioni delle velocities, aggiorna le posizioni dei propri bodies.

# Execution instructions

# Correctness discussion

# Benchmarks with strong and weak scalability.
