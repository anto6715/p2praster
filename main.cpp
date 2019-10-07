#include <iostream>
#include <igraph/igraph.h>
#include <cstring>
#include <random>
#include "raster.h"


void StartTheClock();
double StopTheClock();
void usage(char* cmd);
igraph_t generateGeometricGraph(igraph_integer_t n, igraph_real_t radius);
igraph_t generateBarabasiAlbertGraph(igraph_integer_t n, igraph_real_t power, igraph_integer_t m, igraph_real_t A);
igraph_t generateErdosRenyiGraph(igraph_integer_t n, igraph_erdos_renyi_t type, igraph_real_t param);
igraph_t generateRegularGraph(igraph_integer_t n, igraph_integer_t k);
igraph_t generateRandomGraph(int type, int n);
void printGraphType(int type);

using namespace std;
using namespace std::chrono;

struct Params {
    uint32_t    domainSize;
    long        ni;            // points number of dataset
    int         peers;
    int         p_star;
    string      outputFilename;
    double      convThreshold;
    int         convLimit;
    int         graphType;
    int         fanOut;
    int         roundsToExecute;
    double      delta;
    double      precision;
    int         threshold;
    int         min_size;
    int         radius;
    string      name_file;
};

typedef struct PeerStats {
    int numCandidates;
    int falseNegatives;
    int falsePositives;
    double recall;
    double precision;
    double f1score;
    double f2score;
    double fhalfscore;
    bool theoreticalErrorExceeded;
    double heavyHittersARE;
    double candidatesARE;
    double heavyHittersTAE;
    double candidatesTAE;
    double allItemsARE;
    double allItemsTAE;
} PeerStats;
high_resolution_clock::time_point t1, t2;

int main(int argc, char **argv) {
    /*** Default Parameters***/
    long        ni; // points number
    long        *peerLastItem; // index of a peer last item
    uint32_t    domainSize = 1048575; // number of possible distinct items
    int         peers = 2; // number of peers
    int         fanOut = 3; //fan-out of peers
    int         graphType = 2; // graph distribution: 1 geometric 2 Barabasi-Albert 3 Erdos-Renyi 4 regular (clique)
    double      convThreshold = 0.001; // local convergence tolerance
    int         convLimit = 3; // number of consecutive rounds in which a peer must locally converge
    int         roundsToExecute = -1;
    int         p_star = -1;
    double      delta = 0.04;
    double      precision = 0.0;
    int         threshold = 0;
    int         min_size = 0;
    int         radius = 0;
    string      name_file = "/home/antonio/Scrivania/Datasets/S-sets/s1.csv";

    double          elapsed;
    int             iterations;
    bool            autoseed = false;
    bool            outputOnFile = false;
    string          outputFilename;
    igraph_t        graph;
    Params          params;


    /*** parse command-line parameters ***/
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-di") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing domain size parameter." << endl;
                return -1;
            }
            domainSize = stol(argv[i]);
        } else if (strcmp(argv[i], "-delta") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing delta parameter." << endl;
                return -1;
            }
            delta = stod(argv[i]);
        } else if (strcmp(argv[i], "-p") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of peers parameter." << endl;
                return -1;
            }
            peers = stoi(argv[i]);
        } else if (strcmp(argv[i], "-ps") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing maximum number of peers parameter." << endl;
                return -1;
            }
            p_star = stoi(argv[i]);
        } else if (strcmp(argv[i], "-f") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing fan-out parameter." << endl;
                return -1;
            }
            fanOut = stol(argv[i]);;
        } else if (strcmp(argv[i], "-d") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing graph type parameter" << endl;
                return -1;
            }
            graphType = stoi(argv[i]);
        } else if (strcmp(argv[i], "-ct") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing convergence tolerance parameter." << endl;
                return -1;
            }
            convThreshold = stod(argv[i]);
        } else if (strcmp(argv[i], "-cl") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing # of consecutive rounds in which convergence is satisfied parameter." << endl;
                return -1;
            }
            convLimit = stol(argv[i]);
        } else if (strcmp(argv[i], "-of") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing filename for simulation output." << endl;
                return -1;
            }
            outputFilename = string(argv[i]);
        } else if (strcmp(argv[i], "-r") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of rounds to execute.\n";
                return -1;
            }
            roundsToExecute = stoi(argv[i]);
        } else if (strcmp(argv[i], "-pr") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of precision for raster.\n";
                return -1;
            }
            precision = stod(argv[i]);
        } else if (strcmp(argv[i], "-thr") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of threshold for raster.\n";
                return -1;
            }
            threshold = stoi(argv[i]);
        } else if (strcmp(argv[i], "-ms") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of min size for raster.\n";
                return -1;
            }
            min_size = stoi(argv[i]);
        } else if (strcmp(argv[i], "-dt") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing file name of dataset.\n";
                return -1;
            }
            name_file = argv[i];
        }else if (strcmp(argv[i], "-as") == 0) {
            autoseed = true;
        } else {
            usage(argv[0]);
            return -1;
        }
    }

    /*** read dataset dimensions ***/
    int row, column;
    getDim(name_file, row, column);
    ni = row;

    double* dataset_storage = new double[row*column* sizeof(double)];
    double** dataset = new double*[row * sizeof(double)];
    for (int i = 0; i < row; i++) {
        dataset[i] = &dataset_storage[i*column];
    }
    StartTheClock();
    loadData(dataset, name_file, column);
    double greaddatasettime = StopTheClock();
    if (!outputOnFile) {
        cout <<"Time (seconds) required to load the dataset: " << greaddatasettime << "\n";
    }

    /*** Insert data into an ordered structure  (for vertical partition of dataset) ***/
    set<array<double, 2>> data;
    set<array<double, 2>>::iterator it;
    for (int k = 0; k < row; ++k) {
        data.insert({dataset[k][0], dataset[k][1]});
    }

    /*** Compute last item for each peer***/
    peerLastItem = (long *) calloc(peers, sizeof(long));
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_real_distribution<> distr(-1, 1); // define the range

    for(int i = 0; i < peers - 1; i++){
        float rnd = distr(eng);
        //cerr << "rnd: " << rnd << "\n";
        long last_item = rnd * ((float)ni/(float)peers) * 0.1 + (float) (i+1) * ((float)ni/(float)peers) - 1;
        peerLastItem[i] = last_item;
    }

    peerLastItem[peers - 1] = ni-1;

    /*** check the partitioning correctness ***/
    long sum = peerLastItem[0] + 1;
    //cerr << "peer 0:" << sum << "\n";
    for(int i = 1; i < peers; i++) {
        sum += peerLastItem[i] - peerLastItem[i-1];
        //cerr << "peer " << i << ":" << peerLastItem[i] - peerLastItem[i-1] << "\n";
    }

    if(sum != ni) {
        fprintf(stdout, "ERROR: ni = %ld != sum = %ld\n", ni, sum);
        exit(EXIT_FAILURE);
    }

    /*** assign parameters read from command line ***/
    params.name_file = name_file;
    params.ni = ni;
    params.domainSize = domainSize;
    params.peers = peers;
    params.fanOut = fanOut;
    params.graphType = graphType;
    params.convThreshold = convThreshold;
    params.convLimit = convLimit;
    params.outputFilename = outputFilename;
    params.roundsToExecute = roundsToExecute;
    params.delta = delta;
    params.precision = precision;
    params.threshold = threshold;
    params.min_size = min_size;
    params.radius = radius;
    if (p_star == -1)
        p_star = peers;

    params.p_star = p_star;

    outputOnFile = params.outputFilename.size() > 0;

    if (!outputOnFile) {
        printf("\n\nPARAMETERS:\n");
        //printf("distinct items in the stream = %d\n", params.domainSize);
        cout << "dataset = " << params.name_file << "\n";
        cout << "n° points = " << params.ni << "\n";
        cout << "precision = " << params.precision << "\n";
        cout << "threshold = " << params.threshold << "\n";
        cout << "min size = " << params.min_size << "\n";
        cout << "radius = " << params.radius << "\n";
        cout << "peers = " << params.peers << "\n";
        cout << "fan-out = " << params.fanOut << "\n";
        cout << "graph type = ";
        printGraphType(params.graphType);
        cout << "local convergence tolerance = "<< params.convThreshold << "\n";
        cout << "number of consecutive rounds in which a peer must locally converge = "<< params.convLimit << "\n";
        cout << "\n\n";
    }

    /*** Graph generation ***/
    // turn on attribute handling in igraph
    igraph_i_set_attribute_table(&igraph_cattribute_table);

    // seed igraph PRNG
    igraph_rng_seed(igraph_rng_default(), 42);

    StartTheClock();
    // generate a connected random graph
    graph = generateRandomGraph(params.graphType, params.peers);

    double gengraphtime = StopTheClock();
    if (!outputOnFile) {
        cout <<"Time (seconds) required to generate the random graph: " << gengraphtime << "\n";
    }

    // determine minimum and maximum vertex degree for the graph
    igraph_vector_t result;
    igraph_real_t mindeg;
    igraph_real_t maxdeg;

    igraph_vector_init(&result, 0);
    igraph_degree(&graph, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    igraph_vector_minmax(&result, &mindeg, &maxdeg);
    if (!outputOnFile) {
        cout << "Minimum degree is "  <<  mindeg << ", Maximum degree is " << maxdeg << "\n";
    }

    // add first run of raster algorithm for each peer
    /*** Apply Raster to each peer's dataset ***/
    if (!outputOnFile) {
        printf("\nApply Raster to each peer's dataset...\n");
    }
    vector<unordered_map<array<int, 2>, int, container_hasher>> peers_projection;
    vector<vector<unordered_set<array<int, 2>, container_hasher>>> peers_clusters;
    vector<unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher>> peers_all_points;

    StartTheClock();
    int start = 0;
    for(int peerID = 0; peerID < params.peers; peerID++){
        cout << "peer: " << peerID << endl;
        unordered_map<array<int, 2>, int, container_hasher> projection;
        vector<unordered_set<array<int, 2>, container_hasher>> clusters;
        unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> all_points;


        /*** points agglomeration ***/

        // Raster version
        //mapToTiles(dataset, precision, threshold, projection, start, peerLastItem[peerID]);
        // Raster' version
        mapToTilesPrime(dataset, precision, threshold, projection, all_points, start, peerLastItem[peerID]);

        /*** generate clusters from tiles ***/
        clusteringTiles(projection, min_size, clusters);
        start = peerLastItem[peerID] + 1;

        // save result obtained for each peer for the next steps
        peers_projection.push_back(projection);
        peers_all_points.push_back(all_points);
        peers_clusters.push_back(clusters);
    }

    double rastertime = StopTheClock();
    if (!outputOnFile) {
        printf("Raster done!\n");
        printf("Time (seconds) required by Raster: %f\n", rastertime);
    }
}

void StartTheClock(){
    t1 = high_resolution_clock::now();
}

double StopTheClock() {
    t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    return time_span.count();
}

void usage(char* cmd)
{
    std::cerr
            << "Usage: " << cmd << "\n"
            << "  -ni   number of items\n"
            << "  -di   number of possible distinct items (domain size)\n"
            << "  -p    number of peers\n"
            << "  -f    fan-out of peers\n"
            << "  -d    graph type: 1 geometric 2 Barabasi-Albert 3 Erdos-Renyi 4 regular\n"
            << "  -ct   convergence tolerance\n"
            << "  -cl   number of consecutive rounds in which convergence must be satisfied\n"
            << "  -of   output filename, if specified a file with this name containing all of the peers stats is written\n"
            << "  -pr   number of precision for raster\n"
            << "  -thr  number of threshold for raster\n"
            << "  -ms   number of min size for raster\n"
            << "  -dt   file name containing dataset\n\n";
}

igraph_t generateGeometricGraph(igraph_integer_t n, igraph_real_t radius)
{
    igraph_t G_graph;
    igraph_bool_t connected;

    // generate a connected random graph using the geometric model
    igraph_grg_game(&G_graph, n, radius, 0, 0, 0);

    igraph_is_connected(&G_graph, &connected, IGRAPH_WEAK);
    while(!connected){
        igraph_destroy(&G_graph);
        igraph_grg_game(&G_graph, n, radius, 0, 0, 0);

        igraph_is_connected(&G_graph, &connected, IGRAPH_WEAK);
    }

    return G_graph;
}

igraph_t generateBarabasiAlbertGraph(igraph_integer_t n, igraph_real_t power, igraph_integer_t m, igraph_real_t A)
{

    // n = The number of vertices in the graph
    // power = Power of the preferential attachment. The probability that a vertex is cited is proportional to d^power+A, where d is its degree, power and A are given by arguments. In the classic preferential attachment model power=1
    // m = number of outgoing edges generated for each vertex
    // A = The probability that a vertex is cited is proportional to d^power+A, where d is its degree, power and A are given by arguments

    igraph_t BA_graph;
    igraph_bool_t connected;

    // generate a connected random graph using the Barabasi-Albert model
    igraph_barabasi_game(/* graph=    */ &BA_graph,
            /* n=        */ n,
            /* power=    */ power,
            /* m=        */ m,
            /* outseq=   */ 0,
            /* outpref=  */ 0,
            /* A=        */ A,
            /* directed= */ IGRAPH_UNDIRECTED,
            /* algo=     */ IGRAPH_BARABASI_PSUMTREE,
            /* start_from= */ 0);


    igraph_is_connected(&BA_graph, &connected, IGRAPH_WEAK);
    while(!connected){
        igraph_destroy(&BA_graph);
        igraph_barabasi_game(/* graph=    */ &BA_graph,
                /* n=        */ n,
                /* power=    */ power,
                /* m=        */ m,
                /* outseq=   */ 0,
                /* outpref=  */ 0,
                /* A=        */ A,
                /* directed= */ IGRAPH_UNDIRECTED,
                /* algo=     */ IGRAPH_BARABASI_PSUMTREE,
                /* start_from= */ 0);

        igraph_is_connected(&BA_graph, &connected, IGRAPH_WEAK);
    }

    return BA_graph;
}

igraph_t generateErdosRenyiGraph(igraph_integer_t n, igraph_erdos_renyi_t type, igraph_real_t param)
{
    // n = The number of vertices in the graph
    // type = IGRAPH_ERDOS_RENYI_GNM G(n,m) graph, m edges are selected uniformly randomly in a graph with n vertices.
    //      = IGRAPH_ERDOS_RENYI_GNP G(n,p) graph, every possible edge is included in the graph with probability p

    igraph_t ER_graph;
    igraph_bool_t connected;

    // generate a connected random graph using the Erdos-Renyi model
    igraph_erdos_renyi_game(&ER_graph, type, n, param, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    igraph_is_connected(&ER_graph, &connected, IGRAPH_WEAK);
    while(!connected){
        igraph_destroy(&ER_graph);
        igraph_erdos_renyi_game(&ER_graph, type, n, param, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

        igraph_is_connected(&ER_graph, &connected, IGRAPH_WEAK);
    }

    return ER_graph;
}




igraph_t generateRegularGraph(igraph_integer_t n, igraph_integer_t k)
{
    // n = The number of vertices in the graph
    // k = The degree of each vertex in an undirected graph. For undirected graphs, at least one of k and the number of vertices must be even.


    igraph_t R_graph;
    igraph_bool_t connected;

    // generate a connected regular random graph
    igraph_k_regular_game(&R_graph, n, k, IGRAPH_UNDIRECTED, 0);

    igraph_is_connected(&R_graph, &connected, IGRAPH_WEAK);
    while(!connected){
        igraph_destroy(&R_graph);
        igraph_k_regular_game(&R_graph, n, k, IGRAPH_UNDIRECTED, 0);

        igraph_is_connected(&R_graph, &connected, IGRAPH_WEAK);
    }

    return R_graph;
}



igraph_t generateRandomGraph(int type, int n)
{
    igraph_t random_graph;

    switch (type) {
        case 1:
            random_graph = generateGeometricGraph(n, sqrt(100.0/(float)n));
            break;
        case 2:
            random_graph = generateBarabasiAlbertGraph(n, 1.0, 5, 1.0);
            break;
        case 3:
            random_graph = generateErdosRenyiGraph(n, IGRAPH_ERDOS_RENYI_GNP, 10.0/(float)n);
            // random_graph = generateErdosRenyiGraph(n, IGRAPH_ERDOS_RENYI_GNM, ceil(n^2/3));
            break;
        case 4:
            random_graph = generateRegularGraph(n, n-1);
            break;

        default:
            break;
    }

    return random_graph;

}


void printGraphType(int type)
{

    switch (type) {
        case 1:
            cout << "Geometric random graph\n";
            break;
        case 2:
            cout << "Barabasi-Albert random graph\n";
            break;
        case 3:
            cout << "Erdos-Renyi random graph\n";
            break;
        case 4:
            cout << "Regular random graph\n";
            break;

        default:
            break;
    }

}