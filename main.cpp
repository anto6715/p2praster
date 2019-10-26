#include <iostream>
#include <igraph/igraph.h>
#include <cstring>
#include <random>
#include "raster.h"

void clusterMerge(vector<unordered_set<array<int, 2>, container_hasher>> &clustersIn, vector<unordered_set<array<int, 2>, container_hasher>> &clustersInOut);
void printOrderedProjection(int peerID, int peers, unordered_map<array<int, 2>, double, container_hasher> &projection);
void simultaneousMaxMin(unordered_map<array<int, 2>, double, container_hasher> &projection, int *maxX, int *minX, int *maxY, int *minY);
void getGridSize(int *p, int * q, int peers, int min);
void projectionMerge(unordered_map<array<int, 2>, double, container_hasher> &projectionIn, unordered_map<array<int, 2>, double, container_hasher> &projectionInOut);
void average(double *x, double *y);
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
    int         peers = 10; // number of peers
    int         fanOut = 4; //fan-out of peers
    int         graphType = 3; // graph distribution: 1 geometric 2 Barabasi-Albert 3 Erdos-Renyi 4 regular (clique)
    double      convThreshold = 0.00001; // local convergence tolerance
    int         convLimit = 3; // number of consecutive rounds in which a peer must locally converge
    int         roundsToExecute = -1;
    int         p_star = -1;
    double      delta = 0.04;
    double      precision = -4.2;
    int         threshold = 3;
    int         min_size = 3;
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
            fanOut = stoi(argv[i]);;
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

    auto dataset_storage = new double[row*column];
    auto dataset = new double*[row];
    for (int i = 0; i < row; i++) {
        dataset[i] = &dataset_storage[i*column];
    }
    StartTheClock();
    loadData(dataset, name_file, column);
    double greaddatasettime = StopTheClock();
    if (!outputOnFile) {
        cout <<"Time (seconds) required to load the dataset: " << greaddatasettime << "\n";
    }


    /*** Compute last item for each peer***/
    peerLastItem = new long[peers]();//(long *) calloc(peers, sizeof(long));
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
    int lastNotNull=0;
    long sum = peerLastItem[0] + 1;
    //cerr << "peer 0:" << sum << "\n";
    for(int i = 1; i < peers; i++) {
        if (!peerLastItem[i] && peerLastItem[i] != peerLastItem[i-1]) {
            lastNotNull = i-1;
            //cerr << "peer " << i << ":" << 0 << "\n";
        } else {
            if (peerLastItem[i-1] != 0) {
                sum += peerLastItem[i] - peerLastItem[i-1];
                //cerr << "peer " << i << ":" << peerLastItem[i] - peerLastItem[i-1] << "\n";
            }
            else if (peerLastItem[i]){
                sum += peerLastItem[i] - peerLastItem[lastNotNull];
                //cerr << "peer " << lastNotNull+1 << ":" << peerLastItem[i] - peerLastItem[lastNotNull] << "\n";
            }
        }


    }
    cout<< "sum: " << sum << endl;
    if(sum != ni) {
        cerr << "ERROR: ni = "<< ni << "!= sum = " << sum << endl;
        exit(EXIT_FAILURE);
    }
    cout.flush();

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

    outputOnFile = !params.outputFilename.empty();

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

    /*** Apply first Raster projection to each peer's dataset ***/
    if (!outputOnFile) {
        cout << "\nApply first Raster projection to each peer's dataset...\n";
    }
    auto *projection = new unordered_map<array<int, 2>, double, container_hasher>[params.peers];
    StartTheClock();
    int start = 0;
    for(int peerID = 0; peerID < params.peers; peerID++){
        /*** points agglomeration ***/
        mapToTilesNoThreshold(dataset, precision, threshold, projection[peerID], start, peerLastItem[peerID]);

        start = peerLastItem[peerID] + 1;
    }

    double rastertime = StopTheClock();
    if (!outputOnFile) {
        cout << "Raster done!\n";
        cout << "Time (seconds) required by first (distributed) Raster projection: " << rastertime << endl;
    }
    // this is used to estimate the number of peers
    auto *dimestimate = new double[params.peers]();
    dimestimate[0] = 1;

    auto *prevestimate = new double[params.peers]();

    auto *datasetsizeestimate = new double[params.peers]();
    datasetsizeestimate[0] = peerLastItem[0] + 1;
    for(int i = 1; i < peers; i++)
        datasetsizeestimate[i] = peerLastItem[i] - peerLastItem[i-1];

    int Numberofconverged = params.peers;

    auto *converged = new bool[params.peers]();
    for(int i = 0; i < params.peers; i++)
        converged[i] = false;

    auto *convRounds = new int[params.peers]();

    int rounds = 0;

    if (!outputOnFile) {
        cout <<"\nStarting distributed agglomeration merge..." << endl;
    }

    /*** Merge information about agglomeration ***/
    while( (params.roundsToExecute < 0 && Numberofconverged) || params.roundsToExecute > 0){

        memcpy(prevestimate, dimestimate, params.peers * sizeof(double));

        for(int peerID = 0; peerID < params.peers; peerID++){
            // check peer convergence


            // determine peer neighbors
            igraph_vector_t neighbors;
            igraph_vector_init(&neighbors, 0);
            igraph_neighbors(&graph, &neighbors, peerID, IGRAPH_ALL);
            long neighborsSize = igraph_vector_size(&neighbors);
            if(fanOut < neighborsSize){
                // randomly sample f adjacent vertices
                igraph_vector_shuffle(&neighbors);
                igraph_vector_remove_section(&neighbors, params.fanOut, neighborsSize);
            }

            neighborsSize = igraph_vector_size(&neighbors);
            for(int i = 0; i < neighborsSize; i++){
                int neighborID = (int) VECTOR(neighbors)[i];
                igraph_integer_t edgeID;
                igraph_get_eid(&graph, &edgeID, peerID, neighborID, IGRAPH_UNDIRECTED, 1);

                projectionMerge(projection[neighborID], projection[peerID]);

                double mean = (dimestimate[peerID] + dimestimate[neighborID]) / 2;
                dimestimate[peerID] = mean;
                dimestimate[neighborID] = mean;

                double mean_size = (datasetsizeestimate[peerID] + datasetsizeestimate[neighborID]) / 2;
                datasetsizeestimate[peerID] = mean_size;
                datasetsizeestimate[neighborID] = mean_size;
            }

            igraph_vector_destroy(&neighbors);
        }

        // check local convergence
        if (params.roundsToExecute < 0) {
            for(int peerID = 0; peerID < params.peers; peerID++){

                if(converged[peerID])
                    continue;

                bool dimestimateconv;
                if(prevestimate[peerID])
                    dimestimateconv = fabs((prevestimate[peerID] - dimestimate[peerID]) / prevestimate[peerID]) < params.convThreshold;
                else
                    dimestimateconv = false;

                if(dimestimateconv)
                    convRounds[peerID]++;
                else
                    convRounds[peerID] = 0;

                //printf ("PeerID %d, round %d, convRound %d\n", peerID, rounds, convRounds[peerID]);


                converged[peerID] = (convRounds[peerID] >= params.convLimit);
                if(converged[peerID]){
                    //printf("peer %d rounds before convergence: %d\n", peerID, rounds + 1);
                    Numberofconverged --;
                }
            }
        }
        rounds++;
        cerr << "\r Active peers: " << Numberofconverged << " - Rounds: " << rounds << "          " << endl;

        params.roundsToExecute--;
    }

    /*** Get dimestimate***/
    for(int peerID = 0; peerID < params.peers; peerID++) {
        dimestimate[peerID] = round(1 / dimestimate[peerID]);
        //cout << "dimensione stimata: " << dimestimate[peerID] << endl;
    }

    /*** Remove tiles < threshold***/
    for(int peerID = 0; peerID < params.peers; peerID++){
        unordered_map<array<int, 2>, double, container_hasher>::iterator it;
        it = projection[peerID].begin();
        while (it != projection[peerID].end()) {
            it -> second = round(it->second *dimestimate[peerID]);
            if (it -> second < threshold) {
                projection[peerID].erase(it++);
            } else {
                it++;
            }
        }
    }

    /*** Simultaneous minimum and maximum algorithm***/
    int *minX = (int *) calloc(params.peers, sizeof(int));
    auto minY = (int *) calloc(params.peers, sizeof(int));
    int *maxX = (int *) calloc(params.peers, sizeof(int));
    int *maxY = (int *) calloc(params.peers, sizeof(int));

    for(int peerID = 0; peerID < params.peers; peerID++){
        simultaneousMaxMin(projection[peerID], &maxX[peerID], &minX[peerID], &maxY[peerID], &minY[peerID]);
    }


    /*** per stampare il dataset (ordinato) ottenuto da ogni peer***/
    /*for(int peerID = 0; peerID < params.peers; peerID++) {
        printOrderedProjection(peerID, params.peers, projection[peerID]);
    }*/

    // (x1,y1) is the top-right corner coordinate of a square in checkerboard partition
    // (x2,y2) is the bottom-left corner cordinate of a square in checkerboard partition
    auto y1 = new int*[params.peers];
    auto x1 = new int*[params.peers];
    auto y2 = new int*[params.peers];
    auto x2 = new int*[params.peers];

    // domain of dataset is divide in a grid p*q
    auto q = new int[params.peers];
    auto p = new int[params.peers];
    // how much square for each peer ( if p*q % peers != 0 some peers will have one more square)
    auto squares = new int[params.peers];
    int min =2; // set as input parameter?
    // initial hypothesis
    for (int i = 0; i < params.peers; i++) {
        squares[i] = min;
    }

    /*** Each peer obtain its "square(s)" coordinate(s)***/
    for(int peerID = 0; peerID < params.peers; peerID++) {
        int a, b, restA, restB;
        // find grid p*q in order to assign at least one tile for each peer, with p = int_sup(sqrt(peers))
        getGridSize(&p[peerID], &q[peerID],(int) dimestimate[peerID], squares[peerID]);

        // compute how much tiles for each square on y (a) and x (b) axis
        a =     (maxY[peerID]-minY[peerID]) / q[peerID];
        b =     (maxX[peerID]-minX[peerID]) / p[peerID];
        restA = (maxY[peerID]-minY[peerID]) % q[peerID];
        restB = (maxX[peerID]-minX[peerID]) % p[peerID];
        //cout << "peer: " << peerID << endl;

        // i, j coordinate in grid p*q,
        int *i, *j;
        if (peerID < (p[peerID] * q[peerID]) % (int) dimestimate[peerID] ) {
            squares[peerID]++;
        }
        i = new int[squares[peerID]];
        j = new int[squares[peerID]];

        y1[peerID] = new int[squares[peerID]];
        x1[peerID] = new int[squares[peerID]];
        y2[peerID] = new int[squares[peerID]];
        x2[peerID] = new int[squares[peerID]];


        for (int k = 0; k < squares[peerID]; k++) {
            i[k] = (peerID + k*(int)dimestimate[peerID]) /p[peerID] + 1;
            j[k] = (peerID + k*(int)dimestimate[peerID]) %p[peerID] + 1;

            y1[peerID][k] = i[k] <= restA ? i[k]*(a+1) + minY[peerID] : restA*(a+1) + (i[k]-restA)*a + minY[peerID];
            x1[peerID][k] = j[k] <= restB ? j[k]*(b+1) + minX[peerID] : restB*(b+1) + (j[k]-restB)*b + minX[peerID];
            y2[peerID][k] = (i[k] <= restA && restA != 0) ? y1[peerID][k] - (a+1) : y1[peerID][k] - (a);
            x2[peerID][k] = (j[k] <= restB  && restB != 0) ? x1[peerID][k] - (b+1) : x1[peerID][k] - (b);
            // include the left and bottom edge of grid
            if (y2[peerID][k] == minY[peerID])
                y2[peerID][k] -= 1;

            if (x2[peerID][k] == minX[peerID])
                x2[peerID][k] -= 1;
            //cout << "x1: " << y1[peerID][k]<< " x2: " << x1[peerID][k] << endl;
        }
        delete[] i;
        delete[] j;
    }

    delete[] p;
    delete[] q;
    delete[] squares;
    free(minX);
    free(minY);
    free(maxX);
    free(maxY);

    /*** Each peer find the tiles in own square(s)***/
    auto *squareProjection = new unordered_map<array<int, 2>, double, container_hasher>[params.peers];
    for(int peerID = 0; peerID < params.peers; peerID++) {
        int skip;
        unordered_map<array<int, 2>, double, container_hasher>::iterator it;
        it = projection[peerID].begin();
        while (it != projection[peerID].end()) {
            skip = 0;
            int a = (it -> first)[0];
            int b = (it -> first)[1];
            for (int k = 0; k < squares[peerID]; k++) {
                if (a <= x1[peerID][k] && a > x2[peerID][k] && b <= y1[peerID][k] && b > y2[peerID][k]) {
                    (squareProjection[peerID])[it -> first] = it -> second;
                    projection[peerID].erase(it++);
                    skip = 1;
                    break;
                }
            }
            if (!skip)
                it++;
        }
    }

    for(int peerID = 0; peerID < params.peers; peerID++) {
        delete[] x1[peerID];
        delete[] y1[peerID];
        delete[] x2[peerID];
        delete[] y2[peerID];
    }

    delete[] x1;
    delete[] y1;
    delete[] x2;
    delete[] y2;

    if (!outputOnFile) {
        cout <<"\nStarting local clustering..." << endl;
    }
    StartTheClock();

    /***Each peer clustering its own tiles ***/
    auto *clusters = new vector<unordered_set<array<int, 2>, container_hasher>>[params.peers];
    for(int peerID = 0; peerID < params.peers; peerID++) {
        clusteringTiles(squareProjection[peerID], projection[peerID], params.min_size, clusters[peerID]);
    }


    for(int peerID = 0; peerID < params.peers; peerID++) {
        projection[peerID].clear();
        squareProjection[peerID].clear();
    }
    delete[] projection;
    delete[] squareProjection;

    double clustertime = StopTheClock();
    if (!outputOnFile) {
        cout << "Local clustering done!\n";
        cout << "Time (seconds) required by local Raster clustering: " << clustertime << endl;
    }

    auto clustersestimate = new double[params.peers]();
    auto prevclustersestimate = new double[params.peers]();
    for(int i = 0; i < params.peers; i++)
        clustersestimate[i] = (double )clusters[i].size();

    for(int i = 0; i < params.peers; i++)
        converged[i] = false;

    rounds =0;

    params.roundsToExecute = roundsToExecute;
    Numberofconverged = params.peers;


    if (!outputOnFile) {
        cout <<"\nStarting distributed merge of clusters..." << endl;
    }
    StartTheClock();
    /*** Start distributed cluster merge***/
    while( (params.roundsToExecute < 0 && Numberofconverged) || params.roundsToExecute > 0){

        memcpy(prevclustersestimate, clustersestimate, params.peers * sizeof(double));

        for(int peerID = 0; peerID < params.peers; peerID++){
            // check peer convergence

            // determine peer neighbors
            igraph_vector_t neighbors;
            igraph_vector_init(&neighbors, 0);
            igraph_neighbors(&graph, &neighbors, peerID, IGRAPH_ALL);
            long neighborsSize = igraph_vector_size(&neighbors);
            if(fanOut < neighborsSize){
                // randomly sample f adjacent vertices
                igraph_vector_shuffle(&neighbors);
                igraph_vector_remove_section(&neighbors, params.fanOut, neighborsSize);
            }

            neighborsSize = igraph_vector_size(&neighbors);
            for(int i = 0; i < neighborsSize; i++){
                int neighborID = (int) VECTOR(neighbors)[i];
                igraph_integer_t edgeID;
                igraph_get_eid(&graph, &edgeID, peerID, neighborID, IGRAPH_UNDIRECTED, 1);

                clusterMerge(clusters[neighborID], clusters[peerID]);
                clustersestimate[peerID] = (double) clusters[peerID].size();
                clustersestimate[neighborID] = (double) clusters[neighborID].size();

            }

            igraph_vector_destroy(&neighbors);
        }

        // check local convergence
        if (params.roundsToExecute < 0) {
            for(int peerID = 0; peerID < params.peers; peerID++){

                if(converged[peerID])
                    continue;

                bool clustersestimateconv;
                if(prevclustersestimate[peerID])
                    clustersestimateconv = fabs((prevclustersestimate[peerID] - clustersestimate[peerID]) / prevclustersestimate[peerID]) < params.convThreshold;
                else
                    clustersestimateconv = false;

                if(clustersestimateconv)
                    convRounds[peerID]++;
                else
                    convRounds[peerID] = 0;

                //printf ("PeerID %d, round %d, convRound %d\n", peerID, rounds, convRounds[peerID]);


                converged[peerID] = (convRounds[peerID] >= params.convLimit);
                if(converged[peerID]){
                    //printf("peer %d rounds before convergence: %d\n", peerID, rounds + 1);
                    Numberofconverged --;
                }
            }
        }
        rounds++;
        cerr << "\r Active peers: " << Numberofconverged << " - Rounds: " << rounds << "          " << endl;
        params.roundsToExecute--;
    }

    double mergeclustertime = StopTheClock();
    if (!outputOnFile) {
        cout << "Distributed  merge clusters done!\n";
        cout << "Time (seconds) required by distributed merge clusters: " << mergeclustertime << endl;
    }

    /*** Print info about each peer's clusters***/
    for(int peerID = 0; peerID < params.peers; peerID++) {
        cout << "\npeer: " << peerID << endl;
        printClusters(clusters[peerID], params.precision, peerID);
    }
    
    delete[] dimestimate;
    delete[] peerLastItem;
    delete[] converged;
    delete[] dataset;
    delete[] dataset_storage;
    delete[] prevestimate;
    delete[] prevclustersestimate;
    delete[] convRounds;
    delete[] datasetsizeestimate;

}

void printOrderedProjection(int peerID, int peers, unordered_map<array<int, 2>, double, container_hasher> &projection) {
    ofstream outfile(to_string(peerID) +".csv");
    set<array<int,2>> data;
    unordered_map<array<int, 2>, double, container_hasher>::iterator it;
    set<array<int,2>>::iterator it2;
    it = projection.begin();
    while (it != projection.end()) {
        data.insert(it -> first);
        it++;
    }

    it2 = data.begin();
    while (it2 != data.end()) {
        outfile << "tile: " << (*it2)[0] << "," << (*it2)[1] << "\n";
        it2++;
    }
    outfile.close();

}

void simultaneousMaxMin(unordered_map<array<int, 2>, double, container_hasher> &projection, int *maxX, int *minX, int *maxY, int *minY) {
    int *x1 = new int;
    int *x2 = new int;
    int *y1 = new int;
    int *y2 = new int;
    unordered_map<array<int, 2>, double, container_hasher>::iterator it;
    it = projection.begin();
    if ((projection.size() %2) != 0) {  // odd case
        *minX = (it -> first)[0];  // set as minX the first value
        *minY = (it -> first)[1];  // set as minY the first value
        *maxX = (it -> first)[0];  // set as maxX the first value
        *maxY = (it -> first)[1];  // set as maxY the first value
        it++;
    } else {                            // even case
        *x1 = (it -> first)[0];
        *y1 = (it -> first)[1];
        it++;
        *x2 = (it -> first)[0];
        *y2 = (it -> first)[1];
        it++;
        *minX = (*x1 <= *x2) ? *x1 : *x2;
        *maxX = (*x1 <= *x2) ? *x2 : *x1;
        *minY = (*y1 <= *y2) ? *y1 : *y2;
        *maxY = (*y1 <= *y2) ? *y2 : *y1;
    }

    while (it != projection.end()) {
        *x1 = (it -> first)[0];
        *y1 = (it -> first)[1];
        it++;
        *x2 = (it -> first)[0];
        *y2 = (it -> first)[1];
        it++;
        if (*x1 <= *x2){
            if ( *x1 < *minX) {
                *minX = *x1;
            }

            if ( *x2 > *maxX) {
                *maxX = *x2;
            }
        } else {
            if ( *x2 < *minX) {
                *minX = *x2;
            }

            if ( *x1 > *maxX) {
                *maxX = *x1;
            }
        }

        if (*y1 <= *y2){
            if ( *y1 < *minY) {
                *minY = *y1;
            }

            if ( *y2 > *maxY) {
                *maxY = *y2;
            }
        } else {
            if ( *y2 < *minY) {
                *minY = *y2;
            }

            if ( *y1 > *maxY) {
                *maxY = *y1;
            }
        }
    }

    delete x1;
    delete y1;
    delete x2;
    delete y2;
}

void getGridSize(int *p, int *q, int peers, int min) {
    *q = sqrt(peers*min);
    if ((*q * *q) != peers) {
        *p = *q+1;
        *q = *q * (*q + 1) < peers ? ++*q : *q;
    } else {
        *p = *q;
    }
}

void projectionMerge(unordered_map<array<int, 2>, double, container_hasher> &projectionIn, unordered_map<array<int, 2>, double , container_hasher> &projectionInOut){
    unordered_map<array<int, 2>, double, container_hasher>::iterator itIn;
    unordered_map<array<int, 2>, double , container_hasher>::iterator itInOut;
    itIn = projectionIn.begin();
    while (itIn != projectionIn.end()) {
        itInOut = projectionInOut.find(itIn -> first);
        if (itInOut != projectionInOut.end()) {
            average( &itInOut->second, &itIn->second);
            memcpy(&itIn->second, &itInOut->second, sizeof(double));
        } else {
            projectionInOut[itIn -> first] = (itIn->second)/2.0;
            average(&itIn->second, nullptr);
        }
        itIn++;
    }



}

void clusterMerge(vector<unordered_set<array<int, 2>, container_hasher>> &clustersIn, vector<unordered_set<array<int, 2>, container_hasher>> &clustersInOut) {
    int equals;
    set<int> notCopy; // indexes of common clusters
    unordered_set<array<int, 2>, container_hasher>::iterator it1;
    unordered_set<array<int, 2>, container_hasher>::iterator it2;
    int check = clustersInOut.size();
    for (auto & k : clustersIn) {
        equals = 0;
        for (int l = 0; l < check; l++) {
            if (k.size() == clustersInOut.at(l).size()) {
                it1 = k.begin();
                if ( clustersInOut.at(l).find(*it1) != clustersInOut.at(l).end()) {
                    notCopy.insert(l);
                    equals = 1;
                    break;
                }
            }
        }
        if (!equals) {
            clustersInOut.push_back(k);
        }
    }
    if (notCopy.empty()) {
        for (int l = 0; l < check; l++) {
            clustersIn.push_back(clustersInOut.at(l));
        }
    } else {
        for (int l = 0; l < check; l++) {
            if (notCopy.find(l) == notCopy.end()) {
                clustersIn.push_back(clustersInOut.at(l));
            }
        }

    }
    notCopy.clear();


}

void average(double *x, double *y)  {
    if (x == nullptr) {
        cerr << "x cannot be null!";
        exit(EXIT_FAILURE);
    } else if (y == nullptr) {
        *x = *x/2.0;
    } else {
        *x = (*x+*y)/2.0;
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