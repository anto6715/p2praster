/** @file   */

#include <iostream>
#include <igraph/igraph.h>
#include <cstring>
#include <random>
#include "raster.h"
#include "error.h"


/**
 * This function compute the Manhattan distance between two points p1,p2
 * with coordinate (x1,y1) and (x2,y2):
 * L1(p1,p2) = |x1-x2| + |y1-y2|
 *
 * @param [in] x1 - Coordinate x of P1
 * @param [in] x2 - Coordinate x of P2
 * @param [in] y1 - Coordinate y of P1
 * @param [in] y2 - Coordinate y of P2
 * @return Return Manhattan Distance
 */
int manhattanDistance(int x1, int x2, int y1, int y2);


/**
 * This function compute the centroids for each clusters using a weighted average.
 * Save into centroids[i] the centroid of clusters[i] (for each cluster in clusters)
 *
 * @param [in] clusters - Data structure which contains clusters
 * @param [in,out] centroids - Data structure where insert centroid of respectively cluster
 * @return 0 if success, -3 in case of bad data structures,
 */
int getCentroids(vectorSet3 &clusters, vector<array<int, 2>> &centroids);


/**
 * This function get cluster from peer and its neighbor and do the union.
 * As result we obtain an only cluster that contains all the clusters that
 * are in peer cluster but not in neighbor cluster (and viceversa) and only
 * one copy if a cluster is in both peer and neighbor clusters.
 * Coordinate of centroids is used to establishing if two clusters
 * are the same or not
 *
 * @param [in,out] clustersN - Data structure which contains clusters of peer's neighbor
 * @param [in,out] clustersP - Data structure which contains clusters of peer
 * @param [in,out] centroidsN - Data structure which contains centroids of peer's neighbor
 * @param [in,out] centroidsP - Data structure which contains centroids of peer
 * @param [in] maxDistance - Max distance to consider two cluster similar
 * @return 0 if success, -1 in case of memory error, -25 in case of cluster union error
 */
int clustersMerge(vectorSet3 &clustersN, vectorSet3 &clustersP, vector<array<int, 2>> &centroidsN, vector<array<int, 2>> &centroidsP, int maxDistance);


/**
 * This function merge two cluster that are similar but not equal,
 * it add the tiles that are in common or compute the average for
 * tiles that have the same coordinate but different cardinality
 *
 * @param [in,out] clusterN - Data structure which contains cluster of peer's neighbor
 * @param [in,out] clusterP - Data structure which contains cluster of peer
 * @return 0 if success, -3 in case of bad data structures, -2 in case of insert error, -4 in case of arithmetic error
 */
int clusterUnion(unSet3 &clusterN, unSet3 &clusterP);
// to delete
void printOrderedProjection(int peerID, int peers, unordered_map<array<int, 2>, double, container_hasher> &projection);


/**
 * This function implements the simultaneous Max and Min algorithm on both x and y axes
 *
 * @param [in] projection - Data structure with tiles
 * @param [in,out] maxX - Where save max value for tile x coordinates
 * @param [in,out] minX - Where save min value for tile x coordinates
 * @param [in,out] maxY - Where save max value for tile y coordinates
 * @param [in,out] minY - Where save min value for tile y coordinates
 * @return 0 if success, -3 in case of bad data structure
 */
int simultaneousMaxMin(hashmap &projection, int *maxX, int *minX, int *maxY, int *minY);


/**
 * This function compute the parameter p and q to obtain a checkerboard partition (with grid p x q)
 * of the global projection and guarantees that each peer has at least min "square".
 *
 *
 * @param [in,out] p - How many squares on x axis
 * @param [in,out] q - How many squares on y axis
 * @param [in] peers - Number of peers in the network
 * @param [in] min - Min square for each peer
 * @return 0 if success
 */
int getGridSize(int &p, int &q, int peers, int min);


/**
 * This function merge the tiles of two structures. If a tile is common,
 * it compute the average on the cardinality of each tile and update the value on both structure.
 * If a tiles is present in an only structure, it value is divides by 2 and the tile
 * is copies also in the other structure.
 *
 * @param projectionN - Data structure with peer neighbor tiles
 * @param projectionP - Data structure with  peer tiles
 * @return 0 if success, -2 in case of insert error, -3 in case of bad data structure
 */
int projectionMerge(hashmap &projectionN, hashmap &projectionP);


/**
 * This function compute the average between two double.
 * The result is stored into x variable. If y is null the function
 * simply return x/2, otherwise return the average.
 *
 * @param x - First number that will be update
 * @param y - Second number
 * @return 0 if success, -4 if x is invalid
 */
int average(double *x, double y);


/**
 *
 * @param [in] dataset_storage - Array where store all dataset
 * @param [in] dataset - Array of pointer that pints at dataset_storage
 * @param [in] name_file - Name input file
 * @param [in,out] ni - Nu,ber of points
 * @return 0 if success, -1 in case of memory error, -7 in case of input file error
 */
int readDataset(double **dataset_storage, double ***dataset, string name_file, int &ni);


/**
 * This function compute the last item of the dataset for each peer.
 * Peer with i will obtain all data from peerLastItem[i-1] and peerLastItem[i]
 *
 * @param [in,out] peerLastItem - Data structure where save the result
 * @param [in] peers - Number of peers
 * @param [in] ni - NUmber of points into dataset
 * @return 0 if success
 */
int partitionDataset(long *peerLastItem, int peers, int ni);


/**
 * This function check if previous dataset partition is correct
 *
 * @param peerLastItem - Data structure
 * @param [in] peers - Number of peers
 * @param [in] ni - NUmber of points into dataset
 * @return 0 if success, -9 in case of partition error
 */
int checkPartitioning(long *peerLastItem, int peers, int ni);


/**
 *
 * @param graph - Structure where save graph generated
 * @param peers - Number of peers
 * @param graphType - Type of graph
 * @return 0 if success
 */
int generateGraph(igraph_t &graph, int peers, int graphType);

/**
 * This function print on terminal tha minimum and maximum vertex degree
 *
 * @param graph - Structure that contains graph
 * @return
 */
int graphProperties(igraph_t &graph);


/**
 * This function return the number of neighbor of calling peer,
 * if this fanOut < neighborSize then return fanOut and remove
 * randomly (neighborsSize - fanOut) neighbors from neighbors structure.
 * If fanOut = -1 return neighborsSize.
 *
 * @param [in] graph
 * @param [in,out] neighbors
 * @param [in] peerID
 * @param [in] fanOut
 * @return The dimension of neighbors
 */
int getNeighborsSize(igraph_t &graph, igraph_vector_t &neighbors, int peerID, int fanOut);


int checkFirstConvegence(double prevestimate, double dimestimate, double convThreshold, int &convRounds, bool &converged, int convLimit, int &Numberofconverged);


int restoreCardinality(hashmap &projection, int dimestimate);

int getGroupsTiles( int &nTilesX, int &nTilesY, int &restTilesX, int &restTilesY,
                    int q, int p, int minX, int minY, int maxX, int maxY );

int getSquares(int &squares, int peerID, int p, int q, int dimestimate, int minSquares);

int getSquaresCoordinate(int &x1, int &y1, int &y2, int &x2, int p, int q, int dimestimate,
                         int minX, int minY, int maxX, int maxY, int k, int peerID);


/**
 * When called save the actual time into global variable t1
 */
void StartTheClock();


/**
 * When called save the actual time into gloab variable t2
 * and compute the difference t1-t2
 * @return
 */
double StopTheClock();


/**
 * Print on terminal the input parameter accepted
 * @param cmd - Name of script
 */
void usage(char* cmd);

igraph_t generateGeometricGraph(igraph_integer_t n, igraph_real_t radius);
igraph_t generateBarabasiAlbertGraph(igraph_integer_t n, igraph_real_t power, igraph_integer_t m, igraph_real_t A);
igraph_t generateErdosRenyiGraph(igraph_integer_t n, igraph_erdos_renyi_t type, igraph_real_t param);
igraph_t generateRegularGraph(igraph_integer_t n, igraph_integer_t k);
igraph_t generateRandomGraph(int type, int n);
void printGraphType(int type);

using namespace std;
using namespace std::chrono;

/*!< @struct Params - A structure containing parameters read from command-line.  */
struct Params {
    int        ni;                /*!< The number of point of dataset */
    int         peers;             /*!< The number of peers */
    string      outputFilename;    /*!< The path for the output file */
    double      convThreshold;     /*!< The local convergence tolerance for the consensus algorithm */
    int         convLimit;         /*!< The number of consecutive rounds in which a peer must locally converge */
    int         graphType;         /*!< The graph distribution: 1 geometric, 2 Barabasi-Albert, 3 Erdos-Renyi, 4 regular (clique) */
    int         radius;            /*!< Parameter for graph generation */
    int         fanOut;            /*!< The number of communication that a peer can carry out in a round (-1 ability communication with all neighbors) */
    int         roundsToExecute;   /*!< The number of rounds to carry out in the consensus algorithm */
    double      precision;         /*!< Raster parameters that determines tile dimension */
    int         threshold;         /*!< Raster parameters that establish if maintain or not a tile */
    int         min_size;          /*!< Raster parameters that establish the minimum number of tile in order to form a cluster */
    int         maxDistance;       /*!< Max distance in order to consider two cluster the same */
    int         minSquares;        /*!< minimum squares for each peer from checkerboard partition*/
    string      name_file;         /*!< The path for the input CSV file */
    int         typeAlgorithm;     /*!< 0 if want to make cluster in distributed manner, otherwise each peer after first communication clustering all global projection */
};

void parametersSummary(Params params);

high_resolution_clock::time_point t1, t2;
/**
 *
 * @param argc
 * @param argv
 * @return 0 if success, -1 in case of memory error, -8 in case of dataset error
 */
int main(int argc, char **argv) {
    /*** Default Parameters***/
    int        ni; // points number
    int         peers = 10; // number of peers
    int         fanOut = 3; //fan-out of peers
    int         graphType = 2; // graph distribution: 1 geometric 2 Barabasi-Albert 3 Erdos-Renyi 4 regular (clique)
    double      convThreshold = 0.0001; // local convergence tolerance
    int         convLimit = 3; // number of consecutive rounds in which a peer must locally converge
    int         roundsToExecute = -1;
    double      precision = -4.2;
    int         threshold = 2;
    int         min_size = 3;
    int         radius = 0;
    int         maxDistance = 2;
    int         typeAlgorithm = 0;
    int         minSquares = 1;
    string      name_file = "../Datasets/S-sets/s1.csv";

    bool            outputOnFile = false;
    string          outputFilename;
    igraph_t        graph;
    Params          params;

    int start;

    /*** Array declarations ***/
    long *peerLastItem = nullptr; // index of a peer last item
    double *dataset_storage = nullptr;
    double **dataset = nullptr;
    hashmap *projection = nullptr;
    double *dimestimate = nullptr;
    double *prevestimate = nullptr;
    bool *converged = nullptr;
    int *convRounds = nullptr;
    int *minX = nullptr;
    int *minY = nullptr;
    int *maxX = nullptr;
    int *maxY = nullptr;
    int** y1 = nullptr;
    int** x1 = nullptr;
    int** y2 = nullptr;
    int** x2 = nullptr;
    int* q = nullptr;
    int* p = nullptr;
    int* squares = nullptr;
    hashmap *squareProjection = nullptr;

    /*** parse command-line parameters ***/
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-p") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of peers parameter." << endl;
                return -1;
            }
            peers = stoi(argv[i]);
        }  else if (strcmp(argv[i], "-f") == 0) {
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
        } else {
            usage(argv[0]);
            return -1;
        }
    }

    /*** read dataset dimensions ***/
    StartTheClock();
    if (readDataset(&dataset_storage, &dataset, name_file, ni)) {
        return readDatasetError(__FUNCTION__);

    }
    double greaddatasettime = StopTheClock();
    if (!outputOnFile) {
        cout <<"Time (seconds) required to load the dataset: " << greaddatasettime << "\n";
    }

    /** Compute last item for each peer */
    peerLastItem = new (nothrow) long[peers]();
    if (!peerLastItem)
        return memoryError(__FUNCTION__);
    if (partitionDataset(peerLastItem, peers, ni))
        return partitionError(__FUNCTION__);

    /** check the partitioning correctness */
    if (checkPartitioning(peerLastItem, peers, ni))
        return partitionError(__FUNCTION__);


    /** assign parameters read from command line */
    params.name_file = name_file;
    params.ni = ni;
    params.peers = peers;
    params.fanOut = fanOut;
    params.graphType = graphType;
    params.convThreshold = convThreshold;
    params.convLimit = convLimit;
    params.outputFilename = outputFilename;
    params.roundsToExecute = roundsToExecute;
    params.precision = precision;
    params.threshold = threshold;
    params.min_size = min_size;
    params.radius = radius;
    params.maxDistance = maxDistance;
    params.typeAlgorithm = typeAlgorithm;
    params.minSquares = minSquares;

    outputOnFile = !params.outputFilename.empty();

    if (!outputOnFile) {
        parametersSummary(params);
    }

    /** Graph generation */
    StartTheClock();
    // generate a connected random graph
    generateGraph(graph, peers, params.graphType);

    double gengraphtime = StopTheClock();
    if (!outputOnFile) {
        cout <<"Time (seconds) required to generate the random graph: " << gengraphtime << "\n";
    }

    // determine minimum and maximum vertex degree for the graph
    graphProperties(graph);

    /*** Apply first Raster projection to each peer's dataset ***/
    if (!outputOnFile) {
        cout << "\nApply first Raster projection to each peer's dataset...\n";
    }

    projection = new (nothrow) unordered_map<array<int, 2>, double, container_hasher>[params.peers];
    if(!projection)
        return memoryError(__FUNCTION__);

    StartTheClock();

    start = 0;
    for(int peerID = 0; peerID < params.peers; peerID++){
        /*** points projection ***/
        mapToTiles(dataset, precision, projection[peerID], start, peerLastItem[peerID]);
        start = peerLastItem[peerID] + 1;
    }

    double rastertime = StopTheClock();
    if (!outputOnFile) {
        cout << "Raster done!\n";
        cout << "Time (seconds) required by first (distributed) Raster projection: " << rastertime << endl;
    }

    // this is used to estimate the number of peers
    dimestimate = new (nothrow) double[params.peers]();
    if (!dimestimate)
        return memoryError(__FUNCTION__);
    dimestimate[0] = 1;

    prevestimate = new (nothrow) double[params.peers]();
    if (!prevestimate)
        return memoryError(__FUNCTION__);

    int Numberofconverged = params.peers;

    converged = new (nothrow) bool[params.peers]();
    if (!converged)
        return memoryError(__FUNCTION__);

    for(int i = 0; i < params.peers; i++)
        converged[i] = false;

    convRounds = new (nothrow) int[params.peers]();
    if (!convRounds)
        return memoryError(__FUNCTION__);

    int rounds = 0;

    if (!outputOnFile) {
        cout <<"\nStarting distributed agglomeration merge..." << endl;
    }
    StartTheClock();
    /*** Merge information about agglomeration ***/
    while( (params.roundsToExecute < 0 && Numberofconverged) || params.roundsToExecute > 0){

        memcpy(prevestimate, dimestimate, params.peers * sizeof(double));

        for(int peerID = 0; peerID < params.peers; peerID++){
            // determine peer neighbors
            igraph_vector_t neighbors;
            igraph_vector_init(&neighbors, 0);
            igraph_neighbors(&graph, &neighbors, peerID, IGRAPH_ALL);

            long neighborsSize = getNeighborsSize(graph, neighbors, peerID, fanOut);

            for(int i = 0; i < neighborsSize; i++){
                int neighborID = (int) VECTOR(neighbors)[i];
                igraph_integer_t edgeID;
                igraph_get_eid(&graph, &edgeID, peerID, neighborID, IGRAPH_UNDIRECTED, 1);

                if (projectionMerge(projection[neighborID], projection[peerID])) {
                    return mergeError(__FUNCTION__);
                }

                double mean = (dimestimate[peerID] + dimestimate[neighborID]) / 2;
                dimestimate[peerID] = mean;
                dimestimate[neighborID] = mean;
            }

            igraph_vector_destroy(&neighbors);
        }

        // check local convergence
        if (params.roundsToExecute < 0) {
            for(int peerID = 0; peerID < params.peers; peerID++){
                if(converged[peerID])
                    continue;
                checkFirstConvegence(prevestimate[peerID], dimestimate[peerID], params.convThreshold, convRounds[peerID], converged[peerID], params.convLimit, Numberofconverged);
            }
        }
        rounds++;
        cerr << "\r Active peers: " << Numberofconverged << " - Rounds: " << rounds << "          " << endl;

        params.roundsToExecute--;
    }

    double mergeprojectiontime = StopTheClock();
    if (!outputOnFile) {
        cout << "Merge Raster projection done!\n";
        cout << "Time (seconds) required by merge Raster projection: " << mergeprojectiontime << endl;
    }

    /*** Each peer update dimestimate***/
    for(int peerID = 0; peerID < params.peers; peerID++) {
        dimestimate[peerID] = round(1 / dimestimate[peerID]);
    }
    /*** Each peer restore tiles cardinality***/
    for(int peerID = 0; peerID < params.peers; peerID++){
        restoreCardinality(projection[peerID], dimestimate[peerID]);
    }

    /*** Each peer remove tiles < threshold ***/
    for(int peerID = 0; peerID < params.peers; peerID++){
        projectionThreshold(projection[peerID], threshold);
    }

    if (params.typeAlgorithm == 0) {
        minX = new (nothrow) int[params.peers];
        if (!minX)
            return memoryError(__FUNCTION__);
        minY = new (nothrow) int[params.peers];
        if (!minY)
            return memoryError(__FUNCTION__);
        maxX = new (nothrow) int[params.peers];
        if (!maxX)
            return memoryError(__FUNCTION__);
        maxY = new (nothrow) int[params.peers];
        if (!maxY)
            return memoryError(__FUNCTION__);

        /*** Simultaneous minimum and maximum algorithm***/
        for(int peerID = 0; peerID < params.peers; peerID++){
            simultaneousMaxMin(projection[peerID], &maxX[peerID], &minX[peerID], &maxY[peerID], &minY[peerID]);
        }

        /// (x1,y1) is the top-right corner coordinate of a square in checkerboard partition
        /// (x2,y2) is the bottom-left corner coordinate of a square in checkerboard partition
        /*** x1,x2,y1,y2 are matrices with dimension of n° peers x squares[peerID] (some peer can have an extra square from checkerboard partition) ***/
        y1 = new (nothrow) int*[params.peers];
        if (!y1)
            return memoryError(__FUNCTION__);
        x1 = new (nothrow) int*[params.peers];
        if (!x1)
            return memoryError(__FUNCTION__);
        y2 = new (nothrow) int*[params.peers];
        if (!y2)
            return memoryError(__FUNCTION__);
        x2 = new (nothrow) int*[params.peers];
        if (!x2)
            return memoryError(__FUNCTION__);

        /// domain of dataset is divide in a grid p*q
        q = new (nothrow) int[params.peers];
        if (!q)
            return memoryError(__FUNCTION__);
        p = new (nothrow) int[params.peers];
        if (!p)
            return memoryError(__FUNCTION__);

        /// how much square for each peer ( if p*q % peers != 0 some peers will have one extra square)

        squares = new (nothrow) int[params.peers];
        if (!squares)
            return memoryError(__FUNCTION__);

        /*** Each peer obtain its "square(s)" coordinate(s)***/
        for(int peerID = 0; peerID < params.peers; peerID++) {
            /// find grid p*q in order to assign at least one tile for each peer, with p = int_sup(sqrt(peers))
            getGridSize(p[peerID], q[peerID],(int) dimestimate[peerID], params.minSquares);

            /// peer gets its number of squares
            getSquares(squares[peerID], peerID, p[peerID], q[peerID], dimestimate[peerID], params.minSquares);

            y1[peerID] = new (nothrow) int[squares[peerID]];
            if (!y1[peerID])
                return memoryError(__FUNCTION__);
            x1[peerID] = new (nothrow) int[squares[peerID]];
            if (!x1[peerID])
                return memoryError(__FUNCTION__);
            y2[peerID] = new (nothrow) int[squares[peerID]];
            if (!y2[peerID])
                return memoryError(__FUNCTION__);
            x2[peerID] = new (nothrow) int[squares[peerID]];
            if (!x2[peerID])
                return memoryError(__FUNCTION__);

            for (int k = 0; k < squares[peerID]; k++) {
                /// Compute coordinate for k-th square
                getSquaresCoordinate(x1[peerID][k],y1[peerID][k], y2[peerID][k], x2[peerID][k], p[peerID], q[peerID],
                        dimestimate[peerID], minX[peerID], minY[peerID], maxX[peerID], maxY[peerID], k, peerID);
            }

        }
        delete[] p, p = nullptr;
        delete[] q, q = nullptr;

        delete[] minX, minX = nullptr;
        delete[] minY, minY = nullptr;
        delete[] maxX, maxX = nullptr;
        delete[] maxY, maxY = nullptr;

        /*** Each peer find the tiles in own square(s)***/

        squareProjection = new (nothrow) hashmap[params.peers];
        if (!squareProjection)
            return memoryError(__FUNCTION__);

        for(int peerID = 0; peerID < params.peers; peerID++) {
            int skip;
            hashmap::iterator it;
            it = projection[peerID].begin();
            while (it != projection[peerID].end()) {
                skip = 0;
                int a = (it -> first)[0];
                int b = (it -> first)[1];
                /*** Check if the tile is in one of peer square***/
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

        delete[] squares;

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


        /***Each peer clustering its own tiles and compute centroids of clusters***/
        vectorSet3 *clusters;
        try {
            clusters = new vector<unordered_set<array<int, 3>, container_hasher>>[params.peers];
        } catch (bad_alloc& ba) {
            cerr << "bad_alloc caught: " << ba.what() << '\n';
            exit(-1);
        }

        vector<array<int, 2>> *centroids;
        try {
            centroids = new vector<array<int, 2>>[params.peers];
        } catch (bad_alloc& ba) {
            cerr << "bad_alloc caught: " << ba.what() << '\n';
            exit(-1);
        }
        for(int peerID = 0; peerID < params.peers; peerID++) {
            clusteringTiles(squareProjection[peerID], projection[peerID], params.min_size, clusters[peerID]);
            getCentroids(clusters[peerID], centroids[peerID]);
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

        double *clustersestimate;
        try {
            clustersestimate = new double[params.peers]();
        } catch (bad_alloc& ba) {
            cerr << "bad_alloc caught: " << ba.what() << '\n';
            exit(-1);
        }
        double *prevclustersestimate;
        try {
            prevclustersestimate = new double[params.peers]();
        } catch (bad_alloc& ba) {
            cerr << "bad_alloc caught: " << ba.what() << '\n';
            exit(-1);
        }
        for(int i = 0; i < params.peers; i++)
            clustersestimate[i] = (double )clusters[i].size();

        // Reinitialize some parameters
        for(int i = 0; i < params.peers; i++){
            converged[i] = false;
            convRounds[i] = 0;
        }

        // reinitialize some parameters
        rounds = 0;
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

                // determine peer neighbors
                igraph_vector_t neighbors;
                igraph_vector_init(&neighbors, 0);
                igraph_neighbors(&graph, &neighbors, peerID, IGRAPH_ALL);

                long neighborsSize = getNeighborsSize(graph, neighbors, peerID, fanOut);

                for(int i = 0; i < neighborsSize; i++){
                    int neighborID = (int) VECTOR(neighbors)[i];
                    igraph_integer_t edgeID;
                    igraph_get_eid(&graph, &edgeID, peerID, neighborID, IGRAPH_UNDIRECTED, 1);

                    clustersMerge(clusters[neighborID], clusters[peerID], centroids[neighborID], centroids[peerID], params.maxDistance);

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
                        clustersestimateconv = fabs((prevclustersestimate[peerID] - clustersestimate[peerID])) == 0;
                    else
                        clustersestimateconv = false;

                    if(clustersestimateconv)
                        convRounds[peerID]++;
                    else
                        convRounds[peerID] = 0;

                    converged[peerID] = (convRounds[peerID] >= params.convLimit);
                    if(converged[peerID]){
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
            printClusters(clusters[peerID], peerID);
        }
        cout << "\n\n";

        unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> all_points;
        unordered_map<array<int, 2>, double, container_hasher> proj;
        mapToTilesPrime(dataset, precision, threshold, ni, proj, all_points);
        printAllPointsClustered(clusters[0], all_points);

        delete[] prevclustersestimate;
    } else {

        vector<unordered_set<array<int, 2>, container_hasher>> *clusters;
        try {
            clusters = new vector<unordered_set<array<int, 2>, container_hasher>>[params.peers];
        } catch (bad_alloc& ba) {
            cerr << "bad_alloc caught: " << ba.what() << '\n';
            exit(-1);
        }
        if (!outputOnFile) {
            cout <<"\nEach peer clustering the global projection..." << endl;
        }

        StartTheClock();

        for(int peerID = 0; peerID < params.peers; peerID++){
            clusteringTiles(projection[peerID], params.min_size, clusters[peerID]);
        }

        double clustertime = StopTheClock();
        if (!outputOnFile) {
            cout << "Local clustering done!\n";
            cout << "Time (seconds) required by local clustering: " << clustertime << endl;
        }

        /*** Print info about each peer's clusters***/
        for(int peerID = 0; peerID < params.peers; peerID++) {
            cout << "\npeer: " << peerID << endl;
            printClusters(clusters[peerID], peerID);
        }
        cout << "\n\n";

        unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> all_points;
        unordered_map<array<int, 2>, double, container_hasher> proj;
        mapToTilesPrime(dataset, precision, threshold, ni, proj, all_points);
        printAllPointsClustered(clusters[0], all_points);
    }




    delete[] dimestimate;
    delete[] peerLastItem;
    delete[] converged;
    delete[] dataset;
    delete[] dataset_storage;
    delete[] prevestimate;
    delete[] convRounds;

    return 0;

}

int manhattanDistance(int x1, int x2, int y1, int y2) {
    return abs(x1-x2) + abs(y1-y2);
}

int getCentroids(vectorSet3 &clusters, vector<array<int, 2>> &centroids) {
    if (clusters.size() <= 0) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }
    centroids.clear();
    unSet3::iterator it_tiles;
    double n, x, y;

    /************ for each cluster in clusters ************/
    for (int j = 0; j < clusters.size(); j++) {
        //cout << "Cluster n° " << j << " with size " << cluster.at(j).size() << ": " << endl;
        n = 0, x = 0, y = 0;
        if (!clusters.at(j).size()) {
            cerr << "Bad cluster structure" << endl;
            centroids.clear();
            return dataError(__FUNCTION__);
        }
        it_tiles = clusters.at(j).begin(); // pointer to start of j-th cluster in clusters (cluster = list of tiles, clusters = list of cluster)
        /************ for each tile in cluster j-th ************/
        for (int i = 0; i < clusters.at(j).size(); i++) {  // clusters.at(j).size() represent the number of tiles contained in cluster j-th
            x += (*it_tiles)[0]*(*it_tiles)[2];
            y += (*it_tiles)[1]*(*it_tiles)[2];
            n += (*it_tiles)[2];
            it_tiles++;
        }
        centroids.push_back({(int) (x/n), (int) (y/n)});
    }

    return 0;
}

//remove
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

int simultaneousMaxMin(hashmap &projection, int *maxX, int *minX, int *maxY, int *minY) {
    if (projection.size() <= 0) {
        cerr << "Bad projection data structure" << endl;
        return dataError(__FUNCTION__);
    }
    int x1, x2, y1, y2;

    hashmap::iterator it;
    it = projection.begin();

    if ((projection.size() % 2) != 0) {  // odd case
        *minX = (it -> first)[0];  // set as minX the first value
        *minY = (it -> first)[1];  // set as minY the first value
        *maxX = (it -> first)[0];  // set as maxX the first value
        *maxY = (it -> first)[1];  // set as maxY the first value
        it++;
    } else {                            // even case
        x1 = (it -> first)[0];
        y1 = (it -> first)[1];
        it++;
        x2 = (it -> first)[0];
        y2 = (it -> first)[1];
        it++;
        *minX = (x1 <= x2) ? x1 : x2;
        *maxX = (x1 <= x2) ? x2 : x1;
        *minY = (y1 <= y2) ? y1 : y2;
        *maxY = (y1 <= y2) ? y2 : y1;
    }

    while (it != projection.end()) {
        x1 = (it -> first)[0];
        y1 = (it -> first)[1];
        it++;
        x2 = (it -> first)[0];
        y2 = (it -> first)[1];
        it++;
        // on coordinate x
        if (x1 <= x2){
            if ( x1 < *minX) {
                *minX = x1;
            }

            if ( x2 > *maxX) {
                *maxX = x2;
            }
        } else {
            if ( x2 < *minX) {
                *minX = x2;
            }

            if ( x1 > *maxX) {
                *maxX = x1;
            }
        }
        // on coordinate y
        if (y1 <= y2){
            if ( y1 < *minY) {
                *minY = y1;
            }

            if ( y2 > *maxY) {
                *maxY = y2;
            }
        } else {
            if ( y2 < *minY) {
                *minY = y2;
            }

            if ( y1 > *maxY) {
                *maxY = y1;
            }
        }
    }

    return 0;
}

int getGridSize(int &p, int &q, int peers, int min) {
    q = sqrt(peers*min);
    if ((q * q) != peers*min) { // check if peers*min is a perfect square
        p = q+1;                // if peers+min is not a perfect square chose p as int_sup(sqrt(peers*min))
        q = q * (q + 1) < peers ? ++q : q;
    } else {
        p = q;
    }

    return 0;
}


/**
 * The first while iterates on all neighbor tiles, if a tile is yet present
 * into peer tiles, add the two cardinality, if not insert the tile into
 * peer projection. At the end of the while, peer contains in its hasmap all
 * commons tile with the sums of cardinality, the tiles that was in neighbor
 * but not in peer and the tiles that was in peer but not in neighbor.
 * Now is sufficient with the second while iterate on peer projection
 * and divide by two all tiles cardinality to get the average for each
 * tile in a linear time (include common and uncommon tiles)
 */
int projectionMerge(hashmap &projectionN, hashmap &projectionP){
    if (projectionN.size() <= 0) {
        cerr << "Bad projection data structure" << endl;
        return dataError(__FUNCTION__);
    }

    if (projectionP.size() <= 0) {
        cerr << "Bad projection data structure" << endl;
        return dataError(__FUNCTION__);
    }
    hashmap::iterator itN;
    hashmap::iterator itP;

    itN = projectionN.begin();
    while (itN != projectionN.end()) {
        itP = projectionP.find(itN -> first);
        if (itP != projectionP.end()) { // tile is common
            itP -> second += itN->second;
        } else {
            auto check = projectionP.insert( {itN -> first, itN -> second});
            if (!(check.second)) {
                return insertError(__FUNCTION__);
            }
        }
        itN++;
    }
    // now compute average for each tile ( common or not)
    itP = projectionP.begin();
    while (itP != projectionP.end()) {
        itP -> second /= 2;
        itP++;
    }
    projectionN = projectionP;

    return 0;
}


/**
 * This function use an hashmap to merge in linear time 2 clusters with dimension m and n.
 * A naive method check all m tiles of first cluster with all n tiles of second clusters
 * with a complexity O(n*m) in worst case.
 * If we assume that we can insert and find an element in the hashmap in time O(1) in average case,
 * the first for cycle insert in time O(m) all tiles using as key only the coordinate.
 * The second for cycle try to find the n tiles into the hashmap, it will cost O(n).
 * If is not present add it to the hashmap, also in this case the add operation
 * in the worst case will be O(n).
 */
int clusterUnion(unSet3 &clusterN, unSet3 &clusterP) {
    if (clusterN.size() <= 0) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }
    if (clusterP.size() <= 0) {
        cerr << "Bad clusters data structure" << endl;
        return dataError(__FUNCTION__);
    }
    unSet3::iterator itN = clusterN.begin();
    unSet3::iterator itP = clusterP.begin();
    hashmap::iterator itTmp;

    hashmap tmp;
    unSet3 s;

    array<int, 2> tile;
    int cardinality;

    for (int k = 0; k < clusterN.size(); ++k) {
        tile = { (*itN)[0], (*itN)[1] };
        cardinality =  (*itN)[2];

        auto check = tmp.insert({tile, cardinality});
        if (!(check.second)) {
            return insertError(__FUNCTION__);
        }
        itN++;
    }

    for (int l = 0; l < clusterP.size(); ++l) {
        tile = { (*itP)[0], (*itP)[1] };
        cardinality =  (*itP)[2];

        itTmp = tmp.find(tile);
        if (itTmp != tmp.end()) {
            if (!average(&itTmp -> second, cardinality)) {
                return arithmeticError(__FUNCTION__);
            }
        }

        else {
            auto check = tmp.insert( {tile, cardinality});
            if (!(check.second)) {
                return insertError(__FUNCTION__);
            }
        }
        itP++;
    }

    itTmp = tmp.begin();
    for (int m = 0; m < tmp.size(); m++) {
        auto check = s.insert({ (itTmp ->first)[0], (itTmp ->first)[1], (int) itTmp -> second});
        if (!(check.second)) {
            return insertError(__FUNCTION__);
        }
        itTmp++;
    }
    clusterP.operator=(s);
    clusterN.operator=(s);

    return 0;
}


/**
 * The double for cycle find all common cluster and set to 1 the correspondent index
 * into notCopy structure. For example if the k-th cluster of clusterN is equal at the
 * l-th clusterP, then notCopyN[k] and notCopyP[l] will be set to 1.
 * In this way the nex for cycle using these information will copy only the cluster
 * that are not common.
 */
int clustersMerge(vectorSet3 &clustersN, vectorSet3 &clustersP, vector<array<int, 2>> &centroidsN, vector<array<int, 2>> &centroidsP, int maxDistance) {
    int tmpSize = centroidsP.size(); /// get the actual dimension because it will increase but we are not interested at the clusters that will be added at the end
    int returnValue = -1;
    int *notCopyP = nullptr;  /**!< array of int, if notCopyP[i] is set to 1, means that i-th cluster of clusterP is also present in clusterN */
    int *notCopyN = nullptr;  /**!< array of int, if notCopyP[i] is set to 1, means that i-th cluster of clusterN is also present in clusterP */

    notCopyP = new (nothrow) int[centroidsP.size()]();
    if(!notCopyP)
        return memoryError(__FUNCTION__);


    notCopyN = new (nothrow) int[centroidsN.size()]();
    if(!notCopyN){
        returnValue = memoryError(__FUNCTION__);
        goto ON_EXIT;
    }

    for (int k = 0; k < centroidsN.size(); k++) {
        for (int l = 0; l < tmpSize; l++) {
            if (manhattanDistance(centroidsN.at(k)[0], centroidsP.at(l)[0], centroidsN.at(k)[1], centroidsP.at(l)[1]) <= maxDistance) {
                notCopyP[l] = 1, notCopyN[k] = 1;
                if ( !(clustersP.at(l) == clustersN.at(k))) {
                    if (!clusterUnion(clustersN.at(k), clustersP.at(l))) {
                        cerr << "Cluster Union error" << endl;
                        returnValue = -25;
                        goto ON_EXIT;
                    }

                }
                break;
            }
        }
    }

    for (int k = 0; k < clustersN.size(); k++) {
        if (!notCopyN[k]) {
            clustersP.push_back(clustersN.at(k));
            centroidsP.push_back({centroidsN.at(k)[0], centroidsN.at(k)[1]});
        }
    }
    for (int l = 0; l < tmpSize; l++) {
        if (!notCopyP[l]) {
            clustersN.push_back(clustersP.at(l));
            centroidsN.push_back({centroidsP.at(l)[0], centroidsP.at(l)[1]});
        }
    }

    returnValue = 0;

    ON_EXIT:
    if (notCopyN != nullptr)
        delete[] notCopyN;

    if (notCopyP != nullptr)
        delete[] notCopyP;

    return returnValue;

}

int average(double *x, double y)  {
    if (x == nullptr) {
        cerr << "x cannot be null!" << endl;
        return arithmeticError(__FUNCTION__);
    } else {
        *x = (*x+y)/2.0;
    }

    return 0;
}

int readDataset(double **dataset_storage, double ***dataset, string name_file, int &ni) {
    int row, column;
    if(getDim(name_file, row, column)) {
        return fileError(__FUNCTION__);
    }
    ni = row;

    *dataset_storage = new (nothrow) double[row*column];
    if (!(*dataset_storage)){
        return memoryError(__FUNCTION__);
    }


    (*dataset) = new (nothrow) double*[row];
    if (!(*dataset)){
        delete[] *dataset_storage, *dataset_storage = nullptr;
        return memoryError(__FUNCTION__);
    }

    for (int i = 0; i < row; i++) {
        (*dataset)[i] = &(*dataset_storage)[i*column];
        if (!(*dataset)[i]) {
            cerr << "Not enough memory allocating the " << i << "-th tile_points array" << endl;
            for (int x = 0; x < i; x++) {
                delete[] (*dataset)[x], (*dataset)[x] = nullptr;
            }
            delete[] *dataset, *dataset = nullptr;
            return memoryError(__FUNCTION__);
        }
    }

    if(loadData(*dataset, name_file, column)) {
        delete[] *dataset_storage, *dataset_storage = nullptr;
        for (int i = 0; i < row; i++) {
            delete[] (*dataset)[i], (*dataset)[i] = nullptr;
        }
        delete[] *dataset, *dataset = nullptr;

        return fileError(__FUNCTION__);
    }


    return 0;
}


int partitionDataset(long *peerLastItem, int peers, int ni) {
    std::random_device rd; /** obtain a random number from hardware */
    std::mt19937 eng(rd()); /** seed the generator */
    std::uniform_real_distribution<> distr(-1, 1); /** define the range */

    for(int i = 0; i < peers - 1; i++){
        float rnd = distr(eng);
        //cerr << "rnd: " << rnd << "\n";
        long last_item = rnd * ((float)ni/(float)peers) * 0.1 + (float) (i+1) * ((float)ni/(float)peers) - 1;
        peerLastItem[i] = last_item;
    }

    peerLastItem[peers - 1] = ni-1;
    return 0;
}

int checkPartitioning(long *peerLastItem, int peers, int ni) {
    long sum = peerLastItem[0] + 1;
    //cerr << "peer 0:" << sum << "\n";
    for(int i = 1; i < peers; i++) {
        sum += peerLastItem[i] - peerLastItem[i-1];
        //cerr << "peer " << i << ":" << peerLastItem[i] - peerLastItem[i-1] << "\n";
    }

    if(sum != ni) {
        cerr << "ERROR: ni = "<< ni << "!= sum = " << sum << endl;
        return partitionError(__FUNCTION__);
    }
    return 0;
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
    igraph_grg_game(&G_graph, n, radius, 0, nullptr, nullptr);

    igraph_is_connected(&G_graph, &connected, IGRAPH_WEAK);
    while(!connected){
        igraph_destroy(&G_graph);
        igraph_grg_game(&G_graph, n, radius, 0, nullptr, nullptr);

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

void parametersSummary(Params params) {
    cout << "\n\nPARAMETERS:\n";
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

int generateGraph(igraph_t &graph, int peers, int graphType) {
    // turn on attribute handling in igraph
    igraph_i_set_attribute_table(&igraph_cattribute_table);

    // seed igraph PRNG
    igraph_rng_seed(igraph_rng_default(), 42);

    // generate a connected random graph
    graph = generateRandomGraph(graphType, peers);
    return 0;
}

int graphProperties(igraph_t &graph) {
    igraph_vector_t result;
    igraph_real_t mindeg;
    igraph_real_t maxdeg;

    igraph_vector_init(&result, 0);
    igraph_degree(&graph, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    igraph_vector_minmax(&result, &mindeg, &maxdeg);

    cout << "Minimum degree is "  <<  mindeg << ", Maximum degree is " << maxdeg << "\n";

    return 0;
}

int getNeighborsSize(igraph_t &graph, igraph_vector_t &neighbors, int peerID, int fanOut) {
    long neighborsSize = igraph_vector_size(&neighbors);
    if(fanOut < neighborsSize && fanOut != -1){
        // randomly sample f adjacent vertices
        igraph_vector_shuffle(&neighbors);
        igraph_vector_remove_section(&neighbors, fanOut, neighborsSize);
    }

    neighborsSize = igraph_vector_size(&neighbors);
    return neighborsSize;
}

int checkFirstConvegence(double prevestimate, double dimestimate, double convThreshold, int &convRounds, bool &converged, int convLimit, int &Numberofconverged) {
    bool dimestimateconv;
    if(prevestimate)
        dimestimateconv = fabs((prevestimate - dimestimate) / prevestimate) < convThreshold;
    else
        dimestimateconv = false;

    if(dimestimateconv)
        convRounds++;
    else
        convRounds = 0;

    converged = (convRounds >= convLimit);
    if(converged){
        Numberofconverged --;
    }
}

int restoreCardinality(hashmap &projection, int dimestimate) {
    hashmap::iterator it;
    it = projection.begin();
    while (it != projection.end()) {
        it->second = round(it->second * dimestimate);
        it++;
    }
}

int getGroupsTiles( int &nTilesX, int &nTilesY, int &restTilesX, int &restTilesY,
                    int q, int p, int minX, int minY, int maxX, int maxY ) {
    if (!p)
        return arithmeticError(__FUNCTION__);

    if (!q)
        return arithmeticError(__FUNCTION__);

    nTilesY = (maxY - minY) / q;
    nTilesX = (maxX - minX) / p;
    restTilesY = (maxY - minY) % q;
    restTilesX = (maxX - minX) % p;

    return 0;

}
int getSquares(int &squares, int peerID, int p, int q, int dimestimate, int minSquares) {
    if (peerID < (p * q) % (int) dimestimate ) {
        squares = minSquares+1; // extra square for some peers
    } else
        squares = minSquares;

    return 0;
}

int getSquaresCoordinate(int &x1, int &y1, int &y2, int &x2, int p, int q, int dimestimate,
                         int minX, int minY, int maxX, int maxY, int k, int peerID) {
    int i, j;
    int nTilesY, nTilesX, restTilesY, restTilesX;
    /// compute how group tiles for each square on y and x axis
    getGroupsTiles(nTilesX, nTilesY, restTilesX, restTilesY,
                   q, p, minX, minY, maxX, maxY);
    i = (peerID + k*(int)dimestimate) /p + 1;
    j = (peerID + k*(int)dimestimate) %p + 1;

    y1 = i <= restTilesY ? i * (nTilesY + 1) + minY : restTilesY * (nTilesY + 1) + (i - restTilesY) * nTilesY + minY;
    x1 = j <= restTilesX ? j * (nTilesX + 1) + minX : restTilesX * (nTilesX + 1) + (j - restTilesX) * nTilesX + minX;
    y2 = (i <= restTilesY && restTilesY != 0) ? y1 - (nTilesY + 1) : y1 - (nTilesY);
    x2 = (j <= restTilesX && restTilesX != 0) ? x1 - (nTilesX + 1) : x1 - (nTilesX);
    /// include the left and bottom edge of grid
    if (y2 == minY)
        y2 -= 1;

    if (x2 == minX)
        x2 -= 1;
    //cout << "x1: " << y1[peerID][k]<< " x2: " << x1[peerID][k] << endl;

    return 0;
}


