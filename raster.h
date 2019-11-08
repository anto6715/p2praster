/** @file*/

/**
* @author Antonio Mariani
* @date 05/10/19
*/


#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <cstring>
#include <boost/functional/hash.hpp>
#include <unordered_map>
#include <chrono>
#include <unordered_set>

using namespace std;
#ifndef P2PRASTER_RASTER_H
#define P2PRASTER_RASTER_H

/**
 *
 * @param [in] seed
 * @param [in] value
 */
void hash_combine(std::size_t& seed, std::size_t value);

/**
 * Struct used to obtain hash for data structure with key: array<int, n>
 */
struct container_hasher {
    template<class T>
    std::size_t operator()(const T& c) const {
        std::size_t seed = 0;
        for(const auto& elem : c) {
            hash_combine(seed, std::hash<typename T::value_type>()(elem));
        }
        return seed;
    }
};

/**
 *
 * @param [in] name_file - The input file name
 * @param [in,out] row - Number of csv rows
 * @param [in,out] column - Number of csv columns
 * @return Save in row, column the number of csv rows and columns, return 0 if success, -1 in case of error
 */
int getDim(string name_file, int &row, int &column);

/**
 *
 * @param [in,out] m - Matrix where load dataset
 * @param [in] name_file - Name of dataset file
 * @param [in] column - Number of file column
 * @return Save in m the dataset loaded from name_file, return 0 if success, -1 in case of error
 */
int loadData(double **m, string name_file, int column);

/**
 * This functions uses a simple for cycle in order to obtains the key of the neighbors of coordinate in input,
 * then search this neighbors into projection and return the results into result
 *
 * @param [in] coordinate - Coordinate of tile of which we want to obtain its neighbors
 * @param [in,out] projection - Data structure contains all tiles where search the tile's neighbors
 * @param [in,out] result - Contains all tile's neighbors found into projection
 * @result Move the neighbors found from projection to result
 */
void getNeighbors(array<int, 2> coordinate, unordered_map<array<int, 2>, double, container_hasher> &projection, unordered_set<array<int, 2>, container_hasher> &result);

/**
 * This functions uses a simple for cycle in order to obtains the key of the neighbors of coordinate in input,
 * then search this neighbors into projection and squareProjection and return the results into result
 *
 * @param [in] coordinate - Coordinate of tile of which we want to obtain its neighbors and its cardinality
 * @param [in,out] squareProjection - Data structure contains all peer's tiles where search the tile's neighbors
 * @param [in,out] projection - Data structure contains all remaining tiles where search the tile's neighbors
 * @param [in,out] result - Contains all tile's neighbors found into projection or squareProjection
 * @result Move the neighbors found from projection or squareProjection to result
 */
void getNeighbors(array<int, 3> coordinate, unordered_map<array<int, 2>, double, container_hasher> &squareProjection, unordered_map<array<int, 2>, double, container_hasher> &projection, unordered_set<array<int, 3>, container_hasher> &result);

/**
 * This function performs the first step of raster projection, it simple multiply
 * the points for a scalar and maintains only the integer part of the product,
 * Note: this function don't apply the threshold to the tiles
 *
 * @param [in] m - Dataset
 * @param [in] precision - Parameter for raster algorithm
 * @param [in] threshold - Parameter for raster algorithm
 * @param [in,out] projection - Data structure where save the projection result of raster algorithm
 * @param [in] start - First point in the dataset of peer
 * @param [in] end - Last point in the dataset of peer
 * @result Save into projection the tiles with its cardinality
 */
void mapToTilesNoThreshold(double **m, double precision, unordered_map<array<int, 2>, double, container_hasher> &projection, int start, int end);

/**
 * This function implements the raster' variant of the algorithm.
 * During the projection phase, using an hashmap keeps track of the
 * associations tile-point using the tile as key and a list the relative
 * associated points as value
 *
 * @param [in] m - Dataset
 * @param [in] precision - Parameter for raster algorithm (it is the value for the approximation)
 * @param [in] threshold - Parameter for raster algorithm (all tiles with cardinality < threshold are delete)
 * @param [in] n - Number of points in Dataset
 * @param [in,out] projection - Data structure where save the projection result of raster algorithm
 * @param [in,out] all_points - Data structure where save the mapping tiles-points
 */
void mapToTilesPrime(double **m, double precision, int threshold, int n, unordered_map<array<int, 2>, double, container_hasher> &projection, unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> &all_points);

/**
 * This function take one-by-one each tile in projection, try to find its neighbor
 * (always in projection) and form a cluster with these adjacent tiles.
 * If the cluster has a number of tiles < min_size it is discarded
 *
 *
 * @param [in,out] projection - Data structure with tiles
 * @param [in] min_size - Parameter for raster algorithm (it is the minimum number of tiles in order to form a cluster)
 * @param [in,out] clusters - Data structure where save the clusters found
 * @result Save into clusters the clusters found
 */
void clusteringTiles(unordered_map<array<int, 2>, double, container_hasher> &projection, int min_size, vector<unordered_set<array<int, 2>, container_hasher>> &clusters);

/**
 *
 * @param [in,out] squareProjection - Data structure with peer's tiles
 * @param [in,out] projection - Data structure with remaining tiles
 * @param [in] min_size - Parameter for raster algorithm (it is the minimum number of tiles in order to form a cluster)
 * @param [in,out] clusters - Data structure where save the clusters found
 * @result Save into clusters the clusters found
 */
void clusteringTiles(unordered_map<array<int, 2>, double, container_hasher> &squareProjection, unordered_map<array<int, 2>, double, container_hasher> &projection, int min_size, vector<unordered_set<array<int, 3>, container_hasher>> &clusters);

/**
 * this function print on terminal the total number of tiles clsutered,
 * if outputOnFile = true print also information about cluster and their
 * tiles in a csv file name with peerId.csv
 *
 * @param [in] clusters - Data structure which contains clusters to print
 * @param [in] precision - Parameter for raster algorithm
 * @param [in] peerID - (optional) Peer id that wnat to write on file
 * @param [in] outputOnFile - (optional) True if want to write a csv file
 */
void printClusters(vector<unordered_set<array<int, 2>, container_hasher>> &clusters, int peerID);

/**
 * This function differ from above function in cluster structure,
 * here clusters maintain information about tiles cardinality
 *
 * @param [in] clusters - Data structure which contains clusters to print
 * @param [in] precision - Parameter for raster algorithm
 * @param [in] peerID - (optional) Peer id that wnat to write on file
 * @param [in] outputOnFile - (optional) True if want to write a csv file
 */
void printClusters(vector<unordered_set<array<int, 3>, container_hasher>> &clusters, int peerID);

/**
 * This function using a data structure from raster' variant, print all clusters
 * with all relative point. All points not clustered are assign at the cluster 0.
 *
 * @param [in] clusters - Data structure which contains clusters to print
 * @param all_points - Data structure which contains mapping tile-points for raster prime' variant
 */
void printAllPointsClustered(vector<unordered_set<array<int, 2>, container_hasher>> &clusters, unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> &all_points);

/**
 * This function differ from above function in cluster structure,
 * here clusters maintain information about tiles cardinality
 *
 * @param [in] clusters - Data structure which contains clusters to print
 * @param all_points - Data structure which contains mapping tile-points for raster prime' variant
 */
void printAllPointsClustered(vector<unordered_set<array<int, 3>, container_hasher>> &clusters, unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> &all_points);

// Remove?
/**
 *
 * @param [in] clusters - Data structure which contains clusters
 * @param [in] all_points - Data structure which contains mapping tile-points for raster prime' variant
 * @param [in] precision - Parameter for raster algorithm
 * @Result Print on terminal the mean Shannon entropy and mean Density for each cluster
 */
void analyzeClusters(vector<unordered_set<array<int, 2>, container_hasher>> &clusters, unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> &all_points, double precision);// raster'


#endif //P2PRASTER_RASTER_H