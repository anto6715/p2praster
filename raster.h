//
// Created by antonio on 05/10/19.
//

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

void hash_combine(std::size_t& seed, std::size_t value);

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

int loadData(double **m, const string& name_file, int column);

int getDim(const string& name_file, int &row, int &column);

// unordered
void getNeighbors(array<int, 2> coordinate, unordered_map<array<int, 2>, double, container_hasher> &squareProjection, unordered_map<array<int, 2>, double, container_hasher> &projection, unordered_set<array<int, 2>, container_hasher> &result);

void mapToTiles(double **m, double precision, int threshold, unordered_map<array<int, 2>, int, container_hasher> &projection, int start, int end);

void mapToTilesNoThreshold(double **m, double precision, int threshold, unordered_map<array<int, 2>, double, container_hasher> &projection, int start, int end);

void mapToTilesPrime(double **m, double precision, int threshold, unordered_map<array<int, 2>, int, container_hasher> &projection, unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> &all_points, int start, int end);

void clusteringTiles(unordered_map<array<int, 2>, double, container_hasher> &squareProjection, unordered_map<array<int, 2>, double, container_hasher> &projection, int min_size, vector<unordered_set<array<int, 2>, container_hasher>> &clusters);

void printClusters(vector<unordered_set<array<int, 2>, container_hasher>> &clusters, double precision, int peerID = 0);

void printAllPointsClustered(vector<unordered_set<array<int, 2>, container_hasher>> &clusters, unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> &all_points);

void analyzeClusters(vector<unordered_set<array<int, 2>, container_hasher>> &clusters, unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> &all_points, double precision);// raster'


#endif //P2PRASTER_RASTER_H