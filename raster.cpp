#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include "raster.h"

using namespace std;
using namespace std::chrono;

void hash_combine(std::size_t& seed, std::size_t value) {
    seed ^= value + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

int getDim(string name_file, int &row, int &column) {
    // Initialization
    row = 0;
    column = 0;
    ifstream inputFile(name_file);
    while (inputFile) {
        string s;
        if (!getline(inputFile, s)) break;
        if (row == 0) {
            istringstream ss(s);
            while (ss) {
                string line;
                if (!getline(ss, line, ','))
                    break;
                column++;
            }
        }
        row++;
    }
    if (!inputFile.eof()) {
        cerr << "Could not read file " << name_file << "\n";
        //__throw_invalid_argument("File not found.");
        return -1;
    }
    return 0;
}

int loadData(double **m, string name_file, int column) {
    ifstream inputFile(name_file);
    int row = 0;
    while (inputFile) {
        string s;
        if (!getline(inputFile, s)) break;
        if (s[0] != '#') {
            istringstream ss(s);
            while (ss) {
                for (int i = 0; i < column; i++) {
                    string line;
                    if (!getline(ss, line, ','))
                        break;
                    try {
                        m[row][i] = stod(line);
                    }
                    catch (const std::invalid_argument &e) {
                        cout << "NaN found in file " << name_file << " line " << row
                             << endl;
                        e.what();
                    }
                }
            }
        }
        row++;
    }
    if (!inputFile.eof()) {
        cerr << "Could not read file " << name_file << "\n";
        //__throw_invalid_argument("File not found.");
        return -1;
    }
    return 0;
}

void getNeighbors(array<int, 2> coordinate, unordered_map<array<int, 2>, double, container_hasher> &projection, unordered_set<array<int, 2>, container_hasher> &result) {
    int x = coordinate[0];
    int y = coordinate[1];
    int radius = 1;
    int count;
    unordered_map<array<int, 2>, double, container_hasher>::iterator it;
    unordered_set<array<int, 2>, container_hasher> neighbors;
    unordered_set<array<int, 2>, container_hasher>::iterator it_neighbor;

    // neighboring generation of coordinate
    count = 0;
    for (int i = -radius; i <= radius ; i++) {
        for (int j = -radius; j <= radius ; j++) {
            if (i == 0 && j == 0)
                continue;
            neighbors.insert({x + i, y + j});
            count++;

        }
    }

    // if a neighbor is present in tiles, add it in result
    it_neighbor = neighbors.begin();
    for (int i = 0; i < neighbors.size(); i++) {
        it = projection.find(*it_neighbor);
        it_neighbor++;
        if (it != projection.end()) {
            result.insert(it -> first);
            projection.erase(it++);
        }
    }
}

void getNeighbors(array<int, 3> coordinate, unordered_map<array<int, 2>, double, container_hasher> &squareProjection, unordered_map<array<int, 2>, double, container_hasher> &projection, unordered_set<array<int, 3>, container_hasher> &result) {
    int x = coordinate[0];
    int y = coordinate[1];
    int radius = 1;
    int count;
    unordered_map<array<int, 2>, double, container_hasher>::iterator it;
    unordered_set<array<int, 2>, container_hasher> neighbors;
    unordered_set<array<int, 2>, container_hasher>::iterator it_neighbor;

    // neighboring generation of coordinate
    count = 0;
    for (int i = -radius; i <= radius ; i++) {
        for (int j = -radius; j <= radius ; j++) {
            if (i == 0 && j == 0)
                continue;
            neighbors.insert({x + i, y + j});
            count++;

        }
    }

    // if a neighbor is present in tiles, add it in result
    it_neighbor = neighbors.begin();
    for (int i = 0; i < neighbors.size(); i++) {
        it = squareProjection.find(*it_neighbor);
        if (it != squareProjection.end()) {
            result.insert({(it->first)[0], (it->first)[1],(int) it->second});
            squareProjection.erase(it++);
        } else {
            it = projection.find(*it_neighbor);
            if (it != projection.end()) {
                result.insert({(it->first)[0], (it->first)[1],(int) it->second});
                projection.erase(it++);
            }
        }


        it_neighbor++;
    }
}

void mapToTilesNoThreshold(double **m, double precision, unordered_map<array<int, 2>, double, container_hasher> &projection, int start, int end) {
    double scalar;
    unordered_map<array<int, 2>, double, container_hasher>::iterator it;
    scalar = pow(10, precision);
    int lat, lon;
    for (int i = start; i <= end; i++) {
        lat =(int) (m[i][0] * scalar);
        lon =(int) (m[i][1] * scalar);
        array<int, 2> tile = {lat,lon};
        it = projection.find(tile);
        if (it != projection.end()) {
            it->second++;
        } else {
            projection[tile] = 1.0;
        }
    }
}

void mapToTilesPrime(   double **m,
                        double precision,
                        int threshold,
                        int n,
                        unordered_map<array<int, 2>, double, container_hasher> &projection,
                        unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> &all_points) {
    double scalar;
    unordered_map<array<int, 2>, double, container_hasher>::iterator it;
    unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher>::iterator it_map_all_points;
    scalar = pow(10, precision);
    int lat, lon;
    for (int i = 0; i < n; i++) {
        lat =(int) (m[i][0] * scalar);
        lon =(int) (m[i][1] * scalar);
        array<int, 2> tile = {lat,lon};
        it = projection.find(tile);
        if (it != projection.end()) {
            it->second++;
        } else {
            projection[tile] = 1;
        }

        it_map_all_points = all_points.find(tile);
        if (it_map_all_points != all_points.end()) {
            (it_map_all_points -> second).insert({m[i][0],m[i][1]});
        } else {
            unordered_set<array<double , 2>, container_hasher> point;
            point.insert({m[i][0], m[i][1]});
            all_points[tile] = point;
        }
    }
    // remove tile with count < threshold
    it = projection.begin();
    while (it != projection.end()) {
        if (it -> second < threshold) {
            projection.erase(it++);
        } else {
            it++;
        }
    }

}

void clusteringTiles(unordered_map<array<int, 2>, double, container_hasher> &projection, int min_size, vector<unordered_set<array<int, 2>, container_hasher>> &clusters) {
    unordered_map<array<int, 2>, double, container_hasher>::iterator iterator;

    // read all tiles
    while ((iterator = projection.begin()) != projection.end()) {
        // read and remove first element recursively
        array<int, 2> x = iterator->first;
        projection.erase(iterator++);

        unordered_set<array<int, 2>, container_hasher> visited;
        visited.insert(x);

        // get neighbors of tile in exam
        unordered_set<array<int, 2>, container_hasher> to_check;
        getNeighbors(x, projection, to_check);

        // for each neighbor, try to find respectively neighbor in order to add they to single cluster
        while (!to_check.empty()) {
            array<int, 2> value = *to_check.begin() ;
            to_check.erase((to_check.begin()));
            visited.insert(value);

            unordered_set<array<int, 2>, container_hasher> temp;

            getNeighbors(value, projection, temp);
            while (!temp.empty()) {
                to_check.insert(*temp.begin());
                temp.erase((temp.begin()));
            }
        }
        // validate visited as cluster if is size is >= min size
        if (visited.size() >= min_size) {
            clusters.push_back(visited);
        }

    }
} // unordered

void clusteringTiles(unordered_map<array<int, 2>, double, container_hasher> &squareProjection,unordered_map<array<int, 2>, double , container_hasher> &projection, int min_size, vector<unordered_set<array<int, 3>, container_hasher>> &clusters) {
    unordered_map<array<int, 2>, double, container_hasher>::iterator iterator;

    // read all tiles
    while ((iterator = squareProjection.begin()) != squareProjection.end()) {
        // read and remove first element recursively
        array<int, 3> x;
        x = {(iterator->first)[0], (iterator->first)[1],(int) iterator->second};
        squareProjection.erase(iterator++);

        unordered_set<array<int, 3>, container_hasher> visited;
        visited.insert(x);

        // get neighbors of tile in exam
        unordered_set<array<int, 3>, container_hasher> to_check;
        getNeighbors(x, squareProjection, projection, to_check);

        // for each neighbor, try to find respectively neighbor in order to add they to single cluster
        while (!to_check.empty()) {
            array<int, 3> value = *to_check.begin() ;
            to_check.erase((to_check.begin()));
            visited.insert(value);

            unordered_set<array<int, 3>, container_hasher> temp;

            getNeighbors(value, squareProjection, projection, temp);
            while (!temp.empty()) {
                to_check.insert(*temp.begin());
                temp.erase((temp.begin()));
            }
        }
        // validate visited as cluster if is size is >= min size
        if (visited.size() >= min_size) {
            clusters.push_back(visited);
        }

    }
}

void printClusters(vector<unordered_set<array<int, 2>, container_hasher>> &clusters, int peerID) {
    ofstream outfile(to_string(peerID) +".csv");
    cout <<  "n° cluster: " << clusters.size() << endl;
    unordered_set<array<int, 2>, container_hasher>::iterator it;

    int count_tiles = 0;
    for (int j = 0; j < clusters.size(); j++) {
        outfile << "Cluster n° " << j << " with size " << clusters.at(j).size() << ": " << endl;
        it = clusters.at(j).begin(); // pointer to start of j-th cluster (cluster = list of tiles)
        for (int i = 0; i < clusters.at(j).size(); i++) {
            count_tiles++; // count the total number of tiles clustered
            outfile << (*it)[0] << ",";
            outfile << (*it)[1] << ",";

            outfile << j << endl;
            it++; // next tile of the actual cluster
        }
    }
    outfile.close();
    cout << "Tiles clustered: " << count_tiles << endl;
} // semi-unordered

void printClusters(vector<unordered_set<array<int, 3>, container_hasher>> &clusters, int peerID) {
    ofstream outfile(to_string(peerID) +".csv");
    cout <<  "n° cluster: " << clusters.size() << endl;
    unordered_set<array<int, 3>, container_hasher>::iterator it;

    int count_tiles = 0;
    for (int j = 0; j < clusters.size(); j++) {
        outfile << "Cluster n° " << j << " with size " << clusters.at(j).size() << ": " << endl;
        it = clusters.at(j).begin(); // pointer to start of j-th cluster (cluster = list of tiles)
        for (int i = 0; i < clusters.at(j).size(); i++) {
            count_tiles++; // count the total number of tiles clustered
            outfile << (*it)[0] << ",";
            outfile << (*it)[1] << ",         ";
            //outfile << (*it)[2] << ","; // tile cardinality

            outfile << j << endl;
            it++; // next tile of the actual cluster
        }
    }
    outfile.close();
    cout << "Tiles clustered: " << count_tiles << endl; // on terminal

}

void printAllPointsClustered(   vector<unordered_set<array<int, 3>, container_hasher>> &clusters,
                                unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> &all_points){
    cout.precision(15);
    ofstream outfile("clustered.csv");
    int count_not_clustered = 0;
    int count_tiles = 0;
    int count_points = 0;
    unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher>::iterator it_map_all_points;
    unordered_set<array<double , 2>, container_hasher>::iterator it_set_all_points;
    unordered_set<array<int, 3>, container_hasher>::iterator it_tiles;


    /************ for each cluster in clusters ************/
    for (int j = 0; j < clusters.size(); j++) {
        //cout << "Cluster n° " << j << " with size " << cluster.at(j).size() << ": " << endl;
        it_tiles = clusters.at(j).begin(); // pointer to start of j-th cluster in clusters (cluster = list of tiles, clusters = list of cluster)
        /************ for each tile in cluster j-th ************/
        for (int i = 0; i < clusters.at(j).size(); i++) {
            it_map_all_points = all_points.find({(*it_tiles)[0], (*it_tiles)[1]}); // try to find in all_points the tile (with its list of points) from cluster
            if (it_map_all_points != all_points.end()) {
                it_set_all_points = (it_map_all_points -> second).begin(); // pointer to the first element in the list of points associated to the founded tile
                /************ for each point in the tile ************/
                for (int k = 0; k < (it_map_all_points -> second).size(); k++) {
                    outfile << (*it_set_all_points)[0] << ",";
                    outfile << (*it_set_all_points)[1] << ",";
                    outfile << j + 1 <<endl;
                    //cout << (*it_set_all_points)[0] << ",";
                    //cout << (*it_set_all_points)[1] << ",";
                    //cout << j <<endl;
                    it_set_all_points++;
                    count_points++;
                }
                all_points.erase(it_map_all_points++);
            }
            it_tiles++;
            count_tiles++;
        }
    }
    it_map_all_points = all_points.begin(); // first tile remaining in all_points
    /************ for each tile that are not in the clusters ************/
    for (int i = 0; i < all_points.size(); i++) {
        if (it_map_all_points != all_points.end()) {
            it_set_all_points = (it_map_all_points -> second).begin(); // pointer to the first element in the list of points associated to the founded tile
            /************ for each point in the tiles that are not in the clusters ************/
            for (int k = 0; k < (it_map_all_points -> second).size(); k++) {
                outfile << (*it_set_all_points)[0] << ",";
                outfile << (*it_set_all_points)[1] << ",";
                outfile << 0 <<endl;
                it_set_all_points++;
                count_not_clustered++;
            }
            //all_points.erase(it_map_all_points++);
            it_map_all_points++;
        }
    }
    outfile.close();
    cout << "Total points not clustered: " << count_not_clustered << endl;
    cout << "Tile not clustered: " << all_points.size() << endl;
    cout << "Tiles clustered: " << count_tiles << endl;
    cout << "Clusters: " << clusters.size() << endl;
    cout << "Points clustered: " << count_points << endl;
    cout << "Points analyzed: " << count_points + count_not_clustered << endl;
}

void printAllPointsClustered(   vector<unordered_set<array<int, 2>, container_hasher>> &clusters,
                                unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> &all_points){
    cout.precision(15);
    ofstream outfile("clustered.csv");
    int count_not_clustered = 0;
    int count_tiles = 0;
    int count_points = 0;
    unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher>::iterator it_map_all_points;
    unordered_set<array<double , 2>, container_hasher>::iterator it_set_all_points;
    unordered_set<array<int, 2>, container_hasher>::iterator it_tiles;


    /************ for each cluster in clusters ************/
    for (int j = 0; j < clusters.size(); j++) {
        //cout << "Cluster n° " << j << " with size " << cluster.at(j).size() << ": " << endl;
        it_tiles = clusters.at(j).begin(); // pointer to start of j-th cluster in clusters (cluster = list of tiles, clusters = list of cluster)
        /************ for each tile in cluster j-th ************/
        for (int i = 0; i < clusters.at(j).size(); i++) {
            it_map_all_points = all_points.find((*it_tiles)); // try to find in all_points the tile (with its list of points) from cluster
            if (it_map_all_points != all_points.end()) {
                it_set_all_points = (it_map_all_points -> second).begin(); // pointer to the first element in the list of points associated to the founded tile
                /************ for each point in the tile ************/
                for (int k = 0; k < (it_map_all_points -> second).size(); k++) {
                    outfile << (*it_set_all_points)[0] << ",";
                    outfile << (*it_set_all_points)[1] << ",";
                    outfile << j + 1 <<endl;
                    //cout << (*it_set_all_points)[0] << ",";
                    //cout << (*it_set_all_points)[1] << ",";
                    //cout << j <<endl;
                    it_set_all_points++;
                    count_points++;
                }
                all_points.erase(it_map_all_points++);
            }
            it_tiles++;
            count_tiles++;
        }
    }
    it_map_all_points = all_points.begin(); // first tile remaining in all_points
    /************ for each tile that are not in the clusters ************/
    for (int i = 0; i < all_points.size(); i++) {
        if (it_map_all_points != all_points.end()) {
            it_set_all_points = (it_map_all_points -> second).begin(); // pointer to the first element in the list of points associated to the founded tile
            /************ for each point in the tiles that are not in the clusters ************/
            for (int k = 0; k < (it_map_all_points -> second).size(); k++) {
                outfile << (*it_set_all_points)[0] << ",";
                outfile << (*it_set_all_points)[1] << ",";
                outfile << 0 <<endl;
                it_set_all_points++;
                count_not_clustered++;
            }
            //all_points.erase(it_map_all_points++);
            it_map_all_points++;
        }
    }
    outfile.close();
    cout << "Total points not clustered: " << count_not_clustered << endl;
    cout << "Tile not clustered: " << all_points.size() << endl;
    cout << "Tiles clustered: " << count_tiles << endl;
    cout << "Clusters: " << clusters.size() << endl;
    cout << "Points clustered: " << count_points << endl;
    cout << "Points anallized: " << count_points + count_not_clustered << endl;
}

void analyzeClusters(vector<unordered_set<array<int, 2>, container_hasher>> &clusters,
                    unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher> &all_points,
                    double precision) {
    double scalar = pow(10, precision);
    double area = scalar * scalar;
    double* shannon = new double[clusters.size()];
    double* density = new double[clusters.size()];
    double* size_cluster = new double[clusters.size() * sizeof(double)];
    double pji; // probability that a point is in the i-th tile of the j-th cluster
    unordered_map<array<int, 2>, unordered_set<array<double , 2>, container_hasher>, container_hasher>::iterator it_map_all_points;
    unordered_set<array<double , 2>, container_hasher>::iterator it_set_all_points;
    unordered_set<array<int, 2>, container_hasher>::iterator it_tiles;

    double** tile_points = new double*[clusters.size()];
    // allocating vectors that will contain n° of points for each tile i of respective cluster j
    for (int j = 0; j < clusters.size(); j++) {
        //cout << "Cluster n° " << j << " with size " << clusters.at(j).size() << ": " << endl;
        tile_points[j] = new double[clusters.at(j).size()];
        if (!tile_points[j]) {
            cout << "error at: " << j << endl;
        }
    }

    /************ for each cluster in clusters ************/
    for (int j = 0; j < clusters.size(); j++) {
        //cout << "Cluster n° " << j << " with size " << cluster.at(j).size() << ": " << endl;
        it_tiles = clusters.at(j).begin(); // pointer to start of j-th cluster in clusters (cluster = list of tiles, clusters = list of cluster)
        /************ for each tile in cluster j-th ************/
        for (int i = 0; i < clusters.at(j).size(); i++) {  // clusters.at(j).size() represent the number of tiles contained in cluster j-th
            it_map_all_points = all_points.find((*it_tiles)); // try to find in all_points the tile (with its list of points) from cluster
            size_cluster[j] += (it_map_all_points -> second).size();
            tile_points[j][i] = (it_map_all_points -> second).size();
            it_tiles++;
        }
    }
    /************ for each cluster in clusters ************/
    for (int j = 0; j < clusters.size(); j++) {
        /************ for each tile in cluster j-th ************/
        for (int i = 0; i < clusters.at(j).size(); i++) {  // clusters.at(j).size() represent the number of tiles contained in cluster j-th
            pji = tile_points[j][i]/size_cluster[j];
            shannon[j] += -(pji * log2(pji));
        }
    }

    // calculating density for each cluster
    /************ for each cluster in clusters ************/
    for (int j = 0; j < clusters.size(); j++) {
        density[j] = size_cluster[j]/(clusters.at(j).size() * area);
    }

    double mean_shannon = 0;
    double mean_density = 0;
    for (int j = 0; j < clusters.size(); j++) {
        mean_shannon += shannon[j];
        mean_density += density[j];
    }

    mean_shannon = mean_shannon/clusters.size();
    mean_density = mean_density/clusters.size();

    cout << "Shannon mean: " << mean_shannon << endl;
    cout << "Density mean: " << mean_density << endl;
}

