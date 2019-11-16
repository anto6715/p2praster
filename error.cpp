/** @file*/

/**
* @author Antonio Mariani
* @date 16/11/19
*/

#include <iostream>
#include "error.h"
using namespace std;

int fileError(string nameFunction) {
    cerr << nameFunction << ": Can't open file: " << "clustered.csv" << endl;
    return -7;
}

int readDatasetError(string nameFunction) {
    cerr << nameFunction << ": Can't load dataset: " << "clustered.csv" << endl;
    return -8;
}

int memoryError(string nameFunction) {
    cerr << nameFunction << ": Not enough memory" << endl;
    return -1;
}

int partitionError(string nameFunction) {
    cerr << nameFunction << ": Error partition dataset" << endl;
    return -9;
}

int insertError(string nameFunction) {
    cerr << nameFunction << ": Insert Error" << endl;
    return -2;
}

int findError(string nameFunction) {
    cerr << nameFunction << ": Find Error" << endl;
    return -10;
}

int dataError(string nameFunction) {
    cerr << nameFunction << ": Corrupt data Error" << endl;
    return -3;
}
