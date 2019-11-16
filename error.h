/** @file*/

/**
* @author Antonio Mariani
* @date 16/11/19
*/
#include <iostream>
using namespace std;

#ifndef P2PRASTER_ERROR_H
#define P2PRASTER_ERROR_H

int fileError(string nameFunction);
int readDatasetError(string nameFunction);
int memoryError(string nameFunction);
int partitionError(string nameFunction);
int insertError(string nameFunction);
int findError(string nameFunction);
int dataError(string nameFunction);

#endif //P2PRASTER_ERROR_H
