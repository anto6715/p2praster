/** @file*/

/**
* @author Antonio Mariani
* @date 27/11/19
*/

#include "raster.h"

#ifndef P2PRASTER_STATISTIC_H
#define P2PRASTER_STATISTIC_H

int getPrecision(vectorSet2D &C, vectorSet2D &T, double &Precision);

int getPij(vectorSet2D &C, vectorSet2D &T, double **pij, int ni);

int getPCi(vectorSet2D &C, int ni, double *pCi);

int getNij(unSet2D &Ti, unSet2D &Ci);

int getMax_nij(vectorSet2D &T, unSet2D &Ci);

int getMax_nij_mji(vectorSet2D &T, unSet2D &Ci, int &max_nij, int &mji);

int getRecall(vectorSet2D &C, vectorSet2D &T, double &Recall);

double getFMeasure(double Recall, double Precision);

int getConditionalEntropy(vectorSet2D &C, vectorSet2D &T, int ni, double &Entropy);

#endif //P2PRASTER_STATISTIC_H
