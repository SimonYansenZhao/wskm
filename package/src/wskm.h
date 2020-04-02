#include <stdlib.h>
#include "Utils.h"

void ewkm( // Inputs ----------------------------------------------------------
		double *x,         // Numeric matrix as vector by col (nr*nc)
		int *nr,           // Number of rows
		int *nc,           // Number of columns
		int *k,            // Number of clusters
		double *lambda,    // Learning rate
		int *maxiter,      // Maximum number of iterations
		double *delta,     // Minimum change below which iteration stops
		int *maxrestart,   // Maximum number of restarts
		int *init,         // Initial k prototypes.
		// Outputs ---------------------------------------------------------
		int *iterations,   // Number of iterations
		int *cluster,      // Cluster assignment for each obs (nr)
		double *centers,   // Cluster centers (k*nc)
		double *weights,   // Variable weights (k*nc)
		int *restarts,	   // Number of restarts
		int *totiters,     // Number of iterations including restarts
		double *totss,     // Total sum of squares
		double *withiness  // Vector of sum of square in every cluster
		);

void fgkm(const double *x, const int *nr, const int *nc, const int *k,
		const double *lambda, const double *eta, const int *numGroups,
		const int *groupInfo, const double *delta, const int *maxiter,
		const int *maxrestart, int *init, // const unsigned int *seed,
		int *cluster,
		double *centers, double *featureWeight, double *groupWeight,
		int *iterations, int *restarts, int *totiter, double *totalCost, //
		double *totss, //    total sum of squares
		double *withiness // vector of sum of square in every cluster
		);

void twkm(const double *x, const int *nr, const int *nc, const int *k,
		const double *lambda, const double *eta, const int *numGroups,
		const int *groupInfo, const double *delta, const int *maxiter,
		const int *maxrestart, int *init, // unsigned int *seed,
		int *cluster,
		double *centers, double *featureWeight, double *groupWeight,
		int *iterations, int *restarts, int *totiter, double *totalCost, //
		double *totss, //    total sum of squares
		double *withiness // vector of sum of square in every cluster
		);
