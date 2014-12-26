// Copyright (c) 2011 Shenzhen Institutes of Advanced Technology
// Chinese Academy of Sciences

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <R.h>

#include <math.h>
#include <ctype.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "Utils.h"

void init_featureWeight(double *featureWeight, const int *k, const int *nc,
		const int *numGroups, const int *groupInfo) {
	int l, j, *nums;
	nums = (int *) calloc(*numGroups, sizeof(int));
	for (j = 0; j < *nc; j++) {
		nums[groupInfo[j]]++;
	}
	for (l = 0; l < *k; l++)
		for (j = 0; j < *nc; j++)
			featureWeight[j * (*k) + l] = 1.0 / nums[groupInfo[j]];
	free(nums);
}

void init_groupWeight(double *groupWeight, const int *k, const int *numGroups) {
	int l, t;
	for (l = 0; l < *k; l++)
		for (t = 0; t < *numGroups; t++)
			groupWeight[t * (*k) + l] = 1.0 / *numGroups;
}

void update_cluster(const double *x, const int *nr, const int *nc, const int *k,
		const int *numGroups, const int *groupInfo, int *cluster,
		double *centers, double *featureWeight, double *groupWeight) {
	int i, j, l;
	double min_dist, o_dist;

	for (i = 0; i < *nr; i++) {
		min_dist = 1.79769e+308;
		for (l = 0; l < *k; l++) {

			o_dist = 0.0;
			for (j = 0; j < *nc; j++) {
				o_dist += featureWeight[j * (*k) + l]
						* groupWeight[groupInfo[j] * (*k) + l]
						* eu_distance(centers[j * (*k) + l], x[j * (*nr) + i]);
			}
			if (o_dist <= min_dist) {
				min_dist = o_dist;
				cluster[i] = l;
			}
		}
	}
}

int update_centers(const double *x, const int *nr, const int *nc, const int *k,
		int *cluster, double *centers) {
	int i, j, l, *no_cluster;
	no_cluster = (int *) calloc(*k, sizeof(int));

	if (!no_cluster) {
		error("can not allocate [].\n");
	}

	for (l = 0; l < *k; l++) {
		for (j = 0; j < *nc; j++) {
			centers[j * (*k) + l] = 0.0;
		}
	}

	for (i = 0; i < *nr; i++) {
		no_cluster[cluster[i]]++;
		for (j = 0; j < *nc; j++) {
			centers[j * (*k) + cluster[i]] += x[j * (*nr) + i];
		}
	}

	int flag = 1;
	for (l = 0; l < *k; l++) {
		if (no_cluster[l] == 0) {
			flag = 0;
			break;
		}
		for (j = 0; j < *nc; j++) {
			centers[j * (*k) + l] /= (double) no_cluster[l];
		}
	}
	free(no_cluster);
	return flag;
}

void update_featureWeight(const double *x, const int *nr, const int *nc,
		const int *k, const double *eta, const int *numGroups,
		const int *groupInfo, int *cluster, double *centers,
		double *featureWeight, double *groupWeight) {
	int i, j, l, t;

	int index;
	for (i = 0; i < (*nc) * (*k); i++) {
		featureWeight[i] = 0;
	}

	for (j = 0; j < *nc; j++) {
		for (i = 0; i < *nr; i++) {
			index = j * (*k) + cluster[i];
			featureWeight[index] += groupWeight[groupInfo[j]*(*k)+cluster[i]]
					* eu_distance(x[j * (*nr) + i],
							centers[j * (*k) + cluster[i]]);
		}
	}

	double *sum, *sum2, *max;
	sum = (double*) malloc(*numGroups * sizeof(double));
	sum2 = (double*) malloc(*numGroups * sizeof(double));
	max = (double*) malloc(*numGroups * sizeof(double));

	for (t = 0; t < *numGroups; t++) {
		sum[t] = 0;
		sum2[t] = 0;
		max[t] = -1.79769e+308;
	}

	double minWeight = (double) 0.00001 / (*nc);
	for (l = 0; l < *k; l++) { // every CLUSTER
		for (t = 0; t < *numGroups; t++) {
			sum[t] = 0;
			sum2[t] = 0;
		}
		//find maximum
		for (j = 0; j < *nc; j++) {
			index = j * (*k) + l;
			featureWeight[index] = -featureWeight[index] / (*eta);
			if (max[groupInfo[j]] < featureWeight[index]) {
				max[groupInfo[j]] = featureWeight[index];
			}
		}
		//compute exp()
		for (j = 0; j < *nc; j++) {
			index = j * (*k) + l;
			featureWeight[index] = exp(
					featureWeight[index] - max[groupInfo[j]]);
			sum[groupInfo[j]] += featureWeight[index];
		}
		//normalize
		for (j = 0; j < *nc; j++) {
			index = j * (*k) + l;
			featureWeight[index] /= sum[groupInfo[j]];
			if (featureWeight[index] < minWeight) {
				featureWeight[index] = minWeight;
			}
			sum2[groupInfo[j]] += featureWeight[index];
		}

		for (j = 0; j < *nc; j++) {
			featureWeight[j * (*k) + l] /= sum2[groupInfo[j]];
		}
	}

	free(sum);
	free(sum2);
	free(max);
}

void update_groupWeight(const double *x, const int *nr, const int *nc,
		const int *k, const double *lambda, const int *numGroups,
		const int *groupInfo, int *cluster, double *centers,
		double *featureWeight, double *groupWeight) {
	int i, j, l, t, index;

	for (l = 0; l < (*numGroups) * (*k); l++) {
		groupWeight[l] = 0;
	}

	for (j = 0; j < *nc; j++) {
		index = groupInfo[j] * (*k);
		for (i = 0; i < *nr; i++) {
			groupWeight[index + cluster[i]] += featureWeight[j * (*k)
					+ cluster[i]]
					* eu_distance(centers[j * (*k) + cluster[i]],
							x[j * (*nr) + i]);
		}
	}

	for (l = 0; l < (*numGroups) * (*k); l++) {
		groupWeight[l] = -groupWeight[l] / (*lambda);
	}

	double sum = 0.0, sum2, max = 0.0;

// implement expNormalize()
	double minWeight = (double) 0.00001 / (*numGroups);
	for (l = 0; l < *k; l++) {
		sum = 0.0;
		max = groupWeight[0 * (*k) + l]; // initially assign gw[l][0] to max
		for (t = 1; t < *numGroups; t++) {
			index = t * (*k) + l;
			if (groupWeight[index] >= max)
				max = groupWeight[index];
		}

		for (t = 0; t < *numGroups; t++) {
			index = t * (*k) + l;
			groupWeight[index] = exp(groupWeight[index] - max);
			sum += groupWeight[index];
		}

		sum2 = 0;
		for (t = 0; t < *numGroups; t++) {
			index = t * (*k) + l;
			groupWeight[index] /= sum;
			if (groupWeight[index] < minWeight) {
				groupWeight[index] = minWeight;
			}
			sum2 += groupWeight[index];
		}

		if (sum2 != 1) {
			for (t = 0; t < *numGroups; t++) {
				groupWeight[t * (*k) + l] /= sum2;
			}
		}
	}
}

double calculate_cost(const double *x, const int *nr, const int *nc,
		const int *k, const double *lambda, const double *eta,
		const int *numGroups, const int *groupInfo, int *cluster,
		double *centers, double *featureWeight, double *groupWeight) {

	int i, j, l, t;
	double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, dispersion;

	for (i = 0; i < *nr; i++) {
		for (j = 0; j < *nc; j++) {
			sum0 += groupWeight[groupInfo[j] * (*k) + cluster[i]]
					* featureWeight[j * (*k) + cluster[i]]
					* eu_distance(centers[j * (*k) + cluster[i]],
							x[j * (*nr) + i]);
		}
	}

	for (l = 0; l < *k; l++) {
		for (t = 0; t < *numGroups; t++) {
			sum1 += groupWeight[t * (*k) + l] * log(groupWeight[t * (*k) + l]);
		}

		for (j = 0; j < *nc; j++)
			sum2 += featureWeight[j * (*k) + l]
					* log(featureWeight[j * (*k) + l]);

	}

	dispersion = sum0 + sum1 * (*lambda) + sum2 * (*eta);
	return dispersion;
}

void fgkm(const double *x, const int *nr, const int *nc, const int *k,
		const double *lambda, const double *eta, const int *numGroups,
		const int *groupInfo, const double *delta, const int *maxiter,
		const int *maxrestart, int *init, // const unsigned int *seed,
		int *cluster,
		double *centers, double *featureWeight, double *groupWeight,
		int *iterations, int *restarts, int *totiter, double *totalCost, //
		double *totss, //    total sum of squares
		double *withiness // vector of sum of square in every cluster
		) {

	double dispersion1, dispersion2;
	int flag_not_restart = 1;

	*restarts = 0;
	*iterations = 0;
	*totiter = 0;

	// srand(seed);
	GetRNGstate();
	while ((*restarts) < *maxrestart) {
		if (*init == 0)
		  init_centers(x, nr, nc, k, centers); // assign randomly

		init_featureWeight(featureWeight, k, nc, numGroups, groupInfo); // equal value
		init_groupWeight(groupWeight, k, numGroups); // equal value

		*iterations = 0;
		dispersion2 = 1.79769e+308;

		while ((*iterations) < *maxiter) {
			(*iterations)++;
			(*totiter)++;
			dispersion1 = dispersion2;

			update_cluster(x, nr, nc, k, numGroups, groupInfo, cluster, centers,
					featureWeight, groupWeight);

			flag_not_restart = update_centers(x, nr, nc, k, cluster, centers);

			if (!flag_not_restart) {
				(*restarts)++;
				break;
			}

			update_featureWeight(x, nr, nc, k, eta, numGroups, groupInfo,
					cluster, centers, featureWeight, groupWeight);

			update_groupWeight(x, nr, nc, k, lambda, numGroups, groupInfo,
					cluster, centers, featureWeight, groupWeight);

			dispersion2 = calculate_cost(x, nr, nc, k, lambda, eta, numGroups,
					groupInfo, cluster, centers, featureWeight, groupWeight);

			// if change of dispersion below delta or iterations exceed max iterations,
			//	then, terminate and return.
			//  but if dispersion < 0, then ???
			if ((fabs((dispersion1 - dispersion2) / dispersion1)) <= (*delta)
					|| (*iterations) == *maxiter) {

				*totalCost = dispersion1;

				/**
				 * calculate sum of squares
				 */
				sum_squares(x, nr, nc, k, cluster, centers, totss, withiness);
				// Done.
				PutRNGstate();
				return;
			}
		}
	}

	PutRNGstate();
}
