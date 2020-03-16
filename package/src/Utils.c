// Copyright (c) 2012 Shenzhen Institutes of Advanced Technology
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
#include <stdlib.h>
#include <string.h>
#include <math.h>

void sum_squares(const double *x, // Data matrix
		const int *nr, // number of rows
		const int *nc, // number of columns
		const int *k, // number of clusters
		int *cluster, // A vector of integers (from 1:k) indicating the cluster to which each point is allocated.
		double *centers, // Cluster centers of every cluster
		double *totss, // The total sum of squares.
		double *withiness // Vector of within-cluster sum of squares, one component per cluster.
		) {
	//printf("Running sum of squares\n");
	int i, j;

	// vector to save average of all objects ---------------------------
	double* global_average = (double *) malloc(*nc * sizeof(double));
	for (j = 0; j < *nc; j++) {
		global_average[j] = 0.0;
	}

	for (i = 0; i < *nr; i++) {
		for (j = 0; j < *nc; j++) {
			// global_average[j] += x[i][j]
			global_average[j] += x[j * (*nr) + i];
		}
	}
	for (j = 0; j < *nc; j++) {
		global_average[j] = global_average[j] / (*nr);
	}

	// calculate total sum of square ---------------------------------

	*totss = 0.0;
	double tmp_ss, temp;
	for (i = 0; i < (*nr); i++) {
		tmp_ss = 0.0;
		for (j = 0; j < (*nc); j++) {
			temp = global_average[j] - x[j * (*nr) + i];
			tmp_ss += temp * temp;
		}
		*totss += tmp_ss;
	}

	// calculate withiness --------------------------------------------
	int t;
	for (t = 0; t < (*k); t++) {
		withiness[t] = 0.0;
	}

	for (i = 0; i < (*nr); i++) {

		// calculate the sum of square of object i and the center of the cluster which object i in.
		tmp_ss = 0;
		for (j = 0; j < (*nc); j++) {
			temp = x[j * (*nr) + i] - centers[j * (*k) + cluster[i]];
			tmp_ss += temp * temp;
		}

		// sum the with in class sum of square of the cluster which object i in.
		withiness[cluster[i]] += tmp_ss;
	}

	free(global_average);
	return;
}

void parseGroup(const char **strGroup, // string of group information, formatted as "0,1,2,4;3,5;6,7,8" or "0-2,4;3,4,8;6-8", where ";" defines a group;
		int *numGroups, // no. of groups
		int *groupInfo // group assignment
		) {

	int length, index = -1, end, i = 0, j;
	char *buffer, *p1;
	buffer = (char *)malloc(20*sizeof(char));

	// -1: error, ignore next data until next "," or ";"
	//  0: no data
	//  1: single index
	//  2: range
	short status = -1;

	length = strlen(strGroup[0]);
	*numGroups = 0;

	p1 = &buffer[0];
	status = 0;
	while (i < length) {
		switch (strGroup[0][i]) {
		case '-':
			switch (status) {
			case 1:
			case 2:
				*p1 = '\0';
				index = atoi(buffer);
				p1 = &buffer[0];
				status = 2;
				break;
			case 0:
				status = -1;
				break;
			default:
				break;
			}
			break;
		case ';':
		case ',':
			//generate index
			switch (status) {
			case 1:
				*p1 = '\0';
				index = atoi(buffer);
				groupInfo[index] = *numGroups;
				p1 = &buffer[0];
				status = 0;
				break;
			case 2:
				//range
				*p1 = '\0';
				end = atoi(buffer);
				for (j = index; j <= end; j++) {
					groupInfo[j] = *numGroups;
				}
				p1 = &buffer[0];
				status = 0;
				break;
			default:
				//ignore
				status = -1;
				break;
			}
			if (status == -1) {
				//reset
				status = 0;
			} else if (strGroup[0][i] == ';') {
				//next group
				(*numGroups)++;
			}
			break;
		default:
			switch (status) {
			case 0:
				status = 1;
			case 1:
			case 2:
				*p1 = strGroup[0][i];
				p1++;
				break;
			default:
				//ignore
				break;
			}
			break;
		}

		i++;
	}
	//process left
	//generate index
	switch (status) {
	case 1:
		*p1 = '\0';
		index = atoi(buffer);
		groupInfo[index] = *numGroups;
		p1 = &buffer[0];
		status = 0;
		break;
	case 2:
		//range
		*p1 = '\0';
		end = atoi(buffer);
		for (j = index; j <= end; j++) {
			groupInfo[j] = *numGroups;
		}
		p1 = &buffer[0];
		status = 0;
		break;
	default:
		//ignore
		break;
	}
	// add 1 to the no. of groups
	(*numGroups)++;
	free(buffer);
}

double eu_distance(const double f1, const double f2) {
	return (f1 - f2) * (f1 - f2);
}

void init_centers(const double *x, const int *nr, const int *nc, const int *k,
		double *centers) {
	int i, j, l, flag, index, *random_obj_num;
	random_obj_num = (int *) calloc(*k, sizeof(int));

	if (!random_obj_num) {
		error("can't allocate random_obj_num\n");
	}

	for (l = 0; l < *k; l++) {
		random_obj_num[l] = -1;
	}

	for (l = 0; l < *k; l++) {
		flag = 1;

		while (flag) {
		  // index = (int) (rand() % (*nr));
			index = (int) (*nr-1) * unif_rand();
			flag = 0;
			for (i = 0; i < l; i++) {
				if (random_obj_num[i] == index)
					flag = 1;
			}
		}
		random_obj_num[l] = index;

		for (j = 0; j < *nc; j++)
			centers[j * (*k) + l] = x[j * (*nr) + index];

	} // for l = 1:k
	free(random_obj_num);
}

// ----- Calculate exp^{a}/sum(exp^{a}) -----

void expNormalize(double *a, int length, double minValue) {
	int i;
	double max;
	double sum = 0;

	max = a[0];

	for (i = 0; i < length; i++)
		if (a[i] > max)
			max = a[i];

	for (i = 0; i < length; i++) {
		a[i] = exp(a[i] - max);
		sum += a[i];
	}

	double sum2 = 0;
	for (i = 0; i < length; i++) {
		a[i] /= sum;
		if (a[i] < minValue) {
			a[i] = minValue;
		}
		sum2 += a[i];
	}
	for (i = 0; i < length; i++) {
		a[i] /= sum2;
	}
}
