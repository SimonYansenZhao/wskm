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

/**
 * This function calculate the total sum squares and
 * within class sum squares of every cluster.
 *
 * @author Longfei Xiao <lf.xiao@siat.ac.cn>
 */

#include <stdlib.h>

void sum_squares(const double *x, // Data matrix
		const int *nr, // number of rows
		const int *nc, // number of columns
		const int *k, // number of clusters
		int *cluster, // A vector of integers (from 1:k) indicating the cluster to which each point is allocated.
		double *centers, // Cluster centers of every cluster
		double *totss, // The total sum of squares.
		double *withiness // Vector of within-cluster sum of squares, one component per cluster.
		);

void parseGroup(const char *strGroup, // string of group information, formatted as "0,1,2,4;3,5;6,7,8" or "0-2,4;3,4,8;6-8", where ";" defines a group;
		int *numGroups, // no. of groups
		int *groupInfo // group assignment
		);

double eu_distance(const double f1, const double f2);

void init_centers(const double *x, const int *nr, const int *nc, const int *k,
		double *centers);

// ----- Calculate exp^{a}/sum(exp^{a}) -----

void expNormalize(double *a, int length, double minValue);
