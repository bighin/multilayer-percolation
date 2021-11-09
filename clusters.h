#ifndef __CLUSTERS_H__
#define __CLUSTERS_H__

#include <stdbool.h>

#include <gsl/gsl_rng.h>

#include "bonds.h"
#include "common.h"

struct nclusters_t
{
	int *vals[MAX_NR_OF_LAYERS];
	int lx,ly,nrlayers;

	struct ibond2d_t *bonds[MAX_NR_OF_LAYERS];
	struct ivbond2d_t *ivbonds[MAX_NR_OF_LAYERS];
};

struct nclusters_t *nclusters_init(int x,int y,int nrlayers);
void nclusters_fini(struct nclusters_t *bc);
int nclusters_get_value(struct nclusters_t *nclusters,int x,int y,int layer);
void nclusters_set_value(struct nclusters_t *nclusters,int x,int y,int layer,int value);

struct statistics_t
{
	int cntsingle;
	int cntbilayer;

	int jumps;
	int matches1;
	int matches2;
	int matches1_by_layer[MAX_NR_OF_LAYERS];
	int matches2_by_layer[MAX_NR_OF_LAYERS];

	int nr_percolating1;
	int nr_percolating2;
};

int nclusters_identify_percolation(struct nclusters_t *nclusters,int *jumps,struct statistics_t *stat,int seq,const gsl_rng *rngctx,bool pbcz);

#endif
