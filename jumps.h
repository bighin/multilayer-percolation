#ifndef __JUMPS_H__
#define __JUMPS_H__

#include <stdbool.h>

#include "clusters.h"

int ncluster_evaluate_jumps(struct nclusters_t *nclusters,int id,int spanning,bool pbcz,int *pbins,int *ns);

#endif //__JUMPS_H__
