#ifndef __SITES_H__
#define __SITES_H__

#include "common.h"

/*
	This struct identifies a two-dimensional Ising-like configuration.
*/

struct sites2d_t
{
	int *values;
	int lx,ly;
};

struct sites2d_t *sites2d_init(int x,int y);
void sites2d_fini(struct sites2d_t *s);
int sites2d_get(struct sites2d_t *s,int x,int y);
void sites2d_set(struct sites2d_t *s,int x,int y,int value);

#endif //__SITES_H__
