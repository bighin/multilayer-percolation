#include <assert.h>
#include <stdlib.h>

#include "common.h"
#include "bonds.h"

/*
	A floating point quantity defined on each bond in a two-dimensional lattice
*/

struct bond2d_t *bond2d_init(int x,int y)
{
	struct bond2d_t *ret;

	assert(x>0);
	assert(y>0);

	if(!(ret=malloc(sizeof(struct bond2d_t))))
		return NULL;
	
	ret->vals[0]=malloc(sizeof(double)*x*y);
	ret->vals[1]=malloc(sizeof(double)*x*y);

	if((!ret->vals[0])||(!ret->vals[1]))
	{
		if(ret->vals[0])
			free(ret->vals[0]);

		if(ret->vals[1])
			free(ret->vals[1]);

		if(ret)
			free(ret);
		
		return NULL;
	}
	
	ret->lx=x;
	ret->ly=y;
	
	return ret;
}

void bond2d_fini(struct bond2d_t *b)
{
	if(b)
	{
		if(b->vals[0])
			free(b->vals[0]);

		if(b->vals[1])
			free(b->vals[1]);
		
		free(b);
	}
}

double bond2d_get_value(struct bond2d_t *b,int x,int y,short direction)
{
	assert((x>=0)&&(x<b->lx));
	assert((y>=0)&&(y<b->ly));
	assert((direction==DIR_X)||(direction==DIR_Y));

	return b->vals[direction][MAKE_INDEX(b,x,y)];
}

void bond2d_set_value(struct bond2d_t *b,int x,int y,short direction,double value)
{
	assert((x>=0)&&(x<b->lx));
	assert((y>=0)&&(y<b->ly));
	assert((direction==DIR_X)||(direction==DIR_Y));

	b->vals[direction][MAKE_INDEX(b,x,y)]=value;
}

/*
	An integer quantity defined on each bond in a two-dimensional lattice
*/

struct ibond2d_t *ibond2d_init(int x,int y)
{
	struct ibond2d_t *ret;

	assert(x>0);
	assert(y>0);

	if(!(ret=malloc(sizeof(struct ibond2d_t))))
		return NULL;
	
	ret->vals[0]=malloc(sizeof(int)*x*y);
	ret->vals[1]=malloc(sizeof(int)*x*y);

	if((!ret->vals[0])||(!ret->vals[1]))
	{
		if(ret->vals[0])
			free(ret->vals[0]);

		if(ret->vals[1])
			free(ret->vals[1]);

		if(ret)
			free(ret);
		
		return NULL;
	}
	
	ret->lx=x;
	ret->ly=y;
	
	return ret;
}

void ibond2d_fini(struct ibond2d_t *b)
{
	if(b)
	{
		if(b->vals[0])
			free(b->vals[0]);

		if(b->vals[1])
			free(b->vals[1]);
		
		free(b);
	}
}

int ibond2d_get_value(struct ibond2d_t *b,int x,int y,short direction)
{
	assert((x>=0)&&(x<b->lx));
	assert((y>=0)&&(y<b->ly));
	assert((direction==DIR_X)||(direction==DIR_Y));

	return b->vals[direction][MAKE_INDEX(b,x,y)];
}

void ibond2d_set_value(struct ibond2d_t *b,int x,int y,short direction,int value)
{
	assert((x>=0)&&(x<b->lx));
	assert((y>=0)&&(y<b->ly));
	assert((direction==DIR_X)||(direction==DIR_Y));

	b->vals[direction][MAKE_INDEX(b,x,y)]=value;
}

/*
	The vbond2d_t and ivbond2d_t structures describe a lattice of
	vertical bond variables (defined between two different layers)
	taking double and integer values, respectively.
*/

struct vbond2d_t *vbond2d_init(int x,int y)
{
	struct vbond2d_t *ret;

	assert(x>0);
	assert(y>0);

	if(!(ret=malloc(sizeof(struct vbond2d_t))))
		return NULL;
	
	ret->vals=malloc(sizeof(double)*x*y);

	if(!ret->vals)
	{
		if(ret)
			free(ret);
		
		return NULL;
	}
	
	ret->lx=x;
	ret->ly=y;
	
	return ret;
}

void vbond2d_fini(struct vbond2d_t *vb)
{
	if(vb)
	{
		if(vb->vals)
			free(vb->vals);
		
		free(vb);
	}
}

double vbond2d_get_value(struct vbond2d_t *vb,int x,int y)
{
	assert((x>=0)&&(x<vb->lx));
	assert((y>=0)&&(y<vb->ly));

	return vb->vals[MAKE_INDEX(vb,x,y)];
}

void vbond2d_set_value(struct vbond2d_t *vb,int x,int y,double val)
{
	assert((x>=0)&&(x<vb->lx));
	assert((y>=0)&&(y<vb->ly));

	vb->vals[MAKE_INDEX(vb,x,y)]=val;
}

struct ivbond2d_t *ivbond2d_init(int x,int y)
{
	struct ivbond2d_t *ret;

	assert(x>0);
	assert(y>0);

	if(!(ret=malloc(sizeof(struct ivbond2d_t))))
		return NULL;
	
	ret->vals=malloc(sizeof(int)*x*y);

	if(!ret->vals)
	{
		if(ret)
			free(ret);
		
		return NULL;
	}
	
	ret->lx=x;
	ret->ly=y;
	
	return ret;
}

void ivbond2d_fini(struct ivbond2d_t *vb)
{
	if(vb)
	{
		if(vb->vals)
			free(vb->vals);
		
		free(vb);
	}
}

int ivbond2d_get_value(struct ivbond2d_t *vb,int x,int y)
{
	assert((x>=0)&&(x<vb->lx));
	assert((y>=0)&&(y<vb->ly));

	return vb->vals[MAKE_INDEX(vb,x,y)];
}

void ivbond2d_set_value(struct ivbond2d_t *vb,int x,int y,int val)
{
	assert((x>=0)&&(x<vb->lx));
	assert((y>=0)&&(y<vb->ly));

	vb->vals[MAKE_INDEX(vb,x,y)]=val;
}
