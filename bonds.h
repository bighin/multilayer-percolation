#ifndef __BONDS_H__
#define __BONDS_H__

#define DIR_X	(0)
#define DIR_Y	(1)

struct bond2d_t
{
	double *vals[2];
	int lx,ly;
};

struct bond2d_t *bond2d_init(int x,int y);
void bond2d_fini(struct bond2d_t *b);
double bond2d_get_value(struct bond2d_t *b,int x,int y,short direction);
void bond2d_set_value(struct bond2d_t *b,int x,int y,short direction,double value);

struct ibond2d_t
{
	int *vals[2];
	int lx,ly;
};

struct ibond2d_t *ibond2d_init(int x,int y);
void ibond2d_fini(struct ibond2d_t *b);
int ibond2d_get_value(struct ibond2d_t *b,int x,int y,short direction);
void ibond2d_set_value(struct ibond2d_t *b,int x,int y,short direction,int value);

struct vbond2d_t
{
	double *vals;
	int lx,ly;
};

struct vbond2d_t *vbond2d_init(int x,int y);
void vbond2d_fini(struct vbond2d_t *vb);
double vbond2d_get_value(struct vbond2d_t *vb,int x,int y);
void vbond2d_set_value(struct vbond2d_t *vb,int x,int y,double val);

struct ivbond2d_t
{
	int *vals;
	int lx,ly;
};

struct ivbond2d_t *ivbond2d_init(int x,int y);
void ivbond2d_fini(struct ivbond2d_t *vb);
int ivbond2d_get_value(struct ivbond2d_t *vb,int x,int y);
void ivbond2d_set_value(struct ivbond2d_t *vb,int x,int y,int val);

#endif //__BONDS_H__
