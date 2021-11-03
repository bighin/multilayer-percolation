#ifndef __COMMON_H__
#define __COMMON_H__

#define MAX_NR_OF_LAYERS	(256)

#define MAKE_INDEX(ctx,x,y)	((x)+ctx->lx*(y))

#define MIN(a,b)	(((a)<(b))?(a):(b))
#define MAX(a,b)	(((a)>(b))?(a):(b))

#endif //__COMMON_H__
