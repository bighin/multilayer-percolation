#include <assert.h>
#include <stdlib.h>

#include "sites.h"

struct sites2d_t *sites2d_init(int x,int y)
{
	struct sites2d_t *ret;

	assert(x>0);
	assert(y>0);

	if(!(ret=malloc(sizeof(struct sites2d_t))))
		return NULL;
	
	ret->values=malloc(sizeof(int)*x*y);

	if(!ret->values)
	{
		if(ret)
			free(ret);
		
		return NULL;
	}
	
	ret->lx=x;
	ret->ly=y;
	
	return ret;
}

void sites2d_fini(struct sites2d_t *s)
{
	if(s)
	{
		if(s->values)
			free(s->values);
		
		free(s);
	}
}

int sites2d_get(struct sites2d_t *s,int x,int y)
{
	assert(s!=NULL);
	assert((x>=0)&&(x<s->lx));
	assert((y>=0)&&(y<s->ly));

	return s->values[MAKE_INDEX(s,x,y)];
}

void sites2d_set(struct sites2d_t *s,int x,int y,int value)
{
	assert(s!=NULL);
	assert((x>=0)&&(x<s->lx));
	assert((y>=0)&&(y<s->ly));

	s->values[MAKE_INDEX(s,x,y)]=value;
}
