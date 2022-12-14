#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "common.h"
#include "bonds.h"
#include "clusters.h"
#include "jumps.h"

/*
	Clusters on a nlayer!
*/

struct nclusters_t *nclusters_init(int x,int y,int nrlayers)
{
	struct nclusters_t *ret;
	
	assert(x>0);
	assert(y>0);
	assert(nrlayers>0);
	assert(nrlayers<MAX_NR_OF_LAYERS);

	if(!(ret=malloc(sizeof(struct nclusters_t))))
		return NULL;

	for(int c=0;c<nrlayers;c++)
	{
		if(!(ret->vals[c]=malloc(sizeof(int)*x*y)))
		{
			/*
				If one allocation fails, we free all the
				previous one and return.
			*/

			for(c--;c>=0;c--)
				if(ret->vals[c])
					free(ret->vals[c]);
			
			return NULL;
		}
	}

	ret->lx=x;
	ret->ly=y;
	ret->nrlayers=nrlayers;

	return ret;
}

void nclusters_fini(struct nclusters_t *nc)
{
	if(nc)
	{
		for(int c=0;c<nc->nrlayers;c++)
			if(nc->vals[c])
				free(nc->vals[c]);

		free(nc);
	}
}

int nclusters_get_value(struct nclusters_t *nclusters,int x,int y,int layer)
{
	assert(nclusters!=NULL);
	assert(x>=0);
	assert(y>=0);
	assert(layer>=0);
	assert(x<nclusters->lx);
	assert(y<nclusters->ly);
	assert(layer<nclusters->nrlayers);

	return nclusters->vals[layer][MAKE_INDEX(nclusters,x,y)];
}

void nclusters_set_value(struct nclusters_t *nclusters,int x,int y,int layer,int value)
{
	assert(nclusters!=NULL);
	assert(x>=0);
	assert(y>=0);
	assert(layer>=0);
	assert(x<nclusters->lx);
	assert(y<nclusters->ly);
	assert(layer<nclusters->nrlayers);

	nclusters->vals[layer][MAKE_INDEX(nclusters,x,y)]=value;
}

/*
	Clusters are identified by means of the Hoshen–Kopelman algorithm
*/

#define MAX_NR_OF_CLUSTERS	(2*1024*1024)

int hk_find(int labels[MAX_NR_OF_CLUSTERS+1],int x)
{
	int y=x;

	while(labels[y]!=y)
		y=labels[y];

	/*
		We also collapse the tree of aliases, that's an optional optimization.
	*/

	while (labels[x]!=x)
	{
		int z=labels[x];
		labels[x]=y;
		x=z;
	}

	return y;
}

int hk_union(int labels[MAX_NR_OF_CLUSTERS+1],int x,int y)
{
	return labels[hk_find(labels,x)]=hk_find(labels,y);
}

int find_maximum(const int *vals,int len)
{
	assert(len>0);

	if(len==1)
		return vals[0];

	int maximum=vals[0];

	for(int c=1;c<len;c++)
		if(vals[c]>maximum)
			maximum=vals[c];

	return maximum;
}

int count_non_zeroes(const int *vals,int len)
{
	int cnt=0;

	for(int c=0;c<len;c++)
		if(vals[c]!=0)
			cnt++;

	return cnt;
}

int nclusters_identify_percolation(struct nclusters_t *nclusters,int *jumps,struct statistics_t *stat,int seq,const gsl_rng *rngctx,bool pbcz)
{
	int id=1;

	assert(nclusters);

	for(int x=0;x<nclusters->lx;x++)
		for(int y=0;y<nclusters->ly;y++)
			for(int l=0;l<nclusters->nrlayers;l++)
				nclusters_set_value(nclusters,x,y,l,0);

	int *labels=malloc(sizeof(int)*(MAX_NR_OF_CLUSTERS+1));

	for(int x=0;x<nclusters->lx;x++)
	{
		for(int y=0;y<nclusters->ly;y++)
		{
			for(int l=0;l<nclusters->nrlayers;l++)
			{

#define NR_OF_NEIGHBOURS	(3)

				int neighbours[NR_OF_NEIGHBOURS]={0,0,0};

				if(x!=0)
					if(ibond2d_get_value(nclusters->bonds[l],x-1,y,DIR_X)==1)
						neighbours[0]=nclusters_get_value(nclusters,x-1,y,l);

				if(y!=0)
					if(ibond2d_get_value(nclusters->bonds[l],x,y-1,DIR_Y)==1)
						neighbours[1]=nclusters_get_value(nclusters,x,y-1,l);

				if(l!=0)
					if(ivbond2d_get_value(nclusters->ivbonds[l-1], x, y)==1)
						neighbours[2]=nclusters_get_value(nclusters, x, y, l-1);

				int nr_of_neighbours=count_non_zeroes(neighbours, NR_OF_NEIGHBOURS);

				if(nr_of_neighbours==0)
				{
					labels[id]=id;
					nclusters_set_value(nclusters,x,y,l,id++);

					assert(id<MAX_NR_OF_CLUSTERS);
				}
				else
				{
					int maximum=find_maximum(neighbours,NR_OF_NEIGHBOURS);

					for(int j=0;j<NR_OF_NEIGHBOURS;j++)
						if((neighbours[j]!=0)&&(neighbours[j]!=maximum))
							hk_union(labels,neighbours[j],maximum);

					nclusters_set_value(nclusters,x,y,l,hk_find(labels,maximum));
				}

				assert(nclusters_get_value(nclusters,x,y,l)!=0);
			}
		}
	}

	if(pbcz==true)
		for(int x=0;x<nclusters->lx;x++)
			for(int y=0;y<nclusters->ly;y++)
				if(ivbond2d_get_value(nclusters->ivbonds[nclusters->nrlayers-1], x, y)==1)
					hk_union(labels,nclusters_get_value(nclusters,x,y,0),nclusters_get_value(nclusters,x,y,nclusters->nrlayers-1));

	/*
		Normalization and collection of statistics about the clusters.
	*/

	id=1;

	int *new_labels=malloc(sizeof(int)*(MAX_NR_OF_CLUSTERS+1));

	for(int c=0;c<(MAX_NR_OF_CLUSTERS+1);c++)
		new_labels[c]=0;

	struct cluster_info_t
	{
		int minx,miny,maxx,maxy;
	};

	struct cluster_info_t *info=malloc(sizeof(struct cluster_info_t)*(MAX_NR_OF_CLUSTERS+1));

	for(int x=0;x<nclusters->lx;x++)
	{
		for(int y=0;y<nclusters->ly;y++)
		{
			for(int l=0;l<nclusters->nrlayers;l++)
			{
				int entry=nclusters_get_value(nclusters,x,y,l);
				int r=hk_find(labels, entry);

				if(new_labels[r]==0)
				{
					new_labels[r]=id++;

					info[new_labels[r]].minx=info[new_labels[r]].maxx=x;
					info[new_labels[r]].miny=info[new_labels[r]].maxy=y;
				}
				else
				{
					info[new_labels[r]].minx=MIN(info[new_labels[r]].minx, x);
					info[new_labels[r]].maxx=MAX(info[new_labels[r]].maxx, x);

					info[new_labels[r]].miny=MIN(info[new_labels[r]].miny, y);
					info[new_labels[r]].maxy=MAX(info[new_labels[r]].maxy, y);
				}

				nclusters_set_value(nclusters,x,y,l,new_labels[r]);

				assert(nclusters_get_value(nclusters,x,y,l)!=0);
			}
		}
	}

	if(new_labels)
		free(new_labels);

	/*
		We select a random lattice site to check whether it belongs to the percolating cluster,
		following the criterion in lecture 2 of "An Introduction to Universality" by A. Codello.
	*/

	int rx=gsl_rng_uniform_int(rngctx, nclusters->lx);
	int ry=gsl_rng_uniform_int(rngctx, nclusters->ly);
	int rl=gsl_rng_uniform_int(rngctx, nclusters->nrlayers);

	int rx_by_layer[MAX_NR_OF_LAYERS];
	int ry_by_layer[MAX_NR_OF_LAYERS];
	int rz_by_layer[MAX_NR_OF_LAYERS];

	for(int z=0;z<nclusters->nrlayers;z++)
	{
		rx_by_layer[z]=gsl_rng_uniform_int(rngctx, nclusters->lx);
		ry_by_layer[z]=gsl_rng_uniform_int(rngctx, nclusters->ly);
		rz_by_layer[z]=z;
	}

	/*
		Finally, we count the percolating clusters.

		The percolation criterion corresponds to the extension rule in:
		J. Machta, Y.S. Choi, A. Lucke, T. Schweizer, and L.M. Chayes,
		Phys. Rev. Lett. 75, 2792 (1995);
		Phys. Rev. E 54, 1332 (1996).
	*/

#warning Jumps are calculated only for the first percolating cluster we find. This could be easily extended to all clusters.
#warning The only difference is that then jumps should be divided by nr_percolating

	int maxid=id-1;
	int nr_percolating=0;
	bool first_has_been_found=false;

	for(int thisid=1;thisid<=maxid;thisid++)
	{
		assert(info[thisid].maxx>=info[thisid].minx);
		assert(info[thisid].maxy>=info[thisid].miny);

		int xlength=info[thisid].maxx-info[thisid].minx+1;
		int ylength=info[thisid].maxy-info[thisid].miny+1;

		if(xlength==nclusters->lx)
		{
			if(first_has_been_found==false)
			{
				first_has_been_found=true;

				if(jumps!=NULL)
					*jumps=ncluster_evaluate_jumps(nclusters, thisid, DIR_X, pbcz, stat->pbins, stat->ns);
			}

			if(nclusters_get_value(nclusters, rx, ry, rl)==thisid)
			{
				switch(seq)
				{
					case 1:
					stat->matches1=1;
					break;

					case 2:
					stat->matches2=1;
					break;
				}
			}

			for(int c=0;c<nclusters->nrlayers;c++)
			{
				if(nclusters_get_value(nclusters, rx_by_layer[c], ry_by_layer[c], rz_by_layer[c])==thisid)
				{
					switch(seq)
					{
						case 1:
						stat->matches1_by_layer[c]=1;
						break;

						case 2:
						stat->matches2_by_layer[c]=1;
						break;
					}
				}
			}

			nr_percolating++;
		}
		else if(ylength==nclusters->ly)
		{
			if(first_has_been_found==false)
			{
				first_has_been_found=true;

				if(jumps!=NULL)
					*jumps=ncluster_evaluate_jumps(nclusters, thisid, DIR_Y, pbcz, stat->pbins, stat->ns);
			}

			if(nclusters_get_value(nclusters, rx, ry, rl)==thisid)
			{
				switch(seq)
				{
					case 1:
					stat->matches1=1;
					break;

					case 2:
					stat->matches2=1;
					break;
				}
			}

			for(int c=0;c<nclusters->nrlayers;c++)
			{
				if(nclusters_get_value(nclusters, rx_by_layer[c], ry_by_layer[c], rz_by_layer[c])==thisid)
				{
					switch(seq)
					{
						case 1:
						stat->matches1_by_layer[c]=1;
						break;

						case 2:
						stat->matches2_by_layer[c]=1;
						break;
					}
				}
			}

			nr_percolating++;
		}
	}

	if(info)
		free(info);

	if(labels)
		free(labels);

	return nr_percolating;
}
