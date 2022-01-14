#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <limits.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>

#include "jumps.h"
#include "clusters.h"

/*
	Simply prints a matrix, with a newline after each row.
	It seems like this simple function is missing in GSL.
*/

void gsl_matrix_int_print(gsl_matrix_int *m)
{
	for(size_t i=0;i<m->size1;i++)
	{
		for(size_t j=0;j<m->size2;j++)
			printf("%d ", gsl_matrix_int_get(m, i, j));

		printf("\n");
	}
}

/*
	This is a just a thin wrapper around the standard GSL sparse matrix.

	The reason is that in GSL the default value for element which haven't
	been set yet is 0, whereas here we want INT_MAX as default value.
*/

struct adjacency_t
{
	gsl_spmatrix_int *m;
	int nr_vertices;
};

struct adjacency_t *init_adjacency(int nr_vertices)
{
	struct adjacency_t *ret=malloc(sizeof(struct adjacency_t));

	ret->m=gsl_spmatrix_int_alloc(nr_vertices,nr_vertices);
	ret->nr_vertices=nr_vertices;

	return ret;
}

void fini_adjacency(struct adjacency_t *adj)
{
	gsl_spmatrix_int_free(adj->m);
}

void adjacency_set(struct adjacency_t *adj, int i, int j, int weight)
{
	if(weight==INT_MAX)
		return;

	gsl_spmatrix_int_set(adj->m,i,j,1+weight);
}

int adjacency_get(struct adjacency_t *adj, int i, int j)
{
	int result=gsl_spmatrix_int_get(adj->m,i,j);

	if(result==0)
		return INT_MAX;

	return result-1;
}

/*
	A utility function to find the vertex with minimum distance value,
	from the set of vertices not yet included in shortest path tree.
*/

int minimum_distance(const int *distances, const bool *in_spt, int V)
{
	int min, min_index;

	min=INT_MAX;
	for(int v=0; v<V; v++)
	{
		if((in_spt[v]==false)&&(distances[v]<=min))
		{
			min=distances[v];
			min_index=v;
		}
	}

	return min_index;
}

/*
	Dijkstra algorithm, adapted from: https://gist.github.com/zSANSANz/408c77daeb734cb40936c39b610dc579
*/

int dijkstra_distance(struct adjacency_t *adj, int from, int to)
{
	int V=adj->nr_vertices;

	/*
		The output array. distances[i] will hold the shortest distance from 'from' to i.
	*/

	int *distances=malloc(sizeof(int)*V);

	/*
		in_spt[i] will be true if either vertex i is included in the shortest
		path tree or the shortest distance from 'from' to i is finalized.
	*/

	bool *in_spt=malloc(sizeof(bool)*V);

	/*
		Initialize all distances to infinite and all values in in_spt to false.
	*/

	for (int i=0; i<V; i++)
	{
		distances[i]=INT_MAX;
		in_spt[i]=false;
	}

	/*
		Distance of the starting vertex from itself is always 0
	*/

	distances[from]=0;

	/*
		Find shortest path for all vertices
	*/

	for(int count=0; count<(V-1); count++)
	{
		/*
			Pick the minimum distance vertex from the set of vertices not
			yet processed. u is always equal to src in first iteration.
		*/

		int u=minimum_distance(distances, in_spt, V);

		/*
			Mark the picked vertex as processed
		*/

		in_spt[u]=true;

		/*
			Update distances value of the adjacent vertices of the picked vertex.
		*/

		for(int v=0; v<V; v++)
		{
			/*
				Update distances[v] only if is not in in_spt, there is an edge from
				u to v, and total weight of path from src to v through u is
				smaller than current value of distances[v]
			*/

			int edge_weight=adjacency_get(adj, u, v);

			if((in_spt[v]==false)&&(edge_weight<INT_MAX)&&(distances[u]!=INT_MAX)&&(distances[u]+edge_weight<distances[v]))
				distances[v]=distances[u]+edge_weight;
		}
	}

	int result=distances[to];

	if(distances)
		free(distances);

	if(in_spt)
		free(in_spt);

	return result;
}

void add_edge(struct adjacency_t *adj, int id1, int id2, int weight)
{
	/*
		This could be redundant.
	*/

	adjacency_set(adj,id1,id2,weight);
	adjacency_set(adj,id2,id1,weight);
}

void process_neighbours(struct nclusters_t *vertices, struct adjacency_t *adj, int x, int y, int l, int spanning,bool pbcz)
{
	int id=nclusters_get_value(vertices, x, y, l);

	if((x!=0)&&(ibond2d_get_value(vertices->bonds[l],x-1,y,DIR_X)==1))
	{
		int xprime=x-1;
		int yprime=y;
		int lprime=l;
		int idprime=nclusters_get_value(vertices, xprime, yprime, lprime);

		if(idprime!=-1)
			add_edge(adj,id,idprime,0);
	}

	if((y!=0)&&(ibond2d_get_value(vertices->bonds[l],x,y-1,DIR_Y)==1))
	{
		int xprime=x;
		int yprime=y-1;
		int lprime=l;
		int idprime=nclusters_get_value(vertices, xprime, yprime, lprime);

		if(idprime!=-1)
			add_edge(adj,id,idprime,0);
	}

	if((l!=0)&&(ivbond2d_get_value(vertices->ivbonds[l-1],x,y)==1))
	{
		int xprime=x;
		int yprime=y;
		int lprime=l-1;
		int idprime=nclusters_get_value(vertices, xprime, yprime, lprime);

		if(idprime!=-1)
			add_edge(adj,id,idprime,1);
	}

	if((x!=(vertices->lx-1))&&(ibond2d_get_value(vertices->bonds[l],x,y,DIR_X)==1))
	{
		int xprime=x+1;
		int yprime=y;
		int lprime=l;
		int idprime=nclusters_get_value(vertices, xprime, yprime, lprime);

		if(idprime!=-1)
			add_edge(adj,id,idprime,0);
	}

	if((y!=(vertices->ly-1))&&(ibond2d_get_value(vertices->bonds[l],x,y,DIR_Y)==1))
	{
		int xprime=x;
		int yprime=y+1;
		int lprime=l;
		int idprime=nclusters_get_value(vertices, xprime, yprime, lprime);

		if(idprime!=-1)
			add_edge(adj,id,idprime,0);
	}

	if((l!=(vertices->nrlayers-1))&&(ivbond2d_get_value(vertices->ivbonds[l],x,y)==1))
	{
		int xprime=x;
		int yprime=y;
		int lprime=l+1;
		int idprime=nclusters_get_value(vertices, xprime, yprime, lprime);

		if(idprime!=-1)
			add_edge(adj,id,idprime,1);
	}

	if(pbcz==true)
	{
		if((l==(vertices->nrlayers-1))&&(ivbond2d_get_value(vertices->ivbonds[l],x,y)==1))
		{
			int xprime=x;
			int yprime=y;
			int lprime=0;
			int idprime=nclusters_get_value(vertices, xprime, yprime, lprime);

			if(idprime!=-1)
				add_edge(adj,id,idprime,1);
		}
	}

	if((x==0)&&(spanning==DIR_X))
	{
		add_edge(adj,id,0,0);
	}

	if((x==(vertices->lx-1))&&(spanning==DIR_X))
	{
		add_edge(adj,id,1,0);
	}

	if((y==0)&&(spanning==DIR_Y))
	{
		add_edge(adj,id,0,0);
	}

	if((y==(vertices->ly-1))&&(spanning==DIR_Y))
	{
		add_edge(adj,id,1,0);
	}
}

/*
	A struct and some functions used to find out how many of the sites
	in the percolating cluster belong to each layer.
*/

struct layer_info_t
{
	int id;
	int bin;
};

int compar_layer_info(const void *p, const void *q)
{
	int bin1=((struct layer_info_t *)(p))->bin;
	int bin2=((struct layer_info_t *)(q))->bin;

	return (bin1 - bin2);
}

/*
	From Rosetta Code (http://rosettacode.org/wiki/Permutations/Rank_of_a_permutation#C:_Myrvold_and_Ruskey)
	some functions to determine the rank of a permutation
*/

#define SWAP(a,b) do{t=(a);(a)=(b);(b)=t;}while(0)

int mr_rank1(int n, int *vec, int *inv)
{
	int s, t;

	if (n < 2)
		return 0;

	s = vec[n-1];
	SWAP(vec[n-1], vec[inv[n-1]]);
	SWAP(inv[s], inv[n-1]);

	return s + n*mr_rank1(n-1, vec, inv);
}

int permutation_to_rank(int n, const int *vec)
{
	int i, r, *v, *inv;

	v = malloc(n * sizeof(int));
	inv = malloc(n * sizeof(int));

	for (i = 0; i < n; i++)
	{
		v[i] = vec[i];
		inv[vec[i]] = i;
	}

	r=mr_rank1(n, v, inv);

	free(inv);
	free(v);

	return r;
}
int get_permutation_bin(int nrlayers,const int bins[MAX_NR_OF_LAYERS])
{
#if 0
		/*
			Is this preliminary "rotation" needed?
		*/

		int bins2[MAX_NR_OF_LAYERS];
		int seed=gsl_rng_uniform_int(rngctx, nclusters->nrlayers);

		for(int c=0;c<nclusters->nrlayers;c++)
			bins2[c]=bins[(c+seed)%nclusters->nrlayers];

		for(int c=0;c<nclusters->nrlayers;c++)
			bins[c]=bins2[c];
#endif
	/*
		Actual algorithm: we sort the layers from the one with the most sites
		belong to the percolating clusters, to the one with the least.
	*/

	struct layer_info_t layer_infos[MAX_NR_OF_LAYERS];

	for(int c=0;c<nrlayers;c++)
	{
		layer_infos[c].id=c;
		layer_infos[c].bin=bins[c];
	}

	qsort(layer_infos, nrlayers, sizeof(struct layer_info_t), compar_layer_info);

	/*
		Then we identify the corresponding permutation
	*/

#if 0
	for(int c=0;c<nclusters->nrlayers;c++)
	{
		printf("%d[%d] ", layer_infos[c].id, layer_infos[c].bin);
	}
#endif

	int *ps=malloc(sizeof(int)*nrlayers);

	for(int c=0;c<nrlayers;c++)
		ps[c]=layer_infos[c].id;

	int permutation_rank=permutation_to_rank(nrlayers,ps);

	if(ps)
		free(ps);
#if 0
	printf(" --> %d \n",permutation_rank);
#endif

	return permutation_rank;
}

int ncluster_evaluate_jumps(struct nclusters_t *nclusters,int id,int spanning,bool pbcz,int *pbins)
{
	assert(nclusters!=NULL);
	assert(id!=0);
	assert((spanning==DIR_X)||(spanning==DIR_Y));

	/*
		We create a graph that corresponds to the cluster, augmented by two nodes:
		assuming the cluster is percolating from left to right, a start node (id=0)
		is connected to all the nodes on the left side, whereas an end node (id=1)
		is connected to all the nodes on the right side.

		All the sites in the clusters are also added as vertices of the graph, with
		id starting from 2; the	weight of every edge is 0, except for edges that
		correspond to a jump between different layers, whose weight is 1.

		It is then clear that the minimum number of jumps can be found using Dijkstra
		algorithm from the 'start' node to the 'end' nodes.
	*/

	int nr_vertices=2;

	/*
		We start by assigning a progressive ID to each vertex...
	*/

	struct nclusters_t *vertices=nclusters_init(nclusters->lx,nclusters->ly,nclusters->nrlayers);
	assert(vertices!=NULL);

	int bins[MAX_NR_OF_LAYERS]={0};

	for(int x=0;x<nclusters->lx;x++)
		for(int y=0;y<nclusters->ly;y++)
			for(int l=0;l<nclusters->nrlayers;l++)
				if(nclusters_get_value(nclusters, x, y, l)==id)
				{
					int new_vertex_id;

					/*
						One could give a new id to every site in the cluster, and the
						algorithm would still work OK.

						However, a bit of optimization, uniting already at this level
						sites on the same layer, goes a long way into saving memory and
						computation complexity later.
					*/

					if((x>0)&&(nclusters_get_value(nclusters, x-1, y, l)==id)&&(ibond2d_get_value(nclusters->bonds[l],x-1,y,DIR_X)==1))
						new_vertex_id=nclusters_get_value(vertices, x-1, y, l);
					else if((y>0)&&(nclusters_get_value(nclusters, x, y-1, l)==id)&&(ibond2d_get_value(nclusters->bonds[l],x,y-1,DIR_Y)==1))
						new_vertex_id=nclusters_get_value(vertices, x, y-1, l);
					else
						new_vertex_id=nr_vertices++;

					nclusters_set_value(vertices, x, y, l, new_vertex_id);
					bins[l]++;
				}
				else
				{
					nclusters_set_value(vertices, x, y, l, -1);
				}

	for(int c=0;c<nclusters->nrlayers;c++)
	{
		assert(nclusters->bonds[c]!=NULL);
		vertices->bonds[c]=nclusters->bonds[c];
	}

	for(int c=0;c<nclusters->nrlayers-1;c++)
	{
		assert(nclusters->ivbonds[c]!=NULL);
		vertices->ivbonds[c]=nclusters->ivbonds[c];
	}

	if(pbcz==true)
	{
		int c=nclusters->nrlayers-1;

		assert(nclusters->ivbonds[c]!=NULL);
		vertices->ivbonds[c]=nclusters->ivbonds[c];
	}

	/*
		...and then we create and populate the adjacency matrix.
	*/

	struct adjacency_t *adj=init_adjacency(nr_vertices);

	for(int x=0;x<nclusters->lx;x++)
		for(int y=0;y<nclusters->ly;y++)
			for(int l=0;l<nclusters->nrlayers;l++)
				if(nclusters_get_value(vertices,x,y,l)!=-1)
					process_neighbours(vertices,adj,x,y,l,spanning,pbcz);

	/*
		One could avoid the adjacency matrix, and calculate the weights
		on the flight in dijkstra_distance(), maybe...

		I am not sure how much faster this would be...
	*/

	int jumps=dijkstra_distance(adj, 0, 1);

	pbins[get_permutation_bin(nclusters->nrlayers,bins)]++;

	/*
		Finally, we free the vertices structure along with graph adjacency matrix
		and we return the number of jumps.
	*/

	if(vertices)
		nclusters_fini(vertices);

	if(adj)
		fini_adjacency(adj);

	return jumps;
}
