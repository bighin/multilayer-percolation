#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <limits.h>

#include <gsl/gsl_matrix.h>

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

int dijkstra_distance(gsl_matrix_int *graph, int from, int to)
{
	int V=graph->size1;

	assert(graph->size1==graph->size1);

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

			int edge_weight=gsl_matrix_int_get(graph, u, v);

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

void add_edge(gsl_matrix_int *graph, int id1, int id2, int weight)
{
	/*
		This could be redundant.
	*/

	gsl_matrix_int_set(graph,id1,id2,weight);
	gsl_matrix_int_set(graph,id2,id1,weight);
}

void process_neighbours(struct nclusters_t *vertices, gsl_matrix_int *graph, int x, int y, int l, int spanning,bool pbcz)
{
	int id=nclusters_get_value(vertices, x, y, l);

	if((x!=0)&&(ibond2d_get_value(vertices->bonds[l],x-1,y,DIR_X)==1))
	{
		int xprime=x-1;
		int yprime=y;
		int lprime=l;
		int idprime=nclusters_get_value(vertices, xprime, yprime, lprime);

		if(idprime!=-1)
			add_edge(graph,id,idprime,0);
	}

	if((y!=0)&&(ibond2d_get_value(vertices->bonds[l],x,y-1,DIR_Y)==1))
	{
		int xprime=x;
		int yprime=y-1;
		int lprime=l;
		int idprime=nclusters_get_value(vertices, xprime, yprime, lprime);

		if(idprime!=-1)
			add_edge(graph,id,idprime,0);
	}

	if((l!=0)&&(ivbond2d_get_value(vertices->ivbonds[l-1],x,y)==1))
	{
		int xprime=x;
		int yprime=y;
		int lprime=l-1;
		int idprime=nclusters_get_value(vertices, xprime, yprime, lprime);

		if(idprime!=-1)
			add_edge(graph,id,idprime,1);
	}

	if((x!=(vertices->lx-1))&&(ibond2d_get_value(vertices->bonds[l],x,y,DIR_X)==1))
	{
		int xprime=x+1;
		int yprime=y;
		int lprime=l;
		int idprime=nclusters_get_value(vertices, xprime, yprime, lprime);

		if(idprime!=-1)
			add_edge(graph,id,idprime,0);
	}

	if((y!=(vertices->ly-1))&&(ibond2d_get_value(vertices->bonds[l],x,y,DIR_Y)==1))
	{
		int xprime=x;
		int yprime=y+1;
		int lprime=l;
		int idprime=nclusters_get_value(vertices, xprime, yprime, lprime);

		if(idprime!=-1)
			add_edge(graph,id,idprime,0);
	}

	if((l!=(vertices->nrlayers-1))&&(ivbond2d_get_value(vertices->ivbonds[l],x,y)==1))
	{
		int xprime=x;
		int yprime=y;
		int lprime=l+1;
		int idprime=nclusters_get_value(vertices, xprime, yprime, lprime);

		if(idprime!=-1)
			add_edge(graph,id,idprime,1);
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
				add_edge(graph,id,idprime,1);
		}
	}

	if((x==0)&&(spanning==DIR_X))
	{
		add_edge(graph,id,0,0);
	}

	if((x==(vertices->lx-1))&&(spanning==DIR_X))
	{
		add_edge(graph,id,1,0);
	}

	if((y==0)&&(spanning==DIR_Y))
	{
		add_edge(graph,id,0,0);
	}

	if((y==(vertices->ly-1))&&(spanning==DIR_Y))
	{
		add_edge(graph,id,1,0);
	}
}

int ncluster_evaluate_jumps(struct nclusters_t *nclusters,int id,int spanning,bool pbcz)
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

	for(int x=0;x<nclusters->lx;x++)
		for(int y=0;y<nclusters->ly;y++)
			for(int l=0;l<nclusters->nrlayers;l++)
				if(nclusters_get_value(nclusters, x, y, l)==id)
				{
					int new_vertex_id;

					/*
						One could give a new id to every site in the cluster, and the
						algorithm would still work OK.

						However, a bit of optimization in uniting already at this level
						sites on the same layer goes a long way into saving memory and
						computation complexity later.
					*/

					if((x>0)&&(nclusters_get_value(nclusters, x-1, y, l)==id))
						new_vertex_id=nclusters_get_value(vertices, x-1, y, l);
					else if((y>0)&&(nclusters_get_value(nclusters, x, y-1, l)==id))
						new_vertex_id=nclusters_get_value(vertices, x, y-1, l);
					else
						new_vertex_id=nr_vertices++;

					nclusters_set_value(vertices, x, y, l, new_vertex_id);
				}
				else
					nclusters_set_value(vertices, x, y, l, -1);

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

	gsl_matrix_int *graph=gsl_matrix_int_alloc(nr_vertices,nr_vertices);

	for(int i=0;i<nr_vertices;i++)
		for(int j=0;j<nr_vertices;j++)
			gsl_matrix_int_set(graph,i,j,INT_MAX);

	for(int x=0;x<nclusters->lx;x++)
		for(int y=0;y<nclusters->ly;y++)
			for(int l=0;l<nclusters->nrlayers;l++)
				if(nclusters_get_value(vertices,x,y,l)!=-1)
					process_neighbours(vertices,graph,x,y,l,spanning,pbcz);

	/*
		One could avoid the adjacency matrix, and calculate the weights
		on the flight in dijkstra_distance(), maybe...

		I am not sure how faster this would be...
	*/

	int jumps=dijkstra_distance(graph, 0, 1);

	/*
		Finally, we free the vertices structure along with graph adjacency matrix
		and we return the number of jumps.
	*/

	if(vertices)
		nclusters_fini(vertices);

	if(graph)
		gsl_matrix_int_free(graph);

	return jumps;
}
