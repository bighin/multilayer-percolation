#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

#include <gsl/gsl_rng.h>

#include "bonds.h"
#include "sites.h"
#include "clusters.h"

void seed_rng(gsl_rng *rng)
{
	char *devname="/dev/urandom";
	FILE *dev;

	if((dev=fopen(devname,"r"))!=NULL)
	{
		unsigned long seed;

		fread(&seed,sizeof(unsigned long),1,dev);
		fclose(dev);

		gsl_rng_set(rng,seed);
	}
	else
	{
		printf("Warning: couldn't read from %s to seed the RNG.\n",devname);
	}
}

int get_random_value(double p,gsl_rng *rng_ctx)
{
	if(gsl_rng_uniform(rng_ctx)<p)
		return 1;
	
	return 0;
}

#define TWO_LAYER_PERCOLATION		(1)
#define SINGLE_LAYER_PERCOLATION	(2)

int do_run_bond(int xdim,int ydim,int zdim,double p,double pperp,gsl_rng *rng,int *jumps,int *matches1,int *matches2)
{
	struct nclusters_t *ncs;
	int result=0;
	
	ncs=nclusters_init(xdim,ydim,zdim);
	assert(ncs);
	
	/*
		The random bonds are created...
	*/

	for(int z=0;z<zdim;z++)
	{
		ncs->bonds[z]=ibond2d_init(xdim,ydim);

		for(int x=0;x<xdim;x++)
		{
			for(int y=0;y<ydim;y++)
			{
				ibond2d_set_value(ncs->bonds[z],x,y,DIR_X,get_random_value(p,rng));
				ibond2d_set_value(ncs->bonds[z],x,y,DIR_Y,get_random_value(p,rng));
			}
		}
	}

	for(int z=0;z<(zdim-1);z++)
	{
		ncs->ivbonds[z]=ivbond2d_init(xdim,ydim);

		for(int x=0;x<xdim;x++)
		{
			for(int y=0;y<ydim;y++)
			{
				ivbond2d_set_value(ncs->ivbonds[z],x,y,get_random_value(pperp,rng));
			}
		}
	}

	/*
		...and then the clusters are identified and measured!
	
		First: clusters that can span on more than one layer.
	*/

	if(nclusters_identify_percolation(ncs,jumps,matches1,rng)==true)
			result|=TWO_LAYER_PERCOLATION;

	/*
		Second: the vertical links are removed, so that now we look for
		percolation of clusters living only a single layer.
	*/

	for(int z=0;z<(zdim-1);z++)
		for(int x=0;x<xdim;x++)
			for(int y=0;y<ydim;y++)
				ivbond2d_set_value(ncs->ivbonds[z],x,y,0);

	if(nclusters_identify_percolation(ncs,NULL,matches2,rng)==true)
		result|=SINGLE_LAYER_PERCOLATION;

	/*
		Final cleanup.
	*/

	for(int z=0;z<zdim;z++)
	{
		ibond2d_fini(ncs->bonds[z]);
		
		if(z!=(zdim-1))
			ivbond2d_fini(ncs->ivbonds[z]);
	}

	nclusters_fini(ncs);

	return result;
}

int do_run_site(int xdim,int ydim,int zdim,double p,double pmiddle,gsl_rng *rng,int *jumps,int *matches1,int *matches2)
{
	struct sites2d_t *top,*middle,*bottom;
	
	top=sites2d_init(xdim,ydim);
	middle=sites2d_init(xdim,ydim);
	bottom=sites2d_init(xdim,ydim);
	
	/*
		We start by creating the top, middle and bottom layers.
	*/

	for(int x=0;x<xdim;x++)
	{
		for(int y=0;y<ydim;y++)
		{
			sites2d_set(top,x,y,get_random_value(p,rng));
			sites2d_set(middle,x,y,get_random_value(pmiddle,rng));
			sites2d_set(bottom,x,y,get_random_value(p,rng));
		}
	}
	
	/*
		Then we map the problem to a bond percolation problem
		on a bilayer of dimensions (xdim-1), (ydim-1).
	*/

	struct nclusters_t *ncs;
	int result=0;
	
	ncs=nclusters_init(xdim-1,ydim-1,2);
	assert(ncs);

	ncs->bonds[0]=ibond2d_init(xdim-1,ydim-1);
	ncs->bonds[1]=ibond2d_init(xdim-1,ydim-1);

	for(int x=0;x<xdim-1;x++)
	{
		for(int y=0;y<ydim-1;y++)
		{
			ibond2d_set_value(ncs->bonds[0],x,y,DIR_X,sites2d_get(bottom,x,y)&sites2d_get(bottom,x+1,y));
			ibond2d_set_value(ncs->bonds[0],x,y,DIR_Y,sites2d_get(bottom,x,y)&sites2d_get(bottom,x,y+1));

			ibond2d_set_value(ncs->bonds[1],x,y,DIR_X,sites2d_get(top,x,y)&sites2d_get(top,x+1,y));
			ibond2d_set_value(ncs->bonds[1],x,y,DIR_Y,sites2d_get(top,x,y)&sites2d_get(top,x,y+1));
		}
	}

	ncs->ivbonds[0]=ivbond2d_init(xdim-1,ydim-1);

	for(int x=0;x<xdim-1;x++)
	{
		for(int y=0;y<ydim-1;y++)
		{
			ivbond2d_set_value(ncs->ivbonds[0],x,y,sites2d_get(top,x,y)&sites2d_get(middle,x,y)&sites2d_get(bottom,x,y));
		}
	}

	/*
		At last, the clusters are identified and measured!
	
		First: clusters that can span on more than one layer.
	*/

	if(nclusters_identify_percolation(ncs,jumps,matches1,rng)==true)
			result|=TWO_LAYER_PERCOLATION;

	/*
		Second: the vertical links are removed, so that now we look for
		percolation of clusters living only a single layer.
	*/

	for(int x=0;x<xdim-1;x++)
		for(int y=0;y<ydim-1;y++)
			ivbond2d_set_value(ncs->ivbonds[0],x,y,0);

	if(nclusters_identify_percolation(ncs,NULL,matches2,rng)==true)
		result|=SINGLE_LAYER_PERCOLATION;

	/*
		Final cleanup.
	*/

	for(int z=0;z<=1;z++)
		ibond2d_fini(ncs->bonds[z]);

	ivbond2d_fini(ncs->ivbonds[0]);		

	nclusters_fini(ncs);

	sites2d_fini(top);
	sites2d_fini(middle);
	sites2d_fini(bottom);

	return result;
}

#define BOND_PERCOLATION	(19)
#define SITE_PERCOLATION	(20)

struct config_t
{
	int total_runs;
	int xdim,ydim,nrlayers;

	int type;
	bool measure_jumps;

};

void do_batch(struct config_t *config,char *outfile)
{
	FILE *out=fopen(outfile,"w+");
	assert(out);

	setvbuf(out,(char *)(NULL),_IONBF,0);

#ifdef NDEBUG
#pragma omp parallel for collapse(2) schedule(dynamic)
#endif

	for(int centipperp=0;centipperp<=100;centipperp+=1)
	{
		for(int centip=0;centip<=100;centip+=1)
		{
			double p,pperp;
			gsl_rng *rng_ctx;

			p=0.01*centip;
			pperp=0.01*centipperp;

			rng_ctx=gsl_rng_alloc(gsl_rng_mt19937);
			assert(rng_ctx!=NULL);
			seed_rng(rng_ctx);

			int cntsingle,cntbilayer;
			int totaljumps=0,totalmatches1=0,totalmatches2=0;

			cntsingle=cntbilayer=0;
			for(int c=0;c<config->total_runs;c++)
			{
				int jumps=0,matches1=0,matches2=0;

				int (*do_run_ptr)(int,int,int,double,double,gsl_rng *,int *,int *,int *);

				switch(config->type)
				{
					case BOND_PERCOLATION:
					do_run_ptr=&do_run_bond;
					break;

					case SITE_PERCOLATION:
					do_run_ptr=&do_run_site;
					break;

					default:
					assert(false);
				}

				int *pjumps=(config->measure_jumps==true)?(&jumps):(NULL);
				int *pmatches1=&matches1;
				int *pmatches2=&matches2;

				switch(do_run_ptr(config->xdim,config->ydim,config->nrlayers,p,pperp,rng_ctx,pjumps,pmatches1,pmatches2))
				{
					case 0:
					break;
					
					case TWO_LAYER_PERCOLATION:
					cntbilayer++;
					break;

					case SINGLE_LAYER_PERCOLATION:
					cntsingle++;
					break;

					case SINGLE_LAYER_PERCOLATION|TWO_LAYER_PERCOLATION:
					cntbilayer++;
					cntsingle++;
					break;
				}

				totaljumps+=jumps;
				totalmatches1+=matches1;
				totalmatches2+=matches2;
			}

			gsl_rng_free(rng_ctx);

#pragma omp critical
			{
				fprintf(out,"%f %f %f %f %f %f %f\n",p,pperp,((double)(cntbilayer))/((double)(config->total_runs)),
	                                                                     ((double)(cntsingle))/((double)(config->total_runs)),
	                                                                     ((double)(totaljumps))/((double)(config->total_runs)),
					                                     ((double)(totalmatches1))/((double)(config->total_runs)),
					                                     ((double)(totalmatches2))/((double)(config->total_runs)));
				fflush(out);
			}
		}
	}

	if(out)
		fclose(out);
}

void do_mini_batch(struct config_t *config,char *outfile)
{
	FILE *out=fopen(outfile,"w+");
	assert(out);

	setvbuf(out,(char *)(NULL),_IONBF,0);

#ifdef NDEBUG
#pragma omp parallel for collapse(2) schedule(dynamic)
#endif

	for(int centipperp=50;centipperp<=50;centipperp+=1)
	{
		for(int centip=0;centip<=100;centip+=1)
		{
			double p,pperp;
			gsl_rng *rng_ctx;

			p=0.01*centip;
			pperp=0.01*centipperp;

			rng_ctx=gsl_rng_alloc(gsl_rng_mt19937);
			assert(rng_ctx!=NULL);
			seed_rng(rng_ctx);

			int cntsingle,cntbilayer;
			int totaljumps=0,totalmatches1=0,totalmatches2=0;

			cntsingle=cntbilayer=0;
			for(int c=0;c<config->total_runs;c++)
			{
				int jumps=0,matches1=0,matches2=0;

				int (*do_run_ptr)(int,int,int,double,double,gsl_rng *,int *,int *,int *);

				switch(config->type)
				{
					case BOND_PERCOLATION:
					do_run_ptr=&do_run_bond;
					break;

					case SITE_PERCOLATION:
					do_run_ptr=&do_run_site;
					break;

					default:
					assert(false);
				}

				int *pjumps=(config->measure_jumps==true)?(&jumps):(NULL);
				int *pmatches1=&matches1;
				int *pmatches2=&matches2;

				switch(do_run_ptr(config->xdim,config->ydim,config->nrlayers,p,pperp,rng_ctx,pjumps,pmatches1,pmatches2))
				{
					case 0:
					break;
					
					case TWO_LAYER_PERCOLATION:
					cntbilayer++;
					break;

					case SINGLE_LAYER_PERCOLATION:
					cntsingle++;
					break;

					case SINGLE_LAYER_PERCOLATION|TWO_LAYER_PERCOLATION:
					cntbilayer++;
					cntsingle++;
					break;
				}

				totaljumps+=jumps;
				totalmatches1+=matches1;
				totalmatches2+=matches2;
			}

			gsl_rng_free(rng_ctx);

#pragma omp critical
			{
				fprintf(out,"%f %f %f %f %f %f %f\n",p,pperp,((double)(cntbilayer))/((double)(config->total_runs)),
					                                     ((double)(cntsingle))/((double)(config->total_runs)),
					                                     ((double)(totaljumps))/((double)(config->total_runs)),
					                                     ((double)(totalmatches1))/((double)(config->total_runs)),
					                                     ((double)(totalmatches2))/((double)(config->total_runs)));
				fflush(out);
			}
		}
	}

	if(out)
		fclose(out);
}

int main(int argc,char *argv[])
{
	struct config_t config;

	config.total_runs=100;
	config.type=BOND_PERCOLATION;
	config.measure_jumps=false;

	config.xdim=config.ydim=512;
	config.nrlayers=3;
	do_batch(&config,"trilayer512.dat");

	config.xdim=config.ydim=512;
	config.nrlayers=6;
	do_batch(&config,"esalayer512.dat");

	config.xdim=config.ydim=16;
	config.nrlayers=2;
	do_batch(&config,"bilayer16.dat");

	config.xdim=config.ydim=32;
	config.nrlayers=2;
	do_batch(&config,"bilayer32.dat");

	config.xdim=config.ydim=64;
	config.nrlayers=2;
	do_batch(&config,"bilayer64.dat");

	config.xdim=config.ydim=128;
	config.nrlayers=2;
	do_batch(&config,"bilayer128.dat");

	config.xdim=config.ydim=256;
	config.nrlayers=2;
	do_batch(&config,"bilayer256.dat");

	config.xdim=config.ydim=512;
	config.nrlayers=2;
	do_batch(&config,"bilayer512.dat");

	config.measure_jumps=true;

	config.xdim=config.ydim=16;
	config.nrlayers=2;
	do_mini_batch(&config,"jumps16.dat");

	config.xdim=config.ydim=32;
	config.nrlayers=2;
	do_mini_batch(&config,"jumps32.dat");

	config.xdim=config.ydim=64;
	config.nrlayers=2;
	do_mini_batch(&config,"jumps64.dat");

	config.xdim=config.ydim=128;
	config.nrlayers=2;
	do_mini_batch(&config,"jumps128.dat");

	config.xdim=config.ydim=256;
	config.nrlayers=2;
	do_mini_batch(&config,"jumps256.dat");

	return 0;
}
