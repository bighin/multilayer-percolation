#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

#include <gsl/gsl_rng.h>

#include "bonds.h"
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

int do_run_bond(int xdim,int ydim,int zdim,double p,double pperp,gsl_rng *rng,int *jumps,int *matches1,int *matches2,
		int *nr_percolating1,int *nr_percolating2,bool pbcz)
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

	for(int z=0;z<zdim;z++)
	{
		/*
			If we do not have periodic boundary conditions in the z direction,
			then there are no vertical bonds joining the last and the first layer.
		*/

		if((z==(zdim-1))&&(pbcz==false))
		{
			ncs->ivbonds[z]=NULL;
			continue;
		}

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
		The clusters are identified and measured.

		First: clusters that can span on more than one layer.
	*/

	if((*nr_percolating1=nclusters_identify_percolation(ncs,jumps,matches1,rng,pbcz))>0)
			result|=TWO_LAYER_PERCOLATION;

	/*
		Second: the vertical links are removed, so that now we look for
		percolation of clusters living only a single layer.
	*/

	for(int z=0;z<zdim;z++)
	{
		if((z==(zdim-1))&&(pbcz==false))
			continue;

		for(int x=0;x<xdim;x++)
			for(int y=0;y<ydim;y++)
				ivbond2d_set_value(ncs->ivbonds[z], x, y, 0);
	}

	if((*nr_percolating2=nclusters_identify_percolation(ncs,NULL,matches2,rng,pbcz))>0)
		result|=SINGLE_LAYER_PERCOLATION;

	/*
		Final cleanup.
	*/

	for(int z=0;z<zdim;z++)
	{
		ibond2d_fini(ncs->bonds[z]);
		
		if((z!=(zdim-1))||(pbcz==true))
			ivbond2d_fini(ncs->ivbonds[z]);
	}

	nclusters_fini(ncs);

	return result;
}

struct config_t
{
	int total_runs;
	int xdim,ydim,nrlayers;

	bool measure_jumps;
	bool pbcz;

	int mincentipperp,maxcentipperp;
	int mincentip,maxcentip;
};

void do_batch(struct config_t *config,char *outfile)
{
	FILE *out=fopen(outfile,"w+");
	assert(out);

	setvbuf(out,(char *)(NULL),_IONBF,0);

#ifdef NDEBUG
#pragma omp parallel for collapse(2) schedule(dynamic)
#endif

	for(int centipperp=config->mincentipperp;centipperp<=config->maxcentipperp;centipperp+=1)
	{
		for(int centip=config->mincentip;centip<=config->maxcentip;centip+=1)
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
			int totalpercolating1=0,totalpercolating2=0;

			cntsingle=cntbilayer=0;
			for(int c=0;c<config->total_runs;c++)
			{
				int jumps=0,matches1=0,matches2=0;
				int nr_percolating1=0,nr_percolating2=0;

				int *pjumps=(config->measure_jumps==true)?(&jumps):(NULL);

				switch(do_run_bond(config->xdim,config->ydim,config->nrlayers,p,pperp,rng_ctx,pjumps,&matches1,&matches2,&nr_percolating1,&nr_percolating2,config->pbcz))
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
				totalpercolating1+=nr_percolating1;
				totalpercolating2+=nr_percolating2;
			}

			gsl_rng_free(rng_ctx);

#pragma omp critical
			{
				fprintf(out,"%f %f %f %f %f %f %f %f %f\n",p,pperp,((double)(cntbilayer))/((double)(config->total_runs)),
	                                                                           ((double)(cntsingle))/((double)(config->total_runs)),
	                                                                           ((double)(totaljumps))/((double)(config->total_runs)),
					                                           ((double)(totalmatches1))/((double)(config->total_runs)),
					                                           ((double)(totalmatches2))/((double)(config->total_runs)),
									           ((double)(totalpercolating1))/((double)(config->total_runs)),
									           ((double)(totalpercolating2))/((double)(config->total_runs)));
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
	config.measure_jumps=false;
	config.pbcz=false;
	config.mincentipperp=0;
	config.maxcentipperp=100;
	config.mincentip=0;
	config.maxcentip=100;

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
	config.mincentipperp=50;
	config.maxcentipperp=50;

	config.xdim=config.ydim=16;
	config.nrlayers=2;
	do_batch(&config,"jumps16.dat");

	config.xdim=config.ydim=32;
	config.nrlayers=2;
	do_batch(&config,"jumps32.dat");

	config.xdim=config.ydim=64;
	config.nrlayers=2;
	do_batch(&config,"jumps64.dat");

	config.xdim=config.ydim=128;
	config.nrlayers=2;
	do_batch(&config,"jumps128.dat");

	config.xdim=config.ydim=256;
	config.nrlayers=2;
	do_batch(&config,"jumps256.dat");

	/*
		Same, with periodic boundary conditions along z.
	*/

	config.pbcz=true;

	config.xdim=config.ydim=512;
	config.nrlayers=3;
	do_batch(&config,"trilayer512_pbcz.dat");

	config.xdim=config.ydim=512;
	config.nrlayers=6;
	do_batch(&config,"esalayer512_pbcz.dat");

	config.xdim=config.ydim=16;
	config.nrlayers=2;
	do_batch(&config,"bilayer16_pbcz.dat");

	config.xdim=config.ydim=32;
	config.nrlayers=2;
	do_batch(&config,"bilayer32_pbcz.dat");

	config.xdim=config.ydim=64;
	config.nrlayers=2;
	do_batch(&config,"bilayer64_pbcz.dat");

	config.xdim=config.ydim=128;
	config.nrlayers=2;
	do_batch(&config,"bilayer128_pbcz.dat");

	config.xdim=config.ydim=256;
	config.nrlayers=2;
	do_batch(&config,"bilayer256_pbcz.dat");

	config.xdim=config.ydim=512;
	config.nrlayers=2;
	do_batch(&config,"bilayer512_pbcz.dat");

	config.measure_jumps=true;
	config.mincentipperp=50;
	config.maxcentipperp=50;

	config.xdim=config.ydim=16;
	config.nrlayers=2;
	do_batch(&config,"jumps16_pbcz.dat");

	config.xdim=config.ydim=32;
	config.nrlayers=2;
	do_batch(&config,"jumps32_pbcz.dat");

	config.xdim=config.ydim=64;
	config.nrlayers=2;
	do_batch(&config,"jumps64_pbcz.dat");

	config.xdim=config.ydim=128;
	config.nrlayers=2;
	do_batch(&config,"jumps128_pbcz.dat");

	config.xdim=config.ydim=256;
	config.nrlayers=2;
	do_batch(&config,"jumps256_pbcz.dat");

	return 0;
}
