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

struct config_t
{
	int total_runs;
	int xdim,ydim,nrlayers;

	bool measure_jumps;
	bool pbcz;

	int mincentipperp,maxcentipperp;
	int mincentip,maxcentip;
};

struct statistics_t
{
	int cntsingle;
	int cntbilayer;

	int jumps;
	int matches1;
	int matches2;

	int nr_percolating1;
	int nr_percolating2;
};

void reset_stats(struct statistics_t *st)
{
	st->cntsingle=0;
	st->cntbilayer=0;

	st->jumps=0;
	st->matches1=0;
	st->matches2=0;

	st->nr_percolating1=0;
	st->nr_percolating2=0;
}

void add_stats(struct statistics_t *total,struct statistics_t *st)
{
	total->cntsingle+=st->cntsingle;
	total->cntbilayer+=st->cntbilayer;

	total->jumps+=st->jumps;
	total->matches1+=st->matches1;
	total->matches2+=st->matches2;

	total->nr_percolating1+=st->nr_percolating1;
	total->nr_percolating2+=st->nr_percolating2;
}

#define TWO_LAYER_PERCOLATION		(1)
#define SINGLE_LAYER_PERCOLATION	(2)

int do_run(struct config_t *config,double p,double pperp,gsl_rng *rng,struct statistics_t *stat)
{
	int xdim=config->xdim;
	int ydim=config->ydim;
	int zdim=config->nrlayers;

	struct nclusters_t *ncs=nclusters_init(xdim,ydim,zdim);
	assert(ncs);

	int result=0;

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

		if((z==(zdim-1))&&(config->pbcz==false))
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

		First: clusters that can span more than one layer.
	*/

	int *pjumps=(config->measure_jumps==true)?(&stat->jumps):(NULL);

	if((stat->nr_percolating1=nclusters_identify_percolation(ncs,pjumps,&stat->matches1,rng,config->pbcz))>0)
			result|=TWO_LAYER_PERCOLATION;

	/*
		Second: the vertical links are removed, so that now we look for
		percolation of clusters living only a single layer.
	*/

	for(int z=0;z<zdim;z++)
	{
		if((z==(zdim-1))&&(config->pbcz==false))
			continue;

		for(int x=0;x<xdim;x++)
			for(int y=0;y<ydim;y++)
				ivbond2d_set_value(ncs->ivbonds[z], x, y, 0);
	}

	if((stat->nr_percolating2=nclusters_identify_percolation(ncs,NULL,&stat->matches2,rng,config->pbcz))>0)
		result|=SINGLE_LAYER_PERCOLATION;

	/*
		Final cleanup.
	*/

	for(int z=0;z<zdim;z++)
	{
		ibond2d_fini(ncs->bonds[z]);

		if((z!=(zdim-1))||(config->pbcz==true))
			ivbond2d_fini(ncs->ivbonds[z]);
	}

	nclusters_fini(ncs);

	return result;
}

void do_batch(struct config_t *config,char *outfile)
{
	FILE *out=fopen(outfile,"w+");
	assert(out);

	setvbuf(out,(char *)(NULL),_IONBF,0);

#ifdef NDEBUG
#pragma omp parallel for collapse(2) schedule(dynamic) default(none) shared(config,out,gsl_rng_mt19937)
#endif

	for(int centipperp=config->mincentipperp;centipperp<=config->maxcentipperp;centipperp+=1)
	{
		for(int centip=config->mincentip;centip<=config->maxcentip;centip+=1)
		{
			double p=0.01*centip;
			double pperp=0.01*centipperp;

			gsl_rng *rng_ctx=gsl_rng_alloc(gsl_rng_mt19937);
			assert(rng_ctx!=NULL);
			seed_rng(rng_ctx);

			struct statistics_t total;
			reset_stats(&total);

			for(int c=0;c<config->total_runs;c++)
			{
				struct statistics_t stats;
				reset_stats(&stats);

				switch(do_run(config, p, pperp, rng_ctx, &stats))
				{
					case 0:
					break;
					
					case TWO_LAYER_PERCOLATION:
					stats.cntbilayer++;
					break;

					case SINGLE_LAYER_PERCOLATION:
					stats.cntsingle++;
					break;

					case SINGLE_LAYER_PERCOLATION|TWO_LAYER_PERCOLATION:
					stats.cntbilayer++;
					stats.cntsingle++;
					break;
				}

				add_stats(&total,&stats);
			}

			gsl_rng_free(rng_ctx);

#pragma omp critical
			{
				fprintf(out,"%f %f ",p,pperp);
				fprintf(out,"%f ",((double)(total.cntbilayer))/((double)(config->total_runs)));
				fprintf(out,"%f ",((double)(total.cntsingle))/((double)(config->total_runs)));
				fprintf(out,"%f ",((double)(total.jumps))/((double)(config->total_runs)));
				fprintf(out,"%f ",((double)(total.matches1))/((double)(config->total_runs)));
				fprintf(out,"%f ",((double)(total.matches2))/((double)(config->total_runs)));
				fprintf(out,"%f ",((double)(total.nr_percolating1))/((double)(config->total_runs)));
				fprintf(out,"%f\n",((double)(total.nr_percolating2))/((double)(config->total_runs)));

				fflush(out);
			}
		}
	}

	if(out)
		fclose(out);
}

int go(int id)
{
	struct config_t config;

	config.total_runs=100;
	config.measure_jumps=false;
	config.mincentipperp=0;
	config.maxcentipperp=100;
	config.mincentip=0;
	config.maxcentip=100;

	switch(id)
	{
		case 1:
		config.pbcz=false;
		config.xdim=config.ydim=512;
		config.nrlayers=3;
		do_batch(&config, "trilayer512.dat");
		break;

		case 2:
		config.pbcz=false;
		config.xdim=config.ydim=512;
		config.nrlayers=6;
		do_batch(&config, "esalayer512.dat");
		break;

		case 3:
		config.pbcz=false;
		config.xdim=config.ydim=16;
		config.nrlayers=2;
		do_batch(&config, "bilayer16.dat");
		break;

		case 4:
		config.pbcz=false;
		config.xdim=config.ydim=32;
		config.nrlayers=2;
		do_batch(&config, "bilayer32.dat");
		break;

		case 5:
		config.pbcz=false;
		config.xdim=config.ydim=64;
		config.nrlayers=2;
		do_batch(&config, "bilayer64.dat");
		break;

		case 6:
		config.pbcz=false;
		config.xdim=config.ydim=128;
		config.nrlayers=2;
		do_batch(&config, "bilayer128.dat");
		break;

		case 7:
		config.pbcz=false;
		config.xdim=config.ydim=256;
		config.nrlayers=2;
		do_batch(&config, "bilayer256.dat");
		break;

		case 8:
		config.pbcz=false;
		config.xdim=config.ydim=512;
		config.nrlayers=2;
		do_batch(&config, "bilayer512.dat");
		break;

		case 9:
		config.pbcz=false;
		config.measure_jumps=true;
		config.total_runs=1000;
		config.mincentipperp=50;
		config.maxcentipperp=50;
		config.xdim=config.ydim=16;
		config.nrlayers=2;
		do_batch(&config, "jumps16.dat");
		break;

		case 10:
		config.pbcz=false;
		config.measure_jumps=true;
		config.total_runs=1000;
		config.mincentipperp=50;
		config.maxcentipperp=50;
		config.xdim=config.ydim=32;
		config.nrlayers=2;
		do_batch(&config, "jumps32.dat");
		break;

		case 11:
		config.pbcz=false;
		config.measure_jumps=true;
		config.total_runs=1000;
		config.mincentipperp=50;
		config.maxcentipperp=50;
		config.xdim=config.ydim=64;
		config.nrlayers=2;
		do_batch(&config, "jumps64.dat");
		break;

		case 12:
		config.pbcz=false;
		config.measure_jumps=true;
		config.total_runs=1000;
		config.mincentipperp=50;
		config.maxcentipperp=50;
		config.xdim=config.ydim=128;
		config.nrlayers=2;
		do_batch(&config, "jumps128.dat");
		break;

		case 13:
		config.pbcz=false;
		config.measure_jumps=true;
		config.total_runs=1000;
		config.mincentipperp=50;
		config.maxcentipperp=50;
		config.xdim=config.ydim=256;
		config.nrlayers=2;
		do_batch(&config, "jumps256.dat");
		break;

		/*
			Same, with periodic boundary conditions along z.
		*/

		case 14:
		config.pbcz=true;
		config.xdim=config.ydim=512;
		config.nrlayers=3;
		do_batch(&config, "trilayer512_pbcz.dat");
		break;

		case 15:
		config.pbcz=true;
		config.xdim=config.ydim=512;
		config.nrlayers=6;
		do_batch(&config, "esalayer512_pbcz.dat");
		break;

		case 16:
		config.pbcz=true;
		config.xdim=config.ydim=16;
		config.nrlayers=2;
		do_batch(&config, "bilayer16_pbcz.dat");
		break;

		case 17:
		config.pbcz=true;
		config.xdim=config.ydim=32;
		config.nrlayers=2;
		do_batch(&config, "bilayer32_pbcz.dat");
		break;

		case 18:
		config.pbcz=true;
		config.xdim=config.ydim=64;
		config.nrlayers=2;
		do_batch(&config, "bilayer64_pbcz.dat");
		break;

		case 19:
		config.pbcz=true;
		config.xdim=config.ydim=128;
		config.nrlayers=2;
		do_batch(&config, "bilayer128_pbcz.dat");
		break;

		case 20:
		config.pbcz=true;
		config.xdim=config.ydim=256;
		config.nrlayers=2;
		do_batch(&config, "bilayer256_pbcz.dat");
		break;

		case 21:
		config.pbcz=true;
		config.xdim=config.ydim=512;
		config.nrlayers=2;
		do_batch(&config, "bilayer512_pbcz.dat");
		break;

		case 22:
		config.pbcz=true;
		config.measure_jumps=true;
		config.total_runs=1000;
		config.mincentipperp=50;
		config.maxcentipperp=50;
		config.xdim=config.ydim=16;
		config.nrlayers=2;
		do_batch(&config, "jumps16_pbcz.dat");
		break;

		case 23:
		config.pbcz=true;
		config.measure_jumps=true;
		config.total_runs=1000;
		config.mincentipperp=50;
		config.maxcentipperp=50;
		config.xdim=config.ydim=32;
		config.nrlayers=2;
		do_batch(&config, "jumps32_pbcz.dat");
		break;

		case 24:
		config.pbcz=true;
		config.measure_jumps=true;
		config.total_runs=1000;
		config.mincentipperp=50;
		config.maxcentipperp=50;
		config.xdim=config.ydim=64;
		config.nrlayers=2;
		do_batch(&config, "jumps64_pbcz.dat");
		break;

		case 25:
		config.pbcz=true;
		config.measure_jumps=true;
		config.total_runs=1000;
		config.mincentipperp=50;
		config.maxcentipperp=50;
		config.xdim=config.ydim=128;
		config.nrlayers=2;
		do_batch(&config, "jumps128_pbcz.dat");
		break;

		case 26:
		config.pbcz=true;
		config.measure_jumps=true;
		config.total_runs=1000;
		config.mincentipperp=50;
		config.maxcentipperp=50;
		config.xdim=config.ydim=256;
		config.nrlayers=2;
		do_batch(&config, "jumps256_pbcz.dat");
		break;

		case 101:
		config.pbcz=false;
		config.measure_jumps=true;
		config.mincentipperp=26;
		config.maxcentipperp=26;
		config.mincentipperp=50;
		config.maxcentipperp=50;
		config.xdim=config.ydim=16;
		config.nrlayers=2;
		do_batch(&config, "jumps16debug.dat");
		break;

		default:
		break;
	}

	return 0;
}

int main(int argc,char *argv[])
{
	if(argc!=2)
		return 0;

	int id=atoi(argv[1]);

	return go(id);
}
