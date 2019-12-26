/*	------------------------------------------------------------------------------------

                            * DBC454 *
	unbiased parallel density based clustering of large scale ITS data


    Copyright (C) SIB  - Swiss Institute of Bioinformatics,   2013-2019 Nicolas Guex
    Copyright (C) UNIL - University of Lausanne, Switzerland       2019 Nicolas Guex


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.


	Code:       Nicolas Guex, 2013-2019
	Contact:    Nicolas.Guex@unil.ch
	Repository: https://github.com/sib-swiss/dbc454


	Article:    Density-based hierarchical clustering of pyro-sequences
	            on a large scale—the case of fungal ITS1

                   Bioinformatics. 2013 May 15; 29(10): 1268–1274.
                   https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtt149




	Machine :	Unix
	Language:	C
	Requires:	mpi, pthread

	Version information

	Version:	1.4.3	Jan.  2013	Article published in Bioinformatics (see above).
	Version:	1.4.4	Dec.  2019	Public release of code under GPL2+ license




	Compiling:   (you will need mpi on your system)
	
	mpicc -O3 -o dbc454 dbc454.c -Wall -lpthread -lm



	Testing:

	./unit_test.sh


	------------------------------------------------------------------------------------
*/


	/*------------------------- I N T E R F A C E ----------------------- */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <limits.h>
#include <time.h>
#include "mpi.h"
#include <math.h>
#include <ctype.h>
#include <xlocale.h>
#include <unistd.h>
#include <errno.h>


/* ------------------------------------------------------------------------------------ */

#define USE_THREADS


#ifdef USE_THREADS
#include <pthread.h>
#endif

/* ------------------------------------------------------------------------------------ */


#define kMaxLineBuf 16384
#define kMAXEVENTS 10000000                             /* 10 millions */
#define kMaxCluster 1000000
#define kMaxInputCol 48
#define kMaxConflict 8

#define kLoadingProgressReporting 100000

#define kMaxMergeRequests 786432	/* 6144 Kb */
#define kMaxCPU 64
#define kStartLocalCluster 16000000  /* 16 millions in practice, allows for up to 128 processors but each cpu should not reach more than 16 million distinct clusters */

/* status of computations chunks */
#define kChunkStatusToDo 0
#define kChunkStatusComputing 1

/* MPI messages */
#define kWhichBlocksToCompute 0
#define kClusterMsg1 1
#define kClusterMsg2 2
#define kCPUdoneMsg 4
#define kFinalCntRequest 8
#define kFinalMergeRequestCntRequest 16
#define kSendFinalMergeRequestRequest 32
#define kJoinListRequest 64
#define kJoinListRequestDone 128
#define kRepeatWithNewDistMsg 512
#define kInitialClusterCntMsg 1024

#define kRenameClusterCntRequest 8192
#define kRenameClusterRequest 16384
#define kClosestMsg 32768

/* status of cpu nodes */
#define kCPU_availaible 2147483647
#define kNoMoreBlocks   2147483647

#define kMaxFilename 1024
#define kHashSize 65536

/* ------------------------------------------------------------------------------------ */
typedef	struct	 SEQDATA_struct	SEQDATA;
struct	SEQDATA_struct
{
	unsigned short data[kMaxInputCol];
};

typedef	struct	SEQNAME_struct	SEQNAME;
struct	SEQNAME_struct
{
	unsigned int condition;
};


typedef	struct	CLOSESTSEQ_struct	CLOSESTSEQ;
struct	CLOSESTSEQ_struct
{
	int id;
	unsigned int dist;
};

typedef	struct	REVISED_CLUSTER_ID_struct	REVISED_CLUSTER_ID;
struct	REVISED_CLUSTER_ID_struct
{
	unsigned int pos;
	unsigned int newid;
};

typedef	struct	CONFLICT_struct	CONFLICT;
struct	CONFLICT_struct
{
	unsigned int finalID;
	unsigned int cnt;
};


typedef	struct	CHECK_CLUSTER_ID_struct CHECK_CLUSTER_ID;
struct	CHECK_CLUSTER_ID_struct
{
	unsigned int finalID;
	unsigned int cnt;
	CONFLICT *conflictPtr; 
};

typedef	struct	CHUNK_struct	CHUNK;
struct	CHUNK_struct
{
	unsigned int ii;
	unsigned int jj;
	unsigned int iilast;
	unsigned int jjlast;
	unsigned int status;
};

typedef	struct	MERGECLUSTER_struct	MERGECLUSTER;
struct	MERGECLUSTER_struct
{
	 unsigned int cluster1;
	 unsigned int cluster2;
};

typedef	struct	JOINREQUEST_struct	JOINREQUEST;
struct	JOINREQUEST_struct
{
	unsigned int getFromCPU;
	unsigned int sendToCPU;
};

typedef	struct	CPU_struct	CPU;
struct	CPU_struct
{
	unsigned int ii;
	unsigned int jj;
	unsigned int iilast;
	unsigned int jjlast;
};

typedef	struct	STATS_struct	STATS;
struct	STATS_struct
{
	float dist;
	float pctAssigned;
	int rawClustersCnt;
	int trimmedClustersCnt;
};


typedef	struct	BRANCH_HEAD_struct	BRANCH_HEAD;
struct	BRANCH_HEAD_struct
{
	float dist;
	unsigned int step;
	int srcCluster;
};

typedef	struct	CLUSTERHISTORY_struct	CLUSTERHISTORY;
struct	CLUSTERHISTORY_struct
{
	unsigned int pass;
	float dist;
	int cluster;
	int descendfrom;
	int mergedto;
	char retain;
};

typedef	struct EXECUTIONPLAN_struct  EXECUTIONPLAN;
struct EXECUTIONPLAN_struct
{
	SEQDATA *fastadata;
	unsigned int *clusterid;
	CPU chunk;
};


typedef	struct	DISTLIST_struct	DISTLIST;
struct	DISTLIST_struct
{
	unsigned int squaredDist;
	float dist;
};

/* ------------------------------------------------------------------------------------ */

static unsigned int gTestDist;
#ifdef USE_THREADS
static volatile unsigned int clustercnt;
#else
static unsigned int clustercnt;
#endif
static unsigned int mergerequestcnt;
static MERGECLUSTER mergerequest[kMaxMergeRequests];

/* those constants have to stay like that... */
#define kClusterTooSmall 2
#define kClusterEliminated 1
#define kClusterLargeEnough 0


static int *clustersnum[kMaxCPU];
static 	CLUSTERHISTORY *clusterhistory;
static unsigned int clusterhistorycnt = 0;
static int printwarnmergereq = 1;

/* ------------------------------------------------------------------------------------ */

static void PrintUsage(char *version)
{
	printf("usage:\n\n");
	printf("dbc454 [-i InputFile] [-r ReferenceFile] -n numEventsToKeepCluster  [([-d DistanceCutoff] [-e EndingDistanceCutoff [-s Step]]) | -l SquaredDistanceTestList ] [-t TmpDirectory] [-o OutputFile] [-c ColumnOutputLevel] [ -v level][ -m]\n\n");
	printf("          -i InputFile              : Read Fasta sequences from file instead of stdin\n");
	printf("          -r ReferenceFile          : Reference Fasta sequences used to post-annotate the clusters (will not affect the clustering)\n");
	printf("          -n numEventsToKeepCluster : Keep only clusters with at least this number of events.\n");
	printf("          -t TmpDirectory           : Path of directory where to store tmp files; default is /tmp\n");
	printf("          -o OutputFile             : print results to OutputFile instead of stdout.\n");
	printf("          -c ColumnOutputLevel      : Define the extent of the output. See Output section.\n");
	printf("          -M                        : Report cluster Merging history\n");
	printf("          -v level                  : Specifies the verbose level; default is 0.\n");
	printf("  Euclidian distances to test can be provided as initial,ending and step using -d -e -s arguments\n");
	printf("          -d DistanceCutoff         : Initial Floating point Euclidian distance cutoff value used to place events in the same cluster. Defaults to 3.381\n");
	printf("          -e EndingDistanceCutoff   : Last Euclidian Distance cutoff to test. Defaults to 7.209\n");
	printf("          -s Step                   : Floating point increment of DistanceCutoff to test. Defaults to 0.174\n");
	printf("  Alternately distances to test can be provided as a list of integer values representing the square of the Euclidian distance\n");
	printf("          -l SquaredDistanceTestList: Comma separated list of integer values representing the square of the Euclidian distances to test\n");
	printf("\n%s\n",version);
	printf("Author: Nicolas Guex; 2013-2019\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, released under GPL2+ and you are welcome to redistribute it under certain conditions.\n");
	printf("----------------------------------------------------------------------------------------------------------------\n");
	printf("OUTPUT\n");
	printf("Progress report, warnings and errors are printed on stderr.\n");
	printf("The clustering result is printed on stdout, and consists of several tab-separated columns, with header:\n");
	printf("1.  sequence position in Fasta input file\n");
	printf("2.  final cluster id (zero for unclassified sequences)\n");
	printf("2b. cluster id at each squared distance cutoff level tested (cluster id is zero for unclassified sequences)\n");
	printf("    note that cluster ids are not stable from level to level (e.g. each level has its own cluster id numbering)\n");
	printf("    WARNING: although the structure of the results (e.g number of clusters and membership) will\n");
	printf("    be identical from run to run, the name (cluster id) of each cluster can vary from run to run.\n");
	printf("3.  sequence header (taken from the Fasta input file)\n");
	printf("4.  sequence (taken from the Fasta input file)\n");
	printf("\nBy default, a full output (-c 5) is produced, which can be quite large.\n");
	printf("This can be reduced with the ColumnOutputLevel option -c\n");
	printf("-c 1: output columns 1. and 2.\n");
	printf("-c 2: output columns 1. 2. and 3.\n");
	printf("-c 3: output columns 1. 2. 2b. and 3.\n");
	printf("-c 4: output columns 1. 2. 3. and 4.\n");
	printf("-c 5: output columns 1. 2. 2b. 3. and 4\n");
	printf("----------------------------------------------------------------------------------------------------------------\n");
	printf("EXAMPLES\n");
	printf("Takes fasta sequences from stdin, identify clusters using 8 processors, and output sequences in the same order as input, using default distance parameters\n");
	printf("  cat seq.fas | mpirun -np 8 dbc454 -n 100 -v 2 | sort --key=1 -n  > seq.clustered.txt\n");
	printf("Same as above, but explicitely supplying the list of squared distances to test corresponding to the default parameters\n");
	printf("  cat seq.fas | mpirun -np 8 dbc454 -n 100 -l '11,12,13,15,16,18,19,21,22,24,26,28,29,31,33,35,38,40,42,44,47,49,51' -v 2 | sort --key=1 -n  > seq.clustered.txt\n");
	printf("Takes fasta sequences from a file, identify clusters using default distance parameters and 8 processors, and output sequences by cluster\n");
	printf("  mpirun -np 8 dbc454 -i input.fas -n 100 | sort --key=2 -nr > seq.clustered.txt\n");
	printf("\nEXAMPLES for compute farms (check with your system administrator the specificities of your system)\n");
	printf("Takes fasta sequences from a file, submit job to a compute farm using 8 processors, and output results to a file\n");
	printf("monitor progress \"live\" thanks to -I flag\n");
	printf("bsub -I -n 8 -a mpichp4 -q dee -J jobname \"mpirun.lsf dbc454 -i seq.fas -o seq.clusters -n 100 -e 3.914 -c 4 -v 1\"\n");
	printf("Same as above, but send stde to a file\n");
	printf("bsub -e stdefile -n 8 -a mpichp4 -q dee -J jobname \"mpirun.lsf dbc454 -i seq.fas -o seq.clusters -n 100 -e 3.914 -c 4 -v 1\"\n");
	printf("----------------------------------------------------------------------------------------------------------------\n");

}
/* ------------------------------------------------------------------------------------ */

static unsigned int FillDistListFromRange(DISTLIST *df,float b,float e,float s)
{
	unsigned int i = 0;
	float dist = b;
	
	do {
		df[i].squaredDist = (int)((dist*dist));
		df[i].dist = dist;
		dist += s;
		if ((i==0) || (df[i].squaredDist > df[i-1].squaredDist))  // retain only increasing distance values.
			i++;
	} while(dist < e);
	return(i);
	
} /* FillDistListFromRange */
/* ------------------------------------------------------------------------------------ */

static unsigned int FillDistListFromList(DISTLIST *df,char *s,unsigned int cnt)
{
	unsigned int i = 0;
	int tot = 0;
	int l;
	unsigned int n;
	
	for(n=0;n<=cnt;n++)
	{
		sscanf(&s[tot],"%u,%n",&df[i].squaredDist,&l); tot+=l;
		df[i].dist = sqrt(df[i].squaredDist);
		if ((i==0) || (df[i].squaredDist > df[i-1].squaredDist))  // retain only increasing distance values.
			i++;
	};
	return(i);	

} /* FillDistListFromList */
/* ------------------------------------------------------------------------------------ */

static void PrintDistList(DISTLIST *df,unsigned int cnt)
{
	unsigned int i;
	
	fprintf(stderr,"LOG:The following distances will be used by dbc454\n");				
	fprintf(stderr,"LOG:squared_disances_tested=");				
	for(i=0;i<cnt;i++)
	{
		fprintf(stderr,",%6d",df[i].squaredDist);				
	};
	fprintf(stderr,"\n");				
	fprintf(stderr,"LOG: corresponding_disances=");				
	for(i=0;i<cnt;i++)
	{
		fprintf(stderr,",%6.3f",sqrt(df[i].squaredDist));				
	};
	fprintf(stderr,"\n");				
	fprintf(stderr,"LOG:      supplied_disances=");				
	for(i=0;i<cnt;i++)
	{
		fprintf(stderr,",%6.3f",df[i].dist);				
	};
	fprintf(stderr,"\n");				

} /* PrintDistList */
/* ------------------------------------------------------------------------------------ */

static unsigned int EncodeFasta(FILE *f,FILE *of,FILE *sf,SEQNAME *fastaname,SEQDATA *fastadata,unsigned short *key,unsigned  short *minkeyval,unsigned short *maxkeyval,unsigned int id,int nproc)
{
	SEQDATA *fastap;
	SEQNAME *fastanamep;
	int i;
	char linbuf[kMaxLineBuf];
	unsigned int rowcnt = 0;
	long long sum[kMaxInputCol];
	double score[kMaxInputCol];
	unsigned short cn;
	double  bestscore;
	unsigned int nextprintout = kLoadingProgressReporting;
	unsigned short minval = 65535;
	unsigned short maxval = 0;
    int seqlen = 0;
	char seq[16384];

	
	fastap = &fastadata[0];
	fastanamep = &fastaname[0];

	for (cn = 0; cn < kMaxInputCol; cn++)
		sum[cn] = 0;
		
	fgets(linbuf,kMaxLineBuf,f);
	if (linbuf[0] == '>')
	{
		fprintf(of,"%s",&linbuf[1]);
	}
	else
	{
		fprintf(stderr,"Error: Input File is not in Fasta format (does not start with '>')");
		return(0);
	}
	cn = 0;

	/* read sequences */
	do
	{
		fgets(linbuf,kMaxLineBuf,f);
		if (linbuf[0] == '>')
		{
			fprintf(of,"%s",&linbuf[1]);
		}
		if (feof(f) || (linbuf[0] == '>')) // save sequence.
		{
			int idx;
			seq[seqlen] = 0;
			if (sf)
				fprintf(sf,"%s\n",seq);


			memset(&fastap->data[0],0,kMaxInputCol*sizeof(unsigned short));
			/* adjacent dinucleotides */
			for (i = 0; i < (seqlen-1); i++)
			{
				switch(seq[i])
				{
					case 'A':	idx = 0; break;
					case 'T':	idx = 4; break;
					case 'C':	idx = 8; break;
					case 'G':	idx = 12; break;
					default:	idx = 9999; break;
				}
				switch(seq[i+1])
				{
					case 'A':	idx += 0; break;
					case 'T':	idx += 1; break;
					case 'C':	idx += 2; break;
					case 'G':	idx += 3; break;
					default:	idx = 9999; break;
				}
				if (idx < 9999)
				{
					fastap->data[idx]++;
					sum[idx]++;
				}
			}

			/*  dinucleotides with one spacer */
			for (i = 0; i < (seqlen-2); i++)
			{
				switch(seq[i])
				{
					case 'A':	idx = 16; break;
					case 'T':	idx = 20; break;
					case 'C':	idx = 24; break;
					case 'G':	idx = 28; break;
					default:	idx = 9999; break;
				}
				switch(seq[i+2])
				{
					case 'A':	idx += 0; break;
					case 'T':	idx += 1; break;
					case 'C':	idx += 2; break;
					case 'G':	idx += 3; break;
					default:	idx = 9999; break;
				}
				if (idx < 9999)
				{
					fastap->data[idx]++;
					sum[idx]++;
				}
			}

			/*  dinucleotides with two spacer */
			for (i = 0; i < (seqlen-3); i++)
			{
				switch(seq[i])
				{
					case 'A':	idx = 32; break;
					case 'T':	idx = 36; break;
					case 'C':	idx = 40; break;
					case 'G':	idx = 44; break;
					default:	idx = 9999; break;
				}
				switch(seq[i+3])
				{
					case 'A':	idx += 0; break;
					case 'T':	idx += 1; break;
					case 'C':	idx += 2; break;
					case 'G':	idx += 3; break;
					default:	idx = 9999; break;
				}
				if (idx < 9999)
				{
					fastap->data[idx]++;
					sum[idx]++;
				}
			}
			++rowcnt;
			fastanamep->condition = id++;
			fastanamep++;
			fastap++;
			if (rowcnt >= nextprintout)
			{
				fprintf(stderr,"LOG: %10u sequences loaded so far...\n",rowcnt);
				fflush(stderr);
				nextprintout += kLoadingProgressReporting;
			}

			seqlen = 0;
		}
		else	// accumulate nucleotides in current sequence.
		{
			i = 0;
			while (linbuf[i] != 0)
			{
				linbuf[i] = (char)toupper(linbuf[i]);
				switch(linbuf[i])
				{
					case 'A':
					case 'T':
					case 'C':
					case 'G':
						seq[seqlen++] = linbuf[i]; 
						break;
					case 'R':	// A or G
					case 'Y':	// C or T
					case 'S':	// G or C
					case 'W':	// A or T
					case 'K':	// G or T
					case 'M':	// A or C
					case 'B':	// C or G or T
					case 'D':	// A or G or T
					case 'H':	// A or C or T
					case 'V':	// A or C or G
					case 'N':	
						seq[seqlen++] = 'N';
						break;							
					case 'U':
						seq[seqlen++] = 'T';			// for now, treat U as T when we gather the dinucleotide statistics.
						break;
				}
				i++;
			}
		}
	} while (!feof(f));



	for (cn = 0; cn < kMaxInputCol;cn++)
	{
		/* compute SDDEV for a column */
		int mean = (int)(sum[cn] / rowcnt);
		long long sd = 0;
		fastap = &fastadata[0];
		for (i = 0; i<rowcnt;i++)
		{
			int diff = (int)fastap->data[cn] - mean;
			sd += (long long)(diff * diff); 
			fastap++;
		}
		score[cn] = (double)sd;
	}

	/* select column with largest sddev as sorting key */
	bestscore = 0.0;
	for (cn = 0; cn < kMaxInputCol;cn++)
	{
		if (score[cn] > bestscore)
		{
			bestscore = score[cn];
			*key = cn;
		}
	}

	fastap = &fastadata[0];
	cn = *key;
	for (i = 0; i<rowcnt;i++)
	{
		if (fastap->data[cn] > maxval)
			maxval = fastap->data[cn]; 
		if (fastap->data[cn] < minval)
			minval = fastap->data[cn]; 
		fastap++;
	}
	
	*minkeyval = minval;
	*maxkeyval = maxval;

	return(rowcnt);
	
} /* EncodeFasta */
/* ------------------------------------------------------------------------------------ */

static void mysort(unsigned short *valuestosort,int start, unsigned int last, unsigned int *inputindices, unsigned int *sortedindices,unsigned int minkeyval,unsigned int maxkeyval)
{
	int i;
	unsigned int cnt[65536];

		for (i = 0; i <= 65535; i++)
			cnt[i] = 0;
		for (i = start; i<=last;i++)
		{
			if ((valuestosort[inputindices[i]] < minkeyval) || (valuestosort[inputindices[i]] > maxkeyval))
				fprintf(stderr,"error\n");
			
			cnt[valuestosort[inputindices[i]]]++;
		}
		for (i = minkeyval+1; i <= maxkeyval; i++)
			cnt[i] += cnt[i-1];
		for (i = minkeyval; i <= maxkeyval; i++)
			cnt[i]--;

		for (i = last; i>=start;i--)
		{
			sortedindices[start+cnt[valuestosort[inputindices[i]]]] = inputindices[i];
			cnt[valuestosort[inputindices[i]]]--;
		}
} /* mysort */
/* ------------------------------------------------------------------------------------ */

static int SortEncodedFasta(SEQNAME *fastaname,SEQDATA	*fastadata,SEQNAME *sortedfastaname,SEQDATA	*sortedfastadata,unsigned int rcnt,unsigned short key,unsigned short minkeyval,unsigned short maxkeyval)
{
	unsigned int i;
	SEQDATA	*fastap;
	unsigned int cnt[65536];
	unsigned int *sortedcellnameidx=NULL;

		
		/* write data sorted according to key */
		sortedcellnameidx = malloc(rcnt*sizeof(unsigned int));
		if (sortedcellnameidx)		/* attempt to use fast sort that uses some memory... */
		{
			unsigned int start;
			unsigned int last;
			unsigned short curkey1;
			unsigned int *sortedcellnameidxcopy=NULL;
			unsigned short *hashval;

			for (i = minkeyval; i <= maxkeyval; i++)
				cnt[i] = 0;

			fastap = &fastadata[0];
			for (i = 0; i<rcnt;i++)
			{
				cnt[fastap->data[key]]++;
				fastap++;
			}

			for (i = minkeyval+1; i <= maxkeyval; i++)
				cnt[i] += cnt[i-1];
			for (i = minkeyval; i <= maxkeyval; i++)
				cnt[i]--;

			for (i = rcnt-1; i>=1;i--)
			{
				fastap--;
				sortedcellnameidx[cnt[fastap->data[key]]--] = i;
			}
			fastap--;
			sortedcellnameidx[cnt[fastap->data[key]]--] = 0;


			/* sort by second key */
			sortedcellnameidxcopy = malloc(rcnt*sizeof(unsigned int));
			if (sortedcellnameidxcopy)
			{
				/* make a copy of the sorted by first key, we will read from here and modify the original as we go */
				memcpy(sortedcellnameidxcopy,sortedcellnameidx,rcnt*sizeof(unsigned int));
				hashval = calloc(rcnt,sizeof(unsigned short));
				if (hashval)
				{
					/* compute one hash key to characterize each sequence */
					unsigned int bytecnt = 2*kMaxInputCol;
					for (i = 0; i< rcnt;i++)
					{
						unsigned int j;
						unsigned int hash = 0;
						char *datap = (char*)(&fastadata[i].data[0]);
						for (j=0; j<bytecnt; j++) 
						{
							hash = (hash ^ (*datap++)) * 16777619;
						}
						hash += hash << 13;
						hash ^= hash >> 7;
						hash += hash << 3;
						hash ^= hash >> 17;
						hash += hash << 5;
						hash = (hash>>16) ^ (hash & 0xffff);
						hashval[i] = (unsigned short)hash;
					}
					
					start = 0;
					do
					{
						last = start+1;
						if (last >= rcnt) {break;}
						curkey1 = fastadata[sortedcellnameidxcopy[start]].data[key];
						minkeyval = maxkeyval = hashval[sortedcellnameidxcopy[start]];
						while(fastadata[sortedcellnameidxcopy[last]].data[key] == curkey1)
						{ 
							if ( hashval[sortedcellnameidxcopy[last]] < minkeyval) minkeyval = hashval[sortedcellnameidxcopy[last]];
							if ( hashval[sortedcellnameidxcopy[last]] > maxkeyval) maxkeyval = hashval[sortedcellnameidxcopy[last]];
							last++; 
						}
						last--;
						// now we have position of entries with same value of key1 in [start..last]
						if (last > start) // we need to sort this chunk.
						{
							mysort(hashval,(int)start,last,sortedcellnameidxcopy,sortedcellnameidx,minkeyval,maxkeyval);
						}
						start=last+1;
					} while(start < rcnt);
					
					free(hashval);
				} // if hashval
				free(sortedcellnameidxcopy);
			}

			/* write results */
			for (i = 0; i<rcnt;i++)
			{
				sortedfastaname[i].condition = fastaname[sortedcellnameidx[i]].condition;
				fastap = &fastadata[sortedcellnameidx[i]];
				memcpy(&sortedfastadata[i].data[0],&fastap->data[0],sizeof(unsigned short)*kMaxInputCol);
			}
			free(sortedcellnameidx);
		}
		else	/* failover to slower sort without memory overhead */
		{
			unsigned int idx = 0;
			unsigned short val;
			for (val = minkeyval; val<= maxkeyval;val++)
			{
				fastap = &fastadata[0];
				for (i = 0; i<rcnt;i++)
				{
					if (fastap->data[key] == val)
					{
						sortedfastaname[idx].condition = fastaname[i].condition;
						memcpy(&sortedfastadata[idx].data[0],&fastap->data[0],sizeof(unsigned short)*kMaxInputCol);
						idx++;
					}
					fastap++;
				}
			}
		}
		return(0);
		
} /* SortEncodedFasta */
/* ------------------------------------------------------------------------------------ */

#ifdef USE_THREADS


static pthread_mutex_t clustercntmutex;
static unsigned short sortkey;
static unsigned int thread_mergerequestcnt;
static MERGECLUSTER thread_mergerequest[kMaxMergeRequests];
static void thread_InsertMergeRequest(unsigned int cluster1, unsigned int cluster2)  /* note: cluster2 > cluster1 is always true */
{
	unsigned int ii;
	
		if (thread_mergerequestcnt == 0)
		{
			thread_mergerequest[0].cluster1 = cluster1;
			thread_mergerequest[0].cluster2 = cluster2;
			thread_mergerequestcnt++;
			return;
		}
	
		if (thread_mergerequestcnt < kMaxMergeRequests)
		{
			/* quickly find position of cluster2 - 1 using a binary divide and conquer... */
			int left;
			unsigned int val;
			int right = thread_mergerequestcnt-1;

			left = 0;
			val=cluster2-1;
			do
			{
				ii=(left+right)>>1;
				if (thread_mergerequest[ii].cluster2 > val)
					right=ii-1;
				else
					left=ii+1;
			} while ((thread_mergerequest[ii].cluster2 != val) && (left <= right));
			
			/* move up from this point to do our real check, using cluster2... */
			/* test if and where to insert request */
			for( ;  ii< thread_mergerequestcnt; ii++)
			{
				if 	 (thread_mergerequest[ii].cluster2 > cluster2) 
					break;  /* we can insert and be done */
				if (thread_mergerequest[ii].cluster2 == cluster2) 
				{
					if (thread_mergerequest[ii].cluster1 == cluster1)
						return; /* already known */
					if (thread_mergerequest[ii].cluster1 > cluster1)
					{
						unsigned int tmp = thread_mergerequest[ii].cluster1;
						thread_mergerequest[ii].cluster1 = cluster1;
						thread_InsertMergeRequest(cluster1,tmp);
						return;
					}
					else
					{
						thread_InsertMergeRequest(thread_mergerequest[ii].cluster1,cluster1);
						return;
					}

				}
			}
			
			/* insert request */
			memmove(&thread_mergerequest[ii+1].cluster1, &thread_mergerequest[ii].cluster1, (thread_mergerequestcnt-ii)*sizeof(MERGECLUSTER));

			thread_mergerequest[ii].cluster1 = cluster1;
			thread_mergerequest[ii].cluster2 = cluster2;
			thread_mergerequestcnt++;
		}
		else
		{
			if (printwarnmergereq == 1)
			{
				fprintf(stderr,"LOG: ERROR: Too many merging requests (max = %d)\n",kMaxMergeRequests);
				printwarnmergereq = 0;
			}
		}


} /* thread_InsertMergeRequest */
/* ------------------------------------------------------------------------------------ */

static void *computesim_funcion(void *ep)
{
	SEQDATA *fastapi;
	SEQDATA *fastapj;
	unsigned int *clusterpi;
	unsigned int *clusterpj;
	unsigned int i,j;
	unsigned int starti,startj;
	unsigned int lasti,lastj;
	unsigned  int actualstartj;

	SEQDATA *fasta = ((EXECUTIONPLAN*)ep)->fastadata;
	unsigned int *clusterid = ((EXECUTIONPLAN*)ep)->clusterid;
	CPU *chunk = &((EXECUTIONPLAN*)ep)->chunk;

	starti = chunk->ii;
	lasti  = chunk->iilast;
	startj = chunk->jj;
	lastj  = chunk->jjlast;



		{
			unsigned short valii,valjj;

			fastapi = &fasta[lasti-1];
			fastapj = &fasta[startj];
			valii = fastapi->data[sortkey];
			valjj = fastapj->data[sortkey];
			if (valjj > valii)
			{
				unsigned  int mindist = valjj-valii;
				mindist *= mindist;
				if (mindist > gTestDist)
				{
					goto skipchunk;
				}
			}
		}


		fastapi = &fasta[starti];
		clusterpi = &clusterid[starti];
		for(i = starti;  i< lasti; i++)
		{
			if (startj <= i)
				actualstartj = i+1;
			else 
				actualstartj = startj;
			fastapj = &fasta[actualstartj];
			clusterpj = &clusterid[actualstartj];
			for(j = actualstartj;  j< lastj; j++)
			{
				int diff;
				unsigned int d;

				if ((*clusterpj) && (*clusterpj == *clusterpi))
					goto straight2next;


				diff = (int)fastapj->data[0] - fastapi->data[0];
				d = diff*diff;
				diff = (int)fastapj->data[1] - fastapi->data[1];
				d += diff*diff;
				diff = (int)fastapj->data[2] - fastapi->data[2];
				d += diff*diff;
				diff = (int)fastapj->data[3] - fastapi->data[3];
				d += diff*diff;
				diff = (int)fastapj->data[4] - fastapi->data[4];
				d += diff*diff;
				diff = (int)fastapj->data[5] - fastapi->data[5];
				d += diff*diff;
				diff = (int)fastapj->data[6] - fastapi->data[6];
				d += diff*diff;
				diff = (int)fastapj->data[7] - fastapi->data[7];
				d += diff*diff;
				diff = (int)fastapj->data[8] - fastapi->data[8];
				d += diff*diff;
				diff = (int)fastapj->data[9] - fastapi->data[9];
				d += diff*diff;
				diff = (int)fastapj->data[10] - fastapi->data[10];
				d += diff*diff;
				diff = (int)fastapj->data[11] - fastapi->data[11];
				d += diff*diff;
				diff = (int)fastapj->data[12] - fastapi->data[12];
				d += diff*diff;
				diff = (int)fastapj->data[13] - fastapi->data[13];
				d += diff*diff;
				diff = (int)fastapj->data[14] - fastapi->data[14];
				d += diff*diff;
				diff = (int)fastapj->data[15] - fastapi->data[15];
				d += diff*diff;

				/* computation done for the first 16 columns, time to attempt to gain some time */
				if (d > gTestDist)
					goto straight2next;

				#define kDATAOFFSET 16
				diff = (int)fastapj->data[kDATAOFFSET+0] - fastapi->data[kDATAOFFSET+0];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+1] - fastapi->data[kDATAOFFSET+1];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+2] - fastapi->data[kDATAOFFSET+2];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+3] - fastapi->data[kDATAOFFSET+3];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+4] - fastapi->data[kDATAOFFSET+4];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+5] - fastapi->data[kDATAOFFSET+5];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+6] - fastapi->data[kDATAOFFSET+6];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+7] - fastapi->data[kDATAOFFSET+7];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+8] - fastapi->data[kDATAOFFSET+8];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+9] - fastapi->data[kDATAOFFSET+9];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+10] - fastapi->data[kDATAOFFSET+10];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+11] - fastapi->data[kDATAOFFSET+11];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+12] - fastapi->data[kDATAOFFSET+12];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+13] - fastapi->data[kDATAOFFSET+13];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+14] - fastapi->data[kDATAOFFSET+14];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+15] - fastapi->data[kDATAOFFSET+15];
				d += diff*diff;

				/* computation done for the first 32 columns, time to attempt to gain some time */
				if (d > gTestDist)
					goto straight2next;

				#define kDATAOFFSET32 32
				diff = (int)fastapj->data[kDATAOFFSET32+0] - fastapi->data[kDATAOFFSET32+0];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+1] - fastapi->data[kDATAOFFSET32+1];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+2] - fastapi->data[kDATAOFFSET32+2];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+3] - fastapi->data[kDATAOFFSET32+3];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+4] - fastapi->data[kDATAOFFSET32+4];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+5] - fastapi->data[kDATAOFFSET32+5];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+6] - fastapi->data[kDATAOFFSET32+6];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+7] - fastapi->data[kDATAOFFSET32+7];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+8] - fastapi->data[kDATAOFFSET32+8];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+9] - fastapi->data[kDATAOFFSET32+9];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+10] - fastapi->data[kDATAOFFSET32+10];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+11] - fastapi->data[kDATAOFFSET32+11];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+12] - fastapi->data[kDATAOFFSET32+12];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+13] - fastapi->data[kDATAOFFSET32+13];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+14] - fastapi->data[kDATAOFFSET32+14];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+15] - fastapi->data[kDATAOFFSET32+15];
				d += diff*diff;

				if (d <= gTestDist) 
				{
					if (*clusterpi)
					{
						if (*clusterpj == 0)
						{
							*clusterpj = *clusterpi;
						}
						else /* i=assigned, j=assigned */
						{
							unsigned int cluster1,cluster2;
							if(*clusterpj > *clusterpi)
							{
								cluster1 = *clusterpi;
								cluster2 = *clusterpj;
							}
							else
							{
								cluster1 = *clusterpj;
								cluster2 = *clusterpi;
							}
							thread_InsertMergeRequest(cluster1,cluster2);
						}
						
					}
					else
					{
						if (*clusterpj == 0)   /* i=not yes assigned, j=not yet assigned */
						{
							pthread_mutex_lock(&clustercntmutex);
							*clusterpi = ++clustercnt;
							*clusterpj = clustercnt;
							pthread_mutex_unlock(&clustercntmutex);
						}
						else /* i=not yes assigned, j=assigned */
						{
							*clusterpi = *clusterpj;
						}
					}
				}
		straight2next:		
				fastapj++;
				clusterpj++;
			}
			fastapi++;
			clusterpi++;
		}
skipchunk:
	pthread_exit(NULL);

} /* computesim_funcion */
#endif
/* ------------------------------------------------------------------------------------ */

static void DoMergeLists(int nproc,int verbose)
{
		JOINREQUEST joinRequest;
		int cpu,sendingCPUoffset;
		unsigned int received,joiningprocesses;
		
		sendingCPUoffset = 1;
		do
		{
			cpu = 1;

			if ((cpu + sendingCPUoffset) >= nproc)
				break;

			/* send instruction to join list */
			joiningprocesses = 0;
			do
			{
				if ((cpu + sendingCPUoffset) >= nproc)
					break;
				joinRequest.getFromCPU = cpu + sendingCPUoffset;
				joinRequest.sendToCPU = cpu;				
				MPI_Send(&joinRequest, 2, MPI_INT,  joinRequest.sendToCPU,  kJoinListRequest, MPI_COMM_WORLD);
				MPI_Send(&joinRequest, 2, MPI_INT,  joinRequest.getFromCPU, kJoinListRequest, MPI_COMM_WORLD);
				cpu += sendingCPUoffset+sendingCPUoffset;
				joiningprocesses++;
			} while (1);

			/* wait until done */
			received = 0;
			cpu = 1;
			do
			{
				unsigned int done;
				if ((cpu + sendingCPUoffset) >= nproc)
					break;
				MPI_Recv(&done, 1, MPI_INT,  cpu, kJoinListRequestDone, MPI_COMM_WORLD, MPI_STATUS_IGNORE/*&status*/);
				cpu += sendingCPUoffset+sendingCPUoffset;
			} while (++received < joiningprocesses);

			sendingCPUoffset *= 2;

		} while(1);

		joinRequest.getFromCPU  = 0;
		joinRequest.sendToCPU = 0;
		for (cpu = 1; cpu < nproc; cpu++)
			MPI_Send(&joinRequest, 2, MPI_INT,  cpu,  kJoinListRequest, MPI_COMM_WORLD);

} /* DoMergeLists */
/* ------------------------------------------------------------------------------------ */

static unsigned int DiscardDuplicates(SEQNAME *sortedfastaname, SEQDATA *sortedfastadata,SEQNAME *retainedfastaname, SEQDATA *retainedfasta,unsigned int *clusterid,unsigned int *sameFollowingCnt, unsigned int rowcnt)
{
	SEQDATA *fastap;
	SEQNAME *fastanamep;
	unsigned int cnt;

	SEQDATA mock;
	SEQDATA *prevfastap;
	unsigned int lastuniqueentry = 0;
	unsigned int retained = 0;
	unsigned int clustercnt = (1*kStartLocalCluster); // do as if first computing node did attribute clusters.
	unsigned int duplicateCnt = 0;



	for (cnt =0; cnt < kMaxInputCol; cnt++)
		mock.data[cnt] = 0xffff;
	prevfastap = &mock;

	cnt = 0;
	fastap = &retainedfasta[0];	
	fastanamep = &retainedfastaname[0];
	for (cnt = 0; cnt < rowcnt; cnt++)
	{
		// "read" entry from memory sorted ones. 
		fastanamep->condition = sortedfastaname[cnt].condition;
		memcpy(&fastap->data[0],&sortedfastadata[cnt].data[0],sizeof(short)*kMaxInputCol);
		// check if should keep it.
		if (memcmp(prevfastap,fastap,sizeof(unsigned short)*kMaxInputCol) != 0) // differ from previous entry
		{
			prevfastap = fastap;
			fastap++;
			lastuniqueentry = retained;
			retained++;
		}			
		else
		{
			sameFollowingCnt[lastuniqueentry]++;
			if (sameFollowingCnt[lastuniqueentry] == 1)
				clusterid[lastuniqueentry] = ++clustercnt;
			duplicateCnt++;
		}
		fastanamep++;
	}
	return(rowcnt-duplicateCnt);
	
} /* DiscardDuplicates */
/* ------------------------------------------------------------------------------------ */

static unsigned int DiscardRefSeqDuplicates(SEQNAME *sortedfastaname, SEQDATA *sortedfastadata,SEQNAME *retainedfastaname, SEQDATA *retainedfasta,unsigned int *sameFollowingCnt, unsigned int rowcnt)
{
	SEQDATA *fastap;
	SEQNAME *fastanamep;
	unsigned int cnt;
	unsigned int lastuniqueentry = 0;

	SEQDATA mock;
	SEQDATA *prevfastap;
	unsigned int retained = 0;
	unsigned int duplicateCnt = 0;



	for (cnt =0; cnt < kMaxInputCol; cnt++)
		mock.data[cnt] = 0xffff;
	prevfastap = &mock;

	cnt = 0;
	fastap = &retainedfasta[0];	
	fastanamep = &retainedfastaname[0];
	for (cnt = 0; cnt < rowcnt; cnt++)
	{
		// "read" entry from memory sorted ones. 
		fastanamep->condition = sortedfastaname[cnt].condition;
		memcpy(&fastap->data[0],&sortedfastadata[cnt].data[0],sizeof(short)*kMaxInputCol);
		// check if should keep it.
		if (memcmp(prevfastap,fastap,sizeof(unsigned short)*kMaxInputCol) != 0) // differ from previous entry
		{
			prevfastap = fastap;
			fastap++;
			lastuniqueentry = retained;
			retained++;
		}			
		else
		{
			sameFollowingCnt[lastuniqueentry]++;
			duplicateCnt++;
		}
		fastanamep++;
	}
	return(rowcnt-duplicateCnt);
	
} /* DiscardRefSeqDuplicates */
/* ------------------------------------------------------------------------------------ */

static unsigned int CountAssigned(unsigned int *clusterid,unsigned int *repeatcnt, unsigned int rowcnt,unsigned int maxclusterid)
{
	unsigned int *clusteridp;
	unsigned int i;
	unsigned int assigned = 0;

			clusteridp = &clusterid[0];
			for(i = 0;  i< rowcnt; i++)
			{
					if ((*clusteridp > 0) && (*clusteridp<=maxclusterid)) /* assigned */
						assigned+=(repeatcnt[i]+1);
					clusteridp++;
			}//i

		return(assigned);
} /* CountAssigned */
/* ------------------------------------------------------------------------------------ */

static void WriteClusterIndices(unsigned int *clusterid,unsigned int rowcnt,float distcutoff,char *dir)
{
	unsigned int *clusteridp;
	FILE *af;
	char fn[kMaxFilename];

		sprintf(fn,"%s-%.6f",dir,distcutoff);
		af=fopen(fn,"wb");
		if (af)
		{
			clusteridp = &clusterid[0];	
			fwrite(clusteridp,sizeof(int),rowcnt,af);
			fclose(af);
		}
		else
			fprintf(stderr,"LOG: Error writing %s\n",fn);
	
} /* WriteClusterIndices */
/* ------------------------------------------------------------------------------------ */

static int LoadClusterIndices(unsigned int *clusterid,unsigned int step,unsigned int nbSeqForClustering,unsigned int rowcnt,float distcutoff,char *dir,unsigned int *repeatcnt,unsigned int minevents)
{
	unsigned int *clusteridp;
	FILE *af;
	char fn[kMaxFilename];
	int err = 1;
	unsigned int i;

		unsigned int *tmp = calloc(rowcnt+1,sizeof(unsigned int));
		if (!tmp)
			fprintf(stderr,"LOG: Error allocating memory for reading %s\n",fn);
		else
		{
			sprintf(fn,"%s-%.6f",dir,distcutoff);
			af=fopen(fn,"rb");
			if (af)
			{
				// transfer cluster id onto appropirate chunk of big rslt table.
				clusteridp = &clusterid[step*rowcnt];
				fread(clusteridp,sizeof(int),nbSeqForClustering,af);
				fclose(af);
				// cnt nb of seq by clusters.
				for (i = 0; i < nbSeqForClustering; i++)
				{
					tmp[clusterid[step*rowcnt+i]]+=repeatcnt[i]+1;
				}
				// flag clusters too small.				
				for (i = 1; i< nbSeqForClustering; i++)
				{
					if (tmp[i] < minevents)
					{
						tmp[i] = kClusterTooSmall;
					}
				}	
				for (i = 0; i < nbSeqForClustering; i++)
				{
					if (tmp[clusterid[step*rowcnt+i]] == kClusterTooSmall)
						clusterid[step*rowcnt+i] = 0;
				}
				err = 0;
			}
			else
				fprintf(stderr,"LOG: Error reading %s\n",fn);
			free(tmp);
		}
		return(err);
	
} /* LoadClusterIndices */
/* ------------------------------------------------------------------------------------ */

static unsigned int WriteOutput(FILE *f,FILE *seqf,SEQNAME *fastaname,unsigned int *clusterid,unsigned int *repeatcnt,unsigned int rowcnt,unsigned int loaded,unsigned int nbSeqForClustering,unsigned int maxclusterid,FILE *of,unsigned int OutputLevel,unsigned int *refseqassigned)
{
	unsigned int *clusteridp;
	unsigned int i;
	SEQNAME *fastanamep;
	char linbuf[kMaxLineBuf];
	char seqbuf[16384];
	fpos_t *tmp = NULL;
	fpos_t *seqp = NULL;
	unsigned int assigned = 0;


		if (OutputLevel != 1)
			fprintf(of,"seqnum\tcluster\theader");
		else
			fprintf(of,"seqnum\tcluster");
		if (seqf) // we want sequences
			fprintf(of,"\tsequence");
		fprintf(of,"\n");

		tmp = calloc(rowcnt,sizeof(fpos_t));
		if (!tmp)
			goto bail;

		for(i = 0;  i< rowcnt; i++)
		{
			fgetpos(f,&tmp[i]);
			fgets(linbuf,kMaxLineBuf,f);
		}

		if (seqf) // we want sequences
		{
			seqp = calloc(rowcnt,sizeof(fpos_t));
			if (!seqp)
				goto bail;
			for(i = 0;  i< rowcnt; i++)
			{
				fgetpos(seqf,&seqp[i]);
				fgets(seqbuf,kMaxLineBuf,seqf);
			}
		}


		clusteridp = &clusterid[0];
		fastanamep = &fastaname[0];
		// print non reference sequences
		for(i = 0;  i< nbSeqForClustering; i++)
		{
			unsigned int cnt = 0;
			if ((*clusteridp > 0) && (*clusteridp<=maxclusterid)) /* assigned */
				assigned+=(repeatcnt[i]+1);
			do 
			{
				fsetpos(f,&tmp[fastanamep->condition-1]);
				fgets(linbuf,kMaxLineBuf,f);
				if (*clusteridp<=maxclusterid)
					fprintf(of,"%u\t%u",fastanamep->condition,*clusteridp);
				else
					fprintf(of,"%u\t0",fastanamep->condition);
				if (seqf) // we want sequences
				{
					fsetpos(seqf,&seqp[fastanamep->condition-1]);
					fgets(seqbuf,16384,seqf);
					linbuf[strlen(linbuf)-1]='\t';
					fprintf(of,"\t>%s%s",linbuf,seqbuf);
				}
				else
				{
					if (OutputLevel != 1)
						fprintf(of,"\t>%s",linbuf);
					else
						fputc('\n',of);
				}
				fastanamep++;
			} while(cnt++ < repeatcnt[i]);
			clusteridp++;
		}//i

		// print reference sequences
		for(i = nbSeqForClustering;  i< loaded; i++)
		{
			unsigned int cnt = 0;
			if ((*clusteridp > 0) && (*clusteridp<=maxclusterid)) /* assigned */
				*refseqassigned+=(repeatcnt[i]+1);
			do 
			{
				fsetpos(f,&tmp[fastanamep->condition-1]);
				fgets(linbuf,kMaxLineBuf,f);
				if (*clusteridp<=maxclusterid)
					fprintf(of,"%u\t%u",fastanamep->condition,*clusteridp);
				else
					fprintf(of,"%u\t0",fastanamep->condition);
				if (seqf) // we want sequences
				{
					fsetpos(seqf,&seqp[fastanamep->condition-1]);
					fgets(seqbuf,16384,seqf);
					linbuf[strlen(linbuf)-1]='\t';
					fprintf(of,"\t>%s%s",linbuf,seqbuf);
				}
				else
				{
					if (OutputLevel != 1)
						fprintf(of,"\t>%s",linbuf);
				}
				fastanamep++;
			} while(cnt++ < repeatcnt[i]);
			clusteridp++;
		}//i




bail:
		if (tmp)
			free(tmp);
		if (seqp)
			free(seqp);

		return(assigned);

} /* WriteOutput */
/* ------------------------------------------------------------------------------------ */

static unsigned int WriteFullOutput(FILE *f,FILE *seqf,SEQNAME *fastaname,unsigned int *clusterid,unsigned int *repeatcnt,CLOSESTSEQ *closestSeq,unsigned int rowcnt,unsigned int loaded,unsigned int nbSeqForClustering,unsigned int maxclusterid,FILE *of,unsigned int *rslt, unsigned int nsteps,DISTLIST *distList,unsigned int *refseqassigned)
{
	unsigned int *clusteridp;
	unsigned int i;
	SEQNAME *fastanamep;
	char linbuf[kMaxLineBuf];
	char seqbuf[16384];
	fpos_t *tmp = NULL;
	fpos_t *seqp = NULL;
	int step;
	unsigned int assigned = 0;

		fprintf(of,"seqnum\tcluster");
		for(step=(nsteps-1); step>=0;step--)
			fprintf(of,"\t%u",distList[step].squaredDist);											
		fprintf(of,"\theader");
		if (seqf) // we want sequences
			fprintf(of,"\tsequence");
		fprintf(of,"\n");

		tmp = calloc(rowcnt,sizeof(fpos_t));
		if (!tmp)
			goto bail;

		for(i = 0;  i< rowcnt; i++)
		{
			fgetpos(f,&tmp[i]);
			fgets(linbuf,kMaxLineBuf,f);
		}

		if (seqf) // we want sequences
		{
			seqp = calloc(rowcnt,sizeof(fpos_t));
			if (!seqp)
				goto bail;
			for(i = 0;  i< rowcnt; i++)
			{
				fgetpos(seqf,&seqp[i]);
				fgets(seqbuf,kMaxLineBuf,seqf);
			}
		}

		clusteridp = &clusterid[0];
		fastanamep = &fastaname[0];
		// print non reference sequences
		for(i = 0;  i< nbSeqForClustering; i++)
		{
			unsigned int cnt = 0;
			if ((*clusteridp > 0) && (*clusteridp<=maxclusterid)) /* assigned */
				assigned+=(repeatcnt[i]+1);
			do 
			{
				fsetpos(f,&tmp[fastanamep->condition-1]);
				fgets(linbuf,kMaxLineBuf,f);
				fprintf(of,"%u\t%u",fastanamep->condition,*clusteridp);
				for(step=(nsteps-1); step>=0;step--)
					fprintf(of,"\t%u",rslt[step*loaded+i]);
				if (seqf) // we want sequences
				{
					fsetpos(seqf,&seqp[fastanamep->condition-1]);
					fgets(seqbuf,16384,seqf);
					linbuf[strlen(linbuf)-1]='\t';
					fprintf(of,"\t>%s%s",linbuf,seqbuf);
				}
				else
				{
					fprintf(of,"\t>%s",linbuf);
				}
				fastanamep++;
			} while(cnt++ < repeatcnt[i]);
			clusteridp++;
		}//i

		// print reference sequences
		for(i = nbSeqForClustering;  i< loaded; i++)
		{
			unsigned int cnt = 0;
			if ((*clusteridp > 0) && (*clusteridp<=maxclusterid)) /* assigned */
				*refseqassigned+=(repeatcnt[i]+1);
			do 
			{
				fsetpos(f,&tmp[fastanamep->condition-1]);
				fgets(linbuf,kMaxLineBuf,f);
				fprintf(of,"%u\t%u",fastanamep->condition,*clusteridp);
				if (closestSeq && closestSeq[i-nbSeqForClustering].id>=0)
				{
					for(step=(nsteps-1); step>=0;step--)
					{
						if (closestSeq[i-nbSeqForClustering].dist <= distList[step].squaredDist)
							fprintf(of,"\t%u",rslt[step*loaded+closestSeq[i-nbSeqForClustering].id]);
						else
							fprintf(of,"\t%d",0);
					}
				}
				else
				{
					for(step=(nsteps-1); step>=0;step--)
						fprintf(of,"\t%d",0);
				}
				if (seqf) // we want sequences
				{
					fsetpos(seqf,&seqp[fastanamep->condition-1]);
					fgets(seqbuf,16384,seqf);
					linbuf[strlen(linbuf)-1]='\t';
					fprintf(of,"\t>%s%s",linbuf,seqbuf);
				}
				else
				{
					fprintf(of,"\t>%s",linbuf);
				}
				fastanamep++;
			} while(cnt++ < repeatcnt[i]);
			clusteridp++;
		}//i
bail:
		if (tmp) free(tmp);
		if (seqp) free(seqp);

		fflush(stdout);
		return(assigned);

} /* WriteFullOutput */
/* ------------------------------------------------------------------------------------ */

static void removeclustersnum(int id)
{
	int *cnp;
	int cpu =   id / kStartLocalCluster ;
	int cluster = id - cpu*kStartLocalCluster;
	cnp = clustersnum[cpu];
	cnp[cluster] = kClusterEliminated;  /* flag this specific cluster id as being removed */

} /* removeclustersnum */
/* ------------------------------------------------------------------------------------ */

static void AdjustClustersID(unsigned int *clusters, unsigned int loaded,unsigned int nproc,int trimmedclustercnt,int *firstAvailClusterID)
{
	unsigned int i;
	int j;
	int *cnp;
	unsigned int clusterid = 1; /* first cluster will have id=1 */
	unsigned int tinyClustersId = trimmedclustercnt +1 ;  /* start to pile up number of clusters too small to pass the min size cutoff after "good" clusters */
	unsigned int ii;



	/* loop over each cpu, which has attributed its own ids */
	for (i = 1; i < nproc; i++)
	{
		cnp = clustersnum[i];
		for (j = 1; j<=cnp[0]; j++)
		{
			if (cnp[j]==kClusterLargeEnough)
			{
				/* rename clusterid (attribute clusterid as a new id) */
				cnp[j]=clusterid++;
			}
			else if (cnp[j]==kClusterTooSmall)
			{
				cnp[j]=tinyClustersId++;
			}
			else
				cnp[j] = 0;
		}

		for(ii = 0;  ii< loaded; ii++)
		{
			if ((int)clusters[ii] >= i*kStartLocalCluster)
			{
				int idx = (int)clusters[ii] - i*kStartLocalCluster;
				if (idx < kStartLocalCluster)
				{
						clusters[ii] = cnp[idx];
				}
			}
		}
		
		
	}
	
	*firstAvailClusterID = tinyClustersId-1;
	
} /* AdjustClustersID */
/* ------------------------------------------------------------------------------------ */

static unsigned int RemoveSmallClusters(unsigned int *clusters,unsigned int *repeatcnt,unsigned int loaded,unsigned int nproc,unsigned int minevents)
{
	unsigned int ii;
	int i,j;
	int *cnp;
	unsigned int retainedClusterCnt = 0;

	/* loop over each cpu, which has attributed its own ids */
	for (i = 1; i < nproc; i++)
	{
		unsigned int *cnt;
		cnp = clustersnum[i];
		cnt = calloc((cnp[0]+1),sizeof(unsigned int));
		for(ii = 0;  ii< loaded; ii++)
		{
			int idx = (int)clusters[ii] - i*kStartLocalCluster;
			if ((idx >= 0) && (idx <= cnp[0])) /* cluster belongs to proc i */
			{
				cnt[idx] += repeatcnt[ii]+1; // +1 for actual self, as table is init to zero.
			}
		}
		for (j = 1; j<=cnp[0]; j++)
		{
			if (cnt[j] < minevents)
			{
				if (cnp[j] != kClusterEliminated)
				{
					cnp[j]=kClusterTooSmall;  /* flag this specific cluster id as being removed */
				}
			}
			else
			{
				retainedClusterCnt++;
			}
		}
		free(cnt);
	}
	return(retainedClusterCnt);
	
} /* RemoveSmallClusters */
/* ------------------------------------------------------------------------------------ */

static void PrintClusterStatus(void)
{
	unsigned int ii;
	for (ii = 0; ii<clusterhistorycnt; ii++)
		fprintf(stderr,"LOG:MERGE_HISTORY: %4u %7.3f | %3d ^ %3d -> %3d | %c\n",clusterhistory[ii].pass,clusterhistory[ii].dist,clusterhistory[ii].cluster,clusterhistory[ii].descendfrom,clusterhistory[ii].mergedto,clusterhistory[ii].retain);
		
				
} /* PrintClusterStatus */
/* ------------------------------------------------------------------------------------ */

static void DistributeRefSeqToClosestCluster(SEQDATA *fasta, unsigned int *clusterid,unsigned int usedForClustering,unsigned int starti,unsigned int lasti,CLOSESTSEQ *closestSeq)
{
	unsigned int i,j,col;
	unsigned int cnt = 0;

	for (i = starti; i < lasti; i++)
	{
			SEQDATA *fastapi = &fasta[i];
			unsigned int dmin = 0xffffffff;
			unsigned int jmin = i;
			// compute distance of ref. sequences (fastapi) to each assigned (fastapj) and record closest sequence.
			for (j = 0; j < usedForClustering; j++)
			{
				if (clusterid[j] > 0) /* assigned */
				{
					unsigned int d = 0;
					SEQDATA *fastapj = &fasta[j];
					for (col=0;col < kMaxInputCol;col++)
					{
						int diff = (int)fastapj->data[col] - fastapi->data[col];
						d += diff*diff;
					}
					if (d < dmin)
					{
						dmin = d;
						jmin = j;
						if (d == 0)
						{
							break;
						}
					}
				}
			} // for j
			if (dmin <= gTestDist)
			{
				clusterid[i] = clusterid[jmin]; // ok to attribute clusterid directly, as ref. sequences are not checked against other ref. seq.
				closestSeq[i-starti].id = jmin;
				closestSeq[i-starti].dist = dmin;
				cnt++;
			}
	}

} /* DistributeRefSeqToClosestCluster */
/* ------------------------------------------------------------------------------------ */

static unsigned int DistributeUnassignedToClosestCluster(SEQDATA *fasta, unsigned int *clusterid,unsigned int usedForClustering,unsigned int maxclusterid,unsigned int starti,unsigned int lasti)
{
	unsigned int i,j,col;
	unsigned int dmin;
	unsigned int jmin;
	unsigned int ambiguousFlag;
	unsigned int ambiguousCnt = 0;
	unsigned int reassigned = 0;

	for (i = starti; i < lasti; i++)
	{
		if (clusterid[i] == 9999999) /* flagged for reassignment */
		{
			SEQDATA *fastapi = &fasta[i];
			dmin = 0xffffffff;
			jmin = i; // just to initialize something.
			ambiguousFlag = 0;

			// compute distance of unassigned (fastapi) to each assigned (fastapj) and record closest event of each cluster.
			for (j = 0; j < usedForClustering; j++)
			{
				if ((clusterid[j] > 0) && (clusterid[j] <= maxclusterid)) /* assigned */
				{
					unsigned int d = 0;
					SEQDATA *fastapj = &fasta[j];
					for (col=0;col < kMaxInputCol;col++)
					{
						int diff = (int)fastapj->data[col] - fastapi->data[col];
						d += diff*diff;
					}
					// record closest distance
					if (d < dmin)
					{
						dmin = d;
						jmin = j;
						ambiguousFlag = 0;
					}
					else if ((d == dmin) && (clusterid[j] != clusterid[jmin]))
					{
						ambiguousFlag = 1;
					}
				}
			}
			// set the id of the closest cluster (but with an offset to make sure it will not affect loop above - cluster > maxclusterid for now)
			if (ambiguousFlag)
			{
				clusterid[i] = 0 + maxclusterid + 1;
				ambiguousCnt++;
			}
			else
			{
				clusterid[i] = clusterid[jmin] + maxclusterid + 1;
				reassigned++;
			}
		}
	}
	// adjust clusterid of reassigned sequences.
	for (i = starti; i < lasti; i++)
	{
		if (clusterid[i] > maxclusterid)
			clusterid[i] -= (maxclusterid + 1);
	}
	return(ambiguousCnt);
	
} // DistributeUnassignedToClosestCluster
/* ------------------------------------------------------------------------------------ */

static void UpdateClusterHistory(unsigned int pass, int clustercnt,int prevclustercnt, float dist)
{
	unsigned int ii;
	unsigned int prevcnt = clusterhistorycnt;
	
		if (pass == 0)
		{
			for (ii = 1; ii<=clustercnt; ii++)
			{
				clusterhistory[clusterhistorycnt].pass = pass;
				clusterhistory[clusterhistorycnt].dist = dist;
				clusterhistory[clusterhistorycnt].cluster = ii;
				clusterhistory[clusterhistorycnt].descendfrom = 0;
				clusterhistory[clusterhistorycnt].mergedto = 0;
				clusterhistory[clusterhistorycnt].retain = ' ';
				clusterhistorycnt++;
			}
		}
		else
		{			
			int parent = 1;
			for (ii = 1; ii<=clustercnt; ii++)
			{
				unsigned int x;
				
				if (parent <= prevclustercnt)
				{
					for (x = 0; x<prevcnt; x++)
					{
						if ((clusterhistory[x].pass == (pass-1)) && (clusterhistory[x].cluster == parent))
						{
							if (clusterhistory[x].mergedto != 0)
								parent++;
						}
					}
					if (parent <= prevclustercnt)
						clusterhistory[clusterhistorycnt].descendfrom = parent;
					else
						clusterhistory[clusterhistorycnt].descendfrom = 0;
					parent++;
				}
				else
					clusterhistory[clusterhistorycnt].descendfrom = 0;
				clusterhistory[clusterhistorycnt].pass = pass;
				clusterhistory[clusterhistorycnt].dist = dist;
				clusterhistory[clusterhistorycnt].cluster = ii;
				clusterhistory[clusterhistorycnt].mergedto = 0;
				clusterhistory[clusterhistorycnt].retain = ' ';
				clusterhistorycnt++;
			}
		}
		
		
} /* UpdateClusterHistory */
/* ------------------------------------------------------------------------------------ */
static char CheckClusterNotYetRetained(int from,int cluster)
{
	int x;
	for (x = from-1; x >= 0; x--)
	{
		if (clusterhistory[x].cluster == cluster)
		{
			if (clusterhistory[x].retain != ' ')
				return('n');
			if (clusterhistory[x].descendfrom == 0)
				return('y');
			return(CheckClusterNotYetRetained(x,clusterhistory[x].descendfrom));
		}
	}
//	return('!');	// terminal node.
	return('y');	// return 'y' instead of '!'. this denotes a cluster identified at the last step.
					// cannot be merged to anything, but must retain. Otherwise, its sequences are
					// arbitrarily redistributed to an other cluster, which is obviously very distant...
	
} /* CheckClusterNotYetRetained */
/* ------------------------------------------------------------------------------------ */

static int SelectClusterHistory(unsigned int rowcnt, unsigned int *clusterid,char *dir)
{
	unsigned int ii;
	int x;
	int cnt;
	
	for (ii = 0; ii<clusterhistorycnt; ii++)
	if (clusterhistory[ii].mergedto)
	{
		if (clusterhistory[ii].descendfrom == 0)
			clusterhistory[ii].retain = 'y';
		else
			clusterhistory[ii].retain = CheckClusterNotYetRetained(ii,clusterhistory[ii].descendfrom);
			
		for (x = 0; x<ii; x++)
		{
			if ((clusterhistory[x].pass == (clusterhistory[ii].pass )) && (clusterhistory[x].cluster == clusterhistory[ii].mergedto))
			{
				if (clusterhistory[x].descendfrom == 0)
					clusterhistory[x].retain = 'y';
				else
					clusterhistory[x].retain = CheckClusterNotYetRetained(x,clusterhistory[x].descendfrom);
				break;
			}
		}
		
	}


	for (x = clusterhistorycnt-1; x >= 0; x--)
	{
		if (clusterhistory[x].pass != clusterhistory[clusterhistorycnt-1].pass)
			break;
		if (clusterhistory[x].retain != ' ')
			continue;
		clusterhistory[x].retain = CheckClusterNotYetRetained(x,clusterhistory[x].descendfrom);
	}
	
	

	cnt = 1;  // start new cluster names at 1

	for (ii = 0; ii<clusterhistorycnt; ii++)
	{
		if (clusterhistory[ii].retain == 'y')
		{
			char fn[kMaxFilename];
			FILE *f = NULL;

			sprintf(fn,"%s-%.6f",dir,clusterhistory[ii].dist);
			f=fopen(fn,"rb");
			if (f)
			{
				unsigned int i;
				for(i = 0;  i< rowcnt; i++)
				{
					unsigned int cl;
					fread(&cl,sizeof(int),1,f);
					if (cl == clusterhistory[ii].cluster)
					{
						clusterid[i] = cnt;
					}
				} 
				fclose(f);
			}
			else
				fprintf(stderr,"LOG: Error reading %s\n",fn);

			cnt++;
		}
	}
	
	return(cnt-1);
		
} /* SelectClusterHistory */
/* ------------------------------------------------------------------------------------ */

static 	unsigned int CheckRsltArray(unsigned int *clusterid,unsigned int *rslt,unsigned int loaded,unsigned int maxClusterId)
{
	unsigned int i,j,cl;
	unsigned int fixedCnt = 0;
	CHECK_CLUSTER_ID *check = calloc((maxClusterId+1),sizeof(CHECK_CLUSTER_ID));
	if (check == NULL)
	{
		fprintf(stderr,"Warning: not enough memory to check cluster filiation.\n");
		goto bail;
	}
	
	for(i = 0;  i< loaded; i++)
	{
		if (clusterid[i] == 0)
			continue;
		if (check[clusterid[i]].finalID == 0) // not yet assigned
		{
			check[clusterid[i]].finalID = rslt[i];  // indicate the reported cluster.
			check[clusterid[i]].cnt++;
		}
		else
		{
			if (check[clusterid[i]].finalID == rslt[i]) // consistent.
				check[clusterid[i]].cnt++;
			else // conflict.
			{
				if (check[clusterid[i]].conflictPtr == NULL)
					check[clusterid[i]].conflictPtr = calloc(kMaxConflict,sizeof(CONFLICT));
				if (check[clusterid[i]].conflictPtr==NULL)
					fprintf(stderr,"Warning: not enough memory to fully check cluster filiation.\n");
				else
				{
					for(j = 0;  j< kMaxConflict; j++)
					{
						if (check[clusterid[i]].conflictPtr[j].finalID == rslt[i])  // already recorded, just increase cnt.
						{
							check[clusterid[i]].conflictPtr[j].cnt++;
							break;
						}
						if (check[clusterid[i]].conflictPtr[j].finalID == 0)  // empty slot.
						{
							check[clusterid[i]].conflictPtr[j].finalID = rslt[i];
							check[clusterid[i]].conflictPtr[j].cnt++;
							break;
						}
					}
					if (j >= kMaxConflict) // first time this conflict is encountered.
						fprintf(stderr,"Warning: kMaxConflict reached.\n");
				} // else valid check[clusterid[i]].conflictPtr.
			} // else conflict.
		}
		
	}
bail:

	for(cl = 1;  cl<= maxClusterId; cl++)
	{
		if (check[cl].conflictPtr)
		{
			unsigned int maxcnt = 0;
			unsigned int finalID = 0;
			for(j = 0;  j< kMaxConflict; j++) // identify highest frequency finalID.
			{
				if (check[cl].conflictPtr[j].cnt > maxcnt)
				{
					maxcnt = check[cl].conflictPtr[j].cnt;
					finalID = check[cl].conflictPtr[j].finalID;
				}
			}
			if (check[cl].cnt > maxcnt)
				finalID = check[cl].finalID;
			
			/* now that we have maxcnt, and the associated finalID desired, let's attribute to noise inconsistent final clusterIDs */
			for(i = 0;  i< loaded; i++)
			{
				if ((clusterid[i] == cl) && (rslt[i] != finalID))
				{
					clusterid[i] = 0;
					fixedCnt++;
				}
			}
			free(check[cl].conflictPtr);
		}
	}
	if (check)
		free(check);
	return(fixedCnt);

} /* CheckRsltArray */
/* ------------------------------------------------------------------------------------ */

static 	unsigned int FlagSequencesToReassign(unsigned int *clusterid,unsigned int *rslt,unsigned int loaded)
{
	unsigned int i;
	unsigned int cnt = 0;

	for(i = 0;  i< loaded; i++)
	{
		if (rslt[i] == 0) // make sure that sequence was assigned at last step.
			continue;
		if (clusterid[i]==0) // need to reassign.
		{
			clusterid[i] = 9999999;
			cnt++;
		}
	}
	return(cnt);
	
} /* FlagSequencesToReassign */
/* ------------------------------------------------------------------------------------ */

static unsigned int CheckBranchHistory(SEQDATA *fasta,int maxclusterid,unsigned int loaded,unsigned int *clusterid,unsigned int *rslt,unsigned int nsteps,DISTLIST *distList,unsigned int unassignBranchErrors)
{
	unsigned int ii;
	BRANCH_HEAD *branchHistoryCheck = NULL;
	unsigned int branchErrorCnt;
	unsigned int step;
	unsigned int cnt;


		branchHistoryCheck = malloc((maxclusterid+1)*sizeof(BRANCH_HEAD));
		if (!branchHistoryCheck)
			return(0);
			
		cnt = 1;  // start new cluster names at 1
		for (ii = 0; ii<clusterhistorycnt; ii++)
		{
			if (clusterhistory[ii].retain == 'y')
			{
				for (step = 0; step < nsteps; step++)
				{
					if (distList[step].dist == clusterhistory[ii].dist)
						break;
				}
				branchHistoryCheck[cnt].step = step;
				branchHistoryCheck[cnt].dist = clusterhistory[ii].dist;
				branchHistoryCheck[cnt].srcCluster = clusterhistory[ii].cluster;
				cnt++;
			}
		}
			
		/* check that final cluster of each sequence is valid compared to its cluster id at base. */
		branchErrorCnt = 0;
		for(ii = 0;  ii< loaded; ii++)
		{
			if (clusterid[ii] == 0)
				continue;
			if (rslt[branchHistoryCheck[clusterid[ii]].step*loaded+ii] == 0)
				continue;
			if (branchHistoryCheck[clusterid[ii]].srcCluster != rslt[branchHistoryCheck[clusterid[ii]].step*loaded+ii])
			{
				clusterid[ii] = 99999999;
				branchErrorCnt++;
			}	
		}
		/* for each branch error, identify closest sequence an copy its history */
		for(ii = 0;  ii< loaded; ii++)
		{
			SEQDATA *fastapi;
			unsigned int dmin;
			unsigned int jmin;
			unsigned int j;
			
			if (clusterid[ii] != 99999999)
				continue;
			
			// compute distance of branch error (fastapi) to each assigned (fastapj).
			fastapi = &fasta[ii];
			dmin = 0xffffffff;
			jmin = ii; // just to initialize to make compiler happy.
			for (j = 0; j < loaded; j++)
			{
				if ((clusterid[j] > 0) && (clusterid[j] <= maxclusterid)) /* assigned, and without branch error */
				{
					unsigned int d = 0;
					unsigned int col;
					SEQDATA *fastapj = &fasta[j];
					for (col=0;col < kMaxInputCol;col++)
					{
						int diff = (int)fastapj->data[col] - fastapi->data[col];
						d += diff*diff;
					}
					if (d < dmin)
					{
						jmin = j;
						dmin = d;
						if (d == 0)
						{
							break;
						}
					}
				} // if
			} // j
			if (unassignBranchErrors)
				clusterid[ii] = 0 + maxclusterid + 1;  // update clusterid on the fly with an ID beyond maxclusterid, to make sure not taken into account in the loop above
			else
				clusterid[ii] = clusterid[jmin] + maxclusterid + 1;  // update clusterid on the fly with an ID beyond maxclusterid, to make sure not taken into account in the loop above
				
		}
		for(ii = 0;  ii< loaded; ii++)
		{
			if (clusterid[ii] > maxclusterid) /* was reassigned */
				clusterid[ii] -= (maxclusterid + 1);  
		}
		free(branchHistoryCheck);
		return(branchErrorCnt);

} /* CheckBranchHistory */
/* ------------------------------------------------------------------------------------ */

static unsigned int ProcessMergeRequests(unsigned int *clusterid,unsigned int loaded,unsigned int mergerequestcnt,unsigned int nproc,int previouslyRetainedClusterCnt,unsigned int pass)
{
	unsigned int ii;
	unsigned  int isMergingPreexistingClusters = 0;

	if (previouslyRetainedClusterCnt > -1)
	{
		/* identifies which preexisting clusters might have been merged and should subsequently be split */
		for (ii = 0; ii<mergerequestcnt; ii++)
		{
			if (mergerequest[ii].cluster2 <= (1*kStartLocalCluster+previouslyRetainedClusterCnt))
			{
				unsigned int b;

				isMergingPreexistingClusters = 1;
				for (b = 0; b< clusterhistorycnt; b++)
				{
					if ((clusterhistory[b].pass == (pass-1)) && clusterhistory[b].cluster == (mergerequest[ii].cluster2-(1*kStartLocalCluster)))
						clusterhistory[b].mergedto = (mergerequest[ii].cluster1-(1*kStartLocalCluster));
				}
			}
		}
	}
	
{
	int i,j,mrg;
	// init has table: no change.
	for (i = 1; i < nproc; i++)
	{
		int *cnp = clustersnum[i];
		for (j = 1; j<=cnp[0]; j++)
		{
			cnp[j]=i*kStartLocalCluster+j;
		}
	}
	// update table according to merge request processing.
	for (mrg = mergerequestcnt-1; mrg >= 0 ; mrg--)
	{
		for (i = 1; i < nproc; i++)
		{
			int *cnp = clustersnum[i];
			for (j = 1; j<=cnp[0]; j++)
			{
				if (cnp[j]==mergerequest[mrg].cluster2)
					cnp[j] = mergerequest[mrg].cluster1;
			}
		}
	}
	// update clusterids
	for(ii = 0;  ii< loaded; ii++)
	{
		if (clusterid[ii] > 0)
		{
			int cpuid = (int)clusterid[ii] / kStartLocalCluster;
			int hashidx = (int)clusterid[ii] - cpuid * kStartLocalCluster;
			int *cnp = clustersnum[cpuid];
			clusterid[ii] = cnp[hashidx];
		}
	}
	// get back to calloc state
	for (i = 1; i < nproc; i++)
	{
		int *cnp = clustersnum[i];
		bzero(&cnp[1], cnp[0]*sizeof(int));
	}
}
	
	return(isMergingPreexistingClusters);
	
} /* ProcessMergeRequests */
/* ------------------------------------------------------------------------------------ */

static void InsertMergeRequest(unsigned int cluster1, unsigned int cluster2)  /* note: cluster2 > cluster1 is always true */
{
		if (mergerequestcnt == 0)
		{
			mergerequest[0].cluster1 = cluster1;
			mergerequest[0].cluster2 = cluster2;
			mergerequestcnt++;
			return;
		}
		
		
		if (mergerequestcnt < kMaxMergeRequests)
		{
			/* quickly find position of cluster2 - 1 using a binary divide and conquer... */
			unsigned int ii;
			         int left = 0;
			         int right = mergerequestcnt-1;
			unsigned int val=cluster2-1;
			do
			{
				ii=(left+right)>>1;
				if (mergerequest[ii].cluster2 > val)
					right=ii-1;
				else
					left=ii+1;
			} while ((mergerequest[ii].cluster2 != val) && (left <= right));
			
			/* move up from this point to do our real check, using cluster2... */
			/* test if and where to insert request */
			for( ;  ii< mergerequestcnt; ii++)
			{
				if 	 (mergerequest[ii].cluster2 > cluster2) 
					break;  /* we can insert and be done */
				if (mergerequest[ii].cluster2 == cluster2) 
				{
					if (mergerequest[ii].cluster1 == cluster1)
						return; /* already known */
					if (mergerequest[ii].cluster1 > cluster1)
					{
						unsigned int tmp = mergerequest[ii].cluster1;
						mergerequest[ii].cluster1 = cluster1;
						InsertMergeRequest(cluster1,tmp);
						return;
					}
					else
					{
						InsertMergeRequest(mergerequest[ii].cluster1,cluster1);
						return;
					}

				}
			}
			
			/* insert request */
			memmove(&mergerequest[ii+1].cluster1, &mergerequest[ii].cluster1, (mergerequestcnt-ii)*sizeof(MERGECLUSTER));

			mergerequest[ii].cluster1 = cluster1;
			mergerequest[ii].cluster2 = cluster2;
			mergerequestcnt++;
		}
		else
		{
			if (printwarnmergereq == 1)
			{
				fprintf(stderr,"LOG: ERROR: Too many merging requests (max = %d)\n",kMaxMergeRequests);
				printwarnmergereq = 0;
			}
		}

} /* InsertMergeRequest */
/* ------------------------------------------------------------------------------------ */

static void InsertMergeRequestWhere(unsigned int cluster1, unsigned int cluster2,unsigned int *where)  /* note: cluster2 > cluster1 is always true */
{

		if (mergerequestcnt < kMaxMergeRequests)
		{
			unsigned int ii;
			
			/* test if and where to insert request (start from where we inserted previous one as lists are sorted according to cluster2  */
			for(ii = *where;  ii< mergerequestcnt; ii++)
			{
				if 	 (mergerequest[ii].cluster2 > cluster2) 
					break;  /* we can insert and be done */
				if (mergerequest[ii].cluster2 == cluster2) 
				{
					if (mergerequest[ii].cluster1 == cluster1)
						return; /* already known */
					if (mergerequest[ii].cluster1 > cluster1)
					{
						unsigned int tmp = mergerequest[ii].cluster1;
						mergerequest[ii].cluster1 = cluster1;
						InsertMergeRequest(cluster1,tmp);
						return;
					}
					else
					{
						InsertMergeRequest(mergerequest[ii].cluster1,cluster1);
						return;
					}

				}
			}
			
			/* insert request */
			memmove(&mergerequest[ii+1].cluster1, &mergerequest[ii].cluster1, (mergerequestcnt-ii)*sizeof(MERGECLUSTER));
			mergerequest[ii].cluster1 = cluster1;
			mergerequest[ii].cluster2 = cluster2;
			mergerequestcnt++;
			*where = ii; /* remember where we did insert to save time */
		}
		else
		{
			if (printwarnmergereq == 1)
			{
				fprintf(stderr,"LOG: ERROR: Too many merging requests (max = %d)\n",kMaxMergeRequests);
				printwarnmergereq = 0;
			}
		}


} /* InsertMergeRequestWhere */
/* ------------------------------------------------------------------------------------ */

static void computesim(SEQDATA *fasta, unsigned int *clusterid, CPU *chunk)
{
	SEQDATA *fastapi;
	SEQDATA *fastapj;
	unsigned int *clusterpi;
	unsigned int *clusterpj;
	unsigned int i,j;
	unsigned int starti,startj;
	unsigned int lasti,lastj;
	unsigned  int actualstartj;


	starti = chunk->ii;
	startj = chunk->jj;
	lasti = chunk->iilast;
	lastj = chunk->jjlast;


#ifdef USE_THREADS
		if (chunk->ii != chunk->jj)
		{
				unsigned short valii,valjj;

				fastapi = &fasta[lasti-1];
				fastapj = &fasta[startj];
				valii = fastapi->data[sortkey];
				valjj = fastapj->data[sortkey];
				if (valjj > valii)
				{
					unsigned  int mindist = valjj-valii;
					mindist *= mindist;
					if (mindist > gTestDist)
					{
						goto skipthischunk;
					}
				}
		}		
#endif
		
		fastapi = &fasta[starti];
		clusterpi = &clusterid[starti];
		for(i = starti;  i< lasti; i++)
		{
			if (startj <= i)
				actualstartj = i+1;
			else 
				actualstartj = startj;
			fastapj = &fasta[actualstartj];
			clusterpj = &clusterid[actualstartj];
			for(j = actualstartj;  j< lastj; j++)
			{
				int diff;
				unsigned int d;

				if ((*clusterpj) && (*clusterpj == *clusterpi))
					goto straight2next;


				diff = (int)fastapj->data[0] - fastapi->data[0];
				d = diff*diff;
				diff = (int)fastapj->data[1] - fastapi->data[1];
				d += diff*diff;
				diff = (int)fastapj->data[2] - fastapi->data[2];
				d += diff*diff;
				diff = (int)fastapj->data[3] - fastapi->data[3];
				d += diff*diff;
				diff = (int)fastapj->data[4] - fastapi->data[4];
				d += diff*diff;
				diff = (int)fastapj->data[5] - fastapi->data[5];
				d += diff*diff;
				diff = (int)fastapj->data[6] - fastapi->data[6];
				d += diff*diff;
				diff = (int)fastapj->data[7] - fastapi->data[7];
				d += diff*diff;
				diff = (int)fastapj->data[8] - fastapi->data[8];
				d += diff*diff;
				diff = (int)fastapj->data[9] - fastapi->data[9];
				d += diff*diff;
				diff = (int)fastapj->data[10] - fastapi->data[10];
				d += diff*diff;
				diff = (int)fastapj->data[11] - fastapi->data[11];
				d += diff*diff;
				diff = (int)fastapj->data[12] - fastapi->data[12];
				d += diff*diff;
				diff = (int)fastapj->data[13] - fastapi->data[13];
				d += diff*diff;
				diff = (int)fastapj->data[14] - fastapi->data[14];
				d += diff*diff;
				diff = (int)fastapj->data[15] - fastapi->data[15];
				d += diff*diff;

				/* computation done for the first 16 columns, time to attempt to gain some time */
				if (d > gTestDist)
					goto straight2next;
				#define kDATAOFFSET 16
				diff = (int)fastapj->data[kDATAOFFSET+0] - fastapi->data[kDATAOFFSET+0];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+1] - fastapi->data[kDATAOFFSET+1];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+2] - fastapi->data[kDATAOFFSET+2];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+3] - fastapi->data[kDATAOFFSET+3];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+4] - fastapi->data[kDATAOFFSET+4];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+5] - fastapi->data[kDATAOFFSET+5];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+6] - fastapi->data[kDATAOFFSET+6];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+7] - fastapi->data[kDATAOFFSET+7];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+8] - fastapi->data[kDATAOFFSET+8];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+9] - fastapi->data[kDATAOFFSET+9];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+10] - fastapi->data[kDATAOFFSET+10];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+11] - fastapi->data[kDATAOFFSET+11];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+12] - fastapi->data[kDATAOFFSET+12];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+13] - fastapi->data[kDATAOFFSET+13];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+14] - fastapi->data[kDATAOFFSET+14];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET+15] - fastapi->data[kDATAOFFSET+15];
				d += diff*diff;

				/* computation done for the first 32 columns, time to attempt to gain some time */
				if (d > gTestDist)
					goto straight2next;

				#define kDATAOFFSET32 32
				diff = (int)fastapj->data[kDATAOFFSET32+0] - fastapi->data[kDATAOFFSET32+0];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+1] - fastapi->data[kDATAOFFSET32+1];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+2] - fastapi->data[kDATAOFFSET32+2];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+3] - fastapi->data[kDATAOFFSET32+3];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+4] - fastapi->data[kDATAOFFSET32+4];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+5] - fastapi->data[kDATAOFFSET32+5];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+6] - fastapi->data[kDATAOFFSET32+6];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+7] - fastapi->data[kDATAOFFSET32+7];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+8] - fastapi->data[kDATAOFFSET32+8];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+9] - fastapi->data[kDATAOFFSET32+9];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+10] - fastapi->data[kDATAOFFSET32+10];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+11] - fastapi->data[kDATAOFFSET32+11];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+12] - fastapi->data[kDATAOFFSET32+12];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+13] - fastapi->data[kDATAOFFSET32+13];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+14] - fastapi->data[kDATAOFFSET32+14];
				d += diff*diff;
				diff = (int)fastapj->data[kDATAOFFSET32+15] - fastapi->data[kDATAOFFSET32+15];
				d += diff*diff;

	
				if (d <= gTestDist) 
				{
					if (*clusterpi)
					{
						if (*clusterpj == 0)
						{
							*clusterpj = *clusterpi;
						}
						else /* i=assigned, j=assigned */
						{
							{
								unsigned int cluster1,cluster2;
								if(*clusterpj > *clusterpi)
								{
									cluster1 = *clusterpi;
									cluster2 = *clusterpj;
								}
								else
								{
									cluster1 = *clusterpj;
									cluster2 = *clusterpi;
								}
								InsertMergeRequest(cluster1,cluster2);
							}
						}
						
					}
					else
					{
						if (*clusterpj == 0)   /* i=not yes assigned, j=not yet assigned */
						{
	#ifdef USE_THREADS
							pthread_mutex_lock(&clustercntmutex);
	#endif
							*clusterpi = ++clustercnt;
							*clusterpj = clustercnt;
	#ifdef USE_THREADS
							pthread_mutex_unlock(&clustercntmutex);
	#endif
						}
						else /* i=not yes assigned, j=assigned */
						{
							*clusterpi = *clusterpj;
						}
					}
				}
		straight2next:		
				fastapj++;
				clusterpj++;
			}
			fastapi++;
			clusterpi++;
		}
#ifdef USE_THREADS
skipthischunk: ;
#endif
	
} /* computesim */

/* ------------------------------------------------------------------------------------ */
static void DoComputingSlave(SEQDATA *fastadata,unsigned int *clusterid,int idproc,unsigned int initialClusterCnt)
{
			CPU	  cpudata;
			MPI_Request mpireq;

			/* determine the first value to use to start recording new clusterids for this processor */
			clustercnt = (idproc*kStartLocalCluster+initialClusterCnt);
			
			do
			{
				MPI_Recv(&cpudata, 4, MPI_INT,  0, kWhichBlocksToCompute, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				if (cpudata.ii != kNoMoreBlocks)
				{
					int sndcnt,rcvcnt;

					/* ------ update clusterid for each data of that might be touched by the process */

					rcvcnt =  (cpudata.iilast-cpudata.ii);
					MPI_Recv(&clusterid[cpudata.ii],rcvcnt, MPI_INT,  0, kClusterMsg1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					if (cpudata.jj != cpudata.ii)
					{
						rcvcnt =  (cpudata.jjlast-cpudata.jj);
						MPI_Recv(&clusterid[cpudata.jj],rcvcnt, MPI_INT,  0, kClusterMsg2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					}

					/* ------ do the heavy computation */


#ifdef USE_THREADS
					{
						
						if (cpudata.jj != cpudata.ii)
						{
							pthread_t thread;
							EXECUTIONPLAN ep;
							CPU			cpudatasection;
							unsigned int tmpl;
							unsigned int where;

							ep.fastadata = fastadata;
							ep.clusterid = clusterid;
							
							// block 1 against 1
							ep.chunk.ii = ep.chunk.iilast = cpudata.ii;
							ep.chunk.iilast += ((cpudata.iilast - cpudata.ii) >> 1);
							ep.chunk.jj = ep.chunk.jjlast = cpudata.jj;
							ep.chunk.jjlast += ((cpudata.jjlast - cpudata.jj) >> 1);

							thread_mergerequestcnt = 0;
							thread_mergerequest[0].cluster2 = UINT_MAX;  /* to avoid need for initial test thread_mergerequestcnt == 0 in thread_InsertMergeRequest */
							if (pthread_create (&thread, NULL, &computesim_funcion, &ep)) 
								fprintf(stderr,"Error: Failed creating thread\n");

							// block 2 against 2
							cpudatasection.ii = cpudata.ii + ((cpudata.iilast - cpudata.ii) >> 1) ;
							cpudatasection.iilast = cpudata.iilast;
							cpudatasection.jj = cpudata.jj + ((cpudata.jjlast - cpudata.jj) >> 1) ;
							cpudatasection.jjlast = cpudata.jjlast;
							computesim(fastadata,clusterid,&cpudatasection);

							if (pthread_join (thread, NULL))
								fprintf(stderr,"Error: Failed pthread_join\n");

							where = 0;
							for (tmpl = 0; tmpl < thread_mergerequestcnt; tmpl++)
							{
								InsertMergeRequestWhere(thread_mergerequest[tmpl].cluster1,thread_mergerequest[tmpl].cluster2,&where);
							}

							// block 1 against 2
							ep.chunk.jj = ep.chunk.jjlast;
							ep.chunk.jjlast = cpudata.jjlast;

							thread_mergerequestcnt = 0;
							if (pthread_create (&thread, NULL, &computesim_funcion, &ep)) 
								fprintf(stderr,"Error: Failed creating thread\n");

							// block 2 against 1
							cpudatasection.jjlast = cpudatasection.jj;
							cpudatasection.jj = cpudata.jj;
							computesim(fastadata,clusterid,&cpudatasection);

							if (pthread_join (thread, NULL))
								fprintf(stderr,"Error: Failed pthread_join\n");

							where = 0;
							for (tmpl = 0; tmpl < thread_mergerequestcnt; tmpl++)
							{
								InsertMergeRequestWhere(thread_mergerequest[tmpl].cluster1,thread_mergerequest[tmpl].cluster2,&where);
							}

						}
						else 
							computesim(fastadata,clusterid,&cpudata);

					}
#else
					computesim(fastadata,clusterid,&cpudata);
#endif




					/* ------ send back updated data */

					MPI_Isend(&idproc, 1, MPI_INT,  0, kCPUdoneMsg, MPI_COMM_WORLD,&mpireq);
					sndcnt = cpudata.iilast-cpudata.ii;
					MPI_Send(&clusterid[cpudata.ii], sndcnt, MPI_INT,  0, kClusterMsg1, MPI_COMM_WORLD);
					if (cpudata.jj != cpudata.ii)
					{
						sndcnt = cpudata.jjlast-cpudata.jj;
						MPI_Send(&clusterid[cpudata.jj], sndcnt, MPI_INT,  0, kClusterMsg2, MPI_COMM_WORLD);
					}
				}
				else /* send the final count of "new clusters" allocated by this proc. */
				{
					if (cpudata.jj != kNoMoreBlocks)
					{
						MPI_Send((void*)&clustercnt, 1, MPI_INT,  0, kFinalCntRequest, MPI_COMM_WORLD);
					}
				}
			} while (cpudata.ii != kNoMoreBlocks);

			/* ------- JOIN mergerequest list (cpus obtain the list from an other cpu until everything is merged) ------- */

			do
			{
				JOINREQUEST joinRequest;

				MPI_Recv(&joinRequest, 2, MPI_INT,  0, kJoinListRequest, MPI_COMM_WORLD, MPI_STATUS_IGNORE/*&status*/);

			if (joinRequest.getFromCPU == 0)
				break;
				
				if (joinRequest.getFromCPU == idproc) /*  we are  the one to send */
				{
					MPI_Send(&mergerequestcnt, 1, MPI_INT,  joinRequest.sendToCPU, kFinalMergeRequestCntRequest, MPI_COMM_WORLD);
					if (mergerequestcnt > 0)
						MPI_Send(&mergerequest, mergerequestcnt*2, MPI_INT,  joinRequest.sendToCPU, kSendFinalMergeRequestRequest, MPI_COMM_WORLD);
				}
				else /* we are the one receiving */
				{
					MERGECLUSTER cpumergerequest[kMaxMergeRequests];
					unsigned int cpumergerequestcnt;
					unsigned int where;
					int tmpl;
					MPI_Recv(&cpumergerequestcnt, 1, MPI_INT,  joinRequest.getFromCPU, kFinalMergeRequestCntRequest, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					if (cpumergerequestcnt > 0)
						MPI_Recv(&cpumergerequest, cpumergerequestcnt*2, MPI_INT,  joinRequest.getFromCPU, kSendFinalMergeRequestRequest, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

					where = 0;
					for (tmpl = 0; tmpl < cpumergerequestcnt; tmpl++)
					{
						InsertMergeRequestWhere(cpumergerequest[tmpl].cluster1,cpumergerequest[tmpl].cluster2,&where);
					}
					MPI_Send(&mergerequestcnt,1, MPI_INT,  0, kJoinListRequestDone, MPI_COMM_WORLD);
				}
			} while(1);
			
			if (idproc == 1) /* send final list to master */
			{
				MPI_Send(&mergerequestcnt, 1, MPI_INT,  0, kFinalMergeRequestCntRequest, MPI_COMM_WORLD);
				if (mergerequestcnt > 0)
					MPI_Send(&mergerequest, mergerequestcnt*2, MPI_INT,  0, kSendFinalMergeRequestRequest, MPI_COMM_WORLD);
			}


	
} /* DoComputingSlave */
/* ------------------------------------------------------------------------------------ */


/* function repeatadly called only by the master to identify a suitable computing chunk to asign to an available slave */
static unsigned int AssignChunk(unsigned int *clusterid,CHUNK *chunk, CPU *cpu, unsigned int nproc)
{
	unsigned int i;
	int sndcnt;
	MPI_Request mpireq;
	unsigned int candidate = 0;

	/* test if computing in already in progress somewhere for one of those blocks */
	// start from last to first proc, as according to lsf, first proc will have lowest load and we submit to the last identified avail proc.
	for (i = (nproc-1); i>0; i--)
	{
		if ((cpu[i].ii == chunk->ii) || (cpu[i].jj == chunk->jj) || (cpu[i].ii == chunk->jj) || (cpu[i].jj == chunk->ii))
			return(0); 

		if (cpu[i].ii == kCPU_availaible)
			candidate = i;
	}

	cpu[candidate].ii = chunk->ii;
	cpu[candidate].jj = chunk->jj;
	cpu[candidate].iilast = chunk->iilast;
	cpu[candidate].jjlast = chunk->jjlast;

	MPI_Isend(&cpu[candidate], 4, MPI_INT,  candidate, kWhichBlocksToCompute, MPI_COMM_WORLD,&mpireq);

	sndcnt = cpu[candidate].iilast-cpu[candidate].ii;
	MPI_Isend(&clusterid[cpu[candidate].ii], sndcnt, MPI_INT,  candidate, kClusterMsg1, MPI_COMM_WORLD,&mpireq);
	if (cpu[candidate].ii != cpu[candidate].jj)
	{
		sndcnt = cpu[candidate].jjlast-cpu[candidate].jj;
		MPI_Isend(&clusterid[cpu[candidate].jj], sndcnt, MPI_INT,  candidate, kClusterMsg2, MPI_COMM_WORLD,&mpireq);
	}
	chunk->status = kChunkStatusComputing;

	return(candidate);
	
} /* AssignChunk */
/* ------------------------------------------------------------------------------------ */

static void DoSlave(SEQDATA *fastadata,unsigned int *clusterid,int idproc,int nproc, unsigned int nbSeqForClustering,unsigned int loadedrefseqcnt)
{
			unsigned int initialClusterCnt = 0;

			if (idproc == 1)
			{
				int i;
				for (i = (int)nbSeqForClustering-1; i >= 0; i--)
				{
					if (clusterid[i] != 0)
					{
						initialClusterCnt = clusterid[i]-idproc*kStartLocalCluster;
						break;
					}
				}
			}


			do
			{
				mergerequestcnt = 0;
				mergerequest[0].cluster2 = UINT_MAX;  /* to avoid need for initial test mergerequestcnt == 0 in InsertMergeRequest */

				DoComputingSlave(fastadata,clusterid,idproc,initialClusterCnt);
				memset(clusterid,0,nbSeqForClustering*sizeof(MPI_INT));  // reset clusterid

				MPI_Recv(&gTestDist, 1, MPI_INT,0, kRepeatWithNewDistMsg,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

				if ((idproc == 1) && (gTestDist > 0))
				{
					MPI_Recv(&initialClusterCnt, 1, MPI_INT,0, kInitialClusterCntMsg,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				fflush(stderr);			

				MPI_Barrier(MPI_COMM_WORLD);
				if (gTestDist == 0)
				{
					CLOSESTSEQ *closestSeq = NULL;
					unsigned int trimmedclustercnt;
					unsigned int starti,lasti;
					unsigned int datachunk = 1+(nbSeqForClustering/nproc);
					MPI_Bcast (&gTestDist, 1, MPI_INT, 0, MPI_COMM_WORLD);

					if (gTestDist != 0)
					{
						MPI_Bcast (&trimmedclustercnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
						MPI_Bcast (&clusterid[0], nbSeqForClustering, MPI_INT, 0, MPI_COMM_WORLD);
						starti = idproc*datachunk;
						if (starti < nbSeqForClustering)
						{
							lasti = starti + datachunk;
							if (lasti > nbSeqForClustering)
								lasti = nbSeqForClustering;
							DistributeUnassignedToClosestCluster(fastadata,clusterid,nbSeqForClustering,trimmedclustercnt,starti,lasti);
							MPI_Send(&clusterid[starti], (lasti-starti), MPI_INT,  0, kClusterMsg1, MPI_COMM_WORLD);
						}
					}

					if (loadedrefseqcnt > 0)
					{
						MPI_Bcast (&clusterid[0], nbSeqForClustering, MPI_INT, 0, MPI_COMM_WORLD);

						/* evenly split refseqs onto computing nodes */
						datachunk = 1+(loadedrefseqcnt/nproc);
						starti = idproc*datachunk;
						if (starti < loadedrefseqcnt)
						{
							lasti = starti + datachunk;
							if (lasti > loadedrefseqcnt)
								lasti = loadedrefseqcnt;

							closestSeq = malloc((lasti-starti)*sizeof(CLOSESTSEQ));
							if (closestSeq)
							{
								memset(closestSeq,0xff,(lasti-starti)*sizeof(CLOSESTSEQ));	// initialize closest seq to -1.
								DistributeRefSeqToClosestCluster(fastadata,clusterid,nbSeqForClustering,nbSeqForClustering+starti,nbSeqForClustering+lasti,closestSeq);
								MPI_Send(&clusterid[nbSeqForClustering+starti], (lasti-starti), MPI_INT,  0, kClusterMsg1, MPI_COMM_WORLD);
								MPI_Send(&closestSeq[0], 2*(lasti-starti), MPI_INT,  0, kClosestMsg, MPI_COMM_WORLD);
							}
							else
								fprintf(stderr,"LOG:Warning: proc %d has not enough memory to assign reference sequences to closest cluster.\n",idproc);
						}
					}
					MPI_Barrier(MPI_COMM_WORLD); 
					if (closestSeq)
						free(closestSeq);
		
					break;
				}
			} while(1);

} // DoSlave
/* ------------------------------------------------------------------------------------ */

int main (int argc, char **argv)
{	
	char	fn[kMaxFilename];
	char	ofn[kMaxFilename];
	char	tfn[kMaxFilename];
	char	tmpdir[kMaxFilename];
	char	tmpfn[kMaxFilename];
	char	rfn[kMaxFilename];
	char	distListString[1024];
	char	*version="VERSION 1.4.4; 2019-12-26";
	SEQNAME *fastaname = NULL;
	SEQDATA *fastadata = NULL;
	FILE	 *f=NULL;
	DISTLIST *distList = NULL;
	unsigned int *sameFollowingCnt = NULL;
	unsigned int loaded;
	unsigned int nbSeqForClustering;
	unsigned int nbSeq = 0;
	unsigned int refseqcnt = 0;
	float distcutoff;
	float distcutoffincreasestep;
	float lastdistcutoff;
	unsigned int nbDistToTest = 0;
	unsigned int currentTestingDistanceIndex = 0;
	int readfromSTDIN = 1;
	int idproc, nproc;
	int verbose;
	unsigned int rowcnt = 0;
	unsigned int loadedrefseqcnt = 0;
	unsigned int cntcutoff = 0;
	int c;
	unsigned int desiredBlockSize = 0;
	struct timeval ts;
	unsigned int ii;
	unsigned int jj;
	unsigned int *clusterid = NULL;
	pid_t	pid;
	unsigned int OutputLevel = 5;
	
	int mrg;
	int *cnp;
	struct timeval te;
	unsigned short key = 0;
	unsigned int printClusterStatus = 0;

	/* must be first instruction */
    if (MPI_Init(&argc, &argv))
		return(1);

	pid = getpid();
	gettimeofday(&ts, NULL); 


		/* --------- initialize MPI */
		nproc=0; 
		
		if (MPI_Comm_rank(MPI_COMM_WORLD, &idproc))
			goto abort;

		if (MPI_Comm_size(MPI_COMM_WORLD, &nproc))
			goto abort;

		if ((nproc >= kMaxCPU) || (nproc < 2))
		{
			if (idproc == 0)
			{
				fprintf(stderr,"Software requires at least 2 cpus (and at most %d cpus)\n",(kMaxCPU-1));
				fprintf(stderr,"Please launch with 'mpirun -n cpus'\n");
				PrintUsage(version);
			}
			goto abort;
		}

		/* --------- Let Master process arguments */
		if (idproc == 0)
		{
			rowcnt = 0;
			loaded = 0;
			fn[0] = 0;
			tfn[0] = 0;
			ofn[0] = 0;
			rfn[0] = 0;
			distListString[0] = 0;

			distcutoff = 3.381;
			distcutoffincreasestep = 0.174;
			lastdistcutoff = 7.209;
			verbose = 0;
			strcpy(tmpdir,"/tmp");
			
			opterr = 0;
			while ((c = getopt (argc, argv, "i:r:n:l:d:e:s:t:c:v:o:M")) != -1)
			switch (c)
			{
			  case 'i':
					strcpy(fn,optarg);
					readfromSTDIN = 0;
				break;

			  case 'r':
					strcpy(rfn,optarg);
				break;

			  case 'n':
					sscanf(optarg,"%u",&cntcutoff);
					if (cntcutoff < 2)
						cntcutoff = 2;
				break;

			  case 'l':
					strcpy(distListString,optarg);
				break;

			  case 'd':
					sscanf(optarg,"%f",&distcutoff);
				break;

			  case 'e':
					sscanf(optarg,"%f",&lastdistcutoff);
				break;

			  case 's':
					sscanf(optarg,"%f",&distcutoffincreasestep);
					if (distcutoffincreasestep < 0.0)
						distcutoffincreasestep = -distcutoffincreasestep;
				break;
				
			  case 't':
					strcpy(tmpdir,optarg);
				break;

			  case 'o':
					strcpy(ofn,optarg);
				break;

			  case 'c':
					sscanf(optarg,"%u",&OutputLevel);
				break;

			  case 'M':
					printClusterStatus = 1;
				break;
					
			  case 'v':
					sscanf(optarg,"%d",&verbose);
					if (verbose > 2)
						verbose = 2;
				break;

			}

			if (cntcutoff == 0)
			{
				PrintUsage(version);
				goto abort;
			}
			if (lastdistcutoff < distcutoff)
				lastdistcutoff = distcutoff; 

	

			/* --------- Master: load fasta sequences ------------------------- */
			sprintf(tfn,"%s/dbc454.tmp.%d",tmpdir,pid);
			fprintf(stderr,"LOG:%s\n",version);
			fprintf(stderr,"LOG:using %d processors\n",nproc);
			if (readfromSTDIN)
				fprintf(stderr,"LOG:Input File=<stdin>\n");
			else
				fprintf(stderr,"LOG:Input File=%s\n",fn);
			if (!ofn[0])
				fprintf(stderr,"LOG:Output File=<stdout>\n");
			else
				fprintf(stderr,"LOG:Output File=%s\n",ofn);
			fprintf(stderr,"LOG:Temporary Files=%s...\n",tfn);
			fprintf(stderr,"LOG:Minumum number of sequences to form a cluster=%u\n",cntcutoff);


			// Check if a distance list was supplied, and get how many levels we will have
			{
				int i = 0;
				while(distListString[i] != 0) { 
					if (distListString[i] == ',')
						nbDistToTest++;
					i++;
				};
				if (nbDistToTest == 0)  // distance list not supplied, thus analyze the begin,end,stp parameters
					nbDistToTest = (int)((lastdistcutoff - distcutoff) / distcutoffincreasestep);


				distList = malloc((nbDistToTest+2)*sizeof(DISTLIST));
				if (!distList)
				{
					fprintf(stderr,"LOG:Not enough memory to allocate the list of distances to evaluate\n");
					goto abort;
				}

				if (distListString[0] == 0)
					nbDistToTest = FillDistListFromRange(distList,distcutoff,lastdistcutoff,distcutoffincreasestep);
				else
					nbDistToTest =FillDistListFromList(distList,distListString,nbDistToTest);
				PrintDistList(distList,nbDistToTest);
			}

			gTestDist = distList[0].squaredDist;

			if (readfromSTDIN)
				f = stdin;
			else
				f = fopen(fn,"r");
			if (!f)
				fprintf(stderr,"LOG:Cannot Open InputFile %s\n",fn);
			else
			{
				SEQNAME *sortedfastaname = NULL;
				SEQDATA *sortedfastadata = NULL;
				fastaname = calloc(kMAXEVENTS,sizeof(SEQNAME));
				fastadata = calloc(kMAXEVENTS,sizeof(SEQDATA));
				if (fastaname && fastadata)
				{
					unsigned short minkeyval = 0; // dummy initialization to make compiler stop whining
					unsigned short maxkeyval = 0; // dummy initialization to make compiler stop whining
					FILE *of = NULL;
					FILE *sf = NULL;
					sprintf(tmpfn,"%s.headers.txt",tfn);
					of = fopen(tmpfn,"w");
					if (OutputLevel >= 4)
					{
						sprintf(tmpfn,"%s.seq.txt",tfn);
						sf = fopen(tmpfn,"w");
					}
					if (!of || (!sf && (OutputLevel >= 4)))
						fprintf(stderr,"LOG:Cannot Write TempFile in '%s'\n",tmpfn);
					else
					{

						rowcnt = EncodeFasta(f,of,sf,fastaname,fastadata,&key,&minkeyval,&maxkeyval,1,nproc); //set first id to 1
						nbSeq = rowcnt;
						fprintf(stderr,"LOG:SequencesCnt=%u\n",nbSeq);
						if (of)
							fclose(of);
						if (sf)
							fclose(sf);
					}
					if (nbSeq > cntcutoff)
					{
						sortedfastaname = calloc(nbSeq,sizeof(SEQNAME));
						sortedfastadata = calloc(nbSeq,sizeof(SEQDATA));
						clusterid = calloc(nbSeq,sizeof(MPI_INT));
						sameFollowingCnt = calloc(nbSeq,sizeof(unsigned int));
						if (sortedfastaname && sortedfastadata && sameFollowingCnt && clusterid)
						{
							SortEncodedFasta(fastaname,fastadata,sortedfastaname,sortedfastadata,nbSeq,key,minkeyval,maxkeyval);
							loaded = DiscardDuplicates(sortedfastaname,sortedfastadata,fastaname,fastadata,clusterid,sameFollowingCnt,nbSeq);
						}
					}
					else
						fprintf(stderr,"LOG:Number of sequences < MinSeqCountInCluster. Aborting.\n");
				}
				if (sortedfastaname)
					free(sortedfastaname);
				if (sortedfastadata)
					free(sortedfastadata);

				nbSeqForClustering = loaded;	

				/* read reference sequences */
				if ((loaded > 0) && (rfn[0] != 0))
				{
					FILE *rf = fopen(rfn,"r");
					if (!rf)
					{
						fprintf(stderr,"LOG:Cannot Open ReferenceFile %s\n",rfn);
						loaded = 0; // force abort.
					}
					else  // append reference sequences at end of dataset
					{
						#define kErrCannotAppend 1
						#define kErrNotEnoughMem 2
						#define kErrTooManySeq 3
						int appenderr = 0;
						FILE *of = NULL;
						FILE *sf = NULL;

						sprintf(tmpfn,"%s.headers.txt",tfn);
						of = fopen(tmpfn,"a");
						if (OutputLevel >= 4)
						{
							sprintf(tmpfn,"%s.seq.txt",tfn);
							sf = fopen(tmpfn,"a");
						}
						if (!of || (!sf && (OutputLevel >= 4)))
							appenderr = kErrCannotAppend;
						else
						{
							SEQNAME *refseqname = NULL;
							SEQDATA *refseqdata = NULL;
							SEQNAME *sortedrefseqname = NULL;
							SEQDATA *sortedrefseqdata = NULL;
							unsigned int *refseqclusterid = NULL;
							unsigned int *refseqsameFollowingCnt = NULL;
							refseqname = calloc(kMAXEVENTS,sizeof(SEQNAME));
							refseqdata = calloc(kMAXEVENTS,sizeof(SEQDATA));
							if ((refseqname == NULL) || (refseqdata == NULL))
								appenderr = kErrNotEnoughMem;
							else
							{
								unsigned short refseqkey = 0; // dummy initialization to make compiler stop whining
								unsigned short refseqminkeyval = 0; // dummy initialization to make compiler stop whining
								unsigned short refseqmaxkeyval = 0; // dummy initialization to make compiler stop whining

								refseqcnt = EncodeFasta(rf,of,sf,refseqname,refseqdata,&refseqkey,&refseqminkeyval,&refseqmaxkeyval,nbSeq+1,nproc); // set first id of reference sequence to nbSeq + 1.
								fprintf(stderr,"LOG:ReferenceSequencesCnt=%u\n",refseqcnt);
								if (of)
									fclose(of);
								if (sf)
									fclose(sf);
								sortedrefseqname = calloc(refseqcnt,sizeof(SEQNAME));
								sortedrefseqdata = calloc(refseqcnt,sizeof(SEQDATA));
								refseqclusterid = calloc(refseqcnt,sizeof(MPI_INT));
								refseqsameFollowingCnt = calloc(refseqcnt,sizeof(unsigned int));
								if (sortedrefseqname && sortedrefseqdata && refseqsameFollowingCnt && refseqclusterid)
								{
									SortEncodedFasta(refseqname,refseqdata,sortedrefseqname,sortedrefseqdata,refseqcnt,refseqkey,refseqminkeyval,refseqmaxkeyval);
									loadedrefseqcnt = DiscardRefSeqDuplicates(sortedrefseqname,sortedrefseqdata,refseqname,refseqdata,refseqsameFollowingCnt,refseqcnt);
								}
								if (sortedrefseqname)
									free(sortedrefseqname);
								if (sortedrefseqdata)
									free(sortedrefseqdata);
									
								/* now append refsequences to regular data */	
								if ((loaded+refseqcnt) >= kMAXEVENTS)
								{
									free(refseqname);
									free(refseqdata);
									appenderr = kErrTooManySeq;
								}
								else
								{

									memcpy(&fastaname[nbSeq],refseqname,refseqcnt*sizeof(SEQNAME));
									memcpy(&fastadata[loaded],refseqdata,loadedrefseqcnt*sizeof(SEQDATA));
									free(refseqname);
									free(refseqdata);

									clusterid = realloc(clusterid,(loaded+loadedrefseqcnt)*sizeof(MPI_INT));
									if (errno == ENOMEM) { free(clusterid) ; clusterid = NULL; }
									sameFollowingCnt = realloc(sameFollowingCnt,(loaded+loadedrefseqcnt/*nbSeq+refseqcnt*/)*sizeof(unsigned int));
									if (errno == ENOMEM) { free(sameFollowingCnt) ; sameFollowingCnt = NULL; }
									if ((clusterid == NULL) || (sameFollowingCnt == NULL))
										appenderr = kErrNotEnoughMem;
									memcpy(&clusterid[loaded/*nbSeq*/],refseqclusterid,loadedrefseqcnt/*refseqcnt*/ *sizeof(MPI_INT));
									memcpy(&sameFollowingCnt[loaded/*nbSeq*/],refseqsameFollowingCnt,loadedrefseqcnt/*refseqcnt*/ *sizeof(unsigned int));
									loaded += loadedrefseqcnt;
									rowcnt += refseqcnt;
								}
								if (refseqclusterid)
									free(refseqclusterid);
								if (refseqsameFollowingCnt)
									free(refseqsameFollowingCnt);
							}
							if (appenderr == kErrCannotAppend)
								fprintf(stderr,"LOG:Cannot Append to TempFile in '%s'\n",tmpfn);
							if (appenderr == kErrNotEnoughMem)
								fprintf(stderr,"LOG:Cannot Allocate memory for reference sequences\n");
							if (appenderr == kErrTooManySeq)
								fprintf(stderr,"LOG:Number of sequences exceed the maximum (%d); aborting\n",kMAXEVENTS);
							if (appenderr > 0)
								loaded = 0; // force abort.

						}//else
					}
				}//loaded > 0

				fprintf(stderr,"LOG:**************************************************************\n");
				
			}
		}// idproc==0
		MPI_Bcast (&nbSeqForClustering, 1, MPI_INT, 0, MPI_COMM_WORLD); // only the experimental dataset is used for clustering
		MPI_Bcast (&loaded, 1, MPI_INT, 0, MPI_COMM_WORLD);				// all "unique" sequences loaded in memory (data + reference)
		MPI_Bcast (&key, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast (&gTestDist, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (loaded == 0)
			goto abort;

#ifdef USE_THREADS
		sortkey = key;
#endif

		/* --------- Slave memory allocation --------------- */
		if (idproc != 0)
		{
			fastadata = calloc(loaded,sizeof(SEQDATA));
			if (!fastadata)
				goto abort;
			clusterid = calloc(loaded,sizeof(MPI_INT));
			if (!clusterid)
				goto abort;
		}
		/* here, should add handshaking to make sure memory was allocated properly on each node. */
		
		/* -------------------- Broadcast Data ------------------------------ */
		{
			unsigned int tosend = loaded;
			unsigned int sendfrom = 0;
			unsigned int smallchunk;

			while(tosend > 0)
			{
			    if (tosend <= 1000000)
			    {
					smallchunk = tosend;
					tosend = 0;
			    }
			    else
			    {
					smallchunk = 1000000;
					tosend -= 1000000;
 			    }	
			    MPI_Bcast (&fastadata[sendfrom], smallchunk*sizeof(SEQDATA), MPI_CHAR, 0, MPI_COMM_WORLD);
			    sendfrom += 1000000;
			}
			MPI_Bcast(&clusterid[0], nbSeqForClustering, MPI_INT,  0, MPI_COMM_WORLD); // no need to broadcast clusterid for refseq, they are zero.
		}
		fflush(stderr);
				

		if (idproc != 0)	/* ---------------------- slave node  -------------------- */
		{
			DoSlave(fastadata,clusterid,idproc,nproc,nbSeqForClustering,(loaded-nbSeqForClustering));
			free(fastadata);
			fastadata = NULL;
		}
		else				/* ---------------- master node ------------- */
		{
			CHUNK *chunk;
			CPU	  cpu[kMaxCPU];
			unsigned int whichcpu;
			unsigned int alldone = 1;  // will be initialized, ignore compiler whining.
			unsigned int submitted;
			unsigned int chunkcnt;
			unsigned int chunckcnt;
			unsigned  int unassigned;
			unsigned int processingBlockSize = 131072;
			int trimmedclustercnt;
			int initialClusterCnt;
			unsigned int passcnt = 0;
			STATS stats[2];
			unsigned int nextrefresh = 0;	// dummy initialization to make compiler stop whining
			unsigned int nextrefreshstep = 0; // dummy initialization to make compiler stop whining


			clusterhistory = malloc(kMaxCluster*sizeof(CLUSTERHISTORY));
			if (!clusterhistory)
			{
				fprintf(stderr,"LOG: ERROR: not enough memory [statically alloacted Err03]\n");
				goto abort;
			}

			if (desiredBlockSize == 0)
			{
				
				/* arrange to keep each slave node busy with at least about 100 computations, but do not go below blocksize of 256 events */
				do 
				{
					chunkcnt = ((nbSeqForClustering/processingBlockSize+2)*(nbSeqForClustering/processingBlockSize+2))/2;
					processingBlockSize >>= 1;
				} while( ((chunkcnt / (nproc-1)) < 100) && (processingBlockSize > 128) );
				processingBlockSize <<= 1;
			}
			else
			{
				processingBlockSize = desiredBlockSize;
				chunkcnt = ((nbSeqForClustering/processingBlockSize+2)*(nbSeqForClustering/processingBlockSize+2))/2;
			}

			chunk = calloc(chunkcnt,sizeof(CHUNK));
			if (!chunk)
			{
				fprintf(stderr,"LOG:Not enough memory to allocate %u chunks; recompile with larger processingBlockSize\n",chunkcnt);
				cpu[0].ii = kNoMoreBlocks;
				cpu[0].jj = kNoMoreBlocks;
				for (ii = 1; ii<nproc; ii++)
					MPI_Send(&cpu[0], 4, MPI_INT,  ii, kWhichBlocksToCompute, MPI_COMM_WORLD);
				goto abort;
			}				
			
			stats[0].dist = 0.0;
			stats[0].rawClustersCnt = -1;
			stats[0].trimmedClustersCnt = -1;
			stats[0].pctAssigned = 0.0;
			
			fprintf(stderr,"DBH:totseconds,squaredDist,distcutoff,nbSeqForClustering,assigned,unassigned,pctassigned,pctunassigned,clustercnt,trimclustercnt\n");

repeatWithNewDist:
			distcutoff = distList[currentTestingDistanceIndex].dist;  // get testing distance from distancelist.
			stats[1].dist = distcutoff;
			stats[1].rawClustersCnt = -1;
			stats[1].trimmedClustersCnt = -1;
			stats[1].pctAssigned = 0.0;

			if (verbose > 0)
				fprintf(stderr,"LOG:DistanceCutoff=%.3f\n",distcutoff);


			mergerequestcnt = 0;
			mergerequest[0].cluster2 = UINT_MAX;  /* to avoid need for initial test mergerequestcnt == 0 in InsertMergeRequest */
			trimmedclustercnt = -1;
			for (ii = 1; ii<nproc; ii++)
			{
				cpu[ii].ii = kCPU_availaible;
				cpu[ii].jj = kCPU_availaible;
			}
			chunckcnt = 0;
			for (ii = 0; ii < nbSeqForClustering; ii += processingBlockSize)
			for (jj = ii; jj < nbSeqForClustering; jj += processingBlockSize)
			{
				SEQDATA *fastapi;
				SEQDATA *fastapj;
				unsigned short valii,valjj;
				
				chunk[chunckcnt].ii = ii;
				chunk[chunckcnt].jj = jj;
				if ((ii+processingBlockSize) < nbSeqForClustering)
					chunk[chunckcnt].iilast = ii+processingBlockSize; 
				else
					chunk[chunckcnt].iilast  = nbSeqForClustering;
				if ((jj+processingBlockSize) < nbSeqForClustering)
					chunk[chunckcnt].jjlast = jj+processingBlockSize;
				else
					chunk[chunckcnt].jjlast  = nbSeqForClustering;
					
				if (ii != jj)
				{
					fastapi = &fastadata[chunk[chunckcnt].iilast-1];
					fastapj = &fastadata[jj];
					valii = fastapi->data[key];
					valjj = fastapj->data[key];
					if (valjj > valii)
					{
						unsigned  int mindist = valjj-valii;
						mindist *= mindist;
						if (mindist > gTestDist)
						{
							continue;
						}
					}
				}
				chunk[chunckcnt].status = kChunkStatusToDo;
				chunckcnt++;
			}

			submitted = 0;
			fflush(stderr);
			if (verbose > 1)
			{
				nextrefresh = chunckcnt / 50;
				nextrefreshstep = nextrefresh+1;
				fprintf(stderr,"LOG:");
			}
			do
			{				
				/* if at least one cpu idle, try to submit as many jobs as possible */
				if (submitted < (nproc-1))
				{
					alldone = 1;
					for (ii = 0; ii < chunckcnt; ii++)
					{
						if (chunk[ii].status == kChunkStatusToDo)
						{
							alldone = 0;
							if (AssignChunk(clusterid,&chunk[ii],cpu,nproc))
							{
									submitted++;
									if ((verbose > 1) && (--nextrefresh == 0))
									{
										fputc('.',stderr);
										fflush(stderr);
										nextrefresh = nextrefreshstep;
									}
									if (submitted == (nproc-1))
										break;
							}
						}
					}
				}
				/* wait until one of the slave cpu is done and update available list accordingly */
				if (submitted)
				{
					int rcvcnt;

					MPI_Recv(&whichcpu, 1, MPI_INT,  MPI_ANY_SOURCE, kCPUdoneMsg, MPI_COMM_WORLD, MPI_STATUS_IGNORE/*&status*/);
					rcvcnt =  (cpu[whichcpu].iilast-cpu[whichcpu].ii);
					MPI_Recv(&clusterid[cpu[whichcpu].ii],rcvcnt, MPI_INT,  whichcpu, kClusterMsg1, MPI_COMM_WORLD, MPI_STATUS_IGNORE/*&status*/);
					if (cpu[whichcpu].jj != cpu[whichcpu].ii)
					{
						rcvcnt =  (cpu[whichcpu].jjlast-cpu[whichcpu].jj);
						MPI_Recv(&clusterid[cpu[whichcpu].jj],rcvcnt, MPI_INT,  whichcpu, kClusterMsg2, MPI_COMM_WORLD, MPI_STATUS_IGNORE/*&status*/);
					}
					cpu[whichcpu].ii = cpu[whichcpu].jj = kCPU_availaible;
					submitted--;
				}
			} while ((submitted > 0) || (!alldone));

			if (verbose > 1)
				fputc('\n',stderr);
			
			/* signal to all nodes that they should clean up */
			/* collect cluster number assigned by each proc and adjust clusters from 1..clustercnt */
			cpu[0].ii = kNoMoreBlocks;
			cpu[0].jj = 0;
			for (ii = 1; ii<nproc; ii++)
			{
				int finalCPUcnt;
				MPI_Send(&cpu[0], 4, MPI_INT,  ii, kWhichBlocksToCompute, MPI_COMM_WORLD);
				MPI_Recv(&finalCPUcnt, 1, MPI_INT,  ii, kFinalCntRequest, MPI_COMM_WORLD, MPI_STATUS_IGNORE/*&status*/);
				finalCPUcnt -= (ii*kStartLocalCluster);

				if (verbose > 2)
				{
					fprintf(stderr,"LOG:Final Cluster Cnt For CPU %3u = %6d\n",ii,finalCPUcnt);
					fflush(stderr);
				}

				clustersnum[ii] = calloc((finalCPUcnt+1),sizeof(int));
				if (!clustersnum[ii])
				{
					fprintf(stderr,"LOG:Not enough memory to allocate cluster ID %u\n",ii);
					goto abort;
				}				
				cnp = clustersnum[ii];
				*cnp = finalCPUcnt;

				if (ii == 1)
				{
					clustercnt = finalCPUcnt;
				}
				else
				{
					clustercnt += finalCPUcnt;
				}
			}

			DoMergeLists(nproc,verbose);
			MPI_Recv(&mergerequestcnt, 1, MPI_INT,  1, kFinalMergeRequestCntRequest, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			if (mergerequestcnt > 0)
				MPI_Recv(&mergerequest, mergerequestcnt*2, MPI_INT,  1, kSendFinalMergeRequestRequest, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

			gettimeofday(&te, NULL);
			ProcessMergeRequests(clusterid,nbSeqForClustering,mergerequestcnt,nproc,stats[0].trimmedClustersCnt,passcnt);
			gettimeofday(&te, NULL);

			for (mrg = mergerequestcnt-1; mrg >= 0 ; mrg--)
			{
				removeclustersnum(mergerequest[mrg].cluster2);
			}
			
			
			// Discarding clusters with too few events

			trimmedclustercnt = RemoveSmallClusters(clusterid,sameFollowingCnt,nbSeqForClustering,nproc,cntcutoff);
			stats[1].trimmedClustersCnt = trimmedclustercnt;
			if (verbose > 0)
				fprintf(stderr,"LOG: %12d Clusters retained ( out of %u )\n",trimmedclustercnt,(clustercnt-mergerequestcnt));

			UpdateClusterHistory(passcnt,trimmedclustercnt,stats[0].trimmedClustersCnt,distcutoff);
			AdjustClustersID(clusterid,nbSeqForClustering,nproc,trimmedclustercnt,&initialClusterCnt);
			unassigned = nbSeqForClustering;

			// always write rslt of distance computation.
			WriteClusterIndices(clusterid,nbSeqForClustering,distcutoff,tfn);

			if (trimmedclustercnt == 0)
			{
				fprintf(stderr,"LOG: %12d Assigned        (%5.1f %%)\n",0,0.0);	
				fprintf(stderr,"LOG: %12u Unassigned      (%5.1f %%)\n",unassigned,100.0);	
				fflush(stderr);
			}
			else
			{				
				unassigned = nbSeq-CountAssigned(clusterid,sameFollowingCnt,nbSeqForClustering,trimmedclustercnt);
				stats[1].pctAssigned = (float)100.0*(nbSeq-unassigned)/nbSeq;
				if (verbose > 0)
				{
					fprintf(stderr,"LOG: %12u Assigned        (%5.1f %%)\n",(nbSeq-unassigned),stats[1].pctAssigned);	
					fprintf(stderr,"LOG: %12u Unassigned      (%5.1f %%)\n",unassigned,100.0-stats[1].pctAssigned);	
					fflush(stderr);
				}

			}
						
			for (ii = 1; ii<nproc; ii++)
			{
				if (clustersnum[ii])
					free(clustersnum[ii]);
			}

			gettimeofday(&te, NULL);
		
			fprintf(stderr,"DBV:%d,%u,%.3f,%u,%u,%u,%.3f,%.3f,%u,%d\n",((int)te.tv_sec-(int)ts.tv_sec),distList[currentTestingDistanceIndex].squaredDist,distcutoff,nbSeq,nbSeq-(unassigned),(unassigned),stats[1].pctAssigned,(100.0 - stats[1].pctAssigned),(clustercnt-mergerequestcnt),trimmedclustercnt);
			fflush(stderr);

			/* test if should keep scanning */
			ii = 0;

			if (nbDistToTest == 1)
				ii = 1; // all done.
			else
			{
				currentTestingDistanceIndex++;
				stats[1].rawClustersCnt = (clustercnt-mergerequestcnt);						
				if (currentTestingDistanceIndex >= nbDistToTest) // all done
					ii = 1;
			}
			
			if (ii == 1)  // all done.
			{
				fprintf(stderr,"LOG:**************************************************************\n");

				if (nbDistToTest > 1) // skip this if we did not test multiple distances because there is no history, and all clusters would be reset to zero.
				{
					memset(clusterid,0,nbSeqForClustering*sizeof(MPI_INT));  // reset all clusterid (ref sequences not assigned at this point, still zero, so ignore them)
					trimmedclustercnt = SelectClusterHistory(nbSeqForClustering,clusterid,tfn);
					if (printClusterStatus)
						PrintClusterStatus();
				}
				free(chunk); 
				gTestDist = 0;
				for (ii = 1; ii<nproc; ii++)
					MPI_Send(&gTestDist, 1, MPI_INT,  ii,  kRepeatWithNewDistMsg, MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD); 

				/* assign leftovers */
				{
					FILE *hf=NULL;
					FILE *sf=NULL;
					FILE *of=NULL;
					unsigned int *rslt = NULL;
					CLOSESTSEQ *closestSeq = NULL;
					unsigned int starti,lasti;
					unsigned int datachunk = 1+(nbSeqForClustering/nproc);
					unsigned int step;
					unsigned int refseqassigned = 0;


					if (nbDistToTest == 1) // special case if we did not recurse, skip greedy step of non-ref sequences.
					{
						gTestDist = 0;
						MPI_Bcast (&gTestDist, 1, MPI_INT, 0, MPI_COMM_WORLD);
					}
					else
					{
						rslt = malloc(loaded*nbDistToTest*sizeof(unsigned int));
						if (rslt)
						{
							for (step = 0; step < nbDistToTest; step++)
							{
								if (LoadClusterIndices(rslt,step,nbSeqForClustering,loaded,distList[step].dist,tfn,sameFollowingCnt,cntcutoff))
								{	
									fprintf(stderr,"ERROR: Cannot load clustering results for step %u\n",step);
									break;
								}
							}
							if (step < nbDistToTest) // error.
							{
								free(rslt);
								rslt = NULL;
							}
						}
						if (rslt)
						{
							FlagSequencesToReassign(clusterid,&rslt[loaded*(nbDistToTest-1)],nbSeqForClustering);
						}
						else
							fprintf(stderr,"ERROR: Cannot allocate memory for clustering history.\n");

						gTestDist = distList[nbDistToTest-1].squaredDist;
						gettimeofday(&te, NULL);

						if (verbose > 0)
							fprintf(stderr,"LOG:%d sec; Reassigning assigned sequences to closest cluster high-density regions\n",((int)te.tv_sec-(int)ts.tv_sec));	

						MPI_Bcast (&gTestDist, 1, MPI_INT, 0, MPI_COMM_WORLD);
						MPI_Bcast (&trimmedclustercnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
						MPI_Bcast (&clusterid[0], nbSeqForClustering, MPI_INT, 0, MPI_COMM_WORLD);
						lasti = 0 + datachunk;
						if (lasti > nbSeqForClustering)
							lasti = nbSeqForClustering;
						DistributeUnassignedToClosestCluster(fastadata,clusterid,nbSeqForClustering,trimmedclustercnt,0,lasti);
						gettimeofday(&te, NULL);

						if (verbose > 1)
							fprintf(stderr,"LOG:%d sec; Collecting Results\n",((int)te.tv_sec-(int)ts.tv_sec));	

						fflush(stderr);
						for (ii = 1; ii<nproc; ii++)
						{
							starti = ii*datachunk;
							if (starti < nbSeqForClustering)
							{
								lasti = starti + datachunk;
								if (lasti > nbSeqForClustering)
									lasti = nbSeqForClustering;

								MPI_Recv(&clusterid[starti], (lasti-starti), MPI_INT,  ii, kClusterMsg1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
							}
						}

						if (rslt)
						{
							unsigned int cnt;
							unsigned int adjustpass;
							unsigned int adjusted;
							unsigned int prevcnt = 0;
							unsigned int ambiguous = 0;
							for (adjustpass = 0; adjustpass < 5; adjustpass++)
							{
								cnt = FlagSequencesToReassign(clusterid,&rslt[loaded*(nbDistToTest-1)],nbSeqForClustering);
								if (cnt == 0)
									break;
								ambiguous = DistributeUnassignedToClosestCluster(fastadata,clusterid,nbSeqForClustering,trimmedclustercnt,0,nbSeqForClustering);
								if ((ambiguous == 0) || (cnt == prevcnt))
									break;
								prevcnt = cnt;
							}

							for (adjustpass = 0; adjustpass < 5; adjustpass++)
							{
								adjusted = CheckBranchHistory(fastadata,trimmedclustercnt,loaded,clusterid,rslt,nbDistToTest,distList,0);
								if (adjusted == 0)
										break;
							}
							if (adjusted > 0)
							{
								adjusted = CheckBranchHistory(fastadata,trimmedclustercnt,loaded,clusterid,rslt,nbDistToTest,distList,1);
							}
							ambiguous += CheckRsltArray(clusterid,&rslt[loaded*(nbDistToTest-1)],loaded,trimmedclustercnt);
							if (ambiguous)
							{
								gettimeofday(&te, NULL);
								if (verbose > 0)
								{
									fprintf(stderr,"LOG:%d sec; Note: some sequences previously attributed to a cluster were reattributed to the cluster zero (noise)\n",((int)te.tv_sec-(int)ts.tv_sec)/*,ambiguousCnt*/);
									fprintf(stderr,"LOG:         The final partition will have slightly less sequences than the penultimate partition.\n");
									fflush(stderr);
								}
							}

						}

					}
					if (loadedrefseqcnt > 0)
					{
						
						MPI_Bcast (&clusterid[0], nbSeqForClustering, MPI_INT, 0, MPI_COMM_WORLD);

						/* evenly split refseqs onto computing nodes */
						datachunk = 1+(loadedrefseqcnt/nproc);
						starti = 0;
						lasti = starti + datachunk;
						if (lasti > loadedrefseqcnt)
							lasti = loadedrefseqcnt;

						closestSeq = malloc(loadedrefseqcnt*sizeof(CLOSESTSEQ));
						if (closestSeq)
						{
							memset(closestSeq,0xff,loadedrefseqcnt*sizeof(CLOSESTSEQ));	// initialize closest seq to -1.
							gettimeofday(&te, NULL);
							if (verbose > 0)
								fprintf(stderr,"LOG:%d sec; Processing Reference Sequences with d=%.3f\n",((int)te.tv_sec-(int)ts.tv_sec),distList[nbDistToTest-1].dist);

							DistributeRefSeqToClosestCluster(fastadata,clusterid,nbSeqForClustering,nbSeqForClustering+starti,nbSeqForClustering+lasti,closestSeq);
						}
						else
							fprintf(stderr,"LOG:Warning: proc %d has not enough memory to assign reference sequences to closest cluster.\n",idproc);

						/* collect clusterid assignments made by each node */
						for (ii = 1; ii<nproc; ii++)
						{
							starti = ii*datachunk;
							if (starti < loadedrefseqcnt)
							{
								lasti = starti + datachunk;
								if (lasti > loadedrefseqcnt)
									lasti = loadedrefseqcnt;

								MPI_Recv(&clusterid[nbSeqForClustering+starti], (lasti-starti), MPI_INT,  ii, kClusterMsg1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
							}
						}
						/* collect closestseq assignments made by each node */
						for (ii = 1; ii<nproc; ii++)
						{
							starti = ii*datachunk;
							if (starti < loadedrefseqcnt)
							{
								lasti = starti + datachunk;
								if (lasti > loadedrefseqcnt)
									lasti = loadedrefseqcnt;

								MPI_Recv(&closestSeq[starti], 2*(lasti-starti), MPI_INT,  ii, kClosestMsg, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
							}
						}

					}
					free(fastadata);  // free as much memory as possible
					fastadata = NULL;
					MPI_Barrier(MPI_COMM_WORLD); 
					gettimeofday(&te, NULL);
					fprintf(stderr,"LOG:%d sec; Writing Clustering Results\n",((int)te.tv_sec-(int)ts.tv_sec));

					sprintf(tmpfn,"%s.headers.txt",tfn);
					hf = fopen(tmpfn,"r");
					if (OutputLevel >= 4)
					{
						sprintf(tmpfn,"%s.seq.txt",tfn);
						sf = fopen(tmpfn,"r");
					}
					if (ofn[0])
						of = fopen(ofn,"w");
					else
						of = stdout;
					fflush(stderr);
					if (rslt && ((OutputLevel == 3) || (OutputLevel == 5)))
						unassigned = nbSeq-WriteFullOutput(hf,sf,fastaname,clusterid,sameFollowingCnt,closestSeq,rowcnt,loaded,nbSeqForClustering,trimmedclustercnt,of,rslt,nbDistToTest,distList,&refseqassigned);
					else
						unassigned = nbSeq-WriteOutput(hf,sf,fastaname,clusterid,sameFollowingCnt,rowcnt,loaded,nbSeqForClustering,trimmedclustercnt,of,OutputLevel,&refseqassigned);

					if (ofn[0] && of)
						fclose(of);
					if (rslt)
						free(rslt);
					if (closestSeq)
						free(closestSeq);

					for (step = 0; step < nbDistToTest; step++)
					{
						sprintf(tmpfn,"%s-%.6f",tfn,distList[step].dist);
						remove(tmpfn);
					}
					if (hf)
					{
						fclose(hf);
						sprintf(tmpfn,"%s.headers.txt",tfn);
						remove(tmpfn);
					}
					if (sf)
					{
						fclose(sf);
						sprintf(tmpfn,"%s.seq.txt",tfn);
						remove(tmpfn);
					}

					if (verbose > 0)
						fprintf(stderr,"LOG: %12d Final Clusters\n",trimmedclustercnt);

					fprintf(stderr,"LOG: %12u Assigned        (%5.1f %%)\n",(nbSeq-unassigned),100.0*(nbSeq-unassigned)/nbSeq);
					fprintf(stderr,"LOG: %12u Unassigned      (%5.1f %%)\n",unassigned,100.0*unassigned/nbSeq);	
					if (refseqcnt > 0)
					{
						fprintf(stderr,"LOG: %12u Ref. Assigned   (%5.1f %%)\n",refseqassigned,100.0*refseqassigned/refseqcnt);	
						fprintf(stderr,"LOG: %12u Ref. Unassigned (%5.1f %%)\n",(refseqcnt-refseqassigned),100.0*(refseqcnt-refseqassigned)/refseqcnt);	
					}
					gettimeofday(&te, NULL);
					fprintf(stderr,"DBV:%d,%u,%.3f,%u,%u,%u,%.3f,%.3f,%d,%d\n",((int)te.tv_sec-(int)ts.tv_sec),distList[nbDistToTest-1].squaredDist,distList[nbDistToTest-1].dist,nbSeq,(nbSeq-unassigned),(unassigned),100.0*(nbSeq-unassigned)/nbSeq,100.0*unassigned/nbSeq,trimmedclustercnt,trimmedclustercnt);
					fprintf(stderr,"LOG:**************************************************************\n");
					fprintf(stderr,"LOG:%d sec; Done.\n",((int)te.tv_sec-(int)ts.tv_sec));
					fflush(stderr);
				}

			}
			else
			{
				gTestDist = distList[currentTestingDistanceIndex].squaredDist;
				/* if gDist increases, do as if computing node #1 had discovered valid clusters so keep resuls already valid and reduce the number of merging events */
				for(ii = 0;  ii< nbSeqForClustering; ii++) // (ref sequences not assigned at this point, still zero, so ignore them)
				{
					if ((clusterid[ii] > 0))
					{
						clusterid[ii] += 1*kStartLocalCluster;
					}
					else
						clusterid[ii] = 0;					
				}
				for (ii = 1; ii<nproc; ii++)
					MPI_Send(&gTestDist, 1, MPI_INT,  ii,  kRepeatWithNewDistMsg, MPI_COMM_WORLD);

				MPI_Send(&initialClusterCnt, 1, MPI_INT,  1,  kInitialClusterCntMsg, MPI_COMM_WORLD); // send starting clustercount to slave node #1
				MPI_Barrier(MPI_COMM_WORLD); 
				stats[0] = stats[1];
				passcnt++;

				if (verbose > 0)
					fprintf(stderr,"LOG:**************************************************************\n");
				goto repeatWithNewDist;
			}
		}
abort:
		if (distList)
			free(distList);
		if (fastadata)
			free(fastadata);
		if (fastaname)
			free(fastaname);
		if (clusterid)
			free(clusterid);
		if (sameFollowingCnt)
			free(sameFollowingCnt);
		if (clusterhistory)
			free(clusterhistory);

		MPI_Barrier(MPI_COMM_WORLD); 
		MPI_Finalize();
		return(0);

} /* main */
/* ------------------------------------------------------------------------------------ */
