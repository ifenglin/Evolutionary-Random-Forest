// 
//  Change any of these parameters to match your needs 
//
# define NVARS 15
# define NFEATURES 3
# define PXOVER 0.8
# define PMUTATION 0.15
struct genotype
{
	uint gene[NVARS];
	//gene 0: max trials per node
	//gene 1: minimal node size
	//gene 2: fraction of observations to sample, in one thousandth.
	//gene 3 - 14: select weight of each feature
	double fitness;
	double upper[NVARS];
	double lower[NVARS];
	//double rfitness;
	//double cfitness;
};