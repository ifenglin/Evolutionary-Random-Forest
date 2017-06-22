// 
//  Change any of these parameters to match your needs 
//
# define NVARS 3
# define PXOVER 0.8
# define PMUTATION 0.15
struct genotype
{
	uint gene[NVARS];
	double fitness;
	double upper[NVARS];
	double lower[NVARS];
	double rfitness;
	double cfitness;
};