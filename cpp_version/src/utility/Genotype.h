// 
//  Change any of these parameters to match your needs 
//
# define PXOVER 0.8
# define PMUTATION 0.15
struct genotype
{
	std::vector<uint> gene;
	//gene 0: max trials per node
	//gene 1: minimal node size
	//gene 2: fraction of observations to sample, in one thousandth.
	//gene 3 - n_vars: select weight of each feature
	double fitness;
	std::vector<double> upper;
	std::vector<double> lower;
	genotype(const size_t n_vars) {
		gene = std::vector<uint>(n_vars);
		upper = std::vector<double>(n_vars);
		lower = std::vector<double>(n_vars);
		fitness = 0;
	}
};