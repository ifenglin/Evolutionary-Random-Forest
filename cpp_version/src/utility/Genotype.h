// 
//  Change any of these parameters to match your needs 
//
struct genotype
{
	std::vector<uint> gene;
	//gene 0: max trials per node
	//gene 1: minimal node size
	//gene 2: minimal leaf size, in fraction of minimal node size in %.
	//gene 3: fraction of observations to sample, in %%.
	//gene 4: tree seed
	//gene 5 - n_vars: select weights of each feature
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