## Growing Random Forest with Diverse Trees through Evolution for Semantic Segmentation
### A thesis for a Master of Science at Technical University of Berlin

Random Forest is a machine learning model that has been well researched and has gained established mathematical foundation. An understanding of the underlying mechanisms suggests approaches to improve the model further. By extending the flexibility of the model's technical architecture and designing appropriate criteria for evaluation of its performance, this study progresses from the adjustment of the parameters to a multi-objective optimization of an ensemble classifier that minimizes the correlations between its individual classifiers and maximizes their correctness and precision. Genetic algorithms are used to solve this optimization problem, and through the evolutionary progress, the random forest decision trees achieve higher correctness, greater precision and lower correlation with other trees in the population. Various evaluation techniques and evolutionary strategies are tested to develop the most effective learning model, to which I give the name Evolutionary Random Forest. Finally, the new model is used to solve the semantic segmentation problem with the data in The Cityscapes Dataset. On the contrary to expectations, no significant improvement has been found in the solutions. The new model is helpful on finding local optimums in the solution space, but the conventional model is a fairly effective tool to the problem in this study.

The program is based on Ranger, a fast implementation of Random Forest (see below.) The important introduced implementations are in the following files:

cpp_version\src\main.exe: read settings from forest_init_config.txt,  run a genetic algorithm (GA) that initilize genes randomly, evaluate/reproduce/select members of the population, and produce logs and reports.

cpp_version\src\Genotype: define the genotype of trees.

cpp_version\src\Forest\rf.cpp: load data, case weight, and prior distribution files, and initialize a forest.

cpp_version\src\Forest\Forest.cpp,
cpp_version\src\Forest\ForestClassification.cpp: initilize Evolutionary Random Forest by initializing trees with individual genes, and evaluate the trees and the forest with out-of-bag samples.

cpp_version\src\Forest\Tree.cpp
cpp_version\src\Forest\TreeClassification.cpp: initilize trees differently based on their genes.

## ranger: A Fast Implementation of Random Forests
Marvin N. Wright

### Introduction
ranger is a fast implementation of random forest (Breiman 2001) or recursive partitioning, particularly suited for high dimensional data. Classification, regression, probability estimation and survival forests are supported. Classification and regression forests are implemented as in the original Random Forest (Breiman 2001), survival forests as in Random Survival Forests (Ishwaran et al. 2008). For probability estimation forests see Malley et al. (2012). 

ranger is written in C++, but a version for R is available, too. We recommend to use the R version. It is easy to install and use and the results are readily available for further analysis. The R version is as fast as the standalone C++ version.

### Installation
#### R version
To install the ranger R package from CRAN, just run

```R
install.packages("ranger")
```

R version >= 3.1 is required. With recent R versions, multithreading on Windows platforms should just work. If you compile yourself, the new RTools toolchain is required.

To install the development version from GitHub using `devtools`, run

```R
devtools::install_github("imbs-hl/ranger")
```

#### Standalone C++ version
To install the C++ version of ranger in Linux or Mac OS X you will need a compiler supporting C++11 (i.e. gcc >= 4.7 or Clang >= 3.0) and Cmake. To build start a terminal from the ranger main directory and run the following commands

```bash
cd cpp_version
mkdir build
cd build
cmake ..
make
```

After compilation there should be an executable called "ranger" in the build directory. 

To run the C++ version in Microsoft Windows please cross compile or ask for a binary.

### Usage
#### R version
For usage of the R version see ?ranger in R. Most importantly, see the Examples section. As a first example you could try 

```R  
ranger(Species ~ ., data = iris)
```

#### Standalone C++ version
In the C++ version type 

```bash
ranger --help 
```

for a list of commands. First you need a training dataset in a file. This file should contain one header line with variable names and one line with variable values per sample. Variable names must not contain any whitespace, comma or semicolon. Values can be seperated by whitespace, comma or semicolon but can not be mixed in one file. A typical call of ranger would be for example

```bash
ranger --verbose --file data.dat --depvarname Species --treetype 1 --ntree 1000 --nthreads 4
```

If you find any bugs, or if you experience any crashes, please report to us. If you have any questions just ask, we won't bite. 

Please cite our paper if you use ranger.

### References
* Wright, M. N. & Ziegler, A. (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. Journal of Statistical Software 77:1-17. http://dx.doi.org/10.18637/jss.v077.i01.
* Schmid, M., Wright, M. N. & Ziegler, A. (2016). On the use of Harrell’s C for clinical risk prediction via random survival forests. Expert Systems with Applications 63:450-459. http://dx.doi.org/10.1016/j.eswa.2016.07.018.
* Wright, M. N., Dankowski, T. & Ziegler, A. (2017). Unbiased split variable selection for random survival forests using maximally selected rank statistics. Statistics in Medicine. http://dx.doi.org/10.1002/sim.7212.
* Breiman, L. (2001). Random forests. Machine learning 45:5-32.
* Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. The Annals of Applied Statistics 2:841-860.
* Malley, J. D., Kruppa, J., Dasgupta, A., Malley, K. G., & Ziegler, A. (2012). Probability machines: consistent probability estimation using nonparametric learning machines. Methods of Information in Medicine 51:74-81.
