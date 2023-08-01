# Tasha M. Santiago

# DIABLO code using a reduced HMP2 dataset that includes matched metagenomic, metabolomic and metadata results.
# DIABLO (N-integration) is a framework that is part of the mixOmics package that extends Partial Least Square (PLS) for the integration of multiple data sets.
# The aim of DIABLO is to identify correlated (or co-expressed) variables measured on heterogeneous data sets.
# The biological question should be along the lines of trying to identify a highly correlated multi-omics signature discriminating known groups of samples.

# Install and run mixOmics.To install this package, start R (version "4.0") and enter: 
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

# BiocManager::install("mixOmics")

# Load libraries
library(mixOmics)
library(ggpubr)

# Read needed files
metagenomics <- read.csv("HMP2_metagenomics.csv", row.names = 1, header = TRUE)

metabolomics <- read.csv("HMP2_metabolomics.csv", row.names = 1, header = TRUE)

taxonomy <- read.csv("HMP2_taxonomy_list.csv", row.names = 1, header = TRUE)

dysbiosis <- read.csv("HMP2_dysbiosis_state.csv", row.names = 1, header = TRUE)

# Create a list with the files uploaded above
HMP2 <- list(metagenomics, metabolomics, taxonomy, dysbiosis)

# View(HMP2)

# Rename the variables or categories within the list
names(HMP2)[1]<- "metagenomics"
names(HMP2)[2]<- "metabolomics"
names(HMP2)[3]<- "taxonomy"
names(HMP2)[4]<- "dysbiosis"

# Load the data
  
X <- list(metagenomics = HMP2$metagenomics,
        metabolomics = HMP2$metabolomics)

Y <- HMP2$dysbiosis$diagnosis # Qualitative matrix and is internally recoded as a dummy block matrix that records the membership of each observation

# Set the number of features in each component. This helps in data visualization as opposed to plotting everything
list.keepX <- list(metagenomics = c(5, 5), metabolomics = c(5, 5))

# Typing list.keepX opens the list with the renamed variables or categories with 5 features each

# The PLS regression (now PLS-DA) is ran as if Y was a continuous matrix, unlike the original PLS method. Sparse (s)PLS-DA, a special case of sparse PLS, performs variable selection and classification in a one-step procedure. 

# Run the method. In this case we are running splsda
MyResult.diablo <- block.splsda(X, Y, keepX = list.keepX) 

# Plot the samples. This plot below shows all the labels. Might not be helpful when working with a lot of samples.
plotIndiv(MyResult.diablo) 

# Customize plots. The plot below includes ellipses and no sample IDs.
plotIndiv <- plotIndiv(MyResult.diablo,ind.names = FALSE, legend = TRUE, ellipse = TRUE, size.subtitle = rel(1), size.axis = rel(1), size.xlabel = rel(1.2),
                  legend.title = "Diagnosis", size.ylabel = rel(1.2))

ggexport(plotIndiv, filename = "plotIndiv.png", res = 400, width = 2500, height = 1500)

plotIndiv.2 <- plotIndiv(MyResult.diablo, ind.names = FALSE, legend = TRUE, ellipse = TRUE, size.subtitle = rel(1), size.axis = rel(1), size.xlabel = rel(1.2),
                       legend.title = "Diagnosis", size.ylabel = rel(1.2), rep.space = 'XY-variate', guide = "none")
ggexport(plotIndiv, filename = "plotIndiv.png", res = 400, width = 2500, height = 1500)

# Generate a global overview of the correlation structure between the metabolomic and metagenomic data at the component level:
plotDiablo(MyResult.diablo, ncomp = 1)

# Generate correlations between variables of different types using circos plots.
circos <- circosPlot(MyResult.diablo, cutoff = 0.5, size.labels = 1, size.variables = 0.8)

# Generate plotLoadings. The plotLoadings function visualises the loading weights of each selected variables on each component (default is comp = 1) and each data set. The color indicates the class in which the variable has the maximum level of expression (contrib = "max") or minimum (contrib ="min"), on average (method="mean") or using the median (method ="median"). 
plotLoadings <- plotLoadings(MyResult.diablo, comp = 1, method = "median", contrib = "max", size.name = 0.75, size.subtitle = 1.2, size.legend = 0.8)

# Generate a network to visualize the correlation between the different types of variables. Each colour represents a type of variable. A threshold can also be set using the argument cutoff.
network(MyResult.diablo, blocks = c(1,2),
        color.node = c('blue', 'red'), 
        cutoff = 0.5, save = 'png', name.save = 'DIABLOnetwork')

# Determine classification performance
set.seed(123)
MyPerf.diablo <- perf(MyResult.diablo, validation = 'Mfold', folds = 5, 
                      nrepeat = 10, 
                      dist = 'centroids.dist')

# Generate an AUC plot using the function auroc see (Rohart, Gautier, et al. 2017) for the interpretation of such output as the ROC and AUC criteria are not particularly insightful in relation to the performance evaluation of our methods, but can complement the statistical analysis.
Myauc.diablo <- auroc(MyResult.diablo, roc.component = "metabolomics", roc.comp = 2)

# Generate coordinate plot of both metagenomics and metabolomics datasets
plotVar(MyResult.diablo,cutoff = 0.5, style = "ggplot2", cex = c(3,3), overlap = sT) 
