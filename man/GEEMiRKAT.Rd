\name{GEEMiRKAT}
\alias{GEEMiRKAT}
\title{GEEMiRKAT}
\usage{
GEEMiRKAT(y, id, covs = NULL, Ks, model, n.perm = 5000)
}
\description{
An adaptive microbiome regression-based kernel association test using generalised estimating equations for longitudinal microbiome studies.
}

\arguments{
  \item{y}{
Response variable (i.e., host phenotype of interest). Exponential family of distributions (e.g., Gaussian, Binomial, Poisson) outcomes.
}
  \item{id}{
A vector for identifying the sequence of subjects/clusters of longitudinal data.

}
  \item{covs}{
A data.frame (or matrix/vector) for covariate (e.g., age). Default is covs = NULL.
}
  \item{Ks}{
Kernel matrix based on diverse distances.
}
  \item{model}{
"gaussian" for Gaussian outcomes, "binomial" for Binomial outcomes, "poisson" for Poisson outcomes.
}
  \item{n.perm}{
A number of permutations. Default is n.perm=5000.
}
}

\value{
$Ind.pvs: The p-values for the individual tests based on diverse distances.

$GEEMiRKAT.pvs: The p-values for the adaptive test and omnibus tests.
}

\examples{
# Import requisite R packages
require(phyloseq)
require(vegan)
require(GUniFrac)
require(PGEE)
require(MiRKAT)
require(ACAT)
require(GEEMiRKAT)

# Import example microbiome data
data(CD)
y <- CD@sam_data$CD
covs <- as.factor(CD@sam_data$smoker)
covs <- as.data.frame(covs)
id <- CD@sam_data$id
otu.tab <- CD@otu_table
tree <- CD@phy_tree

# Calculate the similarity matrix between samples
unifracs <- GUniFrac(otu.tab, tree, alpha = c(0.5, 1))$unifracs
D.weighted = unifracs[, , "d_1"]
D.weighted50 = unifracs[, , paste("d_", 0.5, sep = "")]
D.unweighted = unifracs[, , "d_UW"]
D.BC = as.matrix(vegdist(otu.tab , method="bray"))

K.weighted = D2K(D.weighted)
K.weighted50 = D2K(D.weighted50)
K.unweighted = D2K(D.unweighted)
K.BC = D2K(D.BC)
Ks <- list(BC = K.BC, U.UniFrac = K.unweighted, G.UniFrac.50 = K.weighted50, W.UniFrac = K.weighted)
  
# Fit GEEMiHC
set.seed(123)
out <- GEEMiRKAT(y, covs, id = id, Ks = Ks, model = "binomial")
out
}
