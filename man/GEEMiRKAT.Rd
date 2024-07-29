\name{GEEMiRKAT}
\alias{GEEMiRKAT}
\title{GEEMiRKAT}
\usage{
GEEMiRKAT(y, id, covs = NULL, Ks, model, n.perm = 5000)
}
\description{
The microbiome regression-based kernel association test based on generalized estimating equations for longitudinal microbiome data.
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

$Omn.pvs: The p-values for the omnibus test.
}

\examples{
# Import requisite R packages
require(phyloseq)
require(GUniFrac)
require(vegan)
require(PGEE)
require(MiRKAT)
require(ACAT)
require(GEEMiHC)

# Import example microbiome data
otu.tab <- CD_longitudinal@otu_table
tree <- CD_longitudinal@phy_tree
y=CD_longitudinal@sam_data$label
covs <- data.frame(matrix(NA, length(y), 2))
covs[,1] <- as.numeric(sample_data(CD_longitudinal)$age)
covs[,2] <- as.factor(sample_data(CD_longitudinal)$smoker)
id <- sample_data(CD_longitudinal)$id

unifracs <- GUniFrac(otu.tab, tree, alpha = c(0, 0.25, 0.5, 0.75, 1))$unifracs
D.weighted = unifracs[, , "d_1"]
D.weighted00 = unifracs[, , paste("d_", 0, sep = "")]
D.weighted25 = unifracs[, , paste("d_", 0.25, sep = "")]
D.weighted50 = unifracs[, , paste("d_", 0.5, sep = "")]
D.weighted75 = unifracs[, , paste("d_", 0.75, sep = "")]
D.unweighted = unifracs[, , "d_UW"]
D.BC = as.matrix(vegdist(otu.tab , method="bray"))

K.weighted = D2K(D.weighted)
K.weighted00 = D2K(D.weighted00)
K.weighted25 = D2K(D.weighted25)
K.weighted50 = D2K(D.weighted50)
K.weighted75 = D2K(D.weighted75)
K.unweighted = D2K(D.unweighted)
K.BC = D2K(D.BC)
Ks <- list(BC = K.BC, U.UniFrac = K.unweighted, G.UniFrac.50 = K.weighted50, W.UniFrac = K.weighted)
  
# Fit GEEMiHC
set.seed(123)
out <- GEEMiRKAT(y, covs, id = id, Ks = Ks, model = "binomial")
out
}
