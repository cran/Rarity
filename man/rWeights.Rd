\name{rWeights}
\alias{rWeights}
\alias{rWeights-method}
\alias{rWeights,data.frame-method}
\alias{rWeights,matrix-method}
\alias{rWeights,vector-method}
\title{
Calculation of rarity weights at a single or at multiple scales
}
\description{
Calculate rarity weights for a single scale or for multiple scales on the basis of the selected weighting function(s).
}
\usage{
rWeights(occData, Qmax = max(occData), Qmin = min(occData), 
         wMethods = "W", rCutoff = "Gaston", normalized = T, 
         assemblages, extended = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{occData}{
vector, matrix or data.frame. Occurrence data for a single scale (vector) or several scales (matrix or data.frame). 
}
  \item{Qmax}{
integer. Maximum occurrence (see details)
}
  \item{Qmin}{
integer. Minimum occurrence (see details)
}
  \item{wMethods}{
\code{W}, \code{invQ} or \code{oldW}. Chosen function to weight species occurrence (see details)
}
  \item{rCutoff}{
a decimal or a vector of values between 0 and 1, or \code{Gaston} or \code{Leroy}. Indicates the rarity cutoff(s) or the method to use to calculate the rarity cutoff(s) (see details)
}
  \item{normalized}{
TRUE or FALSE. If TRUE, then weights are normalized between 0 and 1.
}
  \item{assemblages}{
matrix or data.frame. Set of assemblages of species to calculate the rarity cutoff point(s) with the \code{Leroy}  method (optional)
}
  \item{extended}{
TRUE or FALSE. Useful in case of multiple scales only. If TRUE, then weights will be given for every input scale in addition to multiscale weights. If FALSE, then only multiscale weights will be provided. 
}
}
\details{
To calculate single-scale weights, simply provide a vector with species occurrences.

To calculate multiscale rarity weights, provide either a matrix or a data.frame where species are in rows, and each column provides occurrence for a particular scale.

The minimum and maximum weights can be set manually, or automatically calculated with the default parameters. Defaults parameters : if \code{occData} is a vector, \code{Qmin = min(Q)} and \code{Qmax = max(Q)}. If \code{occData} is a matrix or a data.frame, \code{Qmin = apply(occData, 2, min)} and \code{Qmax = apply(occData, 2, max)}

NA values are not (yet) allowed (remove them from the dataset before running the function).

Three weighting methods are available (more will become available later):

\bold{1.} \code{W}: This is the method described in Leroy et al. In Press. We recommend using this method for both single and multiscale weight calculations.

\bold{2.} \code{invQ}: This is the inverse of occurrence.

\bold{3.} \code{oldW}: This is the original method described in Leroy et al. 2012. As this method was improved in Leroy et al. In press, we recommend to rather use \code{W}


For methods \code{W} and \code{oldW}, a rarity cutoff point is required. The rarity cutoff point can either be entered manually (a single value for a single scale, a vector of values for multiple scales), or the methods of \code{Gaston} or \code{Leroy} can be used (see references).

 - \code{Gaston} method: the rarity cutoff point is the first quartile of species occurrences, i.e. rare species are the 25 percent species with the lowest occurrence.
 
 - \code{Leroy} method: the rarity cutoff point is the occurrence at which the average proportion of rare species in local assemblages is 25 percent. This method requires \code{assemblages} to calculate the average proportion of rare species in  assemblages.
}
\value{
A data.frame containing the results : species occurrences, rarity statuses and rarity weights.

-  If \code{occData} is a vector (single scale weights): A data.frame with 3 columns : \code{Q} (species occurrence), \code{R} (species rarity status), \code{W}, (species rarity weight)

-	If \code{occData} is matrix or a data.frame (multiscale rarity weights): A data.frame with \emph{n} columns \code{Q} (species occurrences), \emph{n} columns \code{R} (species rarity statuses), one (if \code{extended = F}) or \emph{n} + 1 (if \code{extended = T}) columns \code{W} (species rarity weights) where \emph{n} is the number of scales (number of columns of \code{occData})
}
\references{
Leroy B., Petillon J., Gallon R., Canard A., & Ysnel F. (2012) Improving occurrence-based rarity metrics in conservation studies by including multiple rarity cut-off points. \emph{Insect Conservation and Diversity}, \bold{5}, 159-168.

Leroy B., Canard A., & Ysnel F. In Press. Integrating multiple scales in rarity assessments of invertebrate taxa. \emph{Diversity and Distributions}

}
\author{
Boris Leroy leroy.boris@gmail.com
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Irr}}
}
\examples{
# 1. Single scale rarity weights
# Example with occurrences of species ranging from 1 to 100
# where 25percent of the species have occurrences lower than 7
occurrences <- c(sample(seq(1, 6), 25, replace = TRUE), 
                 sample(seq(7, 100), 75, replace = TRUE))
names(occurrences) <- paste("sp", 1:100)
occurrences
# Calculation of rarity weights
rWeights(occurrences, rCutoff = "Gaston")
rWeights(occurrences, rCutoff = 0.1)
rWeights(occurrences, wMethods = "invQ")
rWeights(occurrences, wMethods = c("W", "invQ"))

# Calculation of rarity weights with the method of Leroy
assemblages.matrix <- cbind(assemblage.1 = sample(c(0, 1), 100, replace = TRUE),
                            assemblage.2 = sample(c(0, 1), 100, replace = TRUE),
                            assemblage.3 = sample(c(0, 1), 100, replace = TRUE),
                            assemblage.4 = sample(c(0, 1), 100, replace = TRUE),
                            assemblage.5 = sample(c(0, 1), 100, replace = TRUE))
rownames(assemblages.matrix) <- paste("sp", 1:100) # Rownames of assemblages.matrix must correspond to rownames in occurrences
rWeights(occurrences, wMethods = "W", rCutoff = "Leroy", assemblages = assemblages.matrix)

# 2. Multiscale rarity weights
# Example with occurrences of species ranging from 1 to 100 in a first dataset
# where 25percent of the species have occurrences lower than 7
# and ranging from 3 to 80 in a second dataset
# where 25 percent of the species have occurrences lower than 20
scale1 <- c(sample(seq(1, 6), 25, replace = TRUE), sample(seq(7, 100), 75, replace = TRUE))
scale2 <- c(sample(seq(3, 20), 25, replace = TRUE), sample(seq(21, 80), 75, replace = TRUE))
dataset <- cbind(S1 = scale1, 
                 S2 = scale2[sample(1:100, 100, replace = FALSE)]) # random attribution of occurrences to species
rownames(dataset) <- paste("sp", 1:100)
head(dataset)

rWeights(dataset, wMethods = "W", rCutoff = "Gaston")
rWeights(dataset, wMethods = "W", rCutoff = "Gaston", extended = TRUE)
rWeights(dataset, wMethods = c("W", "invQ"), rCutoff = "Gaston", extended = TRUE)
rWeights(dataset, wMethods = c("W", "invQ"), rCutoff = "Leroy", assemblages = assemblages.matrix, extended = TRUE) # Provided that you have created "assemblages.matrix" as above
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
