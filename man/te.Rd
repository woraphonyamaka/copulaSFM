\name{te}
\alias{te}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Technical efficiency measure.
}
\description{
Computing and plotting the technical efficiency.
}
\usage{
TE(theta,Y,X,family=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{theta}{ The estimated parameters form the model}
 \item{Y}{ Vector of dependent variable}
\item{X}{ Matrix of independent variable}
\item{family}{ Copula function eg. Gaussain=1, Student-t=2 (see, Vinecopula package)}
}
\details{
Computing and plotting the technical efficiency.
}
\value{

item{Output }{Technical efficiency series}
item{plot }{Plot of Technical efficiency}

}
\references{
Wiboonpongse, A., Liu, J., Sriboonchitta, S., & Denoeux, T. (2015). Modeling dependence between error components of the stochastic frontier model using copula: application to intercrop coffee production in Northern Thailand. International Journal of Approximate Reasoning, 65, 34-44.
}
\author{
Woraphon Yamaka
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Required packages
library(truncnorm)
library(mvtnorm)
library("VineCopula")
library("frontier")

#example simulation data
data=sfa.simu(nob=200, alpha=c(1,2,0.5),sigV=1,sigU=0.5,family=1,rho=0.5)

# Select familty  copula upper and lower bouubd ( look at CDVine package)
# family=1   # 1 is Gaussian, 2 is Student-t, 3 is Clayton and so on....

#Gaussian (-.99, .99)
#Student t (-.99, .99)
#Clayton (0.1, Inf)

model=copSFM(Y=data$Y,X=data$X,family=1,RHO=0.5,LB=-0.99,UB=0.99)

#EX: Plot the technical efficiency
te1=TE1(model$result[,1],Y=data$Y,X=data$X)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
