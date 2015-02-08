\name{NormalizeRUVRand}
\alias{NormalizeRUVRand}
\title{RUV-random normalized data}
\description{Given suitable controls, this function may be used to obtain a normalized metabolomics data matrix for downstream statistical analysis}
\usage{ NormalizeRUVRand(Y,ctl,k=NULL,lambda=NULL,plotk=TRUE) }
\arguments{
  \item{Y}{A metabolomics data matrix with samples in rows and metabolites in columns}
   \item{ctl}{A logical vector indicating which metabolites should be used as negative controls}
    \item{k}{The number of factors of unwanted variation}
    \item{plotk}{A logical indicating whether a bargraph for the proportion of variance explained by the factors of unwanted variation needs to be plotted}
\item{lambda}{The regularization parameter which depends on k. If not entered, it will be estimated. See DeLivera et al, 2014 for details.}
}
\value{The result is an object of class \code{\link[MetNorm:normdata]{normdata}}.}

\references{
De Livera, A. M., Dias, D. A, De Souza, D., Rupasinghe, T., Pyke, J., Tull, D., Roessner, U., McConville, M., and Speed, T. P. (2012). Normalizing and integrating metabolomics data. Analytical chemistry, 84(24), 10768-76.

De Livera, A.M., Aho-Sysi, M., Jacob, L., Gagnon-Bartch, J., Castillo, S., Simpson, J.A., and Speed, T.P. (2014), Statistical methods for handling unwanted variation in metabolomics data
}
\author{
Alysha M De Livera and Laurent Jacob
}
\examples{
data(UV)
Y<-data.matrix(UV[,-c(1:3)]) 

##Empirical controls
IS<-Y[,which(colnames(Y)=="IS")] 
r<-numeric(dim(Y)[2])
for(j in 1:length(r)){
  r[j]<-cor(IS,Y[,j])
}
ctl<-logical(length(r))
ctl[which(r>round(quantile(r,0.7),2))]<-TRUE 

\dontrun{  
ruv<-NormalizeRUVRand(Y=Y,ctl=ctl,k=3) 
ruvclust<-NormalizeRUVRandClust(RUVRand=ruv,
                   maxIter=200, 
                   nUpdate=100,
                   lambdaUpdate=TRUE,
                   p=2)
ruvclustY<-ruvclust$newY
pairs(princomp(ruvclustY,cor=TRUE)$scores[,c(1:3)],
      pch=as.numeric(UV[,2]), col=UV[,3],
      main="RUV random for clustering")
}
}


