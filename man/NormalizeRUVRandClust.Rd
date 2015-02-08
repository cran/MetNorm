\name{NormalizeRUVRandClust}
\alias{NormalizeRUVRandClust}
\title{RUV-random method for clustering}
\description{Given suitable controls and user input, this function may be used to obtain a normalized metabolomics data matrix suitable for clustering}
\usage{NormalizeRUVRandClust(RUVRand,
                              maxIter, 
                              nUpdate=maxIter+1, 
                              lambdaUpdate=TRUE, 
                              p=p,...)
}

\arguments{
  \item{RUVRand}{Output from \code{\link[MetNorm:NormalizeRUVRand]{NormalizeRUVRand}}}
  \item{maxIter}{Maximum number of iterations}
  \item{nUpdate}{Update the unwanted variation component every nUpdate iterations}
  \item{lambdaUpdate}{A logical indicating whether the regularization parameter needs to be updated}
  \item{p}{The number of clusters to be used in the k-means clustering}
  \item{...}{Other arguments for \code{\link[stats:kmeans]{kmeans}}}
  }
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

