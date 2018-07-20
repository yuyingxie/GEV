\name{predict.dsda}
\alias{predict.dsda}
\title{Prediction for direct sparse discriminant analysis}
\description{
Predict the class labels by direct sparse discriminant analysis.
}
\usage{
predict.dsda(obj,x.matrix)
}
\arguments{
\item{obj}{an object returned by dsda or dsda.path}
\item{x.matrix}{an n by p matrix containing the predictors.}
}

\references{
Mai, Q., Zou, H. and Yuan, M. (2013). A direct approach to sparse discriminant analysis in ultra-high dimensions. Biometrika, 99, 29-42.
}
\seealso{
\code{dsda}
\code{dsda.path}
}
\value{
\item{pred}{the the predicted class labels.}
}
\examples{

}