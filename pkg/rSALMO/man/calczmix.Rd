\name{calczmix}
\alias{calczmix}
\title{Calculate Mixing Depth}
\usage{
  calczmix(temp, depths)
}
\arguments{
  \item{temp}{vector of temperatures (from top to bottom
  (?))}

  \item{depths}{vector of depths}
}
\value{
  list with id and depth of zmix
}
\description{
  Calculate mixing depth from temperature profile using a
  simple temperature threshold.
}
