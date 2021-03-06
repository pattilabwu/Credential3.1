% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_knots.R
\name{findzknots}
\alias{findzknots}
\title{Find knots - sets of features with defined spacing and retention time.}
\usage{
findzknots(features, .z = 1, ppmwid = 5, rtwid = 1, cd = 13.00335 - 12)
}
\arguments{
\item{features}{data.frame Rows are features. Columns are "cc"}

\item{.z}{integer Charge state to assume}

\item{ppmwid}{numeric The maximum mass error in ppm.}

\item{rtwid}{numeric The maximum retention time difference in seconds.}

\item{cd}{numeric The mass spacing to search for (defaults to C13 - C12)}
}
\value{
list A list with values "cc_knot" and "knot".  \code{cc_knot} contains the features assignments to a knot. \code{knot} contains aggregate information about each knot.
}
\description{
Find knots - sets of features with defined spacing and retention time.
}
