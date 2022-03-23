# Copyright (C) 2018  Sebastian Sosa, Ivan Puga-Gonzalez, Hu Feng He, Xiaohua Xie, Cédric Sueur
#
# This file is part of Animal Network Toolkit Software (ANTs).
#
# ANT is a free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# ANT is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

#' @title Data stream permutation for association data
#' @description Pre-network permutation on association data for gambit of the group data collection protocol. The data frame must have a column named 'ID'.
#' @param obs a data frame of gambit of the group observations. The data frame must have a column named 'ID'.
#' @param scan  an integer indicating the column of scans of individual associations in obs.
#' @param ctrlf A confounding factor by which to control group associationsin obs.
#' @param nperm number of permutations to perform.
#' @param progress a boolean indicating if the permutation process must be visible.
#' @param index Which type of index of associations to calculate:
#' \itemize{
#' \item 'sri' for Simple ratio index: \eqn{x \div x+yAB+yA+yB}
#' \item 'hw' for Half-weight index: \eqn{x/x+yAB+1/2(yA+yB)}
#' \item 'sr' for Square root index:\eqn{x/sqr((x+yAB+yA)(x+yAB+yB))}
#' }
#' @param measure a character indicating the social network measure to compute (Only those available in ANTs)
#' @param test a character indicating the test to realize to account for the social network measure
#' @param df a data frame of individual characteristics in which store permutations.
#' @param dfid an integer or a string indicating the column with individual ids in argument \emph{df}.
#' @param rf an integer (column id) or a string (column name) indicating the column holding the factor grouping multiple networks in argument \emph{df}.
#' @param ... Additional arguments related to the social network measure to compute (argument measure).
#' @return list of square association index matrices. The first element of the list is the non-permuted association index matrix.
#' @details Data stream permutation is a pre-network permutation approach. It is used on association data based on the gambit of the group.
#' @author Sebastian Sosa, Ivan Puga-Gonzalez.
#' @references Whitehead, H. A. L. (1997). Analysing animal social structure. Animal behaviour, 53(5), 1053-1067.
#' @references Farine, D. R. (2017). A guide to null models for animal social network analysis. Methods in Ecology and Evolution.
#' @references Sosa, S. (2018). Social Network Analysis, \emph{in}: Encyclopedia of Animal Cognition and Behavior. Springer.
#' @examples
#' head(sim.grp)
#' t=perm.double.grp(sim.grp, 'location', 'time', 10, measure = "met.strength")

perm.double.grp <- function(obs, scan, ctrlf = NULL, nperm, progress = TRUE, 
                            index = 'sri', measure, test = "median", df = NULL, dfid = NULL, rf, ...) {
  ## check for the presence of control arguments
  test2 <- check.df(obs)
  if (is.null(test2)) {
    "Argument df is not a data frame or a list of data frames."
  }
  ## argument df is a single dataframe, perform permutations
  if(test2 == "df ok"){
    result <- perm.double.grp.single(obs = obs, scan = scan, ctrlf = ctrlf, nperm = nperm, progress = progress, index = index, measure = measure, test = test, df = df, dfid = dfid, ...)
    attr(result, "ANT") <- "ANT double permutation GoG single matrix"
    attr(result, "scan") <- scan
    attr(result, "ctrlf") <- ctrlf
    attr(result, "method") <- index
    attr(result, "measure") <- measure
    attr(result, "test") <- test
    
    if(is.null(df)){
      result = list(result, lapply(1:nperm, function(i,result){
        vec_sample_all(result)
      }, result = result))
    }else{
      result = perm.net.nl(result, labels =measure, rf = NULL, nperm = nperm, progress = progress)
    }
    return(result)
  }
  ## argument df is a list of dataframes, perform permutations in each element of the list
  if (test2 == "df list ok") {
    result <- lapply(seq(obs), function(i,obs, scan, ctrlf, nperm, progress, index, measure, test, df, dfid, ...){
      perm.double.grp.single(obs[[i]], scan, ctrlf, nperm, progress, index, measure, test, df, dfid, ...)
    }, obs = obs, scan = scan, ctrlf = ctrlf, nperm = nperm, progress = progress, index = index, measure = measure, test = test, df = df, dfid = dfid, ...)
    attr(result, "ANT") <- "ANT double permutation GoG multiple matrices"
    attr(result, "scan") <- scan
    attr(result, "ctrlf") <- ctrlf
    attr(result, "method") <- index
    attr(result, "measure") <- measure
    attr(result, "test") <- test
    
    if(is.null(df)){
      result = lapply(1:(nperm+1), function(i,result){
        if(i == 1){result}
        lapply(result, function(x){vec_sample_all(x)})
      }, result = result)
    }else{
      result = perm.net.nl(result, labels = measure, rf = rf, nperm = nperm, progress = progress)
    }
    return(result)
  }
}

