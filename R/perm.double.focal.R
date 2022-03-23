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

#' @title  Data stream permutation for focal sampling data .
#' @details  Pre-network permutation for focal sampling data, and for symmetric behaviour only.
#' @param obs a data frame of focal observations.
#' @param ego an integer indicating the column of the focal id for the obs.
#' @param alters an integer indicating the column of focal's alters in obs.
#' @param focal a numeric vector indicating the focal number in obs.
#' @param nperm an integer indicating the number of permutations to performed.
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
#' @description Warning, the original function (Farine 2017) uses a control factor, the number of focals and the ids of the focals.
#' @references Farine, D. R. (2017). A guide to null models for animal social network analysis. Methods in ecology and evolution, 8(10), 1309-1320.
#' @references Sosa, S. (2018). Social Network Analysis, \emph{in}: Encyclopedia of Animal Cognition and Behavior. Springer.
#' @examples
#' # Single network without data frame---------------------
#' head(sim.focal.undirected)
#' t=perm.double.focal(obs = sim.focal.undirected, ego = 3, alters = 4, 
#' focal = 1, nperm = 10, progress = FALSE, measure = "met.strength")

#' # Multiple networks with data frames---------------------
#' d1 = data.frame("id" = names(t[[1]]), "period" = 1)
#' d2 = data.frame("id" = names(t[[1]]), "period" = 2)
#' t = list(d1, d2)
#' obs = list(sim.focal.undirected, sim.focal.undirected)
#' t =perm.double.focal(obs = obs, ego = 3, alters = 4, focal = 1, nperm = 10, 
#' measure = "met.strength",  df = t, dfid = "id", rf = "period")

perm.double.focal <- function(obs, ego, alters, focal, nperm, progress = FALSE, 
                              index = 'sri', measure, test = "median", df = NULL, dfid = NULL, rf, ...) {
  ## check for the presence of control arguments
  if (is.null(focal)) {
    stop("Argument focal cannot be empty")
  }
  test2 <- check.df(obs)
  if (is.null(test2)) {
    "Argument df is not a data frame or a list of data frames."
  }
  ## argument df is a single dataframe, perform permutations
  if(test2 == "df ok"){
    result <- perm.double.focal.single(obs, ego, alters, focal, nperm, progress, index, measure, test, df, dfid, ...)
    attr(result, "ANT") <- "ANT double permutation focal sampling single matrix"
    attr(result, "ego") <- ego
    attr(result, "focal") <- focal
    attr(result, "alters") <- alters
    attr(result, "index") <- index
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
    result <- lapply(seq(obs), function(i, obs, ego, alters, focal, nperm, progress, index, measure, test, df, dfid, ...){
      perm.double.focal.single(obs[[i]], ego, alters, focal, nperm, progress, index, measure, test, df[[i]], dfid, ...)
    }, obs = obs, ego = ego, alters = alters, focal = focal, nperm = nperm, progress = progress, index = index, measure = measure, test= test, df = df, dfid = dfid, ... = ...)
    attr(result, "ANT") <- "ANT double permutation focal sampling multiple matrices"
    attr(result, "ego") <- ego
    attr(result, "focal") <- focal
    attr(result, "alters") <- alters
    attr(result, "index") <- index
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
