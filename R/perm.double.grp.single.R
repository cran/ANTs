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

#' @title Double permutation approach for gambit of the group
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
#' @param ... Additional arguments related to the social network measure to compute (argument measure).
#' @details Output need to be incorporated in a data frame and node label permutations with ANTs function perm.net.nl need to be performed before using any ANTs functions "stat.".
#' @return A numeric vector of individuals social measure corrected by double permutation approch (node label permutation can the be perfomed on this output)
#' @references Farine, D. R., & Carter, G. G. (2022). Permutation tests for hypothesis testing with animal social network data: Problems and potential solutions. Methods in Ecology and Evolution, 13, 144- 156. https://doi.org/10.1111/2041-210X.13741

perm.double.grp.single <- function(obs, scan, ctrlf, nperm, progress = TRUE, index = 'sri', measure, test = "median", df = NULL, dfid = NULL,...){
  # Data stream permutations 
  ds = perm.ds.grp(df = obs, scan = scan, ctrlf = ctrlf, nperm = nperm, index = index, progress = progress)
  
  # Compute social network measure
  m = do.call("rbind", lapply(ds, function(x, measure, ...){
    do.call(measure, list(M = x, ...))
  }, measure = measure, ...))
  
  mO = m[1,] 
  m = m[-1,] 
  
  # Double permutation
  result <- mO - apply(m, 2, function(x){(do.call(test, list(x)))})
  if(!is.null(df)){
    if(is.null(dfid)){
      warning("Argument dfid hasn't been declared. obs and df are considered to be ordered exactly in the same way.")
      df[,ncol(df)+1] = result
      colnames(df)[ncol(df)] = measure
      return(df)
    }else{
      col.id <- df.col.findId(df, dfid)
      df <- merge.met(vec = result, names = names(result), df = df, dfid = col.id, met = measure)
      return(df)
    }
  }else{
    return(result)
  }
}



