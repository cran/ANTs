# Copyright (C) 2018  Sebastian Sosa, Ivan Puga-Gonzalez, Hu Feng He,Peng Zhang, Xiaohua Xie, Cédric Sueur
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
#' @title Extracts statistical measures of interest in Generalized Linear Mixed Models in parallel
#' @description Performs Generalized Linear Mixed Models tests in parallel
#' @param ant an output of ANT function \code{perm.net.nl} with random factor stated, or output of ANT 'met' categories functions in which output of ANT functions \code{perm.ds.focal}, \code{perm.ds.grp} or \code{perm.net.lk} where multiple matrices have been used.
#' @param formula two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a ~ operator and the terms, separated by + operators, on the right. Random-effects terms are distinguished by vertical bars (|) separating expressions for design matrices from grouping factors. Two vertical bars (||) can be used to specify multiple uncorrelated random effects for the same grouping variable. (Because of the way it is implemented, the ||-syntax works only for design matrices containing numeric (continuous) predictors; to fit models with independent categorical effects, see dummy or the lmer_alt function from the afex package.).
#' @param family a GLM family, see \code{\link{glm}} and \code{\link{family}}.
#' @param oda the original data frame of associations when argument ant is obtained with perm.ds.grp or perm.ds.focal ANT functions.
#' @param progress a boolean indicating the visualization of the permutation process.
#' @param ncores an integer indicating the number of jobs to create for parallelization.
#' @param ... Extra arguments for \code{lmer} or \code{glmer} function only.
#' @details GLMM with permutation data.
#' @return Returns a list of 3 elements :
#' \itemize{
#' \item An object of class \code{\link{merMod}} (more specifically, an object of subclass lmerMod or glmerMod), for which many methods are available (e.g. methods(class="merMod")).
#' \item A data frame if the estimates of the permuted models.
#' \item A vector of integers indicating the permutations that returned model errors or warnings (e.g. model convergence issues) and for which new permutations were done.
#' }
#' @seealso \code{\link{lmer}} or \code{\link{glmer}}
stat.glmm.parallel <- function(ant, formula, family, oda = NULL, progress = FALSE, ncores = NULL, ...){
    # Arguments ----------------
  argg <- c(as.list(environment()), list(...))
  if(rstudioapi::isAvailable()){}else{stop("Rstudio is not running")}
  
  if (is.null(attributes(ant)$ANT)) {
    stop("Argument ant must be an object returned by perm.ds.grp, per.ds.focal or per.ds.nl functions")
  }
  if (is.character(family)) {
    fam <- family
  }
  else {
    if (attributes(family)$class == "family") {
      fam <- family$family
    }
    else {
      stop("Argument family is not a character or a family function.")
    }
  }
  

  
  # Create a temporary diretory---------------------------------
  tempdir <- tempfile()
  dir.create(tempdir)
  
  # Split data ---------------------------------------------------
  tmp = split(ant,cut(seq_along(ant),ncores,labels = FALSE))
  
  # Run analysis on permutations ---------------------------------------------------
  results = job.info = NULL
  for (a in 1:length(tmp)) {
    # Writing spliting data -------------------------------
    t = tmp[[a]]# data will be in object t
    odf = tmp[[1]][[1]]
    attributes(t) = attributes(ant)
    tmpFile <- file.path(tempdir, paste("tmp", a, ".RData", sep = ""))
    
    # Writing scripts -------------------------------------
    tmpFile2 <- file.path(tempdir, paste("tmp", a, ".R", sep = ""))
    
      if(a == 1){# If a == 1 then first element of the list is the original data 
        save(t, file = tmpFile)
        tmpFile2 <- file.path(tempdir, paste("tmp", a, ".R", sep = ""))
        
        path = gsub("\\", "/", tmpFile, fixed = TRUE)
        cat(paste0("load(",'\'',path,'\'', ")", "\n"),file = tmpFile2, append = TRUE)
        cat(paste0("library(ANTs)", "\n"),file = tmpFile2, append = TRUE)
        cat(paste0("resultANTsJobs",a, " = stat.glmm("), file = tmpFile2, append = TRUE)

      }else{
        save(t,odf, file = tmpFile)
        tmpFile2 <- file.path(tempdir, paste("tmp", a, ".R", sep = ""))
        
        path = gsub("\\", "/", tmpFile, fixed = TRUE)
        cat(paste0("load(",'\'',path,'\'', ")", "\n"),file = tmpFile2, append = TRUE)
        cat(paste0("library(ANTs)", "\n"),file = tmpFile2, append = TRUE)
        cat(paste0("resultANTsJobs",a, " = ANTs:::stat.glmm.no.first.model(odf =", noquote("odf"), ",") ,file = tmpFile2, append = TRUE)

      }
      for (b in 1:length(argg)) {
        if(b == 1){
          if(is.null(argg[[b]])){next}
          if(names(argg)[b] %in% "ant"){
            cat(paste0(names(argg)[b], "=", "t" ),file = tmpFile2, append = TRUE)
            next
          }
          if(names(argg)[b] %in% "formula"){
            cat(paste0(names(argg)[b], "=", Reduce(paste, deparse(argg[[b]]))),file = tmpFile2, append = TRUE)
            next
          }
          if(names(argg)[b] %in% "family"){
            cat(paste0(names(argg)[b], "=", "\'", fam, "\'"),file = tmpFile2, append = TRUE)
            next
          }
          if(names(argg)[b] %in% "ncores"){
            next
          }
          cat(paste0(names(argg)[b], "=", argg[[b]],"," ),file = tmpFile2, append = TRUE)
        }else{
          if(is.null(argg[[b]])){next}
          if(names(argg)[b] %in% "ant"){
            cat(paste0(",", names(argg)[b], "=", "t" ),file = tmpFile2, append = TRUE)
            next
          }
          if(names(argg)[b] %in% "formula"){
            cat(paste0(",", names(argg)[b], "=", Reduce(paste, deparse(argg[[b]]))),file = tmpFile2, append = TRUE)
            next
          }
          if(names(argg)[b] %in% "family"){
            cat(paste0(",",names(argg)[b], "=", "\'",fam, "\'"),file = tmpFile2, append = TRUE)
            next
          }
          if(names(argg)[b] %in% "ncores"){
            next
          }
          cat(paste0(",", names(argg)[b], "=", argg[[b]]),file = tmpFile2, append = TRUE)
        }
      }
    
      cat(paste0(")","\n","rm(t)"),file = tmpFile2, append = TRUE)
      results = c(results, paste("resultANTsJobs",a, sep = ""))
      # Lunch jobs ---------------------------------------------------
      path2 = gsub("\\", "/", tmpFile2, fixed = TRUE)
      job.info[[a]] <-rstudioapi::jobRunScript(path2,
                                             name = paste("Job #", a),
                                             importEnv = FALSE,
                                             exportEnv = "R_GlobalEnv")
    
  }
  r = list(results, unlist(job.info))
  attr( r, "class") <- "ant glmm parallel"
  return(r)
}


