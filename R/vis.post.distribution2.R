#' @title Histogram of posterior distribution
#' @description Create histogram of posterior distribution and compute posterior distribution statistics (p-values and confidence interval).
#' @param x A numeric vector with the first value corresponding to the observed value.
#' @param quantile a numeric vector of length 2 to indicate the lower and upper confidence interval.
#' @param backgroud.color A specification for the default histogram background color. 
#' @param observe.value.color A specification for the default abline line color of the observed value. 
#' @param ci.lower.color A specification for the default abline line color of the lower confidence interval. 
#' @param ci.upper.color A specification for the default abline line color of the lower confidence interval. 
#' @param xlab A string for x axis label.
#' @param main A string for histogram main title.
#' @param legend A bolean to print or not the legend
#' @param legend.position The x co-ordinates to be used to position the legend. They can be specified by keyword or in any way which is accepted by xy.coords:
#' @param record A bolean indicating to return or not the histogram in a R object.
#' @return an histogram of posterior distribution.
#' @export
#' @examples 
#' t=met.strength(sim.m,sim.df,1) # Computing network metric
#' t=perm.net.nl(t,labels='age',rf=NULL,nperm=1000,progress=FALSE) # Node label permutations
#' r.c=stat.cor(t,'age','strength',progress=FALSE) # Permuted correlation test
#' vis.post.distribution(r.c[,1])# Histogram of posterior distribution


vis.post.distribution2 <- function(x, quantile = c(0.05, 0.95),
                 backgroud.color = "gray63", 
                 observe.value.color = "white", 
                 ci.lower.color =  "white",
                 ci.upper.color =  "white",
                 xlab = NULL,
                 main = NULL,
                 legend = TRUE,
                 legend.position = "topright",
                 record = TRUE){
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  if(length(quantile) > 2){stop("Only two bornes are allowed for quantiles")}

  if(is.vector(x)){
    # Stats--------------------------------------
    # Permuted p-values
    p = stat.p(x)
    obs = x[1]
    vec = x[-1]
    
    # Confidence interval
    ci = quantile(vec, quantile) 
    
    # Hist---------------------------------------
    par(bg = backgroud.color)
    h <- suppressWarnings(hist(vec, breaks = length(vec), xaxt = "n", plot = FALSE))
    cuts <- cut(h$breaks, c(obs, Inf))
    cuts <- ifelse(is.na(cuts), "gray10", "gray25")
    plot(h, col = cuts, border = cuts, xlab = xlab, main = main)
    
    
    # Ablines-------------------------------------
    abline(v = obs, col = observe.value.color)
    abline(v = ci[1], col = ci.lower.color,lty = 2)
    abline(v = ci[2], col = ci.upper.color, lty = 2)
    
    # Legend--------------------------------------
    if(legend){
      legend(legend.position, bty = "n",
             legend=c(paste("Observed value: ", round(obs, digits = 3)),
                      paste("Lower ci: ", round(ci[1], digits = 3)),
                      paste("Upper ci: ", round(ci[2], digits = 3)),
                      paste("Left side p-value: ",  round(p[1], digits = 3)), 
                      paste("Right side p-value: ",  round(p[2], digits = 3)), 
                      paste("One side p-value: ",  round(p[3], digits = 3))), 
             col=c(observe.value.color, ci.lower.color, ci.upper.color), 
             lty=c( 1, 2, 2, 0, 0, 0), cex=0.8)
    }
  }
  if(is.data.frame(x)){

    if (ncol(x) <= 4) {
      (
        par(mfrow = c(1, ncol(x)))
      )
    }
    if (ncol(x) > 4 & ncol(x) < 7) {
      par(mfrow = c(2, 3))
    }
    if (ncol(x) > 6 & ncol(x) < 9) {
      par(mfrow = c(2, 4))
    }
    if(ncol(x) > 9){return(warning("Number of permuted factors are higher than 9. Use ANTs function post.dist to plot them one by one."))}
  
    for (a in 1:ncol(x)) {
      # Stats--------------------------------------
      # Permuted p-values
      p = stat.p(x[,a])
      obs = x[1,a]
      vec = x[-1,a]
      
      # Confidence interval
      ci = quantile(vec, quantile) 
      
      # Hist---------------------------------------
      par(bg = backgroud.color)
      h <- suppressWarnings(hist(vec, breaks = length(vec), xaxt = "n", plot = FALSE))
      cuts <- cut(h$breaks, c(obs, Inf))
      cuts <- ifelse(is.na(cuts), "gray10", "gray25")
      plot(h, col = cuts, border = cuts, xlab = xlab, main = colnames(x)[])
      
      
      # Ablines-------------------------------------
      abline(v = obs, col = observe.value.color)
      abline(v = ci[1], col = ci.lower.color,lty = 2)
      abline(v = ci[2], col = ci.upper.color, lty = 2)
      
      # Legend--------------------------------------
      if(legend){
        legend(legend.position, bty = "n",
               legend=c(paste("Observed value: ", round(obs, digits = 3)),
                        paste("Lower ci: ", round(ci[1], digits = 3)),
                        paste("Upper ci: ", round(ci[2], digits = 3)),
                        paste("Left side p-value: ",  round(p[1], digits = 3)), 
                        paste("Right side p-value: ",  round(p[2], digits = 3)), 
                        paste("One side p-value: ",  round(p[3], digits = 3))), 
               col=c(observe.value.color, ci.lower.color, ci.upper.color), 
               lty=c( 1, 2, 2, 0, 0, 0), cex=0.8)
      }
    }
  }
}

