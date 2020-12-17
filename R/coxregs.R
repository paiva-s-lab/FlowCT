#' Univariable Cox Regression
#' @export unicox

unicox <- function(data, covariates, pval.cutoff = 0.05, timevar, censvar, xlim = c(0,6), log.scale = F, hide.nosig = F){
  data <- as.data.frame(data)
  unicox_pfs <- sapply(covariates, function(x) as.formula(paste0("Surv(", timevar, ", ", censvar, ") ~ `", x, "`")))
  unicox_pfs <- lapply(unicox_pfs, function(x) coxph(x, data))
  
  unicox_pfs_res <- data.frame()
  for(i in covariates){
    if(class(data[,i]) != "numeric"){
      x <- summary(unicox_pfs[[i]])
      coef <- x$coefficients[,1]
      p.value <- x$coefficients[,5]
      HR <- x$conf.int[,1]
      HR.confint.lower <- x$conf.int[,3]
      HR.confint.upper <- x$conf.int[,4]
      exp_1 <- cbind(as.data.frame(coef), as.data.frame(p.value), as.data.frame(HR), 
                     as.data.frame(HR.confint.lower), as.data.frame(HR.confint.upper))
      exp_1 <- rbind(data.frame(row.names = "ref", coef = 0, p.value = 0, HR = 1, HR.confint.lower = 1, HR.confint.upper = 1), 
                     exp_1)
      if(length(unique(data[,i])) > 1){
        rownames(exp_1) <- paste0(i, ".", sort(unique(data[,i])))
      }else{
        rownames(exp_1) <- i
      }
    }else{
      x <- summary(unicox_pfs[[i]])
      coef <- x$coefficients[,1]
      p.value <- x$coefficients[,5]
      HR <- x$conf.int[,1]
      HR.confint.lower <- x$conf.int[,3]
      HR.confint.upper <- x$conf.int[,4]
      exp_1 <- cbind(as.data.frame(coef), as.data.frame(p.value), as.data.frame(HR), 
                     as.data.frame(HR.confint.lower), as.data.frame(HR.confint.upper))
      rownames(exp_1) <- i
    }
    unicox_pfs_res <- rbind(unicox_pfs_res, exp_1)
  }
  
  unicox_pfs_res$clinparam <- rownames(unicox_pfs_res)
  
  if(hide.nosig) unicox_pfs_ress <- unicox_pfs_res[unicox_pfs_res$p.value <= pval.cutoff & unicox_pfs_res$p.value != 0,] 
    else unicox_pfs_ress <- unicox_pfs_res
  sigcol <- ifelse(unicox_pfs_ress$HR != 1 & unicox_pfs_ress$p.value <= pval.cutoff, "black", "gray63")
  
   g <- ggplot(unicox_pfs_ress, aes(x = HR, y = clinparam)) +
    geom_point(size = 3, color = sigcol) +
    geom_errorbar(aes(xmax = HR.confint.lower, xmin = HR.confint.upper), width = .25, color = sigcol) +
    scale_y_discrete(limits = rev(unique(unicox_pfs_ress$clinparam))) +
    theme_minimal(base_size = 9) + theme(axis.text.x = element_text(size = 8))
  
  if(log.scale){
    g <- g + geom_vline(xintercept = 1, color = "gray45", linetype = 2) + 
                                labs(x = "HR (log-transformed)", y = "cell clusters") + 
                                scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                              labels = trans_format("log10", math_format(10^.x)))
  }else{
    g <- g + geom_vline(xintercept = 1, color = "gray45", linetype = 2) + labs(x = "HR", y = "cell clusters") + 
            coord_cartesian(xlim = xlim) #to use limits without removing data (bit.ly/2xpNXq7)
  }
  
  print(g)
  return(unicox_pfs_res)
}


#' Multivariable Cox Regression
#' @export multicox

multicox <- function(data, covariates, pval.cutoff = 0.05, timevar, censvar, xlim = c(0,100), log.scale = F, hide.nosig = T){
  mcox <- summary(coxph(as.formula(paste0("Surv(", timevar, ",", censvar, ") ~ ", paste0("`", covariates, "`", collapse = "+"))), 
                        data = data))
  mcox <- as.data.frame(cbind(mcox$coefficients[,c(1,5)], mcox$conf.int[,c(1,3,4)]))
  colnames(mcox) <- c("coef", "pval", "HR", "HR_low", "HR_up")
  
  if(hide.nosig) mcoxs <- mcox[mcox$pval <= pval.cutoff & mcox$pval != 0,] else mcoxs <- mcox
  sigcol <- ifelse(mcoxs$HR != 1 & mcoxs$pval <= pval.cutoff, "black", "gray63")
  # {if (nrow(mcoxs) == 0) {print("None of the variables is significant!");stop()}}

  mcoxs$clinparam <- rownames(mcoxs)
  # if(log.scale) mcoxs[,c("HR", "HR_up", "HR_low")] <- log10(mcoxs[,c("HR", "HR_up", "HR_low")])

  g <- ggplot(mcoxs, aes(x = HR, y = clinparam)) +
    geom_point(size = 3, color = sigcol) +
    geom_errorbar(aes(xmax = HR_up, xmin = HR_low), width = .25, color = sigcol) +
    scale_y_discrete(limits = rev(unique(mcoxs$clinparam))) +
    theme_minimal(base_size = 9) + theme(axis.text.x = element_text(size = 8))
  
  if(log.scale){
    g <- g + geom_vline(xintercept = 1, color = "gray45", linetype = 2) + 
                          labs(x = "HR (log-transformed)", y = "cell clusters") + 
                        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                      labels = trans_format("log10", math_format(10^.x)))
  }else{
    g <- g + geom_vline(xintercept = 1, color = "gray45", linetype = 2) + labs(x = "HR", y = "cell clusters") + 
            coord_cartesian(xlim = xlim) #to use limits without removing data (bit.ly/2xpNXq7)

  }
  
  print(g)
  return(mcox)
}


#' Forest plot for uni/multivariable Cox regressions
#'
#' Plotting (Forest plot) of uni/multivariable Cox regressions.
#' @param pop.cutoff.obj An object generated through \code{\link[FlowCT.v2:pop.cutoff]{FlowCT.v2::pop.cutoff()}}.
#' @param time.var Survival time variable.
#' @param event.var Variable with event censoring.
#' @param cox.type Cox regression type: "unicox" or "multicox". 
#' @param variables Vector with variables for performing the regression. If nothing is detailed (\code{NULL}, default), all immune's cutoffs will be used.
#' @param ref.var Variable for being used as reference. Default, alphabetical order.
#' @param xlim Limits for x-axis plotting. Default = \code{c(0,10)}.
#' @param return.stats Logical indicating if calculated statistics must be returned. Default = \code{TRUE}.
#' @param log.scale Should x-axis be log-transformed (log10). Default = \code{FALSE}.
#' @param hide.nosig Hide non-significant cell populations. Default = \code{FALSE}.
#' @keywords Cox regression forest-plot multivariable univariable
#' @export cox.plot
#' @examples
#' \dontrun{
#' unicox <- cox.plot(pop.cutoff.obj = pop_cuts, time.var, event.var, cox.type = "multicox", ref.var = "low")
#' }

cox.plot <- function(pop.cutoff.obj, time.var, event.var, cox.type, variables, ref.var, xlim = c(0,10), return.stats = T, log.scale = F, hide.nosig = F){
  if(missing(variables)) variables <- grep(".ct", colnames(pop.cutoff.obj), value = T)
  if(!missing(ref.var)) pop.cutoff.obj[,variables] <- lapply(pop.cutoff.obj[,variables], function(x) relevel(x, ref = ref.var))

  if(tolower(cox.type) == "multicox"){
    cox <- multicox(pop.cutoff.obj, covariates = variables, timevar = time.var, censvar = event.var, log.scale = log.scale, xlim = xlim, hide.nosig = hide.nosig)
  }else if(tolower(cox.type) == "unicox"){
    cox <- unicox(pop.cutoff.obj, covariates = variables, timevar = time.var, censvar = event.var, xlim = xlim, log.scale = log.scale, hide.nosig = hide.nosig)
  }else{
    print("Please, specify one valid option: multiCox or uniCox")
  } 

  if(return.stats) return(cox)
}