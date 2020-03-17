# a function to extract posterior samples from jags output
extract.samples <- function(varname, zm){
  dev <- unlist(lapply(zm, FUN = function(x) x[, varname]))
  return(dev)
}

# A function to find 95% CI for predictions
extract.prediction <- function(zm, jags.data){
  if (length(grep("n.transects", names(jags.data))) > 0){
    length.b <- jags.data$n.transects
  } else if (length(grep("n.months", names(jags.data))) > 0){
    length.b <- jags.data$n.months
  }
  
  b0 <- b1 <- b12 <- vector(mode = "list", length = length.b)
  k <- 1
  for (k in 1:length.b){
    b0[[k]] <- extract.samples(paste0("b0[", k, "]"), zm$samples)
    b1[[k]] <- extract.samples(paste0("b1[", k, "]"), zm$samples)
    
    # in case of poly2
    if (length(grep("b12", unlist(dimnames(zm$samples[[1]])))) > 0){
      b12[[k]] <- extract.samples(paste0("b12[", k, "]"), zm$samples)
    } else {
      b12[[k]] <- rep(0, length(b0[[k]]))
    }
  }
  length.year <- jags.data$n.years
  
  b2 <- vector(mode = "list", length = length.year)
  for (k in 1:length.year){
    b2[[k]] <- extract.samples(paste0("b2[", k, "]"), zm$samples)
  }
  
  Dvec <- seq(floor(min(jags.data$Dist_Berm)), 
              (max(jags.data$Dist_Berm)), 
              by = 0.1)  # min and max of Distance from Berm
  
  Yvec <- seq(1, length.year)
  
  pred.Y <- vector(mode = "list", length = (length.year * length.b))
  k <- y <- c <- 1
  for (k in 1:length.b){      # transects/months
    for (y in 1:length.year){    # years
      pred.Y[[c]] <- b0[[k]] + b1[[k]] %*% t(Dvec) + 
        b12[[k]] %*% t(Dvec^2) + b2[[y]] * Yvec[y]
      
      c <- c + 1
    }
    
  }
  
  coefs <- c(b0, b1, b12, b2)
  # get quantiles
  
  qtiles <- lapply(pred.Y, 
                   FUN = function(x) apply(x, MARGIN = 2, 
                                           FUN = quantile, c(0.025, 0.5, 0.975)))
  
  qtiles.mat <- data.frame(do.call(rbind, lapply(qtiles, FUN = t)))
  
  year.vec <- rep(rep(seq(1, length.year), 
                      each = length(Dvec)), 
                  length.b)
  
  b.vec <- rep(seq(1, length.b), each = length(Dvec) * length(Yvec))
  
  if (length(grep("n.transects", names(jags.data))) > 0){
    pred.df <- data.frame(Transect_f = as.factor(b.vec + 1), 
                          Year = year.vec + 2017, 
                          Dist_Berm = rep(Dvec, length.b * length(Yvec)),
                          pred_low = qtiles.mat$X2.5., 
                          pred_med = qtiles.mat$X50., 
                          pred_high = qtiles.mat$X97.5.)
  } else if (length(grep("n.months", names(jags.data))) > 0){
    pred.df <- data.frame(Month2 = b.vec, 
                          Year = year.vec + 2017, 
                          Dist_Berm = rep(Dvec, length.b * length(Yvec)),
                          pred_low = qtiles.mat$X2.5., 
                          pred_med = qtiles.mat$X50., 
                          pred_high = qtiles.mat$X97.5.)
    
  }
  
  out.list <- list(pred.df = pred.df,
                   coefs = coefs)
  
  return(out.list)
}


extract.prediction.poly3 <- function(zm, jags.data){
  #length.b <- jags.data$n.months
  length.year <- jags.data$n.years
  length.month <- jags.data$n.months
  
  Dvec <- seq(floor(min(jags.data$Dist_Berm)), 
              (max(jags.data$Dist_Berm)), 
              by = 0.1)  # min and max of Distance from Berm
  
  b2 <- extract.samples("b2", zm$samples)
  
  b0 <- b1 <- b12 <- b13 <- pred.Y <- vector(mode = "list", 
                                             length = length.month)
  
  k <- k1 <- c <- 1
  for (k1 in 1:length.year){
    for (k in 1:jags.data$n.months){
      b0[[c]] <- extract.samples(paste0("b0[", k, "]"), zm$samples)
      b1[[c]] <- extract.samples(paste0("b1[", k, "]"), zm$samples)
      b12[[c]] <- extract.samples(paste0("b12[", k, "]"), zm$samples)
      b13[[c]] <- extract.samples(paste0("b13[", k, "]"), zm$samples)
      
      pred.Y[[c]] <- b0[[c]] + b1[[c]] %*% t(Dvec) + 
        b12[[c]] %*% t(Dvec^2) + 
        b13[[c]] %*% t(Dvec^3) + b2 * k1
      
      c <- c + 1
    }
  }
  
  coefs <- c(b0, b1, b12, b13, b2)
  # get quantiles
  
  qtiles <- lapply(pred.Y, 
                   FUN = function(x) apply(x, MARGIN = 2, 
                                           FUN = quantile, c(0.025, 0.5, 0.975)))
  
  qtiles.mat <- data.frame(do.call(rbind, lapply(qtiles, FUN = t)))
  
  year.vec <- rep(seq(1, length.year), 
                  each = length(Dvec) * length.month)
  
  month.vec <- rep(rep(seq(1, length.month), 
                       each = length(Dvec)),
                   length.year)
  
  pred.df <- data.frame(Month2 = month.vec, 
                        Year = year.vec + 2017, 
                        Dist_Berm = rep(rep(Dvec, length.month), 
                                        length.year),
                        pred_low = qtiles.mat$X2.5., 
                        pred_med = qtiles.mat$X50., 
                        pred_high = qtiles.mat$X97.5.)
  
  out.list <- list(pred.df = pred.df,
                   coefs = coefs)
  
  return(out.list)
}

# extract posteriors
extract.prediction.poly2 <- function(zm, jags.data){
  length.b <- jags.data$n.months * jags.data$n.years
  length.year <- jags.data$n.years
  
  # using centered and scaled distances
  Dvec <- seq(floor(min(jags.data$Dist_Berm)), 
              (max(jags.data$Dist_Berm)), 
              by = 0.1)  # min and max of Distance from Berm
  #Yvec <- seq(1, length.year)
  b2 <- extract.samples("b2", zm$samples)
  
  b0 <- b1 <- b12 <- pred.Y <- vector(mode = "list", length = length.b)
  k <- k1 <- c <- 1
  for (k1 in 1:length.year){
    for (k in 1:jags.data$n.months){
      b0[[c]] <- extract.samples(paste0("b0[", k, "]"), zm$samples)
      b1[[c]] <- extract.samples(paste0("b1[", k, "]"), zm$samples)
      b12[[c]] <- extract.samples(paste0("b12[", k, "]"), zm$samples)
      
      pred.Y[[c]] <- b0[[c]] + b1[[c]] %*% t(Dvec) + 
        b12[[c]] %*% t(Dvec^2) + b2 * k1
      
      c <- c + 1
    }
  }
  
  # collect all posterior
  coefs <- c(b0, b1, b12, b2)

  # get quantiles
  qtiles <- lapply(pred.Y, 
                   FUN = function(x) apply(x, MARGIN = 2, 
                                           FUN = quantile, c(0.025, 0.5, 0.975)))
  
  qtiles.mat <- data.frame(do.call(rbind, lapply(qtiles, FUN = t)))
  
  year.vec <- rep(seq(1, length.year), 
                  each = length(Dvec) * jags.data$n.months)
  
  month.vec <- rep(rep(seq(1, jags.data$n.months), 
                       each = length(Dvec)),
                   length.year)
  
  pred.df <- data.frame(Month2 = month.vec, 
                        Year = year.vec + 2017, 
                        Dist_Berm = rep(rep(Dvec, jags.data$n.months), 
                                        length.year),
                        pred_low = qtiles.mat$X2.5., 
                        pred_med = qtiles.mat$X50., 
                        pred_high = qtiles.mat$X97.5.)
  
  
  out.list <- list(pred.df = pred.df,
                   coefs = coefs)
  
  return(out.list)
}

# extract posteriors
extract.prediction.1 <- function(zm, jags.data){
  length.b <- jags.data$n.months * jags.data$n.years
  length.year <- jags.data$n.years
  
  # using centered and scaled distances
  Dvec <- seq(floor(min(jags.data$Dist_Berm)), 
              (max(jags.data$Dist_Berm)), 
              by = 0.1)  # min and max of Distance from Berm
  Yvec <- seq(1, length.year)
  b2 <- extract.samples("b2", zm$samples)
  
  b0 <- b1 <- pred.Y <- vector(mode = "list", length = length.b)
  k <- k1 <- c <- 1
  for (k1 in 1:length.year){
    for (k in 1:jags.data$n.months){
      b0[[c]] <- extract.samples(paste0("b0[", k, "]"), zm$samples)
      b1[[c]] <- extract.samples(paste0("b1[", k, "]"), zm$samples)
      
      pred.Y[[c]] <- b0[[c]] + b1[[c]] %*% t(Dvec) + b2 * k1
      
      c <- c + 1
    }
  }
  
  # collect all posterior
  coefs <- c(b0, b1, b2)
  
  # get quantiles
  qtiles <- lapply(pred.Y, 
                   FUN = function(x) apply(x, MARGIN = 2, 
                                           FUN = quantile, c(0.025, 0.5, 0.975)))
  
  qtiles.mat <- data.frame(do.call(rbind, lapply(qtiles, FUN = t)))
  
  year.vec <- rep(seq(1, length.year), 
                  each = length(Dvec) * jags.data$n.months)
  
  month.vec <- rep(rep(seq(1, jags.data$n.months), 
                       each = length(Dvec)),
                   length.year)
  
  pred.df <- data.frame(Month2 = month.vec, 
                        Year = year.vec + 2017, 
                        Dist_Berm = rep(Dvec, length.b),
                        pred_low = qtiles.mat$X2.5., 
                        pred_med = qtiles.mat$X50., 
                        pred_high = qtiles.mat$X97.5.)
  
  
  out.list <- list(pred.df = pred.df,
                   coefs = coefs)
  
  return(out.list)
}


# extract posteriors
extract.prediction.exp <- function(zm, jags.data){
  length.b <- jags.data$n.months * jags.data$n.years
  length.year <- jags.data$n.years
  
  #b2 <- extract.samples("b2", zm$samples)
  # using centered and scaled distances
  Dvec <- seq(floor(min(jags.data$Dist_Berm)), 
              (max(jags.data$Dist_Berm)), 
              by = 0.1)  # min and max of Distance from Berm
  Yvec <- seq(1, length.year)
  
  b0 <- b1 <- pred.Y <- vector(mode = "list", length = length.b)
  k <- k1 <- c <- 1
  for (k1 in 1:length.year){
    for (k in 1:jags.data$n.months){
      b0[[c]] <- extract.samples(paste0("b0[", k, ",", k1, "]"), zm$samples)
      b1[[c]] <- extract.samples(paste0("b1[", k, ",", k1, "]"), zm$samples)
      
      pred.Y[[c]] <- b0[[c]] * exp(-b1[[c]] %*% t(Dvec))
      
      c <- c + 1
    }
  }
  
  # collect all posterior
  #coefs <- c(b0, b1, b2)
  coefs <- c(b0, b1)
  
  # get quantiles
  qtiles <- lapply(pred.Y, 
                   FUN = function(x) apply(x, MARGIN = 2, 
                                           FUN = quantile, c(0.025, 0.5, 0.975)))
  
  qtiles.mat <- data.frame(do.call(rbind, lapply(qtiles, FUN = t)))
  
  year.vec <- rep(seq(1, length.year), 
                  each = length(Dvec) * jags.data$n.months)
  
  month.vec <- rep(rep(seq(1, jags.data$n.months), 
                       each = length(Dvec)),
                   length.year)
  
  pred.df <- data.frame(Month2 = month.vec, 
                        Year = year.vec + 2017, 
                        Dist_Berm = rep(Dvec, length.b),
                        pred_low = qtiles.mat$X2.5., 
                        pred_med = qtiles.mat$X50., 
                        pred_high = qtiles.mat$X97.5.)
  
  
  out.list <- list(pred.df = pred.df,
                   coefs = coefs)
  
  return(out.list)
}

# extract posteriors
extract.prediction.exp2 <- function(zm, jags.data){
  length.b <- jags.data$n.months * jags.data$n.years
  length.year <- jags.data$n.years
  
  # using centered and scaled distances
  Dvec <- seq(floor(min(jags.data$Dist_Berm)), 
              (max(jags.data$Dist_Berm)), 
              by = 0.1)  # min and max of Distance from Berm
  Yvec <- seq(1, length.year)
  b2 <- extract.samples("b2", zm$samples)
  
  b0 <- b1 <- pred.Y <- vector(mode = "list", length = length.b)
  k <- k1 <- c <- 1
  for (k1 in 1:length.year){
    for (k in 1:jags.data$n.months){
      b0[[c]] <- extract.samples(paste0("b0[", k, "]"), zm$samples)
      b1[[c]] <- extract.samples(paste0("b1[", k, "]"), zm$samples)
      
      pred.Y[[c]] <- b0[[c]] * exp(b1[[c]] %*% t(Dvec) + b2 * k1) 
      
      c <- c + 1
    }
  }
  
  # collect all posterior
  coefs <- c(b0, b1, b2)
  
  # get quantiles
  qtiles <- lapply(pred.Y, 
                   FUN = function(x) apply(x, MARGIN = 2, 
                                           FUN = quantile, c(0.025, 0.5, 0.975)))
  
  qtiles.mat <- data.frame(do.call(rbind, lapply(qtiles, FUN = t)))
  
  year.vec <- rep(seq(1, length.year), 
                  each = length(Dvec) * jags.data$n.months)
  
  month.vec <- rep(rep(seq(1, jags.data$n.months), 
                       each = length(Dvec)),
                   length.year)
  
  pred.df <- data.frame(Month2 = month.vec, 
                        Year = year.vec + 2017, 
                        Dist_Berm = rep(Dvec, length.b),
                        pred_low = qtiles.mat$X2.5., 
                        pred_med = qtiles.mat$X50., 
                        pred_high = qtiles.mat$X97.5.)
  
  
  out.list <- list(pred.df = pred.df,
                   coefs = coefs)
  
  return(out.list)
}
