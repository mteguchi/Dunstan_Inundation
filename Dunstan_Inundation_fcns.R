
# computes LOOIC
compute.LOOIC <- function(loglik, data.vector, MCMC.params){
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  
  loglik.vec <- as.vector(loglik)
  
  # each column corresponds to a data point and rows are MCMC samples
  loglik.mat <- matrix(loglik.vec, nrow = n.per.chain * MCMC.params$n.chains)
  
  # take out the columns that correspond to missing data points
  loglik.mat <- loglik.mat[, !is.na(data.vector)]
  # loglik.mat <- matrix(loglik.vec[!is.na(data.vector)], 
  #                      nrow = MCMC.params$n.chains * n.per.chain)
  
  Reff <- relative_eff(exp(loglik.mat),
                       chain_id = rep(1:MCMC.params$n.chains,
                                      each = n.per.chain),
                       cores = 4)
  
  #
  loo.out <- loo(loglik.mat, 
                 r_eff = Reff, 
                 cores = 4, k_threshold = 0.7)
  
  out.list <- list(Reff = Reff,
                   loo.out = loo.out)
  
  return(out.list)  
}


# Extracting posterior samples from jags output:
extract.samples <- function(varname, zm){
  dev <- do.call(rbind, zm[, varname])
  return(dev)
}

# calculating Bayesian R2 values with logit - Gelman et al. 2018
bayes.logit.R2 <- function(X, betas, samples, sigma.name){
  # predicted y (y^pred in Gelman et al. 2018), just Xb
  # but we are using logit(y) so need to do inverse logit
  y.pred <- rstanarm::invlogit(X %*% t(betas))   # y-tilde
  
  # R2 = var.fit/(var.fit + var.res)  -- Gelman et al. 2018
  # var.fit = the variance of the modeled predictive means = the variance among the expectations of the new data
  var.fit <- apply(y.pred, MARGIN = 2, FUN = var)  
  #var.fit <- var(means.y)
  
  # var.res is the modeled residual variance 
  sigma.res <- unlist(samples[,sigma.name])
  var.res <- sigma.res^2
  #var(apply(pred.y - means.y, MARGIN = 1, FUN = mean))
  
  R2 <- var.fit/(var.fit + var.res)
  return(R2)  
}

# calculating Bayesian R2 values - Gelman et al. 2018
bayes.R2 <- function(X, betas, samples, sigma.name){
  # predicted y (y^pred in Gelman et al. 2018), just Xb
  # but we are using logit(y) so need to do inverse logit
  y.pred <- X %*% t(betas)   # y-tilde
  
  # R2 = var.fit/(var.fit + var.res)  -- Gelman et al. 2018
  # var.fit = the variance of the modeled predictive means = the variance among the expectations of the new data
  var.fit <- apply(y.pred, MARGIN = 2, FUN = var)  
  #var.fit <- var(means.y)
  
  # var.res is the modeled residual variance (this is approximate - modeled variance)
  sigma.res <- unlist(samples[,sigma.name])
  var.res <- sigma.res^2
  #var(apply(pred.y - means.y, MARGIN = 1, FUN = mean))
  
  R2 <- var.fit/(var.fit + var.res)
  return(R2)  
  
}



std_dist_berm <- function(x){
  out <- (x - min(x))/sqrt(var(x))
  return(out)
}


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

# poly2 prediction for any distance/yr/month combo
# new.data should have three variables: Month, Year, and scaled_Dist_Berm
# Month should be transformed into nesting month, where October is 1 and 
# April is 7. Year should be transformed so that 2018 is 1 and 2019 is 2.
predict.poly2 <- function(zm, new.data){
  b2 <- extract.samples("b2", zm$samples)
  pred.Y <- vector(mode = "list", length = nrow(new.data))
  k <- 1
  for (k in 1:nrow(new.data)){
    b0 <- extract.samples(paste0("b0[", new.data$Month[k], "]"), zm$samples)
    b1 <- extract.samples(paste0("b1[", new.data$Month[k], "]"), zm$samples)
    b12 <- extract.samples(paste0("b12[", new.data$Month[k], "]"), zm$samples)
    pred.Y[[k]] <- b0 + b1 * new.data$scaled_Dist_Berm[k] + 
      b12 * (new.data$scaled_Dist_Berm[k])^2 + b2 * new.data$Year[k]
    
  }

  # get quantiles
  qtiles <- lapply(pred.Y, 
                   FUN = quantile, c(0.025, 0.5, 0.975),
                   na.rm = T)
  
  qtiles.df <- data.frame(do.call(rbind, lapply(qtiles, FUN = t)))
  
  out.list <- list(pred.Y = pred.Y,
                   qtiles = qtiles.df)
  return(out.list)
  
}

# poly2.1 prediction for any distance/yr/month combo
# new.data should have three variables: Month, Year, and scaled_Dist_Berm
# Month should be transformed into nesting month, where October is 1 and 
# April is 7. Year should be transformed so that 2018 is 1 and 2019 is 2.
predict.poly2.1 <- function(zm, new.data){
  
  pred.Y <- vector(mode = "list", length = nrow(new.data))
  k <- 1
  for (k in 1:nrow(new.data)){
    b0 <- extract.samples(paste0("b0[", new.data$Month[k], 
                                 ",", new.data$Year[k], "]"), zm$samples)
    b1 <- extract.samples(paste0("b1[", new.data$Month[k], "]"), zm$samples)
    b12 <- extract.samples(paste0("b12[", new.data$Month[k], "]"), zm$samples)
    pred.Y[[k]] <- b0 + b1 * new.data$scaled_Dist_Berm[k] + 
      b12 * (new.data$scaled_Dist_Berm[k])^2 
    
  }
  
  # get quantiles
  qtiles <- lapply(pred.Y, 
                   FUN = quantile, c(0.025, 0.5, 0.975),
                   na.rm = T)
  
  qtiles.df <- data.frame(do.call(rbind, lapply(qtiles, FUN = t)))
  
  out.list <- list(pred.Y = pred.Y,
                   qtiles = qtiles.df)
  return(out.list)
  
}


# extract posteriors for 2nd order polynomial with a year term
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
  coefs <- list(b0=b0, b1=b1, b12=b12, b2=b2)

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

# extract posteriors for 2nd order polynomial with a year term and dataset
# specific intercept (not month). Need a lookup table for dataset ID vs year
extract.prediction.poly2.dataset <- function(zm, jags.data, dataset_def){
  length.dataset <- jags.data$n.dataset
  length.year <- jags.data$n.years
  
  # using centered and scaled distances
  Dvec <- seq(floor(min(jags.data$Dist_Berm)), 
              (max(jags.data$Dist_Berm)), 
              by = 0.1)  # min and max of Distance from Berm
  
  #Yvec <- seq(1, length.year)
  b2 <- extract.samples("b2", zm$samples)
  
  b0 <- b1 <- b12 <- pred.Y <- vector(mode = "list", length = length.dataset)
  c <- 1
  for (c in 1:jags.data$n.dataset){
    b0[[c]] <- extract.samples(paste0("b0[", c, "]"), zm$samples)
    b1[[c]] <- extract.samples(paste0("b1[", c, "]"), zm$samples)
    b12[[c]] <- extract.samples(paste0("b12[", c, "]"), zm$samples)
    
    Year <- filter(dataset_def, dataset.ID == c) %>% select(Year2) %>% as.integer()
    
    pred.Y[[c]] <- b0[[c]] + b1[[c]] %*% t(Dvec) + b12[[c]] %*% t(Dvec^2) + b2 * Year
    
  }
  
  # collect all posterior
  coefs <- list(b0=b0, b1=b1, b12=b12, b2=b2)
  
  # get quantiles
  qtiles <- lapply(pred.Y, 
                   FUN = function(x) apply(x, MARGIN = 2, 
                                           FUN = quantile, c(0.025, 0.5, 0.975)))
  
  qtiles.mat <- data.frame(do.call(rbind, lapply(qtiles, FUN = t)))
  
  dataset.vec <- rep(seq(1, jags.data$n.dataset), 
                     each = length(Dvec))
  
  pred.df <- data.frame(dataset.ID = dataset.vec, 
                        Dist_Berm = rep(Dvec, jags.data$n.dataset),
                        pred_low = qtiles.mat$X2.5., 
                        pred_med = qtiles.mat$X50., 
                        pred_high = qtiles.mat$X97.5.)
  
  
  out.list <- list(pred.df = pred.df,
                   coefs = coefs)
  
  return(out.list)
}

# poly2.dataset.month prediction for any distance/yr/month combo
# new.data should have three variables: Month, Year, and Dist_Berm
# Month should be transformed into nesting month, where October is 1 and 
# April is 7. Year should be transformed so that 2018 is 1 and 2019 is 2.
predict.poly2.dataset.month <- function(zm, new.data, dataset_def){
  
  data.idx <- qtiles <- pred.Y <- vector(mode = "list", length = nrow(new.data))
  
  k <- 1
  for (k in 1:nrow(new.data)){
    dataset_def %>% filter(Month2 == new.data$Month[k] & Year2 ==  new.data$Year[k]) %>%
      select(dataset.ID, tide.order, Year2, first.date) %>%
      mutate(tide.in = first.date >= new.data$Date_begin[k] & first.date <= new.data$Date_end[k]) %>%
      arrange(by = tide.order) -> b0.idx 

    data.idx[[k]] <- b0.idx
    b1 <- extract.samples(paste0("b1[", new.data$Month[k], "]"), zm$samples)
    b12 <- extract.samples(paste0("b12[", new.data$Month[k], "]"), zm$samples)
    b2  <- extract.samples(paste0("b2"), zm$samples)
    
    pred.Y.k1 <- vector(mode = "list", length = nrow(b0.idx))
    k1 <- 1
    for (k1 in 1:nrow(b0.idx)){
      b0 <- extract.samples(paste0("b0[", b0.idx[k1,1], "]"), zm$samples)
      pred.Y.k1[[k1]] <- b0 + b1 * new.data$Berm_Dist[k] + 
        b12 * (new.data$Berm_Dist[k])^2 + 
        b2 * as.integer(b0.idx[k1, "Year2"])
    }

    pred.Y[[k]] <- pred.Y.k1
    # get quantiles
    qtiles[[k]] <- lapply(pred.Y[[k]], 
                          FUN = quantile, c(0.025, 0.5, 0.975),
                          na.rm = T)
    
  }
  
  out.list <- list(pred.Y = pred.Y,
                   qtiles = qtiles,
                   data.idx = data.idx)
  return(out.list)
  
}

# extract posteriors for 2nd order polynomial with a year term
extract.prediction.poly2.dataset.month <- function(zm, jags.data, dataset_def){
  length.month <- length(unique(jags.data$month)) 
  length.year <- jags.data$n.years
  length.dataset <- jags.data$n.dataset
  
  # using centered and scaled distances
  Dvec <- seq(floor(min(jags.data$Dist_Berm)), 
              (max(jags.data$Dist_Berm)), 
              by = 0.1)  # min and max of Distance from Berm
  #Yvec <- seq(1, length.year)
  
  b2 <- extract.samples("b2", zm$samples)
  b1 <- b12  <- vector(mode = "list", length = length.month)
  b0 <- pred.Y <- vector(mode = "list", length = length.dataset)
  k <- k1 <- c <- 1

  # month-specific coefficients
  for (k in 1:length.month){
    b1[[k]] <- extract.samples(paste0("b1[", k, "]"), zm$samples)
    b12[[k]] <- extract.samples(paste0("b12[", k, "]"), zm$samples)
  }
  
  k <- 1
  for (k in 1:length.dataset){
    b0[[k]] <- extract.samples(paste0("b0[", k, "]"), zm$samples)
    c.month <- as.integer(dataset_def[k, "Month2"])
    
    pred.Y[[k]] <- b0[[k]] + b1[[c.month]] %*% t(Dvec) + 
      b12[[c.month]] %*% t(Dvec^2) + b2 * as.integer(dataset_def[k, "Year2"])
  }
  
  # collect all posterior samples
  coefs <- list(b0 = b0, b1 = b1, b12 = b12, b2 = b2)
  
  sigma.y <- extract.samples("sigma.y", zm$samples)
  
  # get quantiles
  qtiles <- lapply(pred.Y, 
                   FUN = function(x) apply(x, MARGIN = 2, 
                                           FUN = quantile, c(0.025, 0.5, 0.975)))
  
  qtiles.mat <- data.frame(do.call(rbind, lapply(qtiles, FUN = t)))
  
  dataset.vec <- rep(seq(1, jags.data$n.dataset), 
                     each = length(Dvec))
  
  pred.df <- data.frame(dataset.ID = dataset.vec, 
                        Dist_Berm = rep(Dvec, jags.data$n.dataset),
                        pred_low = qtiles.mat$X2.5., 
                        pred_med = qtiles.mat$X50., 
                        pred_high = qtiles.mat$X97.5.)
  
  out.list <- list(pred.df = pred.df,
                   coefs = coefs,
                   sigma.y = sigma.y)
  
  return(out.list)
}

# extract posteriors for 2nd order polynomial with no year term
# but month/year specific intercepts
extract.prediction.poly2.1 <- function(zm, jags.data){
  length.b <- jags.data$n.months * jags.data$n.years
  length.year <- jags.data$n.years
  
  # using centered and scaled distances
  Dvec <- seq(floor(min(jags.data$Dist_Berm)), 
              (max(jags.data$Dist_Berm)), 
              by = 0.1)  # min and max of Distance from Berm

  b0 <- b1 <- b12 <- pred.Y <- vector(mode = "list", length = length.b)
  k <- k1 <- c <- 1
  for (k1 in 1:length.year){
    for (k in 1:jags.data$n.months){
      b0[[c]] <- extract.samples(paste0("b0[", k,",", k1, "]"), zm$samples)
      b1[[c]] <- extract.samples(paste0("b1[", k, "]"), zm$samples)
      b12[[c]] <- extract.samples(paste0("b12[", k, "]"), zm$samples)
      
      pred.Y[[c]] <- b0[[c]] + b1[[c]] %*% t(Dvec) + 
        b12[[c]] %*% t(Dvec^2) 
      
      c <- c + 1
    }
  }
  
  sigma.y <- extract.samples("sigma.y", zm$samples)
  
  # collect all posterior
  coefs <- list(b0=b0, b1=b1, b12=b12)
  
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
                   coefs = coefs,
                   sigma.y = sigma.y)
  
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
