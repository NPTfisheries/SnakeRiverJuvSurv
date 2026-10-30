
#code used to develop 
#m-array model in RTMB
#for a single release group
#time-varying survival and detection

library(RTMB)

head(ch_long)

ch_long <- ch_long %>%
  mutate(srr = '12W')

ch_long_list <- ch_long %>%
  group_by(srr, release_season, release_site) %>%
  nest() %>%
  left_join(obs_tbl, by = "release_site")

ch_long_list <- ch_long_list %>%
  mutate(
    ch = map2(
      data,
      obs_loc,
      ~ build_capture_histories(
        tag_history = .x,
        locs_def = .y,
        site_col = "obs_site",
        tag_col = "tag_code",
        time_col = "first_obs",
        censor_col = "last_censored",
        covariate_cols = NULL
      )
    )
  )

marray_list <- ch_long_list %>%
  mutate(m_array = map(
    ch, ~(marray = .x$m_array))) %>%
  select(srr, release_site, m_array)
marray_list$m_array





#M-array for a single release group
#results match the m-array model in PITmodelR

SALRSF <- marray_list$m_array[[14]]
SECESR <- marray_list$m_array[[8]]

#count_cols <- c(grep("^m_", names(SALRSF), value = TRUE), "never")
count_cols <- c(grep("^m_", names(SECESR), value = TRUE), "never")

#m-array formatted with only necessary data
#m_counts <- SALRSF[, count_cols, drop = FALSE]
m_counts <- SECESR[, count_cols, drop = FALSE]

#number released by site
rel <- SALRSF$released
rel <- SECESR$released

#number of sites
n_sites <- ncol(m_counts)

dat = list(mar = as.matrix(m_counts), rel = rel, n_sites = n_sites)
par = list(logit_phi = structure(c(rep(0, n_sites - 2))), logit_p = structure(c(rep(0, n_sites - 2))), logit_lambda = 0)


m_array_fun <- function(par){
  getAll(dat, par)

  lambda <- numeric(1)
  lambda = plogis(logit_lambda)
  
  phi <- p <- q <- rep(0, n_sites - 2)
  
  for (t in 1:(n_sites-2)) {
      phi[t] = plogis(logit_phi[t]) 
      p[t] = plogis(logit_p[t])
      q[t] = 1 - p[t]
  }
  
  pi <- matrix(0, nrow = n_sites - 1, ncol = n_sites)  
  
  #until the final site
  for (t in 1:(n_sites - 1)) {
    
    if (t < (n_sites - 1)) {
    
    #vector of probability of non-recapture
    #q[t] <- 1 - p[t]
    
    #main diagonal: probability of survival and recapture
    pi[t,t] <- phi[t] * p[t]
    
    
    #above the diagonal
    if (t < n_sites - 2) {
     for (j in (t+1):(n_sites - 2)) {
       pi[t, j] = prod(phi[t:j]) * prod(q[t:(j-1)]) * p[j]  
     }
    }
    
    pi[t, n_sites - 1] = prod(phi[t:(n_sites - 2)]) * prod(q[t:(n_sites - 2)]) * lambda
       
    } else {
  
      pi[t, t] <- lambda
    }
  }

  #last column: probability of non-recapture
  for (t in 1:(n_sites - 1)) {
    pi[t, n_sites] <- 1 - sum(pi[t,1:(n_sites-1)])
    
  } 
  
  nll <- 0
  
  for (t in 1:(n_sites - 1)) {
    nll <- nll - dmultinom(mar[t, t:n_sites], prob = pi[t, t:n_sites], log = TRUE)
  }
  
  #need to add confidence interval calculation

  REPORT(phi)
  REPORT(p)
  REPORT(pi)
  REPORT(lambda)
  
  ADREPORT(phi)
  ADREPORT(p)
  ADREPORT(lambda)
  
  nll
         
}

obj <- RTMB::MakeADFun(m_array_fun, par)
obj$fn()
obj$gr()

opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he,control = list(eval.max = 1e4, iter.max = 1e4))
sdr <- sdreport(obj)
opt
sdr

obj$report()
rep <- obj$report()

st_err <- as.list(sdr, "Std", report = T)





#### model for two groups with different survival and same detection
# test using Secesh and Grouse, which have the same detection sites

#eventually wrap the data formatting into a function
#results are slightly different than CJS multiple release group

GROUSC <- marray_list$m_array[[4]]
count_cols <- c(grep("^m_", names(GROUSC), value = TRUE), "never")
mar_gr <- as.matrix(GROUSC[, count_cols, drop = FALSE])

SECESR <- marray_list$m_array[[8]]
count_cols <- c(grep("^m_", names(SECESR), value = TRUE), "never")
mar_se <- as.matrix(SECESR[, count_cols, drop = FALSE])

n_sites <- ncol(mar_se)

n_groups <- 2

#probably will want to list or array data and parameters for different groups

dat <- list(mar = array(c(mar_se, mar_gr), dim = c(n_sites - 1, n_sites, n_groups)), 
            n_sites = n_sites)

par <- list(logit_phi = matrix(0, nrow = n_sites - 2, ncol = n_groups), #unique phi for each group
           logit_p = matrix(0, nrow = n_sites - 2), #shared p
           logit_lambda = structure(rep(0, n_groups))) #unique lambda (?)

m_array_group_fun <- function(par){
  getAll(dat, par)

  lambda = plogis(logit_lambda)
  phi <- matrix(0, nrow  = n_sites - 2, ncol = n_groups)
  
  p <- q <- numeric(n_sites - 2)
  
  for (t in 1:(n_sites-2)) {
    
    p[t] = plogis(logit_p[t])
    
    q[t] <- 1 - p[t]
    
    for (g in 1:n_groups) {
      phi[t, g] <- plogis(logit_phi[t, g])
    }
    
  }
  
  #match m-array dimensions
  pi <- array(0, dim = c(n_sites - 1, n_sites, n_groups))
  
  
  for (g in 1:n_groups) { #loop over release groups
    for (t in 1:(n_sites - 1)) { #loop over release occasions
    
    if (t < (n_sites - 1)) {
      
      #vector of probability of non-recapture
      #q[t] <- 1 - p[t]
      
      #main diagonal: probability of survival and recapture
      pi[t,t,g] <- phi[t,g] * p[t]
      
      
      #above the diagonal
      if (t < n_sites - 2) {
        for (j in (t+1):(n_sites - 2)) {
          pi[t, j, g] = prod(phi[t:j, g]) * prod(q[t:(j-1)]) * p[j]  
        }
      }
      
      pi[t, n_sites - 1, g] = prod(phi[t:(n_sites - 2), g]) * prod(q[t:(n_sites - 2)]) * lambda[g]
      
    } else {
      
      pi[t, t, g] <- lambda[g]
      }
    }
  }
  
  #last column: probability of non-recapture
for (g in 1:n_groups) {
  for (t in 1:(n_sites - 1)) {
    pi[t, n_sites, g] <- 1 - sum(pi[t, 1:(n_sites-1), g])
    
  } 
}
  
  nll <- 0
  
for (g in 1:n_groups) {
  for (t in 1:(n_sites - 1)) {
    nll <- nll - dmultinom(mar[t, t:n_sites, g], prob = pi[t, t:n_sites, g], log = TRUE)
  }
}
  
  #need to add confidence interval calculation
  
  REPORT(phi)
  REPORT(p)
  REPORT(pi)
  REPORT(lambda)
  
  ADREPORT(phi)
  ADREPORT(p)
  ADREPORT(lambda)
  
  nll
  
}

obj <- RTMB::MakeADFun(m_array_group_fun, par)
obj$fn()
obj$gr()

opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he,control = list(eval.max = 1e4, iter.max = 1e4))
sdr <- sdreport(obj)
opt
sdr

obj$report()
rep <- obj$report()

st_err <- as.list(sdr, "Std", report = T)




### model for two groups with a subset of shared detection sites
# test with SECESR and SALRSF
# SECESR - ZEN - SFG - LGR
# SALRSF - ESS - SFG - LGR

#these have the same number of sites
#probably will need to add leading 0s to shorter m-arrays when groups have different numbers of sites

#need one additional p parameter
# same number of phi parameters

n_groups <- 2

MARSHC <- marray_list$m_array[[13]]
count_cols <- c(grep("^m_", names(MARSHC), value = TRUE), "never")
mar_ma <- as.matrix(MARSHC[, count_cols, drop = FALSE])
site_ma <- MARSHC[, 1:2]

SALRSF <- marray_list$m_array[[14]]
count_cols <- c(grep("^m_", names(SALRSF), value = TRUE), "never")
mar_sf <- as.matrix(SALRSF[, count_cols, drop = FALSE])
site_sf <- SALRSF[, 1:2]

SECESR <- marray_list$m_array[[8]]
count_cols <- c(grep("^m_", names(SECESR), value = TRUE), "never")
mar_se <- as.matrix(SECESR[, count_cols, drop = FALSE])
site_se <- SECESR[, 1:2]

SAEFSF <- marray_list$m_array[[7]]
count_cols <- c(grep("^m_", names(SAEFSF), value = TRUE), "never")
mar_ef <- as.matrix(SAEFSF[, count_cols, drop = FALSE])
site_ef <- SAEFSF[, 1:2]

#matrix with detection sites by group
#site_mat <- matrix(cbind(SALRSF$site, SECESR$site, SAEFSF$site), ncol = n_groups)[-1, ]
site_mat <- matrix(cbind(SALRSF$site[-1], site_ma$site[-1]), ncol = 2)
unique_sites <- unique(c(SALRSF$site[-1], SECESR$site[-1], SAEFSF$site[-1]))

#doing this manually to check model. will fix later
site_mat <- matrix(c("KRS", "SFG", "LGR", "MAR", "LGR", NA), ncol = 2)
unique_sites <- c("KRS", "SFG", "LGR", "MAR")

#number of detection sites, also the number of p parameters to estimate
n_sites <- length(unique(c(SALRSF$site[-1], SECESR$site[-1], SAEFSF$site[-1])))
n_sites <- length(unique_sites)

#changed from site to occasion to avoid confusion with detection sites
n_occ <- c(ncol(mar_sf), ncol(mar_ma))
max_n_occ <- max(n_occ)

#COME BACK TO THIS
#sites are by alphabetical order
p_index_mat <- matrix(as.numeric(as.factor(site_mat)), nrow = n_occ - 2, byrow = F)
p_index <- list(p_index_mat[, 1][!is.na(p_index_mat[, 1])], 
                p_index_mat[, 2][!is.na(p_index_mat[, 2])])

#doing this manually for now, fix later if it works
phi_index_mat <- matrix(c(1,2,3,4,5,NA), nrow = n_occ - 2, byrow = F)

mar_list = list(mar_sf, mar_ma)

target_dim <- dim(mar_list[[1]])

pad_to_dim <- function(mat, target_dim) {
  padded <- matrix(0, nrow = target_dim[1], ncol = target_dim[2])
  # place original values into matching column names, matching rows (top-aligned)
  padded[1:nrow(mat), 1:ncol(mat)] <- mat
  padded
}

mar_list_padded <- lapply(mar_list, pad_to_dim,
                          target_dim = target_dim)

mar_array <- simplify2array(mar_list_padded)

mar_array



dat <- list(mar = mar_array,
            p_index = p_index_mat, 
            phi_index = phi_index_mat,
            n_occ = n_occ, # a vector
            n_groups = n_groups, 
            max_n_occ = max_n_occ)

#define phi lists
#phi_sf <- rep(0, dim(mar_sf)[1] - 1)
#phi_ma <- rep(0, dim(mar_ma)[1] - 1)

par <- list(logit_phi = matrix(0, nrow = sum(n_occ) - n_groups*2),
            #logit_phi = matrix(0, nrow = max_n_occ - 2, ncol = n_groups), #unique phi for each group, set to max dimensions
            logit_p = matrix(0, nrow = n_sites), #one p for each detection site
            #logit_p = matrix(0, nrow = n_occ - 2, ncol = n_groups),
            logit_lambda = structure(rep(0, n_groups))) #unique lambda for each group

m_array_multiple_fun <- function(par){
  getAll(dat, par)
  
  lambda = plogis(logit_lambda)
  
  phi <- numeric(sum(n_occ) - n_groups*2) #n_occ - 2 for each group
  p <- numeric(n_sites)
  
  p <- plogis(logit_p)
  phi <- plogis(logit_phi)
  
  nll <- 0
  
  for (g in 1:n_groups) { #loop over release groups
    
    n_occ_g <- n_occ[g]
    
    q_g <- numeric(n_occ_g - 2)
    
    pi_g <- matrix(0, nrow = n_occ_g - 1, ncol = n_occ_g)
  
    
    for (t in 1:(n_occ_g - 1)) { #loop over release occasions
      
      if (t < (n_occ_g - 1)) {
        
        #define q (probably of non-detection)
        p_idx <- p_index_mat[t, g]
        
        q_g[t] <- 1 - p[p_idx]
      }
    }
        
    for (t in 1:(n_occ_g - 1)) { #loop over release occasions
      
      if (t < (n_occ_g - 1)) {
        
        #index with p and phi parameters to use based on occasion and group
        phi_idx <- phi_index_mat[t, g]
        
        p_idx <- p_index_mat[t, g]
        
        #main diagonal: probability of survival and recapture
        pi_g[t,t] <- phi[phi_idx] * p[p_idx]
        
        #above the diagonal
        if (t < n_occ_g - 2) {
          for (j in (t+1):(n_occ_g - 2)) {
            pj <- p[p_index[j, g]]
            phij <- phi[phi_index[, g]]
            pi_g[t, j] = prod(phij[t:j]) * prod(q_g[t:(j-1)]) * pj 
          }
        }
        
        pi_g[t, n_occ_g - 1] = prod(phij[t:(n_occ_g - 2)]) * prod(q_g[t:(n_occ_g - 2)]) * lambda[g]
        
      } else {
        
        pi_g[t, t] <- lambda[g]
      }
    }

  #last column: probability of non-recapture
    for (t in 1:(n_occ_g - 1)) {
      pi_g[t, n_occ_g] <- 1 - sum(pi_g[t, 1:(n_occ_g-1)])
    } 

    #likelihood
    for (t in 1:(n_occ_g - 1)) {
      nll <- nll - dmultinom(mar[t, t:n_occ_g, g], prob = pi_g[t, t:n_occ_g], log = TRUE)
    }
  }
  
  REPORT(phi)
  REPORT(p)
  REPORT(pi_g)
  REPORT(lambda)
  
  ADREPORT(phi)
  ADREPORT(p)
  ADREPORT(lambda)
  
  nll
  
}



obj <- RTMB::MakeADFun(m_array_multiple_fun, par)
obj$fn()
obj$gr()

opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he,control = list(eval.max = 1e4, iter.max = 1e4))
sdr <- sdreport(obj)
opt
sdr

obj$report()
rep <- obj$report()

st_err <- as.list(sdr, "Std", report = T)


### wrap the multiple release group model into a function

#m_array_list is a tibble with column release_group and column m_array. m_array column contains a df with the m-array
#still missing 95% CIs and formatting the output

fit_marray_cjs_multiple <- function(marray_list) {
  
  #compile the m-array data
  mar_list <- lapply(marray_list$m_array, function(df) {
    count_cols <- c(grep("^m_", names(df), value = TRUE), "never")
    as.matrix(df[, count_cols, drop = FALSE])
  })
    
  n_occ <- sapply(mar_list, ncol)
  max_n_occ <- max(n_occ)
  
  max_group <- which.max(n_occ)
  target_dim <- dim(mar_list[[max_group]])
  
  pad_to_dim <- function(mat, target_dim) {
    padded <- matrix(0, nrow = target_dim[1], ncol = target_dim[2])
    # place original values into matching column names, matching rows (top-aligned)
    padded[1:nrow(mat), 1:ncol(mat)] <- mat
    padded
  }
  
  mar_list_padded <- lapply(mar_list, pad_to_dim,
                            target_dim = target_dim)
  
  mar_array <- simplify2array(mar_list_padded) 
  
  n_groups = dim(mar_array)[3]
  
  #p index matrix
  site_lists <- lapply(marray_list$m_array, function(df) df$site[-1])
  max_len <- max(sapply(site_lists, length))
  
  site_mat <- sapply(site_lists, function(x) {
    length(x) <- max_len  
    x
  })
  
  n_sites = length(unique(unlist(site_lists)))
  
  #each detection site is assigned a number
  #they are alphabetical
  p_index_mat <- matrix(as.numeric(as.factor(site_mat)), nrow = max(n_occ) - 2, byrow = F)
  
  #phi index matrix
  #using the entries in the p_index_mat as a template
  phi_index_mat <- matrix(NA_integer_, nrow = nrow(site_mat), ncol = ncol(site_mat))
  
  phi_index_mat[!is.na(site_mat)] <- matrix(seq_len(sum(!is.na(site_mat))))
  
  dimnames(phi_index_mat) <- dimnames(site_mat)
  
  #compile data and parameters
  dat <- list(mar = mar_array,
              p_index = p_index_mat, 
              phi_index = phi_index_mat,
              n_occ = n_occ, # a vector
              n_groups = n_groups, #number of release groups 
              max_n_occ = max_n_occ, 
              n_sites = n_sites)
  
  par <- list(logit_phi = matrix(0, nrow = sum(n_occ) - n_groups*2),
              logit_p = matrix(0, nrow = n_sites), #one p for each detection site
              logit_lambda = structure(rep(0, n_groups))) #unique lambda for each group

  m_array_multiple_fun <- function(par){
  getAll(dat, par)
  "[<-" <- ADoverload("[<-")
  
  lambda = plogis(logit_lambda)
  
  phi <- numeric(sum(n_occ) - n_groups*2) #n_occ - 2 for each group
  p <- numeric(n_sites)
  
  p <- plogis(logit_p)
  phi <- plogis(logit_phi)
  
  nll <- 0
  
  for (g in 1:n_groups) { #loop over release groups
    
    n_occ_g <- n_occ[g]
    
    q_g <- numeric(n_occ_g - 2)
    
    pi_g <- matrix(0, nrow = n_occ_g - 1, ncol = n_occ_g)
    
    
    for (t in 1:(n_occ_g - 1)) { #loop over release occasions
      
      if (t < (n_occ_g - 1)) {
        
        #define q (probably of non-detection)
        p_idx <- p_index[t, g]
        
        q_g[t] <- 1 - p[p_idx]
      }
    }
    
    for (t in 1:(n_occ_g - 1)) { #loop over release occasions
      
      if (t < (n_occ_g - 1)) {
        
        #index with p and phi parameters to use based on occasion and group
        phi_idx <- phi_index[t, g]
        
        p_idx <- p_index[t, g]
        
        #main diagonal: probability of survival and recapture
        pi_g[t,t] <- phi[phi_idx] * p[p_idx]
        
        #above the diagonal
        if (t < n_occ_g - 2) {
          for (j in (t+1):(n_occ_g - 2)) {
            pj <- p[p_index[j, g]]
            phij <- phi[phi_index[, g]]
            pi_g[t, j] = prod(phij[t:j]) * prod(q_g[t:(j-1)]) * pj 
          }
        }
        
        pi_g[t, n_occ_g - 1] = prod(phij[t:(n_occ_g - 2)]) * prod(q_g[t:(n_occ_g - 2)]) * lambda[g]
        
      } else {
        
        pi_g[t, t] <- lambda[g]
      }
    }
    
    #last column: probability of non-recapture
    for (t in 1:(n_occ_g - 1)) {
      pi_g[t, n_occ_g] <- 1 - sum(pi_g[t, 1:(n_occ_g-1)])
    } 
    
    #likelihood
    for (t in 1:(n_occ_g - 1)) {
      nll <- nll - dmultinom(mar[t, t:n_occ_g, g], prob = pi_g[t, t:n_occ_g], log = TRUE)
    }
  }
  
  REPORT(phi)
  REPORT(p)
  REPORT(lambda)
  REPORT(logit_phi)
  REPORT(logit_p)
  
  ADREPORT(phi)
  ADREPORT(p)
  ADREPORT(lambda)
  ADREPORT(logit_phi)
  ADREPORT(logit_p)
  
  nll
  
}

  obj <- RTMB::MakeADFun(m_array_multiple_fun, par)
  obj$fn()
  obj$gr()

  opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he,control = list(eval.max = 1e4, iter.max = 1e4))
  sdr <- sdreport(obj)

  obj$report()
  rep <- obj$report()

  st_err <- as.list(sdr, "Std", report = T)

  fit <- list(rep, st_err)
  fit
  
  

#format output

  #format p
  site_id <- as.data.frame(cbind(unique(unlist(site_lists)), as.factor(unique(unlist(site_lists)))))
  colnames(site_id) <- c("site", "num")
  
  p <- as.data.frame(cbind(fit[[1]]$p, fit[[1]]$logit_p, fit[[2]]$logit_p))
  p <- p %>%
    rownames_to_column(var = "num") %>%
    rename(p = V1, logit_p = V2, se_logit_p = V3) %>%
    mutate(logit_lcl = logit_p - 1.96*se_logit_p, 
           logit_ucl = logit_p + 1.96*se_logit_p, 
           lcl = plogis(logit_lcl), 
           ucl = plogis(logit_ucl)) %>%
    full_join(site_id, by = "num") %>%
    select(p, lcl, ucl, site)
  
  #format phi
  fit[[1]]$phi
  
  from_site <- lapply(marray_list$m_array, function(df) df$site[-length(df$site)])
  
  phi_id <- as.data.frame(site_mat)
  colnames(phi_id) <- marray_list$release_group
  
  phi_id <- phi_id %>% 
    pivot_longer(cols = 1:dim(phi_id)[2], names_to = "release_group", values_to = "to_site", cols_vary = "slowest") %>%
    filter(!is.na(to_site))
  
  phi <- cbind(phi_id, fit[[1]]$phi, fit[[1]]$logit_phi, fit[[2]]$logit_phi) %>%
    rename(phi = 'fit[[1]]$phi', 
           logit_phi = 'fit[[1]]$logit_phi', 
           se_logit_phi = 'fit[[2]]$logit_phi') %>%
    mutate(logit_lcl = logit_phi - 1.96*se_logit_phi, 
           logit_ucl = logit_phi + 1.96*se_logit_phi, 
           lcl = plogis(logit_lcl), 
           ucl = plogis(logit_ucl)) %>%
    select(release_group, to_site, phi, lcl, ucl) %>%
    mutate(from_site = unlist(from_site))
  
  marray_output <- list(
    p = p, 
    phi = phi, 
    fit = rep, 
    st_err = st_err)
  marray_output

}


m_array_fit <- fit_marray_cjs_multiple(marray_list)







#this will go in a wrapper function for consistency
ch_long_list <- ch_long_list %>%
      mutate(
        ch = map2(
          data,
          obs_loc,
          ~ build_capture_histories(
            tag_history = .x,
            locs_def = .y,
            site_col = "obs_site",
            tag_col = "tag_code",
            time_col = "first_obs",
            censor_col = "last_censored",
            covariate_cols = NULL
          )
        )
      )
    
    #individual m array for each release group
    marray_list <- ch_long_list %>%
      mutate(m_array = map(
        ch, ~(marray = .x$m_array))) %>%
      select(srr, release_site, m_array) %>%
      rename(release_group = release_site)
    