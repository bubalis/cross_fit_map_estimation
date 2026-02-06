library(dplyr)
library(ranger)
library(rmeta)



cluster_bootstrap_for_rf <- function(data, group_col, R = 100) {
 
  data$id <- 1:nrow(data)
  
  data_samp <- data[, c('id', group_col)]
  
  # Get the list of unique group IDs (names of the split list)
  unique_groups <- unique(data_samp[[group_col]])
  num_groups <- length(unique_groups)
  
  # 2. Initialize a list to store the R bootstrap samples
  boot_samples <- vector("list", R)
  
  # 3. Perform the Bootstrap
  for (i in 1:R) {
    # Sample group names WITH replacement
    # If we have 10 groups, we pick 10 groups, but some might be duplicates.
    
    sampled_group_names <- sample(unique_groups, num_groups, replace = TRUE)
    
    #Create a vector of the number of times each datapoint is sampled.
    samples_vector <- table(sampled_group_names) %>% 
      data.frame() %>%
      rename(cluster = sampled_group_names ) %>% 
      mutate(cluster = as.integer(as.character(cluster))) %>% 
      right_join(data_samp, by = 'cluster') %>% 
      mutate(Freq = ifelse(is.na(Freq), 0, Freq)) %>% 
      arrange(id) %>% 
      pull(Freq)
    
    boot_samples[[i]] <-  samples_vector
  }
  
  return(boot_samples)
}

#' Simulate a Sample from Clustered Data
#'
#' This function selects a specific number of clusters from a dataset with probability 
#' proportional to their size (as defined by `cluster_counts`), and then draws a 
#' fixed number of random samples from within each selected cluster.
#'
#' @param data A data frame containing the population data. It must contain a column 
#'   named `cluster` which corresponds to the names in `cluster_counts`.
#' @param n.clusters Integer. The number of distinct clusters to select.
#' @param n.per.cluster Integer. The number of observations to sample from each 
#'   selected cluster.
#' @param cluster_counts A named integer or numeric vector. The values represent the 
#'   size (or weight) of the clusters, and the names must correspond to the cluster 
#'   IDs found in `data$cluster`.
#'
#' @return A data frame containing the combined rows of the sampled clusters.
#'
#' @importFrom dplyr filter sample_n bind_rows %>%
#'
#' @examples
#' \dontrun{
#' # Assuming 'my_data' has a 'cluster' column and 'my_counts' is a named vector
#' sampled_df <- sim.sample(data = my_data, 
#'                          n.clusters = 5, 
#'                          n.per.cluster = 10, 
#'                          cluster_counts = my_counts)
#' }
#'
#' @export
sim.sample <- function(data, n.clusters, n.per.cluster, cluster_counts){
  avail <- as.integer(names(cluster_counts)) %in% data$cluster
  clusters_chosen <- sample(names(cluster_counts)[avail], 
                            n.clusters, 
                            replace = F, 
                            prob = cluster_counts[avail] /nrow(data))
  pts <- data.frame()
  for (c in clusters_chosen){
    new <- data %>% 
      filter(cluster == c) %>%
      sample_n(n.per.cluster)
    pts <- bind_rows(pts, new)
  }
  pts
}



#' Limit and Recycle Cross-Validation Folds
#'
#' Checks if the number of unique folds in the data exceeds a maximum limit. 
#' If it does, it reassigns the fold identifiers by cycling through the range 
#' `1:max_folds` using modulo arithmetic.
#'
#' @param train.data A data frame containing the training data. Must include a 
#'   column named `fold`.
#' @param max_folds Integer. The maximum allowable number of folds.
#'
#' @return A data frame. If the original number of folds is less than or equal to 
#'   `max_folds`, the original data is returned. Otherwise, the `fold` column is 
#'   updated.
#'
#' @importFrom dplyr mutate %>%
#'
#' @export
set_max_folds <- function(train.data, max_folds) {
  n <- length(unique(train.data$fold))
  if (n <= max_folds){
    return(train.data)
  }else{
    return(train.data %>% 
             mutate(
               fold =  ceiling(fold %% max_folds + 1)
             ))  
  }
}

#' Estimate Variance for Clustered Sampling Without Replacement
#'
#' Calculates a variance estimate that combines between-cluster and within-cluster 
#' variance, weighted by cluster areas and adjusted for the fraction of the total 
#' area sampled.
#'
#' @param data A data frame containing the sampled observations. Must contain 
#'   columns `cluster` and the variable specified by `target`.
#' @param target Character string. The name of the column in `data` for which 
#'   the variance is being calculated.
#' @param cluster_areas A data frame containing the total areas (or weights) for 
#'   the clusters. Must contain columns `cluster` and `area`.
#'
#' @return A single numeric value representing the estimated variance.
#'
#' @details 
#' The variance is calculated as a weighted combination of:
#' 1. **Between-cluster variance**: Variance of the cluster means.
#' 2. **Within-cluster variance**: Variance of observations within clusters, 
#'    weighted by cluster area.
#' 
#' These are combined using the fraction of the total area sampled (`frac_sampled`) 
#' as a weighting factor: 
#' \code{between_var * (1 - frac_sampled) + within_var * frac_sampled}
#'
#' @importFrom dplyr group_by summarize ungroup pull left_join filter %>%
#' @importFrom stats var weighted.mean
#'
#' @export
cluster_wo_replacement_variance <- function(data, 
                                            target, 
                                            cluster_areas){
  # Calculate variance between cluster means
  between_var <- data %>% 
    group_by(cluster) %>% 
    summarize(m = mean(get(target))) %>% 
    ungroup() %>%
    summarize(x = var(m) / n()) %>%
    pull(x)
  
  # Calculate variance within clusters, weighted by area
  within_var <- data  %>% 
    group_by(cluster) %>% 
    summarize(v = var(get(target)), n = n()) %>%
    left_join(cluster_areas, by = 'cluster') %>% 
    ungroup() %>% 
    summarize(x = weighted.mean(v, area) / sum(n)) %>% 
    pull(x)
  
  frac_sampled <- (cluster_areas %>% 
                     filter(cluster %in% unique(data$cluster)) %>% 
                     summarize(x = sum(area)) %>% pull(x)
                   
  )/sum(cluster_areas$area)
  
  
  return(between_var*(1-frac_sampled) + within_var*frac_sampled )
}


coverage_vals <- c(seq(.2, .8, .2), .9, .95, .99)
coverage_cols <- paste('covered_', as.character(coverage_vals), sep = '')


add.coverage.cols <- function(data){
  data$p_val <- pt(data$t,  data$df)
  for (c in coverage_vals){
    name <- paste('covered_', as.character(c), sep = '')
    data[,name] <- (data$p_val > (.5 - c/2)) &  (data$p_val < (.5 + c/2))
  }
  return(data)
}




#' Estimate the mean of a population using the difference estimator and a 
#' random forest model, using a clustered without replacement design
#'
#' Calculates a variance estimate that combines between-cluster and within-cluster 
#' variance, weighted by cluster areas and adjusted for the fraction of the total 
#' area sampled.
#'
#' @param full_grd A data frame containing the covariates for the full mapped area. 
#' Must contain columns `cluster` and the variable specified by `target`.
#'   @param formula A formula object for the random forest model
#'   
#'     @param formula A formula object for the random forest model
#'     
#' @param train_data A data frame containing the sampled observations used for 
#' fitting the machine learning model. 
#' Must contain  columns `cluster` and the variable specified by `target`.
#' 
#' @param test_data A data frame containing the sampled observations used for 
#' estimating the model's error
#' Must contain  columns `cluster` and the variable specified by `target`.

#'   
#' @param target Character string. The name of the column in `data` for which 
#'   the variance is being calculated.
#' @param cluster_areas A data frame containing the total areas (or weights) for 
#'   the clusters. Must contain columns `cluster` and `area`.
#'
#' @return A dataframe of the estiamted target value, its variance and the estimated
#' model bias
#'
#' @details 
#' The variance is calculated as a weighted combination of:
#' 1. **Between-cluster variance**: Variance of the cluster means.
#' 2. **Within-cluster variance**: Variance of observations within clusters, 
#'    weighted by cluster area.
#' 
#' These are combined using the fraction of the total area sampled (`frac_sampled`) 
#' as a weighting factor: 
#' \code{between_var * (1 - frac_sampled) + within_var * frac_sampled}
#'
#' @importFrom dplyr group_by summarize ungroup pull left_join filter %>%
#' @importFrom stats var weighted.mean
#'
#' @export
estimate_gdif <- function(full_grd, formula, train_data, test_data, 
                          cluster_areas,
                          target = 'y' 
                          ){
  rf.mod <- ranger::ranger(formula, train_data)
  test_data$pred <- predict(rf.mod, test_data)$predictions
  test_data$err <- test_data$pred - test_data[, target]
  
  me <- mean(test_data$err)
  mean_ml_pred <- mean(predict(rf.mod, full_grd)$predictions)
  full_pred <-  mean_ml_pred- me
  
  between_var <- test_data %>% 
    group_by(cluster) %>% 
  summarize(m = mean(err)) %>% 
    ungroup() %>%
    summarize(x = var(m) / n()) %>%
    pull(x)
  
  within_var <- test_data  %>% 
    group_by(cluster) %>% 
    summarize(v = var(err), n = n()) %>%
    left_join(cluster_areas, by = 'cluster') %>% 
    ungroup() %>% 
    summarize(x = weighted.mean(v, area) / sum(n)) %>% 
    pull(x)
  
  
  frac_sampled <- (cluster_areas %>% 
    filter(cluster %in% unique(test_data$cluster)) %>% summarize(x = sum(area)) %>% pull(x)
    )/sum(cluster_areas$area)
  
  if (is.na(between_var)){
    var_est <- within_var
  }else{
    var_est <- between_var*(1-frac_sampled) + within_var*frac_sampled
  }
  
  return(data.frame(pred = full_pred, 
                    var_est = var_est,
                    mean_ml_pred = mean_ml_pred,
                    me = me))
}



#'Estimate the generalized difference estimator by using cross-fitting
#' @param full_grd A data frame containing the covariates for the full mapped area. 
#' Must contain columns `cluster` and the variable specified by `target`.
#' @param data A data frame containing the sample data covariates for the full mapped area. 
#' Must contain columns `cluster` and the variable specified by `target`.
#' 
#'   @param formula A formula object for the random forest model
#'   #' @param target Character string. The name of the column in `data` for which 
#'   the variance is being calculated.
#' @param cluster_areas A data frame containing the total areas (or weights) for 
#'   the clusters. Must contain columns `cluster` and `area`.
#'
#'   
x.fit.gdif <- function(full_grd, data, formula, n.folds, cluster_areas,
                       target=  'y'){
  data <-data %>% 
    mutate(fold = dense_rank(cluster)) %>% set_max_folds(n.folds)
  
  gdif_ests <- data.frame()
  
  for (f in unique(data$fold)){
    train_data <- data %>% filter(fold != f)
    test_data <- data %>% filter(fold == f)
    .r <- estimate_gdif(full_grd, formula, train_data, test_data, 
                        cluster_areas)
    gdif_ests <- bind_rows(gdif_ests, .r)
    gdif_ests$df <- length(unique(test_data$cluster)) -1
    
  }
  
  
  fe_sum <- meta.summaries(gdif_ests$pred, sqrt(gdif_ests$var_est), method = 'fixed')
  re_sum <- meta.summaries(gdif_ests$pred, sqrt(gdif_ests$var_est), method = 'random')
  
  
    
  full_est <- data.frame(pred = c(fe_sum$summary, 
                                  re_sum$summary, 
                                  mean(gdif_ests$pred)
                                  ),
                         se = c(fe_sum$se.summary, 
                                re_sum$se.summary, 
                                sqrt(var(gdif_ests$pred)  / nrow(gdif_ests))
                               ),
                         
                               
                         method = c('fixed.effects', 
                                    'random.effects', 
                                    'simple_avg'
                                    )
                         ) %>%
    bind_rows(
      data.frame(pred = gdif_ests$pred, se = sqrt(gdif_ests$var_est), 
                         method = 'partial'))
  
  if (min(gdif_ests$df >2)){
    re_df_adj <- meta.summaries(gdif_ests$pred, 
                                sqrt(gdif_ests$var_est * (gdif_ests$df/(gdif_ests$df-2)) ), method = 'random')
    
    df_adj_est <- data.frame(pred =mean(gdif_ests$pred),
                              se =re_df_adj$se.summary,
                              method =  're.degree.of.freedom.adj'
                                          )
    
    full_est <- bind_rows(full_est, df_adj_est)
  }
  
  return(full_est %>% mutate(n.folds = n.folds))
}



#' Run and experiment on the gdif estimator using a random forest model,
#' cross-fitting, and a 2-stage cluster without-replacement design.
run.experiment <- function(dataset, formula, n.clusters.sampled,
                           n.sampled.per.cluster, n.folds, 
                           n.reps = 100){

counts <-  table(dataset$cluster)
cluster_areas <- data.frame(counts) %>% 
  rename(cluster = Var1, area = Freq) %>% 
  mutate(cluster = as.integer(cluster))






gdif.sim.res <- data.frame()
for (i in 1:n.reps){
  print(i)
  samps <- sim.sample(dataset, n.clusters.sampled, 
                      n.sampled.per.cluster, counts)
  
  
  db_var_est <- cluster_wo_replacement_variance(samps, 'y', cluster_areas)
  db_pred <- mean(samps$y)
    
  sim.fun <- function(n.fold){
    x.fit.gdif(dataset, samps, formula, n.fold, cluster_areas)}
  
  new_data <- do.call('rbind',
                      lapply(n.folds, sim.fun)
                      )
  
  
  
  clustered_inbag <- cluster_bootstrap_for_rf(samps, 'cluster', 500)
  
  rf.mod <- ranger(formula, samps, inbag =clustered_inbag)
  pred_full <- mean(predict(rf.mod, dataset)$predictions)
  
  samps <- samps %>%
    mutate(pred = rf.mod$predictions, err_oob = pred - y)
  
  pred_oob_gdif <- pred_full - mean(samps$err_oob)
  var_oob_gdif <- cluster_wo_replacement_variance(samps, 'err_oob', cluster_areas)
  
  new_data <-   bind_rows(new_data,
 data.frame(pred = c(db_pred, pred_oob_gdif),
                     se = c(sqrt(db_var_est), sqrt(var_oob_gdif)), 
                     method = c('db', 'oob_err_gdif'), 
                     n.folds = n.clusters.sampled
            )
 ) %>% mutate(run_i = i)
  
  gdif.sim.res <- gdif.sim.res %>% 
    bind_rows(new_data)
}

gdif.sim.res %>%
  mutate(true.mean = mean(dataset$y),
         err = pred - true.mean, 
         t = err / se,
         df = n.folds - 1, 
         n.clusters = n.clusters.sampled) %>% 
  add.coverage.cols()

}


#This experiment takes about 1/2 an hour to run on my machine, simpler example

amazonia_vars <- c('SWIR2', 'Terra_PP', 'Prec_dm', 'Elevation', 'Clay', 'x1', 'x2')
#grdAmazonia <- grdAmazonia %>% mutate(y = AGB)
formula.amaz <- as.formula(paste('y ~ ', paste(amazonia_vars, collapse = ' + '), sep = '' ) )
experiment <- 'amazonia.agb'
load('data/grd_resampled.rda')


test <- run.experiment(grd_resampled, 
                       formula.amaz, 
                       n.clusters.sampled = 8,
                       n.sampled.per.cluster = 12, 
                       n.folds = c(2,4,8), 
                       n.reps = 2)


res_amaz <- run.experiment(grd_resampled, 
                           formula.amaz, 
                           n.clusters.sampled = 20,
                           n.sampled.per.cluster = 12, 
                           n.folds = c(2,5,10,20))




print('Coverage summary')
res_amaz %>% 
  filter(method %in% c('db', 'simple_avg', 'oob_err_gdif' )) %>%
  group_by(method, n.folds) %>% 
  summarize_at(vars(coverage_cols), mean)
  

print('Performance summary')
res_amaz %>%  
  filter(method %in% c('db', 'simple_avg', 'oob_err_gdif' )) %>%
  group_by(method, n.folds) %>%
  summarize(bias = mean(err), variance_of_errs = var(err), mean_variance_of_estimator = mean(se**2),
            var_t = var(t), expected_var_t = mean(df /(df-2))) 

write.csv(res_amaz, 'results/amazon_agb_exp.csv')

#This experiment takes longer to run... more datapoints in the prediction set
#and a larger 
#number of covariates in the dataset

load('data/grd_ocs_eu.rda')
grd.ocs.eu$y <- grd.ocs.eu$ocs
ocs.vars <- colnames(grd.ocs.eu)[4:22]
formula.ocs.eu <- as.formula(paste('y ~', paste(ocs.vars, collapse = ' + '), sep = ' '))

res_ocs_eu <- run.experiment(grd.ocs.eu, 
                             formula.ocs.eu, 
                             20,
                             12, 
                             n.folds = c(2,5,10,20))


print('Coverage summary')
res_ocs_eu %>% 
  filter(method %in% c('db', 'simple_avg', 'oob_err_gdif' )) %>%
  group_by(method, n.folds) %>% 
  summarize_at(vars(coverage_cols), mean)

print('Performance summary')
res_ocs_eu %>%  
  filter(method %in% c('db', 'simple_avg', 'oob_err_gdif' )) %>%
  group_by(method, n.folds) %>%
  summarize(bias = mean(err), variance_of_errs = var(err), mean_variance_of_estimator = mean(se**2),
            var_t = var(t), expected_var_t = mean(df /(df-2))) 

write.csv(res_ocs_eu, 'results/ocs_eu_exp.csv')
