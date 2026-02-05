library(dplyr)
library(ranger)
library(rmeta)


sim.sample <- function(data, n.clusters, n.per.cluster, counts){
  avail <- as.integer(names(counts)) %in% data$cluster
  clusters_chosen <- sample(names(counts)[avail], n.clusters, replace = F, prob = counts[avail] /nrow(data))
  pts <- data.frame()
  for (c in clusters_chosen){
    new <- data %>% 
      filter(cluster == c) %>%
      sample_n(n.per.cluster)
    pts <- bind_rows(pts, new)
  }
  pts
}


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

cluster_wo_replacement_variance <- function(data, target, cluster_areas){
  between_var <- data %>% 
    group_by(cluster) %>% 
    summarize(m = mean(get(target))) %>% 
    ungroup() %>%
    summarize(x = var(m) / n()) %>%
    pull(x)
  
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


estimate_gdif <- function(grd, formula, train_data, test_data, 
                          cluster_areas,
                          type=  'cluster', target = 'y' 
                          ){
  rf.mod <- ranger::ranger(formula, train_data)
  test_data$pred <- predict(rf.mod, test_data)$predictions
  test_data$err <- test_data$pred - test_data[, target]
  
  me <- mean(test_data$err)
  mean_ml_pred <- mean(predict(rf.mod, grd)$predictions)
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

x.fit.gdif <- function(grd, data, formula, n.folds, cluster_areas){
  data <-data %>% 
    mutate(fold = dense_rank(cluster)) %>% set_max_folds(n.folds)
  
  gdif_ests <- data.frame()
  
  for (f in unique(data$fold)){
    train_data <- data %>% filter(fold != f)
    test_data <- data %>% filter(fold == f)
    .r <- estimate_gdif(grd, formula, train_data, test_data, 
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


run.experiment <- function(dataset, formula, n.clusters.sampled,
                           n.samples.per.cluster, n.folds, 
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
    
  sim.fun <- function(n.folds){
    x.fit.gdif(grd_resampled, samps, formula.amaz, 2, cluster_areas)}
  
  new_data <- do.call('rbind',
                      lapply(n.folds, sim.fun)
                      )
  
  new_data <-   bind_rows(new_data,
 data.frame(pred = db_pred, se = sqrt(db_var_est), method = 'db', n.folds = n.clusters.sampled)
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

amazonia_vars <- c('SWIR2', 'Terra_PP', 'Prec_dm', 'Elevation', 'Clay', 'x1', 'x2')
#grdAmazonia <- grdAmazonia %>% mutate(y = AGB)
formula.amaz <- as.formula(paste('y ~ ', paste(amazonia_vars, collapse = ' + '), sep = '' ) )
experiment <- 'amazonia.agb'
load('data/grd_resampled.rda')


test <- run.experiment(grd_resampled, 
                       formula.amaz, 
                       n.clusters.sampled = 8,
                       n.samples.per.cluster = 12, 
                       n.folds = c(2,4,8), n.reps = 2)


res_amaz <- run.experiment(grd_resampled, 
                           formula.amaz, 
                           n.clusters.sampled = 20,
                           n.samples.per.cluster = 12, 
                           n.folds = c(2,5,10,20))





res_amaz %>% 
  group_by(method, n.folds) %>% 
  summarize_at(vars(coverage_cols), mean) %>%
  print(n = 25)

res_amaz %>%  group_by(method, n.folds) %>%
  summarize(bias = mean(err), v_err = var(err), var_est = mean(se**2),
            v_t = var(t)) %>% 
  print(n = 50)





load('data/grd_ocs_eu.rda')
grd.ocs.eu$y <- grd.ocs.eu$ocs
ocs.vars <- colnames(grd.ocs.eu)[4:22]
formula.ocs.eu <- as.formula(paste('y ~', paste(ocs.vars, collapse = ' + '), sep = ' '))

res_ocs_eu <- run.experiment(grd.ocs.eu, 
                             formula.ocs.eu, 
                             20,
                             12, 
                             n.folds = c(2,5,10,20))
