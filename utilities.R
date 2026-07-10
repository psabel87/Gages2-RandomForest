require(tidyverse)
require(zoo)
require(modifiedmk)

#negate %in%
`%nin%` <- negate(`%in%`)

#read info for selected gages. Headwaters selected in script a_GageSelection_2,
#downstream gages and headwater-downstream connections selected in script a_GageSelection_3,
#matched downstream gages selected in script a_GageSelection_4,
#all gage info compiled in a_GageSelection_5
read_gage_info <- function(type = 'headwaters'){
  #type = c('headwaters', 'downstream', 'connections')
  files <- c('headwaters' = 'hw_gage_info_exp.csv',
             'downstream' = 'ds_gage_info_exp.csv',
             'connections' = 'hw_ds_connections_exp.csv',
             'selected' = 'selected_sites.csv',
             'downstream_matched' = 'ds_matched_gage_info_exp.csv',
             'all' = 'all_gage_info_exp.csv')
  
  info_path <- file.path('data', 'gages', files[names(files) == type])
  
  gage_info <- read_csv(info_path)
  
  return(gage_info)
}

#quickly save the most recent plot
qsave <- function(w = 7,h = 3.5){
  return(ggsave('figures/temp.png',width = w, height = h,dpi = 700))
}

#gather individual .csv's for each gage into one. Assumes gages are identified
#by 'site_no'
gather_loose <- function(path){
  files <- list.files(path, full.names = T)
  out <- map(files, ~read_csv(.x, col_types = cols(site_no = 'c'))) %>%
    list_rbind()
  return(out)
}

#calculate modified trend statistics using the 'modifiedmk' library
trendinator <- function(x, length_thresh = 5){
  x_len <- sum(!is.na(x) & !is.nan(x))
  
  if(x_len < length_thresh){
    return(data.frame(tau = NA, sen = NA, p = NA, p0 = NA,
                      nn = NA, s = NA, s0 = NA))
  }
  
  trends <- tryCatch(mmkh3lag(x), error = function(e) NA) %>%
    round(., 5)
  
  if(length(trends) == 1){
    return(data.frame(tau = NA, sen = NA, p = NA, p0 = NA,
                      nn = NA, s = NA, s0 = NA))
  }

  out <- data.frame(t(trends)) %>%
    select(tau = Tau, sen = Sen.s.slope, p = new.P.value, p0 = old.P.value,
           nn = N.N., s = new.variance, s0 = old.variance)
  return(out)
}

#Classify trends based on significance and direction
trend_classifier <- function(x, p, alpha = 0.05){
  classes <- case_when(
    is.na(x) ~ NA,
    p > alpha ~ 'none',
    x > 0 & p <= alpha ~ 'pos',
    x < 0 & p <= alpha ~ 'neg',
    x == 0 ~ 'none',
    .default = NA)
  return(classes)
}

#false discovery rate for multiple hypothesis testing
#Wilks 2016
get_fdr_p <- function(p, fdr_a = 0.1){
  p = p[!is.na(p)]
  p_rank = sort(p)
  thresh = (1:length(p))/length(p)*fdr_a
  fdr_p = max(p_rank[p_rank <= thresh])
  return(fdr_p)
}

#classifier performance metrics
class_performance <- function(pred, obs, wide = F){
  m <- table(pred, obs)
  
  n <- sum(m)
  prevalence <- colSums(m)
  correct <- diag(m)
  
  acc <- sum(correct)/n
  recalls <- correct/prevalence
  macro_recall <- mean(recalls, na.rm = T)
  
  out <- data.frame(
    var = c('acc', 'mrc', paste0('rc_',names(recalls))),
    val = c(acc, macro_recall, recalls)
  )
  if(wide == T) out <- pivot_wider(out, names_from = var, values_from = val)
  
  return(out)
}

#KGE
kge <- function(x,y, modified = T){
  r = cor(x,y,use = 'pairwise.complete')
  ux = mean(x, na.rm = T)
  uy = mean(y, na.rm = T)
  sx = sqrt(var(x, na.rm = T))
  sy = sqrt(var(y, na.rm = T))
  
  b = uy/ux
  g = sy/sx
  if(modified == T){g = (sy/uy)/(sx/ux)}
  
  t1 = (r-1)^2
  t2 = (b-1)^2
  t3 = (g-1)^2
  kge = 1-sqrt(t1+t2+t3)
  return(kge)
}

#R2
r2 <- function(pred, obs){
  pred_nas = is.na(pred)
  obs_nas = is.na(obs)
  nas = pred_nas | obs_nas
  pred = pred[!nas]
  obs = obs[!nas]
  ssr <- sum((obs-pred)^2, na.rm = T)
  obs_mean <- mean(obs, na.rm = T)
  sst <- sum((obs - obs_mean)^2, na.rm = T)
  
  r2 <- 1 - (ssr/sst)
  return(r2)
}

rmse <- function(pred, obs){
  se <- (pred-obs)^2
  rmse <- sqrt(mean(se, na.rm = T))
  return(rmse)
}

mae <- function(pred, obs){
  ae <- abs(pred-obs)
  mae <- mean(ae)
  return(mae)
}

nse <- function(pred, obs){
  if(length(unique(obs)) < 2) return(NA_real_)
  1 - sum((pred - obs)^2, na.rm = T) /
    sum((obs - mean(obs, na.rm = T))^2, na.rm = T)
}

perf_summary <- function(data, lev = NULL, model = NULL){
  pred = data[,'pred']
  obs = data[,'obs']
  
  rmse <- rmse(pred, obs)
  mae <- mae(pred, obs)
  r2 <- r2(pred, obs)
  kge <- kge(pred, obs)
  
  out <- c(rmse, mae, r2, kge)
  names(out) <- c('RMSE', 'MAE', 'R2', 'KGE')
  return(out)
}

#regression performance metrics
regress_performance <- function(pred, obs, wide = F){
  r2 <- r2(pred, obs)
  rmse <- rmse(pred, obs)
  pcor <- cor(pred, obs, use = 'pairwise.complete')
  
  out <- data.frame(
    var = c('r2', 'rmse', 'p_cor'),
    val = c(r2, rmse, pcor)
  )
  if(wide == T) out <- pivot_wider(out, names_from = var, values_from = val)
  
  return(out)
}

region_recoder <- function(eco2){
  recodes <- c('nf' = 5.2,
               'nf' = 5.3,
               'mw' = 6.2,
               'mw' = 7.1,
               'ef' = 8.1,
               'cp' = 8.2,
               'ef' = 8.3,
               'ap' = 8.4,
               'ef' = 8.5,
               'gp' = 9.2,
               'gp' = 9.3,
               'gp' = 9.4,
               'gp' = 9.5,
               'sw' = 10.1,
               'sw' = 10.2,
               'sw' = 11.1,
               'sw' = 12.1,
               'sw' = 13.1)
  match_order <- match(eco2, recodes)
  recode <- names(recodes)[match_order]
  return(recode)
}
theme_CMH <- function(...){
  theme_bw(base_size = 8,
           base_line_size = 0.25) +
    theme(
      legend.text = element_text(size = rel(0.9)),
      axis.text = element_text(size = rel(0.9)),
      plot.title=element_text(face="bold", size=rel(1)),
      panel.grid=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin=unit(c(1,1,1,1), "mm"),
      panel.background = element_rect(fill = fill_alpha('white',0)),
      plot.background = element_rect(fill = fill_alpha('white',0), color = alpha('white',0)),
      legend.background = element_rect(fill = fill_alpha('white',0), color = alpha('white',0)),
      strip.background=element_blank()
    )
}
theme_set(theme_CMH())
