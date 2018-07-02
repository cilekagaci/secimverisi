eixrc <- function(df, row_total_name, col_total_name, row_names, col_names, start_list = NULL,  ...) {
  # Remove dots from column names (it messes up with ei.MD.bayes library)
  col_total_name = col_total_name %>% map_chr(str_replace_all, '\\.', '_')
  row_total_name = row_total_name %>% map_chr(str_replace_all, '\\.', '_')
  row_names = row_names %>% map_chr(str_replace_all, '\\.', '_')
  col_names = col_names %>% map_chr(str_replace_all, '\\.', '_')
  
  names(df) = names(df) %>% map(str_replace_all, '\\.', '_')
  df =  df %>% as.data.frame
  
  slack <- df[, col_total_name] - df[, row_total_name]
  df$row_slack <- pmax(0, slack)
  df$col_slack <- pmax(0, -slack)
  
  df$row_others <- df[, row_total_name] - rowSums(df[, row_names])
  df$col_others <- df[, col_total_name] - rowSums(df[, col_names])
  if (any(df$row_others < 0)) {stop(paste0('Row columns sum to more than ', row_total_name))}
  if (any(df$col_others < 0)) {stop(paste0('Col columns sum to more than ', col_total_name))}
  
  if(any(rowSums(df[, row_names]) + df$row_slack + df$row_others != rowSums(df[, col_names]) + df$col_slack + df$col_others)) {stop('Error :/')}
  
  cols = paste(c(col_names, 'col_slack', 'col_others'), collapse = ',')
  rows = paste(c(row_names, 'row_slack', 'row_others'), collapse = ',')
  
  f = as.formula(paste0('cbind(', cols, ')~cbind(', rows, ')'))
  
  model <- ei.MD.bayes(f, data=df, verbose=10000, start_list, ...)
  return(model)
}


do.tune <- function(df, row_total_name, col_total_name, row_names, col_names, start_list = NULL, lambda1 = 4, lambda2 = 2) {
  # Remove dots from column names (it messes up with ei.MD.bayes library)
  col_total_name = col_total_name %>% map_chr(str_replace_all, '\\.', '_')
  row_total_name = row_total_name %>% map_chr(str_replace_all, '\\.', '_')
  row_names = row_names %>% map_chr(str_replace_all, '\\.', '_')
  col_names = col_names %>% map_chr(str_replace_all, '\\.', '_')
  
  names(df) = names(df) %>% map(str_replace_all, '\\.', '_')
  df =  df %>% as.data.frame
  
  slack <- df[, col_total_name] - df[, row_total_name]
  df$row_slack <- pmax(0, slack)
  df$col_slack <- pmax(0, -slack)
  
  df$row_others <- df[, row_total_name] - rowSums(df[, row_names])
  df$col_others <- df[, col_total_name] - rowSums(df[, col_names])
  if (any(df$row_others < 0)) {stop(paste0('Row columns sum to more than ', row_total_name))}
  if (any(df$col_others < 0)) {stop(paste0('Col columns sum to more than ', col_total_name))}
  
  if(any(rowSums(df[, row_names]) + df$row_slack + df$row_others != rowSums(df[, col_names]) + df$col_slack + df$col_others)) {stop('Error :/')}
  
  cols = paste(c(col_names, 'col_slack', 'col_others'), collapse = ',')
  rows = paste(c(row_names, 'row_slack', 'row_others'), collapse = ',')
  
  f = as.formula(paste0('cbind(', cols, ')~cbind(', rows, ')'))
  tune_nocov <- tuneMD(f, data = df, ntunes = 5, totaldraws = 5000, verbose = 1000, lambda1=lambda1, lambda2=lambda2)
  return (tune_nocov)
}


mcmc2matrix <- function(cc, prefix) {
  get2 <- function(x) x[2]
  tnames <- strsplit(colnames(cc), paste0(prefix, "."))
  idx <- strsplit(sapply(tnames, get2), ".", fixed = TRUE)
  idx <- as.list(as.data.frame(matrix(unlist(idx), byrow = TRUE, 
                                      nrow = length(idx), ncol = length(idx[[1]]))))
  idx <- lapply(idx, as.character)
  idx <- lapply(idx, unique)
  mcpar.cc <- mcpar(cc)
  if (prefix == 'beta') {
    sims <- nrow(cc)
    cc <- array(t(cc), c(sapply(idx, length), sims), dimnames = list(idx[[1]], idx[[2]], idx[[3]], 1:sims))
  } else if(prefix == 'alpha') {
    cc <- array(t(cc), sapply(idx, length), dimnames = list(idx[[1]], idx[[2]]))
  } else if(prefix == 'ccount') {
    sims <- nrow(cc)
    cc <- array(t(cc), c(sapply(idx, length), sims), dimnames = list(idx[[1]], 
                                                                     idx[[2]], 1:sims))
  }
  cc
}


get_start_list <- function(model) {
  b = mcmc2matrix(model$draws$Beta, 'beta')
  dims = dim(b)[1:3]
  start.betas = b[,,,dim(b)[4]] %>% as.matrix
  #start.betas = b[,,,1] %>% as.matrix
  start.alphas = mcmc2matrix(model$draws$Alpha, 'alpha')
  dim(start.betas) = dims
  names(dim(start.betas)) = NULL
  list(start.alphas, start.betas)
}


run_steps <- function(df, burnin, thin, sample, start_list = NULL, row_names, col_names, row_total_name, col_total_name, ...) {
  model = eixrc(df, 
                row_total_name = row_total_name,
                col_total_name = col_total_name,
                row_names = row_names,
                col_names = col_names,
                burnin = burnin,
                thin = thin,
                sample = sample,
                start.list = start_list,
                ret.mcmc = TRUE,
                ...
  )
  return(model)
}


lambda.MD2 <- function (object, columns, rows, ret.mcmc = TRUE) 
{
  rows = rows %>% map_chr(str_replace_all, '\\.', '_')
  columns = col_names %>% map_chr(str_replace_all, '\\.', '_')
  idx = list(rows = rows, columns = columns)
  NG <- length(idx[[1]])
  NP <- length(idx[[2]])
  NI <- length(columns)
  cc = object$draws$Cell.counts
  sims = sims <- dim(cc)[3]
  lambda.out <- array(NA, c(NG, NI, sims), dimnames = list(idx[[1]], 
                                                           columns, 1:sims))
  for (ii in idx[[1]]) {
    lambda.out[ii, , ] <- t(t(cc[ii, columns, ])/apply(cc[ii, 
                                                          columns, ], 2, sum))
  }
  if (ret.mcmc) {
    lambda.out <- t(matrix(lambda.out, NG * NI, sims))
    colnames(lambda.out) <- paste("lambda", matrix(rep(idx[[1]], 
                                                       NI), NG, NI), matrix(rep(columns, NG), NG, NI, byrow = TRUE), 
                                  sep = ".")
    lambda.out <- coda::mcmc(lambda.out)
    attr(lambda.out, "mcpar") <- mcpar.cc
  }
  class(lambda.out) <- c("lambdaMD", class(lambda.out))
  lambda.out
}


run_simulation <- function(df, prefix,
                           row_names, col_names,
                           row_total_name, col_total_name,
                           tune = FALSE,
                           starting_step = NULL, burnin = 5000, thin = 10, sample = 1000, save_every = 10, major_steps = 20,
                           lambda1 = 4,
                           lambda2 = 2) {
  model = NULL
  steps = 0
  start_list = NULL
  draws_list = NULL
  if (!is.null(starting_step)) {
    start_list_file_name = sprintf('%s_step%d.Rdata', prefix, starting_step)
    load(start_list_file_name)
    steps = starting_step
  }
  tune_list = NULL
  if (tune == TRUE) {
    print('tuning: yes')
    tune_list = do.tune(df, start_list = start_list,
                        row_names = row_names,
                        col_names = col_names,
                        row_total_name = row_total_name,
                        col_total_name = col_total_name,
                        lambda1 = lambda1, lambda2 = lambda2)
  }
  for (major_step in 1:major_steps) {
    model = run_steps(df, burnin = burnin, thin = thin, sample = sample, start_list = start_list,
                      row_names = row_names,
                      col_names = col_names,
                      row_total_name = row_total_name,
                      col_total_name = col_total_name,
                      tune.list = tune_list,
                      lambda1 = lambda1,
                      lambda2 = lambda2)
    steps = steps + burnin + thin * sample
    fname = sprintf('%s_step%d.Rdata', prefix, steps)
    
    start_list = get_start_list(model)
    draws = lambda.MD(model, columns = c('col_others', 'col_slack', col_names %>% map_chr(str_replace_all, '\\.', '_')))
    if (is.null(draws_list)) {
      draws_list = draws
    } else {
      draws_list = rbind(draws_list, draws)
    }
    
    save(start_list, draws_list, file = fname)
    print(draws)
    print(fname)
    if ((major_step %% save_every) == 0) {
      fname = sprintf('%s_model%d.Rdata', prefix, steps)
      save(model, file = fname)
    }
  }
  draws_list
}


# Analysis and convergence checks

# Plot MCMC parameter estimates as a function of step number. Used to check if
# MCMC simulation has converged or not.
plot_district_level <- function(pattern) {
  fname = (Sys.glob(pattern) %>% sort(decreasing = TRUE))[1]
  print(fname)
  load(fname)
  
  dfx = draws_list %>% as_tibble %>%
    mutate(n = row_number()) %>%
    gather(key = metric, value = value, matches('lambda')) %>%
    separate(metric, into = c('dummy', 'row', 'col'), sep = '\\.') %>%
    filter(!col %in% c('col_slack')) %>%
    filter(!col %in% c('invalid_vote_count_2018meclis')) %>%
    filter(!row %in% c('row_slack', 'invalid_vote_count_Kasim'))
  dfx %>%
    #filter(col != 'col_others') %>%
    ggplot(aes(x = n, y = value, 
               color = col, group = col)) +
    geom_line(size = 1) +
    facet_wrap(~row, nrow = 1) +
    theme_bw() +
    scale_y_continuous(minor_breaks = seq(0, 1, 0.05)) +
    scale_color_brewer(type = 'qual', palette = 'Paired')
}