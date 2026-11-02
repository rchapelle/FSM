library('opencesp')
library('randomForestSRC')
library('AppliedPredictiveModeling')
library('dplyr')
library('tidyr')
library('fastDummies')
library('here')
library('mice')
library('MLDataR')
library('kSamples')

n_sims <- 100

n_rand_orders <- 100
n_synth <- 100
n_rand_ru <- 100
n_synth_ru <- 100
n_lvls_ru <- 4
n_lvls_ru_fcs <- 5

synth_mb <- 50
pMSE_md <- 8

def_ctrl <- rpart::rpart.control(minbucket = 1, minsplit = 2, cp = -1, xval = 0, maxsurrogate = 0, maxcompete = 0)
get_def_ctrl <- function(mb, md = 30) {
  ctrl <- def_ctrl
  ctrl$minbucket <- mb
  ctrl$maxdepth <- md
  ctrl
}

bn_restart <- 100
bn_perturb <- 100

n_cores <- 10

is_categorical <- function(vec) {
  return(length(unique(vec)) <= 4)
}

is_discrete <- function(vec) {
  return(!is_categorical(vec) && length(unique(vec)) <= length(vec)/10)
}

is_continuous <- function(vec) {
  return(!is_categorical(vec) && !is_discrete(vec))
}

syn_cart <- function(df, visit_sequence, mb, md, m, phs = TRUE) { #phs = post hoc smoothing
  df <- df[, visit_sequence]

  res_l <- list()
  mods_l <- list()
  bw_vec <- numeric(ncol(df))
  prec_vec <- numeric(ncol(df))
  if(is_continuous(df[, 1])) {
    prec_vec[1] <- get_prec(df[, 1])
    bw_vec[1] <- bw.SJ(df[, 1])
  }
  for(j in 2:ncol(df)) {
    if(md > 0) mods_l[[j]] <- rpart::rpart(as.formula(paste0(colnames(df)[j],"~",".")), data = df[, 1:j], control = get_def_ctrl(mb, md))
    else mods_l[j] <- list(NULL)
    if(is_continuous(df[, j])) {
      prec_vec[j] <- get_prec(df[, j])
      bw_vec[j] <- bw.SJ(df[, j])
    }
  }

  for(i in 1:m) {
    synth <- df[sample(nrow(df), nrow(df), TRUE), ]
    synth[, 1] <- csample(df[, 1], size = nrow(df), replace = TRUE)
    if(is_continuous(df[, 1]) && !phs) {
      synth[, 1] <- round(synth[, 1] + rnorm(nrow(df), sd = bw_vec[1]), digits = prec_vec[1])
    }
    for(j in 2:ncol(df)) {
      if(!is.null(mods_l[[j]])) {
        cart_mod <- mods_l[[j]]
        preds_nodes <- treeClust::rpart.predict.leaves(cart_mod, newdata = synth[, 1:j], type = "where")
        synth[, j] <- sapply(preds_nodes, function(x) csample(df[[j]][which(cart_mod$where == x)], 1))
      }
      else synth[, j] <- csample(df[, j], size = nrow(df), replace = TRUE)

      if(is_continuous(df[, j]) && !phs) {
        synth[, j] <- round(synth[, j] + rnorm(nrow(df), sd = bw_vec[j]), digits = prec_vec[j])
      }
    }
    if(phs) {
      for(j in 1:ncol(df)) {
        if(is_continuous(df[, j])) {
          synth[, j] <- round(synth[, j] + rnorm(nrow(df), sd = bw_vec[j]), digits = prec_vec[j])
        }
      }
    }
    res_l[[i]] <- synth[, order(visit_sequence)]
  }
  return(res_l)
}

set_type <- function(x, keep_integer = FALSE) {
  if(is_categorical(x)) as.factor(x)
  else if (keep_integer && is_discrete(x)) as.integer(x)
  else as.numeric(x)
}

set_types <- function(df, keep_integer = FALSE) {
  orig_names <- colnames(df)
  df <- data.frame(lapply(df, function(x) set_type(x, keep_integer)))
  colnames(df) <- orig_names
  df
}

set_num <- function(df) {
  orig_names <- colnames(df)
  df <- data.frame(lapply(df, function(x) {
    if(is.factor(x) && length(levels(x)) == 2) {
      levels(x) <- c(0, 1)
      as.numeric(levels(x))[x]
    }
    else as.numeric(x)
  }))
  colnames(df) <- orig_names
  df
}

score_cart <- function(y, x, mb) {
  if(is_categorical(y)) {
    if (is.null(x)) {
      freq <- table(y)
      p <- freq / length(y)
      log_p <- log(p[as.character(y)])
      return(sum(log_p))
    }
    data <- data.frame(y, x)
    fit <- rpart::rpart(y ~ ., data = data, control = get_def_ctrl(mb))
    leaves <- fit$where
    tab <- table(leaves, y)
    leaves_f <- factor(leaves, levels = rownames(tab))
    y_factor <- factor(y, levels = colnames(tab))
    nsim <- tab[cbind(as.numeric(leaves_f), as.numeric(y_factor))]
    tot_leaves <- rowSums(tab)[as.numeric(leaves_f)]
    return(sum(log(nsim) - log(tot_leaves)))
  }
  else {
    if(is_continuous(y)) densbw <- bw.SJ(y)
    else densbw <- bw.nrd0(y)
    prec <- get_prec(y)
    if (is.null(x)) {
      return(sum(log(sapply(y, function(val) {
        mean(dnorm(val, y, sd = densbw))
      }))))
    }
    data <- data.frame(y, x)
    fit <- rpart::rpart(formula = y ~ ., data = data, control = get_def_ctrl(mb))
    leaves <- factor(fit$where)
    donors <- split(y, leaves)
    unique_donors <- lapply(donors, unique)
    freqs_list <- lapply(donors, function(d) tabulate(match(d, unique(d))))
    return(sum(log(sapply(seq_along(y), function(i) {
      leaf_name <- leaves[i]
      ud <- unique_donors[[leaf_name]]
      ff <- freqs_list[[leaf_name]]
      sum((ff / sum(ff)) * dnorm(y[i], ud, sd = densbw))
    }))))
  }
}


process_order_sim <- function(n_vars_vec, n_iter, n_cores) {
  df_order <- data.frame(n_var = n_vars_vec)
  df_order$cor_fit <- NA
  df_order$cor_fitcart <- NA
  df_order$cor_bn <- NA
  df_order$cor_vines <- NA
  df_order$time_fit <- NA
  df_order$time_fitcart <- NA
  df_order$time_bn <- NA
  df_order$time_vines <- NA

  score_lin <- function(y, x) {
    if(is.null(x)) resids <- y-mean(y)
    else {
      df <- cbind.data.frame(y, x)
      fit <- lm(y ~ ., data = df)
      resids <- fit$residuals
    }
    sd_hat <- sd(resids)
    return(sum(dnorm(resids, sd = sd_hat, log = TRUE)))
  }

  future::plan(future_strategy(), workers = n_cores)
  n_vec <- get_chunks(n_iter, n_cores)

  part_vines <- FALSE
  part_vines_step <- 0
  thr_time <- 100

  for(n_var in df_order$n_var) {
    print(paste0(n_var, "..."))

    process_iter <- function(n) {
      dist_mat_loc <- matrix(nrow = n, ncol = 8)
      for(j in 1:n) {
        sim <- opencesp::get_sim_data(500, n_var, func = rnorm,  sd = 0.1)
        target_order <- sample(1:n_var, n_var)
        sim <- sim[, target_order]
        target_order <- order(target_order)

        blocks <- NULL
        time_struct <- 0

        if(part_vines) {
          time_struct_vines <- system.time({
            blocks_vines <- opencesp::get_ind_blocks(sim, max_size = part_vines_step, n_cores = 1)
          })[3]
        }

        dist_mat_loc[j, 5] <- system.time({
          order_fit <- opencesp::get_dep_order(sim, blocks = NULL, method = "fit-hclust",  fit.n_cores = 1, fit.bf_first = FALSE, score_func = score_lin)
        })[3]

        dist_mat_loc[j, 6] <- system.time({
          order_fitcart <- opencesp::get_dep_order(sim, blocks = NULL, method = "fit-hclust",  fit.n_cores = 1, fit.bf_first = FALSE, score_func = score_cart, mb = synth_mb)
        })[3]

        dist_mat_loc[j, 7] <- system.time({
          order_bn <- opencesp::get_dep_order(sim, blocks = NULL, method = "bn",  bn.restart = bn_restart, bn.perturb = bn_perturb)
        })[3]

        dist_mat_loc[j, 8] <- system.time({
          order_vines <- opencesp::get_dep_order(sim, blocks = switch(2 - part_vines, blocks_vines, NULL), method = "vines", vines.n_cores = 1)
        })[3] + ifelse(part_vines, time_struct_vines, 0)

        dist_mat_loc[j, 1] <- cor(order_fit, target_order, method="kendall")
        dist_mat_loc[j, 2] <- cor(order_fitcart, target_order, method="kendall")
        dist_mat_loc[j, 3] <- cor(order_bn, target_order, method="kendall")
        dist_mat_loc[j, 4] <- cor(order_vines, target_order, method="kendall")
      }
      return(dist_mat_loc)
    }

    dist_mat <- do.call(rbind, future.apply::future_lapply(n_vec, process_iter, future.seed = TRUE))

    df_order[df_order$n_var == n_var, "cor_fit"] <- colMeans(dist_mat)[1]
    df_order[df_order$n_var == n_var, "cor_fitcart"] <- colMeans(dist_mat)[2]
    df_order[df_order$n_var == n_var, "cor_bn"] <- colMeans(dist_mat)[3]
    df_order[df_order$n_var == n_var, "cor_vines"] <- colMeans(dist_mat)[4]
    df_order[df_order$n_var == n_var, "time_fit"] <- colMeans(dist_mat)[5]
    df_order[df_order$n_var == n_var, "time_fitcart"] <- colMeans(dist_mat)[6]
    df_order[df_order$n_var == n_var, "time_bn"] <- colMeans(dist_mat)[7]
    df_order[df_order$n_var == n_var, "time_vines"] <- colMeans(dist_mat)[8]

    if(!part_vines && df_order[df_order$n_var == n_var, "time_vines"] > thr_time) {
      part_vines <- TRUE
      part_vines_step <- n_var
    }
  }

  future::plan(future::sequential)
  return(list(df_order = df_order, vines_step = part_vines_step))
}

process_dep_struct <- function(df, max_sizes, n_cores) {
  res_l <- list()
  for(max_size in max_sizes) {
    new_el <- list(max_size = max_size)
    new_el$blocks <- get_ind_blocks(df, maxs_size = max_size, n_cores = n_cores)
    res_l[[length(res_l) + 1]] <- new_el
  }
  return(res_l)
}

parallel_syn_orders <- function(df, orders, n_cores, metrics = c("PCD", "pMSE", "Hellinger", "PMR"), ...) {
  orig_names <- colnames(df)
  df <- set_types(df)
  colnames(df) <- paste0("var.", colnames(df))
  ellipsis_args <- list(...)

  orig_num <- set_num(df)

  f <- function(orders_subset) {
    res_synth <- list()
    for (i in seq_along(orders_subset)) {
      order <- orders_subset[[i]]

      synth <- syn_cart(df, order, ...)

      df_res <- setNames(data.frame(matrix(ncol = length(metrics), nrow = length(synth))), metrics)

      for(j in 1:length(synth)) {
        synth_num <- set_num(synth[[j]])
        if("Hellinger" %in% metrics) df_res$Hellinger[j] <- mean_hellinger(df, synth[[j]])
        if("PCD" %in% metrics) df_res$PCD[j] <- get_PCD_num(orig_num, synth_num)
        if("pMSE" %in% metrics) df_res$pMSE[j] <- pMSE(df, synth[[j]], control = get_def_ctrl(mb = 1, md = pMSE_md))
        if("PMR" %in% metrics) df_res$PMR[j] <- get_matches_prop(orig_num, synth_num, method = "euclidean", thr= 0.02)
      }

      res_synth[[length(res_synth) + 1]] <- df_res
    }
    return(res_synth)
  }

  future::plan(future_strategy(), workers = n_cores)
  split_orders <- split(orders, cut(seq_along(orders), n_cores, labels = FALSE))
  lists <- future.apply::future_lapply(split_orders, f, future.seed = TRUE)
  processed <- unlist(lists, recursive = FALSE)
  future::plan(future::sequential)
  return(processed)
}

process_util_var <- function(df_list, n_md) {
  res_list <- list()
  for(df in df_list) {
    p_vec <- numeric(n_md)
    for(i in 1:n_md) {
      syn_res <- parallel_syn_orders(df, get_rand_orders(ncol(df), n_rand_orders), n_cores = n_cores, metrics = c("pMSE"), m = n_synth, mb = 1, md = i-1)
      p_vec[i] <- ad.test(lapply(syn_res, function(x) x$pMSE))$ad[1,3]
    }
    res_list[[length(res_list) + 1]] <- data.frame(md = 0:(n_md-1), p = p_vec)
  }
  return(res_list)
}

process_dataset_syn <- function(df, n_rand, n_synth, max_struct_cons, n_cores, metrics = c("PCD", "pMSE", "Hellinger", "PMR"), struct_method = "fit") {
  df <- set_types(df)
  names(df) <- seq_along(df)

  print("Analyzing dataset structure...")

  res <- list()

  if(max_struct_cons < ncol(df)) {
    blocks <- get_ind_blocks(df, max_size = max_struct_cons, n_cores = n_cores)
  }
  else blocks <- NULL

  if(struct_method == "fit") res[[1]] <- list(rand_order = FALSE, order = get_dep_order(df, blocks = blocks, method = "fit-hclust",  fit.n_cores = 1, fit.bf_first = ncol(df) < 15, n_bf = 1, score_func = score_cart, mb = synth_mb))
  else if(struct_method == "bn") res[[1]] <- list(rand_order = FALSE, order = get_dep_order(df, blocks = blocks, method = "bn",  bn.restart = bn_restart, bn.perturb = bn_perturb))
  else if(struct_method == "vines") res[[1]] <- list(rand_order = FALSE, order = get_dep_order(df, blocks = blocks, method = "vines", vines.n_cores = n_cores))

  print("Generating random visit sequences...")
  rand_orders <- get_rand_orders(ncol(df), n_rand)

  for(order in rand_orders) {
    res[[length(res) + 1]] <- list(order = order, rand_order = TRUE)
  }

  all_orders <- lapply(res, function(x) x$order)

  print("Starting CART synthesis...")

  cart_res <- parallel_syn_orders(df, all_orders, n_cores = n_cores, metrics = metrics, mb = synth_mb, md = 30, m = n_synth)

  for(i in seq_along(res)) {
    res[[i]]$cart_res <- cart_res[[i]]
  }

  print("Done.")
  return(res)
}

process_profile <- function(orig, res_list, n_rand, n_synth, n_lvls, n_cores, metrics = c("pMSE", "PMR")) {
  orig <- set_types(orig)
  names(orig) <- seq_along(orig)

  n_erand <- sum(sapply(res_list, function(x) x$rand_order))
  if(n_rand > n_erand) stop("Not enough random dataframes.")

  indexes <- c(1, tail(seq_along(res_list), n_rand))

  for(id in indexes) {
    res_list[[id]]$profile_mb <- list()
    res_list[[id]]$profile_res <- list()
  }

  orders <- lapply(res_list[indexes], function(x) x$order)
  buckets <- c(25, 50, 100, 200, 400)[1:n_lvls]
  if(nrow(orig) > 900) buckets <- buckets * 2

  for(i in seq_along(buckets)) {
    res_md <- parallel_syn_orders(orig, orders, n_cores = n_cores, metrics = metrics, m = n_synth, mb = buckets[i], md = 30)
    for(j in seq_along(indexes)) {
      res_list[[indexes[j]]]$profile_mb[[i]] <- buckets[i]
      res_list[[indexes[j]]]$profile_res[[i]] <- res_md[[j]]
    }
  }

  return(res_list)
}

process_all <- function(orig, vine_mod, filename) {
  orig <- set_types(orig)
  orig_num <- set_num(orig)
  res_list <- process_dataset_syn(orig, n_rand = n_rand_orders, n_synth = n_synth, max_struct_cons = Inf, n_cores = n_cores)
  res_list <- process_profile(orig, res_list, n_rand = n_rand_ru, n_synth = n_synth_ru, n_lvls = n_lvls_ru, n_cores = n_cores)

  if(vine_mod) {
    for(i in seq_along(res_list)) {
      res_list[[i]]$score_fit <- get_fit_score(orig, res_list[[i]]$order, score_cart, mb = synth_mb)
      res_list[[i]]$score_bn <- try({bnlearn::score(bn_from_order(orig_num, res_list[[i]]$order), data = orig_num, type = "loglik-g")}, silent = TRUE)
      if(inherits(res_list[[i]]$score_bn, "try-error")) res_list[[i]]$score_bn <- NA
      if(ncol(orig) <= 10) res_list[[i]]$score_vines <- vine_from_order(orig_num, res_list[[i]]$order, n_cores)$logLik
    }
  }
  saveRDS(res_list, here(paste0("examples/visit_order/", filename)))
}

process_hyperparams <- function(df, mb_vec, mb_ref) {
  df <- set_types(df)
  res_l <- list()
  t_vec <- numeric(length(mb_vec))
  p_vec <- numeric(length(mb_vec))
  for(i in seq_along(mb_vec)) {
    res_l[[i]] <- get_dep_order(df, score_func = score_cart, mb = mb_vec[i])
  }
  rand_orders <- get_rand_orders(ncol(df), n_rand_orders)
  ref_id <- which(mb_vec == mb_ref)
  for(i in seq_along(t_vec)) {
    t_vec[i] <- cor(res_l[[i]], res_l[[ref_id]], method="kendall")
    scores_vec <- sapply(rand_orders, function(x) get_fit_score(df, x, score_func = score_cart, mb = mb_vec[i]))
    ref_fit <- get_fit_score(df, res_l[[ref_id]], score_func = score_cart, mb = mb_vec[i])
    p_vec[i] <- sum(scores_vec > ref_fit)
  }
  return(list(mbs = mb_vec, orders = res_l, taus = t_vec, p_vals = p_vec))
}

process_targeted <- function(df, n_rand, subsets_l, filename) {
  df <- set_types(df)
  main_order <- get_dep_order(df, blocks = NULL, method = "fit-hclust",  fit.n_cores = n_cores, score_func = score_cart, mb = synth_mb)
  orders_l <- list()
  rand_subsets <- list()
  for(l in subsets_l) {
    for(i in 1:n_rand) {
      candidate <- sample(1:ncol(df), l)
      while(any(sapply(rand_subsets, function(x) identical(sort(x), sort(candidate))))) candidate <- sample(1:ncol(df), l)
      rand_subsets[[length(rand_subsets) + 1]] <- candidate
      orders_l[[length(orders_l) + 1]] <- get_dep_order(df, blocks = NULL, method = "fit-hclust",  fit.n_cores = n_cores, priority = candidate, score_func = score_cart, mb = synth_mb)
    }
  }

  res_pMSE_main <- lapply(1:length(rand_subsets), function(x) rep(NA, n_synth))
  res_pMSE_target <- lapply(1:length(rand_subsets), function(x) rep(NA, n_synth))

  synth_l_main <- syn_cart(df, visit_sequence = main_order, mb = synth_mb, md = 30, m = n_synth)

  for(i in seq_along(rand_subsets)) {
    synth_l <- syn_cart(df, visit_sequence = orders_l[[i]], mb = synth_mb, md = 30, m = n_synth)
    orig_rd <- df[, rand_subsets[[i]]]

    for(j in 1:n_synth) {
      main_rd <- synth_l_main[[j]][, rand_subsets[[i]]]
      target_rd <- synth_l[[j]][, rand_subsets[[i]]]
      res_pMSE_main[[i]][j] <- pMSE(orig_rd, main_rd, control = get_def_ctrl(mb = 1, md = pMSE_md))
      res_pMSE_target[[i]][j] <- pMSE(orig_rd, target_rd, control = get_def_ctrl(mb = 1, md = pMSE_md))
    }
  }

  res <- list(subsets = rand_subsets, main_order = main_order, orders_l = orders_l, pMSE_main = res_pMSE_main, pMSE_target = res_pMSE_target)
  saveRDS(res, here(paste0("examples/visit_order/", filename)))
}

process_errprop <- function(df, l_vec, n_subs, n_rand, filename) {
  df <- set_types(df)

  get_rand_perm_const <- function(d, fixed_p, fixed_vals) {
    not_fixed <- sample(setdiff(1:d, fixed_vals), d - length(fixed_p))
    res <- numeric(d)
    res[fixed_p] <- sample(fixed_vals)
    res[fixed_p] <- fixed_vals
    res[setdiff(1:d, fixed_p)] <- not_fixed
    return(res)
  }

  process_l <- function(l) {
    res <- list()
    for(i in 1:n_subs) {
      subset <- sample(1:ncol(df), l)
      possible_pos <- 1:(ncol(df)-length(subset) + 1)
      sub_list <- list()
      for(pos in possible_pos) {
        pMSE_vec <- numeric(n_rand)
        for(j in 1:n_rand) {
          rand_ord <- get_rand_perm_const(ncol(df), (1:length(subset))+(pos-1), subset)
          synth_l <- syn_cart(df, rand_ord, md = 30, mb = synth_mb, m = n_synth)
          pMSE_vec[j] <- mean(sapply(synth_l, function(x) pMSE(df[, subset], x[, subset], control = get_def_ctrl(1, pMSE_md))))
        }
        sub_list[[pos]] <- pMSE_vec
      }
      res[[i]] <- list(subset = subset, res = sub_list)
    }
    return(res)
  }

  res_l <- list()
  for(i in seq_along(l_vec)) {
    l <- l_vec[i]
    res_l[[i]] <- list(l = l, res = process_l(l))
  }

  saveRDS(res_l, here(paste0("examples/visit_order/", filename)))
}

process_errprop_phs <- function(df, l_vec, n_subs, n_rand, filename) {
  df <- set_types(df)

  get_rand_perm_const <- function(d, fixed_p, fixed_vals) {
    not_fixed <- sample(setdiff(1:d, fixed_vals), d - length(fixed_p))
    res <- numeric(d)
    res[fixed_p] <- sample(fixed_vals)
    res[fixed_p] <- fixed_vals
    res[setdiff(1:d, fixed_p)] <- not_fixed
    return(res)
  }

  process_l <- function(l) {
    res <- list()
    for(i in 1:n_subs) {
      subset <- sample(1:ncol(df), l)
      possible_pos <- c(1, (ncol(df)-length(subset) + 1))
      sub_list <- list()
      for(p in 1:(length(possible_pos) + 1)) {
        pMSE_vec <- numeric(n_rand)
        pos <- possible_pos[min(p, length(possible_pos))]
        for(j in 1:n_rand) {
          rand_ord <- get_rand_perm_const(ncol(df), (1:length(subset))+(pos-1), subset)
          synth_l <- syn_cart(df, rand_ord, md = 30, mb = synth_mb, m = n_synth, phs = p <= length(possible_pos))
          pMSE_vec[j] <- mean(sapply(synth_l, function(x) pMSE(df[, subset], x[, subset], control = get_def_ctrl(1, pMSE_md))))
        }
        sub_list[[p]] <- pMSE_vec
      }
      res[[i]] <- list(subset = subset, res = sub_list)
    }
    return(res)
  }

  res_l <- list()
  for(i in seq_along(l_vec)) {
    l <- l_vec[i]
    res_l[[i]] <- list(l = l, res = process_l(l))
  }

  saveRDS(res_l, here(paste0("examples/visit_order/", filename)))
}

process_blocks <- function(n_obs, n_vars, max_size_vec, filename) {
  get_sim_blocks_data <- function(max_block_size, func = rnorm, ...) {
    df_list <- list()
    blocks <- numeric(n_vars)
    n_sim <- 0
    block_id <- 1
    while(n_sim < n_vars) {
      n_add <- sample(1:min(c(max_block_size, n_vars - n_sim)), 1)
      df_list[[length(df_list) + 1]] <- get_sim_data(n_obs, n_add, func = func,  ...)
      blocks[(n_sim + 1):(n_sim + n_add)] <- block_id
      n_sim <- n_sim + n_add
      block_id <- block_id + 1
    }
    return(list(df = do.call(cbind, df_list), blocks = blocks))
  }

  get_block_prec <- function(b_ref, b_est) {
    ids <- expand.grid(i = seq_along(b_ref), j = seq_along(b_est))
    prec_vec <- numeric(nrow(ids))
    for(k in 1:nrow(ids)) {
      id <- ids[k, ]
      if(b_ref[id$i] != b_ref[id$j]) prec_vec[k] <- NA
      else {
        if(b_est[id$i] == b_est[id$j]) prec_vec[k] <- 1
        else prec_vec[k] <- 0
      }
    }
    return(mean(prec_vec, na.rm = TRUE))
  }

  process_sizes <- function(max_sizes) {
    res_vec <- numeric(length(max_sizes))
    for(i in 1:length(max_sizes)) {
      id_res <- numeric(n_sims)
      for(j in 1:n_sims) {
        sim_blocks <- get_sim_blocks_data(max_sizes[i], sd = 2)
        retrieved <- get_ind_blocks(sim_blocks$df, max_size = max_sizes[i], n_cores = 1)
        id_res[j] <- get_block_prec(sim_blocks$blocks, retrieved)
      }
      res_vec[i] <- mean(id_res)
    }
    return(res_vec)
  }

  future::plan(future_strategy(), workers = n_cores)
  split_sizes <- split(max_size_vec, cut(seq_along(max_size_vec), n_cores, labels = FALSE))
  lists <- future.apply::future_lapply(split_sizes, process_sizes, future.seed = TRUE)
  res_vec <- unlist(lists, recursive = FALSE)
  future::plan(future::sequential)

  saveRDS(list(max_sizes = max_size_vec, res = res_vec), here(paste0("examples/visit_order/", filename)))
}

process_FCS <- function(df, filename) {
  df <- set_types(df)
  df_num <- set_num(df)

  syn_FCS <- function(orig, m = 1, n_iter = 50, mb) {
    mods_l <- list()
    for(j in 1:ncol(orig)) {
      mods_l[[j]] <- rpart::rpart(as.formula(paste(colnames(orig)[j], "~ .")), data = orig, control = get_def_ctrl(mb))
    }
    res_l <- list()
    for(s_id in 1:m) {
      synth <- orig
      for(j in 1:ncol(orig)) {
        synth[[j]] <- csample(orig[[j]], nrow(orig), replace = TRUE)
      }
      for(it in 1:n_iter) {
        for(j in sample(1:ncol(orig))) {
          fit <- mods_l[[j]]
          preds_nodes <- treeClust::rpart.predict.leaves(fit, synth, type = "where")
          new <- sapply(preds_nodes, function(x) csample(orig[[j]][which(fit$where == x)], 1))
          synth[[j]] <- new
        }
      }
      res_l[[s_id]] <- synth
    }
    return(res_l)
  }

  buckets <- c(10, 25, 50, 100, 200)[1:n_lvls_ru_fcs]
  pMSE_list_seq <- list()
  pMSE_list_fcs <- list()
  PMR_list_seq <- list()
  PMR_list_fcs <- list()
  time_vec_seq <- numeric(length(buckets))
  time_vec_fcs <- numeric(length(buckets))
  for(i in seq_along(buckets)) {
    syn_order <- get_dep_order(df, fit.n_cores = 1, score_func = score_cart, mb = buckets[i])
    time_vec_seq[i] <- system.time({
      synth_l_seq <- syn_cart(df, visit_sequence = syn_order, mb = buckets[i], md = 30, m = n_synth)
    })[3]
    time_vec_fcs[i] <- system.time({
      synth_l_fcs <- syn_FCS(df, m = n_synth, mb = buckets[i])
    })[3]
    pMSE_vec_seq <- numeric(n_synth)
    pMSE_vec_fcs <- numeric(n_synth)
    PMR_vec_seq <- numeric(n_synth)
    PMR_vec_fcs <- numeric(n_synth)
    for(j in 1:n_synth) {
      pMSE_vec_seq[j] <- pMSE(df, synth_l_seq[[j]], control = get_def_ctrl(pMSE_md))
      pMSE_vec_fcs[j] <- pMSE(df, synth_l_fcs[[j]], control = get_def_ctrl(pMSE_md))
      PMR_vec_seq[j] <- get_matches_prop(df_num, set_num(synth_l_seq[[j]]), method = "euclidean", thr = 0.02)
      PMR_vec_fcs[j] <- get_matches_prop(df_num, set_num(synth_l_fcs[[j]]), method = "euclidean", thr = 0.02)
    }
    pMSE_list_seq[[i]] <- pMSE_vec_seq
    pMSE_list_fcs[[i]] <- pMSE_vec_fcs
    PMR_list_seq[[i]] <- PMR_vec_seq
    PMR_list_fcs[[i]] <- PMR_vec_fcs
  }

  saveRDS(list(mb = buckets,
               pMSE_seq = pMSE_list_seq,
               pMSE_fcs = pMSE_list_fcs,
               PMR_seq = PMR_list_seq,
               PMR_fcs = PMR_list_fcs,
               time_seq = time_vec_seq,
               time_fcs = time_vec_fcs), here(paste0("examples/visit_order/", filename)))
}

set.seed(1)

data("PreDiabetes")
pd <- PreDiabetes
pd$Time_Pre_To_Diabetes <- NULL
pd <- pd[round(pd$BMI, 5) != 33.30411, ]

load_data("wbcd")
wbcd <- impute_rf(wbcd)

data("heartdisease")
hd <- heartdisease
hd$RestingECG[hd$RestingECG %in% c("LVH", "ST")] <- "Abnormal"

data(AlzheimerDisease, package = "AppliedPredictiveModeling")
ad <- predictors
rm(predictors)
ad$diag <- diagnosis
rm(diagnosis)
orig_ad_modalities <- levels(ad$Genotype)
ad <- fastDummies::dummy_cols(ad, remove_first_dummy = TRUE, remove_selected_columns = TRUE)

pd <- set_types(pd)
wbcd <- set_types(wbcd)
hd <- set_types(hd)
ad <- set_types(ad)

set.seed(2)
res_sim_orders_norm <- process_order_sim(c(5:10, (6:12)*2, (3:5)*10), n_sims, n_cores)
saveRDS(res_sim_orders_norm, here("examples/visit_order/res_sim_orders_norm.rds"))

res_wbcd_struct <- process_dep_struct(wbcd, (1:4)*2, n_cores = n_cores)
saveRDS(res_wbcd_struct, here("examples/visit_order/res_wbcd_struct.rds"))

res_wbcd_ddr <- hclust_fit_order(wbcd, score_func = score_cart, mb = synth_mb, n_bf = 1)$hc
saveRDS(res_wbcd_ddr, here("examples/visit_order/res_wbcd_ddr.rds"))

set.seed(3)
res_util_var <- process_util_var(list(pd, wbcd, hd), 6)
saveRDS(res_util_var, here("examples/visit_order/res_util_var.rds"))

set.seed(4)
process_all(pd, TRUE, "res_pd.rds")

set.seed(5)
process_all(wbcd, TRUE, "res_wbcd.rds")

set.seed(6)
process_all(hd, TRUE, "res_hd.rds")

set.seed(7)
time_p <- system.time({
  res_ad <- process_dataset_syn(ad, n_rand = n_rand_orders, n_synth = n_synth, max_struct_cons = Inf, n_cores = n_cores, struct_method = "fit")
})[3]
res_ad[[1]]$time <- time_p
for(i in 2:4) {
  time_p <- system.time({
    res_ad_red <- process_dataset_syn(ad, n_rand = 0, n_synth = n_synth, max_struct_cons = round(ncol(ad) / i), n_cores = n_cores, struct_method = "fit")
  })[3]
  res_ad_red[[1]]$max_cons <- round(ncol(ad) / i)
  res_ad_red[[1]]$time <- time_p
  res_ad[[length(res_ad) + 1]] <- res_ad_red[[1]]
}
ad_num <- set_num(ad)
for(i in seq_along(res_ad)) {
  res_ad[[i]]$score_fit <- get_fit_score(ad, res_ad[[i]]$order, score_cart, mb = synth_mb)
  res_ad[[i]]$score_bn <- try({bnlearn::score(bn_from_order(ad_num, res_ad[[i]]$order), data = ad_num, type = "loglik-g")}, silent = TRUE)
  if(inherits(res_ad[[i]]$score_bn, "try-error")) res_ad[[i]]$score_bn <- NA
}
saveRDS(res_ad, here("examples/visit_order/res_ad.rds"))

set.seed(8)
res_ad_bn <- process_dataset_syn(ad, n_rand = 0, n_synth = n_synth, max_struct_cons = Inf, n_cores = n_cores, struct_method = "bn")
res_ad_bn[[1]]$score_bn <- bnlearn::score(bn_from_order(ad_num, res_ad_bn[[1]]$order), data = ad_num, type = "loglik-g")
saveRDS(res_ad_bn, here("examples/visit_order/res_ad_bn.rds"))

set.seed(9)
res_pd_bn <- process_dataset_syn(pd, n_rand = 0, n_synth = n_synth, max_struct_cons = Inf, n_cores = n_cores, struct_method = "bn")
pd_num <- set_num(pd)
res_pd_bn[[1]]$score_bn <- bnlearn::score(bn_from_order(pd_num, res_pd_bn[[1]]$order), data = pd_num, type = "loglik-g")
saveRDS(res_pd_bn, here("examples/visit_order/res_pd_bn.rds"))

set.seed(10)
res_pd_vines <- process_dataset_syn(pd, n_rand = 0, n_synth = n_synth, max_struct_cons = Inf, n_cores = n_cores, struct_method = "vines")
res_pd_vines[[1]]$score_vines <- vine_from_order(pd, res_pd_vines[[1]]$order, n_cores)$logLik
saveRDS(res_pd_vines, here("examples/visit_order/res_pd_vines.rds"))

set.seed(11)
res_wbcd_bn <- process_dataset_syn(wbcd, n_rand = 0, n_synth = n_synth, max_struct_cons = Inf, n_cores = n_cores, struct_method = "bn")
wbcd_num <- set_num(wbcd)
res_wbcd_bn[[1]]$score_bn <- bnlearn::score(bn_from_order(wbcd_num, res_wbcd_bn[[1]]$order), data = wbcd_num, type = "loglik-g")
saveRDS(res_wbcd_bn, here("examples/visit_order/res_wbcd_bn.rds"))

set.seed(12)
res_wbcd_vines <- process_dataset_syn(wbcd, n_rand = 0, n_synth = n_synth, max_struct_cons = Inf, n_cores = n_cores, struct_method = "vines")
res_wbcd_vines[[1]]$score_vines <- vine_from_order(wbcd, res_wbcd_vines[[1]]$order, n_cores)$logLik
saveRDS(res_wbcd_vines, here("examples/visit_order/res_wbcd_vines.rds"))

set.seed(13)
res_hd_bn <- process_dataset_syn(hd, n_rand = 0, n_synth = n_synth, max_struct_cons = Inf, n_cores = n_cores, struct_method = "bn")
hd_num <- set_num(hd)
res_hd_bn[[1]]$score_bn <- bnlearn::score(bn_from_order(hd_num, res_hd_bn[[1]]$order), data = hd_num, type = "loglik-g")
saveRDS(res_hd_bn, here("examples/visit_order/res_hd_bn.rds"))

set.seed(14)
res_hd_vines <- process_dataset_syn(hd, n_rand = 0, n_synth = n_synth, max_struct_cons = Inf, n_cores = n_cores, struct_method = "vines")
res_hd_vines[[1]]$score_vines <- vine_from_order(hd, res_hd_vines[[1]]$order, n_cores)$logLik
saveRDS(res_hd_vines, here("examples/visit_order/res_hd_vines.rds"))

set.seed(15)
process_errprop(wbcd, 3:5, 3, 20, "res_errprop.rds")

set.seed(16)
process_targeted(wbcd, 3, 3:5, "res_targeted.rds")

set.seed(17)
process_errprop_phs(ad, 3:5, 3, 20, "res_errprop_phs.rds")

set.seed(18)
process_blocks(100, 30, (1:15)*2, "res_blocks.rds")

set.seed(19)
process_FCS(wbcd, "res_fcs.rds")

set.seed(20)
res_hp <- process_hyperparams(wbcd, (1:10)*10, synth_mb)
saveRDS(res_hp, here("examples/visit_order/res_hp.rds"))
