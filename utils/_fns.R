
# theme -------------------------------------------------------------------

theme_vt <- function(){
  theme_classic(base_size=14) %+replace%
    theme(axis.ticks.length = unit(0.3, "cm"),
          axis.ticks = element_line(color = 'black'),
          axis.text = element_text(color="black"),
          axis.title = element_text(color="black"),
          legend.text = element_text(color="black"),
          legend.title = element_text(color="black"))
}

# HDP - without priors -------------------------------------------------------
run_hdp_unsupervised <- function(df_spectra){
  
  require(hdp)
  # df_spectra = sbs96_spectra %>% dplyr::select(1:8)
  mat_spectra <- df_spectra %>% 
    tibble::column_to_rownames(var = "MutationType") %>% 
    t() %>% as.data.frame()
  
  hdp_mut <- hdp_init(ppindex = c(0, rep(1, 3), rep(2:4, each=100)), # index of parental nodes
                      cpindex = c(1, rep(2, 3), rep(3:5, each=100)), # index of concentration param
                      hh=rep(1, ncol(mat_spectra)), # prior is uniform over the 96 mutation categories
                      alphaa=rep(1, 5), # shape hyperparams for five different CPs
                      alphab=rep(1, 5)) # rate hyperparams for five different CPs
  # add data to leaf nodes (one per cancer sample, in row order of mut_count)
  req_num <- numdp(hdp_mut)-dim(mat_spectra)[1]+1
  
  
  hdp_mut <- hdp_setdata(hdp_mut, 
                         dpindex = req_num:numdp(hdp_mut), # index of nodes to add data to
                         mat_spectra) # input data (mutation counts, sample rows match up with specified dpindex)
  #hdp_mut
  chlist <- vector("list", 10)
  for (i in 1:10){
    
    # activate DPs, 10 initial components
    hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10, seed=i*200)
    
    chlist[[i]] <- hdp_posterior(hdp_activated,
                                 burnin=5000,
                                 n=200,
                                 space=200,
                                 cpiter=3,
                                 seed=i*1e3)
  }
  
  # example multi object
  mut_example_multi <- hdp_multi_chain(chlist)
  mut_example_multi <- hdp_extract_components(mut_example_multi)
  
  
  comp_distn <- comp_categ_distn(mut_example_multi)
  
  ncomp <- nrow(comp_distn$mean)-1
  comp_to_plot <- rownames(comp_distn$mean)
  
  processes_df <- data.frame()
  for (ii in seq_along(comp_to_plot)){
    
    cname <- comp_to_plot[ii]
    
    # mean categorical distribution (sig), and credibility interval
    sig <- comp_distn$mean[cname,]
    processes_df <- rbind(processes_df,sig)
  }
  #colnames(processes_df) <- 1:96
  
  ## keep process-0, however, this potentially is not necessary
  ## check while performing decomposition with nnls
  
  ## signature-names
  if(dim(mat_spectra)[2] == 96){
    
    sig_names <- paste0("SBS96", LETTERS[1:nrow(processes_df)])
  }else{
    
    sig_names <- paste0("ID83", LETTERS[1:nrow(processes_df)])
  }
  
  
  colnames(processes_df) <- colnames(mat_spectra)
  rownames(processes_df) <- sig_names
  
  # signatures matrix
  processes_mat <- processes_df %>% 
    t() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "MutationType")
  
  ## exposures matrix
  dp_distn <- comp_dp_distn(mut_example_multi)
  ndp <- nrow(dp_distn$mean)
  ncomp <- ncol(dp_distn$mean)
  
  dpindices <- req_num:numdp(hdp_mut)
  
  dps <- dp(final_hdpState(chains(mut_example_multi)[[1]]))[dpindices]
  pps <- ppindex(final_hdpState(chains(mut_example_multi)[[1]]))[dpindices]
  
  numdata <- sapply(dps, function(x) x@numdata)
  dp_order <- order(numdata, decreasing = TRUE)
  
  exposures <- t(dp_distn$mean[dpindices, , drop = FALSE])
  
  rownames(exposures) <- colnames(processes_mat)[-1]
  colnames(exposures) <- rownames(mat_spectra)
  
  exposures_mat <- exposures %>% 
    t() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "PID")
  
  ## final res
  f_res <- NULL
  f_res[['signatures']] <- processes_mat
  f_res[['exposures']] <- exposures_mat
  
  return(f_res)
  
}

# HDP - with priors -------------------------------------------------------
run_hdp_supervised <- function(df_spectra, sp_ref, priors){
  
  require(hdp)
  # df_spectra = inf_sbs96
  # sp_ref = sp_sbs96_ref
  # priors = c("SBS1", "SBS3", "SBS5", "SBS8", "SBS18", "SBS25", "SBS31", "SBS32", "SBS35")
  mat_spectra <- df_spectra %>% 
    tibble::column_to_rownames(var = "MutationType") %>% 
    t() 
  
  ## priors
  if(is.null(priors)){
    
    prior_sigs <- sp_ref
  }else{
    
    prior_sigs <- sp_ref %>% 
      dplyr::select(one_of(c("Type", priors)))
  }
  prior_sigs <- prior_sigs %>% 
    tibble::column_to_rownames(var = "Type") %>% 
    as.matrix()
  
  ## hdp
  nps <- ncol(prior_sigs)
  nsamples <- nrow(mat_spectra)
  
  ttype_prior <- hdp_prior_init(prior_distn = prior_sigs,
                                prior_pseudoc = rep(1000, nps),
                                hh=rep(1, ncol(mat_spectra)),
                                alphaa=c(1, 1),
                                alphab=c(1, 1))
  
  ttype_prior <- hdp_addconparam(ttype_prior,
                                 alphaa = c(1,1),
                                 alphab = c(1,1))
  
  ttype_prior <- hdp_adddp(ttype_prior,
                           numdp = nsamples+1,
                           ppindex = c(1, rep(1+nps+1, nsamples)),
                           cpindex = c(3, rep(4, nsamples)))
  
  ttype_prior <- hdp_setdata(ttype_prior,
                             dpindex = (1+nps+1)+1:nsamples,
                             mat_spectra[1:nsamples,])
  
  
  chlist <- vector("list", 10)
  
  for(i in 1:10){
    
    ttype_activated <- dp_activate(ttype_prior,
                                   dpindex = (1+nps+1)+0:nsamples,
                                   initcc = nps+10,
                                   seed = i*1000)
    
    chlist[[i]] <- hdp_posterior(ttype_activated,
                                 burnin = 5000,
                                 n = 100,
                                 space = 200,
                                 cpiter = 3,
                                 seed = i*1e6,
                                 verbosity = 1)
  }
  ttype_multi <- hdp_multi_chain(chlist)
  ttype_multi <- hdp_extract_components(ttype_multi)
  
  ## extract signatures & exposures
  
  # activities ##
  tnt <- mat_spectra
  toutmeans <- (comp_dp_distn(ttype_multi)$mean)
  mymeans <- (toutmeans[(nrow(toutmeans)-nrow(tnt) + 1):nrow(toutmeans),])
  rownames(mymeans) <- rownames(tnt)
  
  tochange <- colnames(mymeans)[grep("P", colnames(mymeans))]
  colnames(mymeans)[grep("P", colnames(mymeans))] <- colnames(prior_sigs)[as.numeric(gsub("P", "", tochange))]
  
  exposures_df <- mymeans * rowSums(mat_spectra)
  exposures_df <- exposures_df %>% as.data.frame() %>% 
    tibble::rownames_to_column(var = "PID")
  
  ## signature profiles ##
  processes_df <- t(comp_categ_distn(ttype_multi)$mean) 
  tochange <- colnames(processes_df)[grep("P", colnames(processes_df))]
  
  colnames(processes_df)[grep("P", colnames(processes_df))] <- colnames(prior_sigs)[as.numeric(gsub("P", "", tochange))]
  
  ## add rownames
  rownames(processes_df) <- rownames(prior_sigs)
  processes_df <- processes_df %>% as.data.frame() %>% 
    tibble::rownames_to_column(var = "MutationType")
  
  f_res <- list()
  f_res[['signatures']] <- processes_df
  f_res[['exposures']] <- exposures_df
  
  return(f_res)
}


# PCAWG style overview plots ----------------------------------------------

filter_mat <- function(df, min_mut){
  
  # df = sp_sbs96
  df$total <- df %>% 
    dplyr::select(-c(PID, CANCER_TYPE)) %>% 
    rowSums()
  
  df %>% 
    dplyr::filter(total > min_mut) %>% 
    dplyr::select(-c(total))
}

process_exposure_table <- function(df, context = 'SBS'){
  
  ## functions process exposures table & prepares required format for plotting
  ## with all filterings & stats applied
  # df = dat
  
  ## need data in the format
  ## PID CANCER_TYPE, SBS1/ID1 SBS2/ID2 .... SBSn/IDn
  
  all_ctypes <- df$CANCER_TYPE %>% as.character() %>% unique()
  
  perc_comps <- function(x){
    
    # x = 'ACC'
    x_dat <- df %>% dplyr::filter(CANCER_TYPE == get("x")) %>% 
      dplyr::select(-c(CANCER_TYPE))
    x_total_pid <- dim(x_dat)[1]
    
    
    suppressMessages(x_mlt <- x_dat %>% 
                       tibble::column_to_rownames(var = "PID") %>% 
                       dplyr::select_if(colSums(.) != 0) %>% 
                       tibble::rownames_to_column(var = "PID") %>% 
                       reshape2::melt() %>% 
                       dplyr::filter(value != 0))
    
    
    ## fractions
    x_fraction_active <- x_mlt %>% 
      dplyr::count(variable) %>% 
      dplyr::mutate(FRACTION = n/x_total_pid) %>% 
      dplyr::select(SIGNATURE = variable, FRACTION)
    
    ## medians
    suppressMessages(x_frac_medians <- x_mlt %>% 
                       dplyr::select(-c(PID)) %>% 
                       dplyr::group_by(variable) %>% 
                       dplyr::summarise(MEDIAN = median(value)) %>% 
                       dplyr::mutate(MEDIAN = MEDIAN/2800) %>% 
                       dplyr::rename(SIGNATURE = variable) %>% 
                       merge(x_fraction_active, by = "SIGNATURE") %>% 
                       dplyr::mutate(CANCER_TYPE = x) %>% 
                       dplyr::select(CANCER_TYPE, SIGNATURE, MEDIAN, FRACTION))
    
    return(x_frac_medians)
  }
  
  all_ctype_res <- lapply(all_ctypes, perc_comps) %>% 
    plyr::ldply()
  
  ## add undetetcted signatures, MEDIAN=0, FRACTION=0
  sbs_ref <- c("SBS1","SBS2","SBS3","SBS4","SBS5","SBS6","SBS7a","SBS7b","SBS7c",
               "SBS7d","SBS8","SBS9","SBS10a","SBS10b","SBS11","SBS12","SBS13",
               "SBS14","SBS15","SBS16","SBS17a","SBS17b","SBS18","SBS19","SBS20",
               "SBS21","SBS22","SBS23","SBS24","SBS25","SBS26","SBS27","SBS28",
               "SBS29","SBS30","SBS31","SBS32","SBS33","SBS34","SBS35","SBS36",
               "SBS37","SBS38","SBS39","SBS40","SBS41","SBS42","SBS44")
  
  id_ref <- c("ID1", "ID2", "ID3", "ID4", "ID5", "ID6", "ID7", "ID8", "ID9", "ID10", "ID11", "ID12", "ID13", "ID14", "ID15", "ID16", "ID17")
  
  if(context == 'SBS'){
    
    undetected_sigs <- setdiff(sbs_ref, all_ctype_res$SIGNATURE %>% unique()) %>% 
      as.data.frame() %>% dplyr::rename(SIGNATURE = ".")
  }else{
    
    undetected_sigs <- setdiff(id_ref, all_ctype_res$SIGNATURE %>% unique()) %>% 
      as.data.frame() %>% dplyr::rename(SIGNATURE = ".")
    
  }
  
  all_sigs_median <- all_ctype_res %>% 
    dplyr::select(SIGNATURE, CANCER_TYPE, MEDIAN) %>% 
    reshape2::dcast(SIGNATURE~CANCER_TYPE, value.var = "MEDIAN") %>% 
    plyr::rbind.fill(undetected_sigs) %>% 
    reshape2::melt() %>% 
    dplyr::rename(CANCER_TYPE = variable, MEDIAN = value)
  
  all_sigs_meds_fracs <- all_ctype_res %>% 
    dplyr::select(SIGNATURE, CANCER_TYPE, FRACTION) %>% 
    reshape2::dcast(SIGNATURE~CANCER_TYPE, value.var = "FRACTION") %>% 
    plyr::rbind.fill(undetected_sigs) %>% 
    reshape2::melt() %>% 
    dplyr::rename(CANCER_TYPE = variable, FRACTION = value) %>% 
    merge(all_sigs_median, by = c("CANCER_TYPE", "SIGNATURE"))
  
  
  ## add total number of tumors per indication
  all_sigs_meds_fracs <- df %>% 
    dplyr::count(CANCER_TYPE) %>% 
    dplyr::mutate(SAMPLE_N = paste0(CANCER_TYPE, " (n=", n, ")")) %>% 
    dplyr::select(-c(n)) %>% 
    merge(all_sigs_meds_fracs, by = "CANCER_TYPE", all.y = TRUE)
  
  return(all_sigs_meds_fracs)
}

fmt_sigtable_forplotting <- function(df){
  
  df <- df %>% 
    dplyr::rename(ENTITY = CANCER_TYPE, 
                  MedianSignature = MEDIAN,
                  SAMPLES_ACTIVE = FRACTION,
                  ENTITY_N = SAMPLE_N)
  
  df_x <- df %>% 
    dplyr::select(ENTITY, SIGNATURE, MedianSignature) %>% 
    reshape2::dcast(ENTITY~SIGNATURE, value.var = "MedianSignature")
  
  tile_nums <- rep(c(1,0), 100000)
  
  df_x$tile_col <- tile_nums[1:dim(df_x)[1]]
  df_x$tile_col <- factor(df_x$tile_col, levels = c(0, 1))
  
  df_x <- df_x %>% 
    reshape2::melt() %>% 
    dplyr::mutate(PID = paste0(ENTITY, " ", variable)) 
  
  df_final <- df %>% 
    dplyr::mutate(PID = paste0(ENTITY, " ", SIGNATURE)) %>% 
    dplyr::select(-c(MedianSignature, ENTITY, SIGNATURE)) %>% 
    merge(df_x, by = "PID") %>% 
    dplyr::select(-c(PID)) %>% 
    dplyr::rename(MedianSignature = value,
                  SIGNATURE = variable) 
  
  df_final$ENTITY_N <- factor(df_final$ENTITY_N, levels = unique(df_final$ENTITY_N))
  
  ## order signatures
  ## starts with SBS & ends with [1-9]
  all_sigs <- df_final$SIGNATURE %>% as.character() %>% unique()
  
  all_sigs_order <- c(all_sigs[grepl("^SBS", all_sigs)] %>% gtools::mixedsort(),
                      all_sigs[!grepl("^SBS", all_sigs)] %>% gtools::mixedsort())
  
  df_final$SIGNATURE <- factor(df_final$SIGNATURE, 
                               levels = rev(all_sigs_order))
  
  return(df_final)
  
}

plot_signature_overview <- function(df, cl_limits, cl_breaks){
  
  ## colors
  colPal <- colorRampPalette(c('#FFC651', "#FF6B6B", "#6767F7", "#003D7C"))
  plt_colors <- colPal(length(cl_breaks))
  
  ## final plot
  ggplot(df, aes(x = ENTITY_N, y = SIGNATURE, fill = tile_col, size = SAMPLES_ACTIVE, color = MedianSignature))+
    geom_tile(size = 0.5, color = "white")+
    scale_fill_manual(values = c("1" = "gray95", "0" = "gray90"))+
    geom_point()+
    theme_minimal(base_size = 14)+
    theme(axis.text = element_text(color = "black", size = 14),
          axis.text.x = element_text(color = "black", size = 14),
          axis.title = element_text(size = 14),
          panel.grid = element_blank(),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14))+
    xlab("")+ggpubr::rotate_x_text(angle = 90)+
    ylab("SBS Signatures")+
    guides(fill = 'none')+
    binned_scale(aesthetics = "color",
                 scale_name = "test", 
                 palette = function(x) plt_colors,
                 breaks = cl_breaks,
                 limits = cl_limits,
                 n.breaks = 5,
                 oob = scales::squish,
                 show.limits = FALSE, 
                 guide = guide_coloursteps(barwidth=1.3, barheight = 15,
                                           label.position="right",
                                           title.position = "right",
                                           frame.colour = 'black',
                                           frame.linewidth = 1.5,
                                           title.theme = element_text(angle = 90, size = 18)))+
    scale_size_area(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
}


# SBS96 profile plot ------------------------------------------------------

# SBS-96 profile plot -------------------------------------------------------------
format_sbs_profile <- function(x){
  
  # mutationType, SBS8
  # x <- mb_g3_sbs8
  
  # plot order
  x_ordder <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T",
                "C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A",
                "G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C",
                "T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G",
                "A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T",
                "G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A",
                "T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C",
                "A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G",
                "C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T",
                "T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A",
                "A[T>A]C","A[T>A]G","A[T>A]T",'C[T>A]A',"C[T>A]C",
                "C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G",
                "G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T",
                "A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A",
                "C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C",
                "G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G",
                "T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T",
                "C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A",
                "G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C",
                "T[T>G]G","T[T>G]T")
  
  # add 6 subs. types as a new column
  x <- x %>%
    dplyr::mutate(SUB = MutationType) %>%
    dplyr::mutate(SUB = gsub("].", "", SUB),
                  SUB = gsub(".\\[", "", SUB))
  
  ## order according to x_order
  x <- x[match(x_ordder, x$MutationType), ]
  
  x$MutationsType <- factor(x$MutationType,
                            levels = x_ordder)
  
  return(x)
  
}

plot_sbs_profile_old <- function(x, sig_name){
  
  # sig_name <- "SBS8"
  
  # colors
  subs_colpal <- c("#03BCEE", "#010101", "#E32926",
                   "#CAC9C9", "#A1CE63", "#EBC6C4")
  
  names(subs_colpal) <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  
  
  plt <- ggplot(x, aes_string(x = "MutationsType", y = sig_name, fill = "SUB"))+
    geom_bar(stat = "identity")+
    facet_grid(.~SUB, scales = "free_x")+
    scale_fill_manual(values = subs_colpal)+
    theme_classic(base_size = 12)+
    theme(axis.text = element_text(color = 'black'),
          axis.text.x = element_text(angle=90, size = 6),
          panel.grid = element_blank(),
          legend.position = "none",
          axis.line.x = element_blank(),
          panel.spacing = unit(0.001, "lines"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.ticks.x = element_blank(),
          strip.background.x = element_rect(fill = c("gray"), color = NA))+
    scale_y_continuous(expand = expansion(mult = c(0, .1)))+
    xlab("")+ylab("Relative contribution")
  
  return(plt)
  
}

plot_sbs_profile_notop <- function(x, sig_name){
  
  # sig_name <- "SBS8"
  
  # colors
  subs_colpal <- c("#03BCEE", "#010101", "#E32926",
                   "#CAC9C9", "#A1CE63", "#EBC6C4")
  
  names(subs_colpal) <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  
  plt <- ggplot(x, aes_string(x = "MutationsType", y = sig_name, fill = "SUB"))+
    geom_bar(stat = "identity")+
    facet_grid(.~SUB, scales = "free_x")+
    scale_fill_manual(values = subs_colpal)+
    theme_classic(base_size = 12)+
    theme(axis.text = element_text(color = 'black'),
          axis.text.x = element_text(angle=90, size = 6),
          panel.grid = element_blank(),
          legend.position = "none",
          axis.line.x = element_blank(),
          panel.spacing = unit(0.001, "lines"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.ticks.x = element_blank(),
          strip.background.x = element_blank(),
          strip.text = element_blank())+
    scale_y_continuous(expand = expansion(mult = c(0, .1)))+
    xlab("")+ylab("Relative contribution")
  
  return(plt)
  
}

plot_sbs_profile_topcolor <- function(x, sig_name){
  
  # sig_name <- "SBS8"
  
  # colors
  subs_colpal <- c("#03BCEE", "#010101", "#E32926",
                   "#CAC9C9", "#A1CE63", "#EBC6C4")
  
  names(subs_colpal) <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  
  
  p <- ggplot(x, aes_string(x = "MutationsType", y = sig_name, fill = "SUB"))+
    geom_bar(stat = "identity")+
    facet_grid(.~SUB, scales = "free_x")+
    scale_fill_manual(values = subs_colpal)+
    theme_classic(base_size = 12)+
    theme(axis.text = element_text(color = 'black'),
          axis.text.x = element_text(angle=90, size = 6),
          panel.grid = element_blank(),
          legend.position = "none",
          axis.line.x = element_blank(),
          panel.spacing = unit(0.001, "lines"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.ticks.x = element_blank(),
          strip.background.x = element_rect(fill = c("gray"), color = NA),
          strip.text.x = element_text(color = 'white'))+
    scale_y_continuous(expand = expansion(mult = c(0, .1)))+
    xlab("")+ylab("Relative contribution")+
    ggtitle(sig_name)
  
  ### edit strip-background colors
  
  require(grid)
  g <- grid.force(ggplotGrob(p))
  
  grobs_df <- do.call(cbind.data.frame, grid.ls(g, print = FALSE))
  
  grobs_df$gPath_full <- paste(grobs_df$gPath, grobs_df$name, sep = "::")
  grobs_df$gPath_full <- gsub(pattern = "layout::", 
                              replacement = "", 
                              x = grobs_df$gPath_full, 
                              fixed = TRUE)
  
  ##
  strip_bg_gpath <- grobs_df$gPath_full[grepl(pattern = ".*strip\\.background.*", 
                                              x = grobs_df$gPath_full)]
  
  ## coloring
  ## subs_colpal
  fills <- subs_colpal
  
  for (i in 1:length(strip_bg_gpath)){
    g <- editGrob(grob = g, gPath = strip_bg_gpath[i], gp = gpar(fill = fills[i]))
  }
  
  grid.draw(g)
  
}

plot_sbs_profile_topcolor_minimal <- function(x, sig_name){
  
  # sig_name <- "SBS8"
  
  # colors
  subs_colpal <- c("#03BCEE", "#010101", "#E32926",
                   "#CAC9C9", "#A1CE63", "#EBC6C4")
  
  names(subs_colpal) <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  
  
  p <- ggplot(x, aes_string(x = "MutationsType", y = sig_name, fill = "SUB"))+
    geom_bar(stat = "identity")+
    facet_grid(.~SUB, scales = "free_x")+
    scale_fill_manual(values = subs_colpal)+
    theme_classic(base_size = 12)+
    theme(axis.text = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none",
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          panel.spacing = unit(0.001, "lines"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background.x = element_rect(fill = c("gray"), color = NA),
          strip.text.x = element_text(color = 'white'))+
    scale_y_continuous(expand = expansion(mult = c(0, .1)))+
    xlab("")+ylab("")+
    ggtitle(sig_name)
  
  ### edit strip-background colors
  
  require(grid)
  g <- grid.force(ggplotGrob(p))
  
  grobs_df <- do.call(cbind.data.frame, grid.ls(g, print = FALSE))
  
  grobs_df$gPath_full <- paste(grobs_df$gPath, grobs_df$name, sep = "::")
  grobs_df$gPath_full <- gsub(pattern = "layout::", 
                              replacement = "", 
                              x = grobs_df$gPath_full, 
                              fixed = TRUE)
  
  ##
  strip_bg_gpath <- grobs_df$gPath_full[grepl(pattern = ".*strip\\.background.*", 
                                              x = grobs_df$gPath_full)]
  
  ## coloring
  ## subs_colpal
  fills <- subs_colpal
  
  for (i in 1:length(strip_bg_gpath)){
    g <- editGrob(grob = g, gPath = strip_bg_gpath[i], gp = gpar(fill = fills[i]))
  }
  
  grid.draw(g)
  
}


# ID-83 profile plot ------------------------------------------------------

format_id_profile <- function(x){
  
  # x <- mb_g3_id1
  
  x_context_column <- x %>%
    dplyr::mutate(SUB = MutationsType) %>%
    dplyr::mutate(SUB = sub(":[^:]+$", '', SUB))
  
  # x-axis
  x_order_id <- c(rep(c(1:5, "6+"), 2), rep(c(0:4, "5+"), 2), rep(c(1:5, "6+"), 4),
                  rep(c(0:4, "5+"), 4), c("1"), c("1", "2"), c("1", "2", "3"), c(1:4, "5+"))
  
  # # context-name
  # x_context <- c(rep("C", 6), rep(T, "6"), rep("C", 6), rep("T", 6),
  # rep("2", 6), rep("3", 6), rep("4", 6), rep("5+", 6),
  # rep("2", 6), rep("3", 6), rep("4", 6), rep("5+", 6),
  # rep("2", 1), rep("3", 2), rep("4", 3), rep("5+", 5))
  
  x_context_column$x_Length <- x_order_id
  # x$x_Context <- x_context
  
  # factorize orders
  # x$x_Length <- factor(x$x_Length, levels = unique(x_order_id))
  # x$x_Context <- factor(x$x_Context, levels = unique(x_context))
  x_context_column_uniq <- x_context_column$SUB %>% as.character() %>% unique()
  
  x_context_column$SUB <- factor(x_context_column$SUB,
                                 levels = x_context_column_uniq)
  
  return(x_context_column)
  
}

plot_id_profile_old <- function(x, sig_name){
  
  # x <- mb_g3_id1_format
  # sig_name = "ID5"
  
  # colors
  id_subs_colpal <- c("#FCBD6F", "#FE8002", "#AFDC8A", "#36A02E",
                      "#FCC9B4", "#FB896A", "#F04432", "#BB191A",
                      "#CFE0F1", "#93C3DE", "#4A97C8", "#1764AA",
                      "#E1E1EE", "#B5B5D7", "#8582BC", "#62409A")
  
  names(id_subs_colpal) <- x$SUB %>% levels()
  
  
  plt <- ggplot(x, aes_string(x = "x_Length", y = sig_name, fill = "SUB"))+
    geom_bar(stat = "identity")+
    facet_grid(.~SUB, scales = "free_x")+
    scale_fill_manual(values = id_subs_colpal)+
    theme_classic(base_size = 13)+
    theme(axis.text = element_text(color = 'black'),
          axis.text.x = element_text(angle=90, size = 8),
          panel.grid = element_blank(),
          legend.position = "none",
          axis.line.x = element_blank(),
          panel.spacing = unit(0.001, "lines"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.ticks.x = element_blank(),
          strip.background.x = element_rect(fill = c("gray"), color = NA))+
    scale_y_continuous(expand = expansion(mult = c(0, .1)))+
    xlab("")+ylab("Relative contribution")
  
  return(plt)
  
}

plot_id_profile_topcolor <- function(x, sig_name){
  
  # x <- mb_g3_id1_format
  # sig_name = "ID5"
  
  # colors
  id_subs_colpal <- c("#FCBD6F", "#FE8002", "#AFDC8A", "#36A02E",
                      "#FCC9B4", "#FB896A", "#F04432", "#BB191A",
                      "#CFE0F1", "#93C3DE", "#4A97C8", "#1764AA",
                      "#E1E1EE", "#B5B5D7", "#8582BC", "#62409A")
  
  names(id_subs_colpal) <- x$SUB %>% levels()
  
  
  p <- ggplot(x, aes_string(x = "x_Length", y = sig_name, fill = "SUB"))+
    geom_bar(stat = "identity")+
    facet_grid(.~SUB, scales = "free_x")+
    scale_fill_manual(values = id_subs_colpal)+
    theme_classic(base_size = 13)+
    theme(axis.text = element_text(color = 'black'),
          axis.text.x = element_text(angle=90, size = 8),
          panel.grid = element_blank(),
          legend.position = "none",
          axis.line.x = element_blank(),
          panel.spacing = unit(0.001, "lines"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.ticks.x = element_blank(),
          strip.text = element_text(color = 'white'),
          strip.background.x = element_rect(fill = c("gray"), color = NA))+
    scale_y_continuous(expand = expansion(mult = c(0, .1)))+
    xlab("")+ylab("Relative contribution")+
    ggtitle(sig_name)
  
  
  require(grid)
  g <- grid.force(ggplotGrob(p))
  
  grobs_df <- do.call(cbind.data.frame, grid.ls(g, print = FALSE))
  
  grobs_df$gPath_full <- paste(grobs_df$gPath, grobs_df$name, sep = "::")
  grobs_df$gPath_full <- gsub(pattern = "layout::", 
                              replacement = "", 
                              x = grobs_df$gPath_full, 
                              fixed = TRUE)
  
  ##
  strip_bg_gpath <- grobs_df$gPath_full[grepl(pattern = ".*strip\\.background.*", 
                                              x = grobs_df$gPath_full)]
  
  ## coloring
  ## subs_colpal
  fills <- id_subs_colpal
  
  for (i in 1:length(strip_bg_gpath)){
    g <- editGrob(grob = g, gPath = strip_bg_gpath[i], gp = gpar(fill = fills[i]))
  }
  
  grid.draw(g)
  
}



