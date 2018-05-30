context('score_sapphire')

test_that('scoring sapphire plots', {
  silent <- F
  plt_stff <- !silent
  if(!silent) {require(testthat); require(CampaRi)} 
  
  stdd <- 3; n_dim <- 10; n_snap <- 3000; n_tot <- n_dim*n_snap/3; if(!silent) print(n_tot)
  ann <- c(rep(1, n_snap/6), rep(2, n_snap/3), rep(1, n_snap/6), rep(3, n_snap/3))
  file.pi <- system.file("extdata", "REPIX_000000000021.dat", package = "CampaRi")
  if(plt_stff) sapphire_plot(sap_file = file.pi, ann_trace = ann, only_timeline=F, timeline = T)

  optimal_bas <- CampaRi::basin_optimization(the_sap = file.pi,                                    # PI data
                                             how_fine_search = 5,                                  # number of bins searched between nbins_x_min and nbins_x_max
                                             nbins_x_min = 20,                                     # minimum number of bins (range)
                                             nbins_x_max = 30,                                     # maximum number of bins (range)
                                             basin_optimization_method = "MI_barrier_weighting",   # ranking method
                                             cl.stat.MI_comb = 'kin_MI',                           # Final score combination 
                                             cl.stat.nUni = c(5,10,15,20,25,30,40,50,60),          # MI curve splits 
                                             force_matching = T,                                   # match kin and tempann for the barriers
                                             number_of_clusters = 3,                               # number of cluster approach. NULL in this case it must be
                                             denat_opt = 'poly_interpolation',                     # kin ann is corrected for parabolic artefacts
                                             cl.stat.denat.MI = 7,                                 # if a number also the MI curve is corrected
                                             plot_basin_identification = plt_stff,                 # final plot?
                                             # dbg_basin_optimization = F,                           # debug?
                                             silent = silent)                                      # silent?  
  
  
  expect_error(qres <- score_sapphire(the_sap = file.pi, ann = ann, 
                                      basin_obj = optimal_bas$bas, plot_pred_true_resume = T,
                                      silent = silent), NA)

  
  # new tests using the nSBR function
  # ---------------------------------------------- nSBR manual_barriers option
  expect_error(a1 <- CampaRi::nSBR(data = file.pi, n.cluster = 3, 
                                  comb_met = c('MIC'),
                                  unif.splits = seq(5, 80, 4),  
                                  pk_span = 500, ny = 40, plot = T, 
                                  silent = silent, dbg_nSBR = F, return_plot = T), NA)
  expect_error(ahahscore <- CampaRi::score_sapphire(the_sap = file.pi, ann = ann, plot_pred_true_resume = plt_stff, multi_cluster_policy = 'popup',
                                                    manual_barriers = a1$barriers[1:2], silent = silent), NA)
  
  
  
  
  # test if the number of barriers != ncl (firstly < and secondly >)
  expect_error(ahahscore <- CampaRi::score_sapphire(the_sap = file.pi, ann = ann, plot_pred_true_resume = plt_stff,
                                                    multi_cluster_policy = 'popup', # not relevant for this case
                                                    manual_barriers = a1$barriers[1], silent = silent), NA)
  
  # more cluster than annotation
  expect_error(a1 <- CampaRi::nSBR(data = file.pi, n.cluster = 10, 
                                   comb_met = c('MIC'),
                                   unif.splits = seq(5, 80, 2),  
                                   pk_span = 100, ny = 40, plot = T, 
                                   silent = silent, dbg_nSBR = F, return_plot = T), NA)
  
  expect_error(ahahscore <- CampaRi::score_sapphire(the_sap = file.pi, ann = ann,  plot_pred_true_resume = plt_stff, dbg_score_sapphire = F, 
                                                    multi_cluster_policy = 'popup',
                                                    manual_barriers = a1$barriers[1:9], silent = silent), NA)
  expect_error(ahahscore <- CampaRi::score_sapphire(the_sap = file.pi, ann = ann,  plot_pred_true_resume = plt_stff, dbg_score_sapphire = F, 
                                                    multi_cluster_policy = 'keep',
                                                    manual_barriers = a1$barriers[1:9], silent = silent), NA)
  expect_error(ahahscore <- CampaRi::score_sapphire(the_sap = file.pi, ann = ann,  plot_pred_true_resume = plt_stff, dbg_score_sapphire = F, 
                                                    multi_cluster_policy = 'merge_previous',
                                                    manual_barriers = a1$barriers[1:9], silent = silent), NA)
  
  
  # =====================================================================================================================================================
  # Testing different metrics for multiclass interclustering similarity.
  # --------------------------------------------------------------------
  #
  # Tests used (first the chance properties):
  # 1. Cross-similarity. Create n_cluing clusterings randomly using a number of growing clusters 2____Kmax and compute similarity between pairs
  # 2. Similarity to Ktrue. Generate a random clustering with Ktrue = Kmax/2 clusters and run random clusterings from 2____Kmax (this is prob. redundant)
  # 
  # METHODS USED:
  # A. ARI
  # B. NMI
  # C. 
  #
  
  if(F){
    library(ClusterR); library(ggplot2)
    # init
    n_cluing <- 100
    n_points <- 150
    ks <- 2:50
    cat('Tot numbers to generate are:', max(ks)*n_cluing*n_points, '\n')
    tcl <- list()
    
    # cluster creation - tcl contains n_cluing runs for each ks
    tcl <- lapply(ks, function(ki){
      sapply(1:n_cluing, FUN = function(x) sample.int(ki, n_points, replace = T))
    })
    
    # all(sapply(tcl, max) == ks) # check for effective correct construction - it can be commented
    v_combinations <- combn(1:n_cluing, 2, simplify = T) # no diag and only one tri
    
    all_the_methods <- c('overlap.coef', 'rand_index', 'adjusted_rand_index', 'jaccard_index', 'fowlkes_mallows_index', 'mirkin_metric', 'purity', 'entropy', 'nmi') 
    # mirking and nmi not nec
    
    #, 'morisitas.index', 'horn.index') # bugs ;  'cosine.similarity' is not normalized
    # tcR package
    # tcR::tversky.index(x, y, .a = 0.5, .b = 0.5)
    # tcR::cosine.similarity(.alpha, .beta, .do.norm = NA, .laplace = 0)
    # tcR::overlap.coef(.alpha, .beta)
    # tcR::morisitas.index(.alpha, .beta, .do.unique = T)
    # tcR::horn.index(.alpha, .beta, .do.unique = T)
    # else if(which(i == all_the_methods2) == 2) tcR::cosine.similarity(.alpha = a, .beta = b, .do.norm = NA, .laplace = 0) # horrid warnings + not normalized
    # else if(which(i == all_the_methods2) == 4) print('ohoh') #tcR::morisitas.index(.alpha = a, .beta = b, .do.unique = T) # bug
    # else if(which(i == all_the_methods2) == 5) tcR::horn.index(.alpha = a, .beta = b, .do.unique = T) # bug
    out2 <- list()
    for(i in all_the_methods[1]){
      out2[[i]] <- lapply(1:length(tcl), FUN = function(x){
        cat('Doing', i,'part', x, 'out of', length(tcl), '\r')
        if(x == length(tcl)) cat('\n')
        sapply(1:ncol(v_combinations), function(y) {
          a <- tcl[[x]][,v_combinations[1,y]]
          b <- tcl[[x]][,v_combinations[2,y]]
          if(which(i == all_the_methods2) == 1) tcR::tversky.index(x = a, y = b, .a = 1, .b = 0)
          else if(which(i == all_the_methods2) == 2) tcR::overlap.coef(.alpha = a, .beta = b)
          else ClusterR::external_validation(true_labels = a, clusters = b, method = i) 
        })
      })
    }
    # save(out, file = 'cluster_varlidation.Rdata')
    # load(file = 'cluster_varlidation.Rdata')
    sc.m <- lapply(out, function(x) sapply(x, mean))
    sc.sd <- lapply(out, function(x) sapply(x, sd))
    dfpl <- data.frame('msc' = NULL, 'sdsc' = NULL, 'ks' = NULL, 'method' = NULL)
    for(i in 1:length(all_the_methods)){
      if(all_the_methods[i] == 'mirkin_metric') next
      dfpl <- rbind(dfpl, data.frame('msc' = sc.m[[i]], 'sdsc' = sc.sd[[i]], 'ks' = ks, 'method' = all_the_methods[[i]]))
    }
    ggplot(data = dfpl) + theme_classic() + 
      geom_line(mapping = aes(x=ks, y=msc, group = method, col = method)) + 
      geom_errorbar(mapping = aes(x = ks, ymin = msc - sdsc, ymax = msc + sdsc, group = method, col = method)) +
      ggtitle('Analysis of different clustering comparison measures', subtitle = 'Random generated vectors of 150 points. Std.dev made with 100 runs.') +
      xlab('# clusters') + ylab('Normalized score')
    
    ggsave('analysis_different_clu_comparisons.png', width = 8, height = 6)
    # other methods for ARI (and maybe something else)
    
    
    # Timings:
    a <- tcl[[40]][,1]; b <- tcl[[40]][,2]
    mo <- microbenchmark::microbenchmark(
      rand_index = external_validation(a, b, method = 'rand_index'), 
      rand_index21 = clues::adjustedRand(cl1 = a, cl2 = b, randMethod = 'Rand'), # this is a simple rand
      adjusted_rand_index1 = external_validation(a, b, method = 'adjusted_rand_index'),
      adjusted_rand_index22 = clues::adjustedRand(cl1 = a, cl2 = b, randMethod = 'HA'), 
      adjusted_rand_index23NO = clues::adjustedRand(cl1 = a, cl2 = b, randMethod = 'MA'), # diff number
      adjusted_rand_index24NO = clues::adjustedRand(cl1 = a, cl2 = b, randMethod = 'FM'), # diff number
      adjusted_rand_index3NO = fossil::adj.rand.index(group1 = a, group2 = b),            # diff number
      adjusted_rand_index4 = shallot::adj.rand.index(c1 = a, c2 = b),
      adjusted_rand_index5 = mclust::adjustedRandIndex(x = a, y = b),
      adjusted_rand_index6 = phyclust::RRand(trcl = a, prcl = b),                      
      jaccard_index = external_validation(a, b, method = 'jaccard_index'), 
      fowlkes_Mallows_index = external_validation(a, b, method = 'fowlkes_mallows_index'), 
      mirkin_metric = external_validation(a, b, method = 'mirkin_metric'), 
      purity = external_validation(a, b, method = 'purity'), 
      entropy = external_validation(a, b, method = 'entropy'), 
      # nmi = external_validation(a, b, method = 'nmi'), # superslow and bad
      times = 50
    )
    autoplot(mo) + theme_minimal() +
      ggtitle('Timings for all the method I found', subtitle = 'There are some flaws in the approx for approx of ari')
    ggsave('timings_of_methods.png')
    
    # lets check if they are simular
    fin <- list(adjusted_rand_index1 = external_validation(a, b, method = 'adjusted_rand_index'),
              adjusted_rand_index21 = clues::adjustedRand(cl1 = a, cl2 = b, randMethod = 'Rand'), # 
              adjusted_rand_index22 = clues::adjustedRand(cl1 = a, cl2 = b, randMethod = 'HA'), # 
              adjusted_rand_index23 = clues::adjustedRand(cl1 = a, cl2 = b, randMethod = 'MA'), # 
              adjusted_rand_index24 = clues::adjustedRand(cl1 = a, cl2 = b, randMethod = 'FM'), # 
              adjusted_rand_index3 = fossil::adj.rand.index(group1 = a, group2 = b),
              adjusted_rand_index4 = shallot::adj.rand.index(c1 = a, c2 = b),
              adjusted_rand_index5 = mclust::adjustedRandIndex(x = a, y = b),
              adjusted_rand_index6 = phyclust::RRand(trcl = a, prcl = b))
    fin
    
    
    # Only best ARI
    a <- tcl[[40]][,1]; b <- tcl[[40]][,2]
    mo <- microbenchmark::microbenchmark(
      ClusterR = ClusterR::external_validation(a, b, method = 'adjusted_rand_index'),
      clues = clues::adjustedRand(cl1 = a, cl2 = b, randMethod = 'HA'), 
      shallot = shallot::adj.rand.index(c1 = a, c2 = b),
      mclust = mclust::adjustedRandIndex(x = a, y = b),
      phyclust = phyclust::RRand(trcl = a, prcl = b),                      
      times = 50
    )
    autoplot(mo) + theme_minimal() + ggtitle('Timings for equally scoring adjusted random indexes')
    ggsave('Final_timings_ARI.png')
    # some notes on other measures
    # ~~~~~~~~
    # clues 
    # Available indices are: 
#     "Rand", 
#     "HA" (Hubert and Arabie's adjusted Rand index), 
#     "MA" (Morey and Agresti's adjusted Rand index), 
#     "FM" (Fowlkes and Mallows's index), 
#     "Jaccard" (Jaccard index). By default, all 5 indices will be output.
    
    
  }
})
