context('basins_recognition')

test_that('Test for basin recognition with ext files', {
  
  silent <- T
  plt_stff <- !silent
  if(!silent) {require(testthat); require(CampaRi)}
  ## dir.cur <- getwd()
  ## str.st <- rev(gregexpr('CampaRi', dir.cur)[[1]])[1]
  ## dir.root <- substr(dir.cur, 1, str.st + 7)
  file.pi <- system.file("extdata", "PROGIDX_000000000001_sbrtest.dat", package = "CampaRi")
  ## data.pi <- data.table::fread(paste0(dir.root, "inst/extdata/PROGIDX_000000000001_sbrtest.dat"), data.table=FALSE)[,c(1,3,4,6,7,5)]
  data.pi<- data.table::fread(file.pi, data.table=FALSE)[,c(1,3,4,6,7,5)]
  
  ## basins_recognition(data.pi, nx=50, ny=50, match=F, plot=T, out.file=F, silent=F)
  nxb <- 50
  nyb <- 50
  lcs <- c(F,T)
  mat <- c(T,F)
  avg <- c('movav', 'SG')
  
  parspace <- expand.grid(nxb, nyb, lcs, mat, avg)
  colnames(parspace) <- c("nxbb", "nybb", "lcsb", "matb", "avgb")
  
  set.seed(342)
  if(!silent) cat("Testing", nrow(parspace), "parameter combinations on basins_recognition\n")
  attach(parspace)
  for(ii in seq(nrow(parspace))) {
    expect_error(tmp <- basins_recognition(data=file.pi, nx=nxbb[ii], ny=nybb[ii], ny.aut=F, local.cut=lcsb[ii], 
                                           match=matb[ii], dyn.check=2, avg.opt=as.character(avgb[ii]), plot=F, out.file=(ii==1), silent=T, pol.degree=3), NA)
    if(ii==1) {
      expect_error(out <- data.table::fread("BASINS_000000000001_sbrtest.dat", data.table=FALSE), NA)
      expect_true(all(seq(nrow(tmp$tab.st)) == sort(unique(out[,3]))))
    } 
  }
  detach(parspace)
  
  ## Testing other less relevant parameter
  expect_error(tmp <- basins_recognition(data=file.pi, nx=10, ny=10, 
                                         ny.aut=T, local.cut=F, match=T, dyn.check=2, avg.opt="SG", plot=T, out.file=T, silent=F, pol.degree=3), NA)
  
  expect_error(tmp <- basins_recognition(data=file.pi, nx=2, ny=2, ny.aut=T, local.cut=F, match=T))
  expect_error(tmp <- basins_recognition(data=data.table::fread(file.pi)[1:100,], nx=7, ny=7, out.file = T, plot = T), NA)
  if(file.exists("./BASINS_21.dat")) file.remove("./BASINS_21.dat")
  # list.files()
  expect_error(out <- data.table::fread("./BASINS_000000000001_sbrtest.dat", data.table=FALSE), NA)
  if(file.exists("./BASINS_000000000001_sbrtest.dat")) file.remove("./BASINS_000000000001_sbrtest.dat")
  expect_true(all(seq(nrow(tmp$tab.st)) == sort(unique(out[,3]))))
  expect_warning(tmp <- basins_recognition(data=file.pi, nx=parspace$nxbb[1], ny=parspace$nybb[1], only.kin=T, 
                                           local.cut=F, match=T, dyn.check=1, avg.opt="SG", plot=F, out.file=F, silent=F, pol.degree=1.5))
  
  # ------------------------------------------------------------------------------------------------------------------------------------- dgarol  
  # testing for minimum number of bins without crashing
  pi.table <- data.table::fread(file.pi, verbose = F, data.table = F)
  expect_error(CampaRi::basins_recognition(pi.table[1:50,], nx = 7, new.dev = F, plot = T, out.file = F, silent = T))
  expect_error(CampaRi::basins_recognition(pi.table[1:80,], nx = 6, new.dev = F, plot = T, out.file = F, silent = T))
  
  # testing various analysis of the clusters
  file.pi <- system.file("extdata", "REPIX_000000000021.dat", package = "CampaRi")
  # file.pi <- system.file("extdata", "PROGIDX_000000000001_sbrtest.dat", package = "CampaRi")
  
  nbin <- round(sqrt(3000*10)); if(!silent) print(nbin) # hard wired
  expect_error(bas <- basins_recognition(data = file.pi, nx = nbin, match = F, plot = plt_stff, out.file = F, new.dev = F, 
                                         cl.stat = T,
                                         cl.stat.KL = T,
                                         cl.stat.TE = T,
                                         cl.stat.wMI = T,
                                         cl.stat.stft = T,
                                         cl.stat.entropy = T, cl.stat.denat = 'process_subtraction', dbg_basins_recognition = F,
                                         plot.cl.stat = T), NA)
  
  expect_error(bas <- basins_recognition(data = file.pi, nx = 7, match = F, plot = plt_stff, out.file = F, plot.cl.stat = plt_stff,
                                         new.dev = F, cl.stat.weight.barriers = T), NA)
  
  
  # ------------------------------------------------------- neuro tests
  # the following tests have been executed only 
  # extensive testing on real cases. It needs not to be
  # executed at all
  do_it <- FALSE
  # fpi_empty  
  expect_error(bas <- basins_recognition(data = file.pi, nx = 150, match = F, plot = T, out.file = F, new.dev = F, 
                                         cl.stat = T, cl.stat.denat = 'poly_interpolation', dbg_basins_recognition = F,
                                         plot.cl.stat = T, cl.stat.MI_comb = 'MI',
                                         cl.stat.nUni = c(5,10,15,20,25,30,40,50,60)), NA)
  if(do_it){
    # fpi_best  
    bas <- basins_recognition(data = fpi_best, nx = 150, match = T, plot = T, out.file = F, new.dev = F, 
                              cl.stat = T, cl.stat.denat = 'process_subtraction', cl.stat.nUni = 20,
                              plot.cl.stat = T)
    bout1 <- CampaRi::basins_recognition(fpi_best, nx = 67, plot = T, match = T, out.file = F, new.dev = F, silent = F,
                                         cl.stat.weight.barriers = T, cl.stat.denat = T, plot.cl.stat = T, # true is like process_subtraction
                                         cl.stat.nUni = c(5,10,15,20,25,30,40,50,60), dbg = F, cl.stat.MI_comb = 'kin_MI', dbg_basins_recognition = F) #hard wired
    bout2 <- CampaRi::basins_recognition(fpi_best, nx = 133, plot = T, match = T, out.file = F, new.dev = F, silent = F,
                                         cl.stat.weight.barriers = T, cl.stat.denat = T, plot.cl.stat = T, cl.stat.MI_comb = 'kin_MI', cl.stat.denat.MI = -1,
                                         cl.stat.nUni = seq(50,200,5)) #hard wired
    # fpi_ave  
    bas <- basins_recognition(data = fpi_ave, nx = 150, match = T, plot = T, out.file = F, new.dev = F, 
                              cl.stat = T, cl.stat.denat = 'process_subtraction', cl.stat.nUni = 20,
                              plot.cl.stat = T, dbg_basins_recognition = F)
    # fpi_worst  
    bas <- basins_recognition(data = fpi_worst, nx = 150, match = T, plot = T, out.file = F, new.dev = F, 
                              cl.stat = T, cl.stat.denat = 'process_subtraction', cl.stat.nUni = c(5,10,15,20,25,30,40,50,60),
                              plot.cl.stat = T)
    bas <- basins_recognition(data = fpi_worst, nx = 150, match = T, plot = T, out.file = F, new.dev = F, 
                              cl.stat = T, cl.stat.denat = 'poly_interpolation', cl.stat.nUni = c(5,10,15,20,25,30,40,50,60),
                              plot.cl.stat = T)
  }
  
  
  # x <- rnorm(10000, mean = 0, sd = 10)
  # require(microbenchmark); require(ggplot2)
  # a <- microbenchmark(bbmisc = BBmisc::normalize(x, method = 'standardize', range = c(0,1)),
  #                     classic =  (x-min(x))/(max(x)-min(x)), `` = 50)
  # autoplot(a)
  
  cat("Test on basins_recognition done\n\n")
  if(file.exists('Rplots.pdf')) system('rm Rplot*')
  if(file.exists('BASINS_000000000001_sbrtest.dat')) file.remove('BASINS_000000000001_sbrtest.dat')
  if(file.exists('BASINS_1.dat')) file.remove('BASINS_1.dat')
})
