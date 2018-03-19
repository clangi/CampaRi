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
  nxb <- seq(10,80,by=20)
  nyb <- seq(10,80,by=20)
  lcs <- c(F,T)
  mat <- c(T,F)
  avg <- c('movav', 'SG')
  
  parspace <- expand.grid(nxb, nyb, lcs, mat, avg)
  colnames(parspace) <- c("nxbb", "nybb", "lcsb", "matb", "avgb")
  
  set.seed(342)
  cat("Testing", nrow(parspace), "parameter combinations on basins_recognition")
  attach(parspace)
  for(ii in seq(nrow(parspace))) {
    expect_error(tmp <- basins_recognition(data=file.pi, nx=nxbb[ii], ny=nybb[ii], ny.aut=F, local.cut=lcsb[ii], out.file = F, 
                                           match=matb[ii], dyn.check=2, avg.opt=as.character(avgb[ii]), plot=F, silent=T, pol.degree=3), NA)
    if(ii==1) {
      expect_error(out <- data.table::fread("BASINS_000000000001_sbrtest.dat", data.table=FALSE), NA)
      expect_true(all(seq(nrow(tmp$tab.st)) == sort(unique(out[,3]))))
    } 
  }
  detach(parspace)
  
  ## Testing other less relevant parameter
  expect_error(tmp <- basins_recognition(data=file.pi, nx=parspace$nxbb[1], ny=parspace$nybb[1], 
                                         ny.aut=T, local.cut=F, match=T, dyn.check=2, avg.opt="SG", plot=T, out.file=F, silent=F, pol.degree=3), NA)
  # list.files()
  expect_error(out <- data.table::fread("./BASINS_000000000001_sbrtest.dat", data.table=FALSE), NA)
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
                                         cl.stat.entropy = T, #cl.stat.denat = T,
                                         plot.cl.stat = T), NA)
  expect_error(bas <- basins_recognition(data = file.pi, nx = 7, match = F, plot = plt_stff, out.file = F, plot.cl.stat = plt_stff,
                                         new.dev = F, cl.stat.weight.barriers = T), NA)
  
  
  # ------------------------------------------------------- neuro tests
  # the following tests have been executed only 
  # extensive testing on real cases. It needs not to be
  # executed at all
  do_it <- FALSE
  if(do_it){
    fpi_empty <- data.table::fread('/disk2a/dgarolini/neuro/simulated_data/analysis/ISI_networks_end17/PROGIDX_000000000100.dat', data.table = F)
    fpi_best <- data.table::fread('/disk2a/dgarolini/neuro/simulated_data/analysis/ISI_networks_end17/PROGIDX_000000000102.dat', data.table = F)
    fpi_ave <- data.table::fread('/disk2a/dgarolini/neuro/simulated_data/analysis/ISI_networks_end17/PROGIDX_000000001102.dat', data.table = F)
    fpi_worst <- data.table::fread('/disk2a/dgarolini/neuro/simulated_data/analysis/ISI_networks_end17/PROGIDX_000000001281.dat', data.table = F)
  }else{
    fpi_empty <- fpi_best <- fpi_ave <- fpi_worst <- file.pi
  }
  
  # fpi_empty  
  expect_error(bas <- basins_recognition(data = fpi_empty, nx = 150, match = F, plot = T, out.file = F, new.dev = F, 
                                         cl.stat = T, cl.stat.denat = T, dbg = F,
                                         plot.cl.stat = T,
                                         cl.stat.nUni = c(5,10,15,20,25,30,40,50,60)), NA)
  if(do_it){
    # fpi_best  
    bas <- basins_recognition(data = fpi_best, nx = 150, match = T, plot = T, out.file = F, new.dev = F, 
                              cl.stat = T, cl.stat.denat = T, cl.stat.nUni = 20,
                              plot.cl.stat = T)
    # fpi_ave  
    bas <- basins_recognition(data = fpi_ave, nx = 150, match = T, plot = T, out.file = F, new.dev = F, 
                              cl.stat = T, cl.stat.denat = T, cl.stat.nUni = 20,
                              plot.cl.stat = T)
    # fpi_worst  
    bas <- basins_recognition(data = fpi_worst, nx = 150, match = T, plot = T, out.file = F, new.dev = F, 
                              cl.stat = T, cl.stat.denat = T, cl.stat.nUni = c(5,10,15,20,25,30,40,50,60),
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
