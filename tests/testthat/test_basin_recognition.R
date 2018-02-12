context('basins_recognition')

test_that('Test for basin recognition with ext files', {

    ## dir.cur <- getwd()
    ## str.st <- rev(gregexpr('CampaRi', dir.cur)[[1]])[1]
    ## dir.root <- substr(dir.cur, 1, str.st + 7)
    file.pi <- system.file("extdata", "PROGIDX_000000000001_sbrtest.dat", package = "CampaRi")
    ## data.pi <- data.table::fread(paste0(dir.root, "inst/extdata/PROGIDX_000000000001_sbrtest.dat"), data.table=FALSE)[,c(1,3,4,6,7,5)]
    data.pi <- data.table::fread(file.pi, data.table=FALSE)[,c(1,3,4,6,7,5)]
    

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
        expect_error(tmp <- basins_recognition(data=data.pi, nx=nxbb[ii], ny=nybb[ii], ny.aut=F, local.cut=lcsb[ii], 
                                               match=matb[ii], dyn.check=2, avg.opt=as.character(avgb[ii]), plot=F, out.file=(ii==1), silent=T, pol.degree=3), NA)
        if(ii==1) {
            expect_error(out <- data.table::fread("BASINS_1.dat", data.table=FALSE), NA)
            expect_true(all(seq(nrow(tmp$tab.st)) == sort(unique(out[,3]))))
        } 
    }
    detach(parspace)

    ## Testing other less relevant parameter
    expect_error(tmp <- basins_recognition(data=file.pi, nx=parspace$nxbb[1], ny=parspace$nybb[1], 
                                           ny.aut=T, local.cut=F, match=T, dyn.check=2, avg.opt="SG", plot=T, out.file=T, silent=F, pol.degree=3), NA)
    ## list.files()
    expect_error(out <- data.table::fread("./BASINS_000000000001_sbrtest.dat", data.table=FALSE), NA)
    expect_true(all(seq(nrow(tmp$tab.st)) == sort(unique(out[,3]))))
    expect_warning(tmp <- basins_recognition(data=file.pi, nx=parspace$nxbb[1], ny=parspace$nybb[1], only.kin=T, local.cut=F, match=T, dyn.check=1, avg.opt="SG", plot=F, out.file=F, silent=F, pol.degree=1.5))
    
    # testing for minimum number of bins without crashing
    pi.table <- data.table::fread(file.pi, verbose = F, data.table = F)
    expect_error(CampaRi::basins_recognition(pi.table[1:50,], nx = 7, new.dev = F, plot = T, silent = T))
    expect_error(CampaRi::basins_recognition(pi.table[1:80,], nx = 6, new.dev = F, plot = T, silent = T))
    
    cat("Test on basins_recognition done\n\n")
    if(file.exists('BASINS_000000000001_sbrtest.dat')) file.remove('BASINS_000000000001_sbrtest.dat')
    if(file.exists('BASINS_1.dat')) file.remove('BASINS_1.dat')
    
})
