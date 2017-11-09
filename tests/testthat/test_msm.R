context('basins_recognition')

test_that('Test for basin recognition with ext files', {

    cat("TESTING MSM function\n")
    
    file.pi <- system.file("extdata", "PROGIDX_000000000001_sbrtest.dat", package = "CampaRi")
    data.pi <- data.table::fread(file.pi, data.table=FALSE)[,c(1,3,4,6,7,5)]
    expect_error(seqq <- basins_recognition(data=data.pi, nx=50, ny=20, local.cut=F, match=F, dyn.check=1, avg.opt="SG", plot=F, out.file=F, silent=F)$seq.st, NA)


    ## basins_recognition(data.pi, nx=50, ny=50, match=F, plot=T, out.file=F, silent=F)
    out1 <- expect_error(msm(seqq, 10, tm.opt="mle", CK.test=T, CK.lags=c(1:10), silent=F, eig.plot=T, dev.new.eig=T), NA)
    out2 <- expect_error(msm(seqq, 1, tm.opt="symm", CK.test=T, setA=c(1,2), silent=T, dev.new.CK=T), NA)

    cat("Test on msm_recognition done\n\n")
})
