context('run_campari')

test_that('Test run_campari from installation', {
  
  silent <- F
  if(!silent) {require(testthat); require(CampaRi)} 
  # if(F) setwd('projects/CampaR/garbage/')
  
  # installing (fastly) campari SOMEWHERE
  expect_error(CampaRi::install_campari(install_ncminer = T, silent_built = F, no_optimization = T), NA) # to do it usually
  
  # setting the binary exe
  main_dir <- system.file('extdata/for_campari/', package = "CampaRi")
  bin_dir <- paste0(main_dir, '/', dir(main_dir)[1], '/bin/')
  ca_exe <- paste0(bin_dir, '/', dir(bin_dir)[1])
  print(ca_exe)
  print(list.files(ca_exe))
  expect_true(is.character(ca_exe))
  expect_true(file.exists(ca_exe))
  expect_true(ca_exe != "")
  
  # NBU.IN definition
  data.table::fwrite(list('NBU'), file = 'nbu.in', row.names = F, col.names = F)
  data.table::fwrite(list('END'), file = 'nbu.in', append = T, row.names = F, col.names = F)
  
  # --------------------------------------------------------------------------- nbu simulation
  expect_error(CampaRi::run_campari(FMCSC_SEQFILE="nbu.in", campari_exe = paste0(ca_exe, '/campari'), # you must have it defined according to CAMPARI's rules
                          # FMCSC_BASENAME="NBU", # lets try the base_name option
                          base_name = "NBU", print_status = T, # it will take 55 s in background ~
                          PARAMETERS="oplsaal.prm", # if this variable it is not supplied will be assigned to <full path to folder>/campari/params/abs3.2_opls.prm
                          FMCSC_SC_IPP=0.0,  FMCSC_PDBANALYZE = 0,
                          FMCSC_SC_BONDED_T=1.0,
                          FMCSC_DYNAMICS=3,
                          FMCSC_FRICTION=3.0,
                          FMCSC_TIMESTEP=0.005,
                          FMCSC_TEMP=400.0,
                          FMCSC_NRSTEPS=1000,
                          FMCSC_EQUIL=0,
                          FMCSC_XYZOUT=1,
                          FMCSC_XYZPDB=3,
                          FMCSC_TOROUT=1,
                          FMCSC_COVCALC=20000000,
                          FMCSC_SAVCALC=20000000,
                          FMCSC_POLCALC=20000000,
                          FMCSC_RHCALC=20000000,
                          FMCSC_INTCALC=20000000,
                          FMCSC_POLOUT=20000000,
                          FMCSC_ENSOUT=20000000,
                          FMCSC_ENOUT=20000000,
                          FMCSC_RSTOUT=20000000
  ), NA)
  
  # taking the dihedral angles
  trj <- data.table::fread("FYC.dat", header = F, skip = 1, data.table = FALSE)[,-1]
  trj <- sapply(trj, as.numeric) # always be sure that it is numeric!
  trj <- matrix(trj, nrow = 1000, ncol =3) # always be sure that it is numeric!
  
  # --------------------------------------------------------------------------- direct trj insertion
  expect_error(run_campari(trj = trj, base_name = "ascii_based_analysis", campari_exe = ca_exe,
                          FMCSC_CPROGINDMODE=1, #mst
                          FMCSC_CCOLLECT=1, print_status = T,
                          FMCSC_CMODE=4,
                          FMCSC_CDISTANCE=7, #rmsd without alignment 7 - dihedral distances need a complete analysis (pdb_format dcd pdb etc...) 
                          FMCSC_CPROGINDSTART=21, #starting snapshot 
                          # FMCSC_CPROGINDRMAX=1000, #search att
                          # FMCSC_BIRCHHEIGHT=2, #birch height
                          FMCSC_CMAXRAD=10880, #clustering
                          FMCSC_CRADIUS=10880,
                          FMCSC_CCUTOFF=10880,
                          FMCSC_CPROGINDWIDTH=1000), NA) 
  
  # --------------------------------------------------------------------------- tsv file
  expect_error(run_campari(data_file = 'ascii_based_analysis.tsv', base_name = "ascii_based_analysis", campari_exe = ca_exe,
                           FMCSC_CPROGINDMODE=1, #mst
                           FMCSC_CCOLLECT=1, print_status = T,
                           FMCSC_CMODE=4,
                           FMCSC_CDISTANCE=7, #rmsd without alignment 7 - dihedral distances need a complete analysis (pdb_format dcd pdb etc...) 
                           FMCSC_CPROGINDSTART=21, #starting snapshot 
                           # FMCSC_CPROGINDRMAX=1000, #search att
                           # FMCSC_BIRCHHEIGHT=2, #birch height
                           FMCSC_CMAXRAD=10880, #clustering
                           FMCSC_CRADIUS=10880,
                           FMCSC_CCUTOFF=10880,
                           FMCSC_CPROGINDWIDTH=1000), NA) 
  expect_error(run_campari(base_name = "ascii_based_analysis", campari_exe = ca_exe, FMCSC_NCDM_ASFILE = 'ascii_based_analysis.tsv', 
                           FMCSC_CPROGINDMODE=1, #mst
                           FMCSC_CCOLLECT=1, print_status = T,
                           FMCSC_CMODE=4,
                           FMCSC_CDISTANCE=7, #rmsd without alignment 7 - dihedral distances need a complete analysis (pdb_format dcd pdb etc...) 
                           FMCSC_CPROGINDSTART=21, #starting snapshot 
                           # FMCSC_CPROGINDRMAX=1000, #search att
                           # FMCSC_BIRCHHEIGHT=2, #birch height
                           FMCSC_CMAXRAD=10880, #clustering
                           FMCSC_CRADIUS=10880,
                           FMCSC_CCUTOFF=10880,
                           FMCSC_CPROGINDWIDTH=1000), NA)
  
  # --------------------------------------------------------------------------- dcd file
  expect_error(run_campari(data_file = 'NBU_traj.dcd', base_name = "ascii_based_analysis", campari_exe = ca_exe,
                           FMCSC_CPROGINDMODE=1, #mst
                           FMCSC_CCOLLECT=1, print_status = T,
                           FMCSC_CMODE=4,
                           FMCSC_CDISTANCE=5, #rmsd without alignment 7 - dihedral distances need a complete analysis (pdb_format dcd pdb etc...) 
                           FMCSC_CPROGINDSTART=21, #starting snapshot 
                           FMCSC_CMAXRAD=10880, #clustering
                           FMCSC_CRADIUS=10880,
                           FMCSC_CCUTOFF=10880,
                           FMCSC_CPROGINDWIDTH=1000))  
  expect_error(run_campari(data_file = 'NBU_traj.dcd', base_name = "ascii_based_analysis", campari_exe = ca_exe,
                           FMCSC_CPROGINDMODE=1, #mst
                           FMCSC_NCDM_NRFRMS= 14,
                           FMCSC_CCOLLECT=1, print_status = T,
                           FMCSC_CMODE=4,
                           FMCSC_CDISTANCE=5, #rmsd without alignment 7 - dihedral distances need a complete analysis (pdb_format dcd pdb etc...) 
                           FMCSC_CPROGINDSTART=21, #starting snapshot 
                           FMCSC_CMAXRAD=10880, #clustering
                           FMCSC_CRADIUS=10880,
                           FMCSC_CCUTOFF=10880,
                           FMCSC_CPROGINDWIDTH=1000))  
  expect_error(run_campari(data_file = 'NBU_traj.dcd', base_name = "ascii_based_analysis", seq_in = 'nbu.in',
                           campari_exe = ca_exe,
                           FMCSC_CPROGINDMODE=2, #mst
                           FMCSC_NRSTEPS = 1000,
                           FMCSC_CCOLLECT=1, print_status = T,
                           FMCSC_CMODE=4,
                           FMCSC_CDISTANCE=5, #rmsd without alignment 7 - dihedral distances need a complete analysis (pdb_format dcd pdb etc...) 
                           FMCSC_CPROGINDSTART=21, #starting snapshot 
                           FMCSC_CMAXRAD=10880, #clustering
                           FMCSC_CRADIUS=1,
                           FMCSC_CPROGINDRMAX=10, #search att
                           FMCSC_CCUTOFF=1,
                           FMCSC_CPROGINDWIDTH=1000), NA) 
  expect_error(run_campari(base_name = "ascii_based_analysis", seq_in = 'nbu.in',
                           campari_exe = ca_exe, FMCSC_DCDFILE = 'NBU_traj.dcd', FMCSC_PDBANALYZE = 1,
                           FMCSC_CPROGINDMODE=2, #mst
                           FMCSC_NRSTEPS = 1000,
                           FMCSC_CCOLLECT=1, print_status = T,
                           FMCSC_CMODE=4,
                           FMCSC_CDISTANCE=5, #rmsd without alignment 7 - dihedral distances need a complete analysis (pdb_format dcd pdb etc...) 
                           FMCSC_CPROGINDSTART=21, #starting snapshot 
                           FMCSC_CMAXRAD=10880, #clustering
                           FMCSC_CRADIUS=1,
                           FMCSC_CPROGINDRMAX=10, #search att
                           FMCSC_CCUTOFF=1,
                           FMCSC_CPROGINDWIDTH=1000), NA)
  
  # --------------------------------------------------------------------------- bash_rc append for further search of ca_exe
  expect_error(if(!file.copy(from = '~/.bashrc', to = '~/.bashrc_tmp', overwrite = T)) stop('error'), NA)
  system(paste0('echo "export PATH=$PATH:', ca_exe, '" >> ~/.bashrc')) # just not to mess up with my bashrc
  expect_error(run_campari(trj = trj, base_name = "ascii_based_analysis", 
                           FMCSC_CPROGINDMODE=1, #mst
                           FMCSC_CCOLLECT=1, print_status = T,
                           FMCSC_CMODE=4,
                           FMCSC_CDISTANCE=7, #rmsd without alignment 7 - dihedral distances need a complete analysis (pdb_format dcd pdb etc...) 
                           FMCSC_CPROGINDSTART=21, #starting snapshot 
                           # FMCSC_CPROGINDRMAX=1000, #search att
                           # FMCSC_BIRCHHEIGHT=2, #birch height
                           FMCSC_CMAXRAD=10880, #clustering
                           FMCSC_CRADIUS=10880,
                           FMCSC_CCUTOFF=10880,
                           FMCSC_CPROGINDWIDTH=1000), NA) 
  expect_error(if(!file.copy(from = '~/.bashrc_tmp', to = '~/.bashrc', overwrite = T)) stop('error'), NA)
  file.copy(from = '~/.bashrc_tmp', to = '~/.bashrc', overwrite = T)
  
  
  
  if(file.exists('Makefile')) file.remove('Makefile')
  if(file.exists('VERSION')) file.remove('VERSION')
  if(file.exists('FYC.dat')) file.remove('FYC.dat')
  if(file.exists('nbu.in')) file.remove('nbu.in')
  if(file.exists('config.log')) file.remove('config.log')
  if(file.exists('config.status')) file.remove('config.status')
  if(file.exists('FRAMES_NBL.nc')) file.remove('FRAMES_NBL.nc')
  if(file.exists('PROGIDX_000000000021.dat')) file.remove('PROGIDX_000000000021.dat')
  file_to_d <- c(list.files(pattern = 'STRUCT_CLUSTERING*'), list.files(pattern = 'ascii_based_analysis.*'), list.files(pattern = 'NBU.*'))
  for(i in file_to_d) {if(file.exists(i)) file.remove(i)}
  
  if(F){
    # some extra tests
    ddir <- '../CampaRi/'
    a <- paste0(ddir, '/', list.files(ddir, recursive = T))
    a <- a[!grepl(pattern = '/doc/', x = a)]; a
    a <- a[!grepl(pattern = '/lib/', x = a)]; a  
    a <- a[!grepl(pattern = '/bin/', x = a)]; a  
    a <- a[!grepl(pattern = '/to_d/', x = a)]; a  
    # a <- a[!grepl(pattern = '/src/', x = a)]  # to check the src!!
    b <- lapply(a, function(x) capture.output(tools::showNonASCII(readLines(x)), type = "message"))
    b <- lapply(a, function(x) grep(pattern = 'utf', x = readLines(x), fixed = T))
    lb <- sapply(b, length)
    findf <- data.frame('n.nonA' = lb[lb != 0], 'file.name' = a[lb != 0])
    findf[order(findf$n.nonA),]
    
    
    # remove the mod etc...
    a[grepl(pattern = '/bin/', x = a)]
    a[grepl(pattern = '/lib/', x = a)]
    system('rm -rf inst/extdata/for_campari/bin')
    system('rm -rf inst/extdata/for_campari/lib')
    a <- list.files('inst/extdata/for_campari', recursive = T)
    a[grepl(pattern = '.mod$', x = a)]
    file.remove(paste0('inst/extdata/for_campari/', a[grepl(pattern = '.mod$', x = a)]))
    file.remove(paste0('inst/extdata/for_campari/source/DEPENDENCIES'))
    file.remove(paste0('inst/extdata/for_campari/source/config.guess'))
    file.remove(paste0('inst/extdata/for_campari/source/config.sub'))
    file.remove(paste0('inst/extdata/for_campari/source/install-sh'))
  }
  
})