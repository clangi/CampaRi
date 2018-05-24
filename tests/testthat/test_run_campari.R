context('run_campari')

test_that('Test run_campari from installation', {
  if(F) setwd('projects/CampaR/garbage/')
  expect_error(CampaRi::install_campari(install_ncminer = T, silent_built = T), NA) # to do it usually
  bin_dir <- system.file('extdata/for_campari/bin/', package = "CampaRi")
  ca_exe <- paste0(bin_dir, '/', dir(bin_dir)[1])
  print(ca_exe)
  expect_true(is.character(ca_exe))
  expect_true(file.exists(ca_exe))
  expect_true(ca_exe != "")
  # system('printf "NBU\nEND" &> nbu.in')
  data.table::fwrite(list('NBU'), file = 'nbu.in', row.names = F, col.names = F)
  data.table::fwrite(list('END'), file = 'nbu.in', append = T, row.names = F, col.names = F)
  expect_error(CampaRi::run_campari(FMCSC_SEQFILE="nbu.in", campari_exe = ca_exe, # you must have it defined according to CAMPARI's rules
                          # FMCSC_BASENAME="NBU", # lets try the base_name option
                          base_name = "NBU", print_status = T, # it will take 55 s in background ~
                          PARAMETERS="oplsaal.prm", # if this variable it is not supplied will be automatically assigned to <full path to folder>/campari/params/abs3.2_opls.prm
                          FMCSC_SC_IPP=0.0,
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
  
  trj <- data.table::fread("FYC.dat", header = F, skip = 1, data.table = FALSE)[,-1]
  trj <- sapply(trj, as.numeric) # always be sure that it is numeric!
  trj <- matrix(trj, nrow = 1000, ncol =3) # always be sure that it is numeric!
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
                           FMCSC_CPROGINDWIDTH=1000), NA) #local cut is automatically adjusted to 1/10 if it is too big (as here)
  # debugonce(run_campari)
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
  
  
  
  if(file.exists('Makefile')) file.remove('Makefile')
  if(file.exists('VERSION')) file.remove('VERSION')
  if(file.exists('FYC.dat')) file.remove('FYC.dat')
  if(file.exists('NBU.key')) file.remove('NBU.key')
  if(file.exists('NBU.log')) file.remove('NBU.log')
  if(file.exists('nbu.in')) file.remove('nbu.in')
  if(file.exists('NBU_END.int')) file.remove('NBU_END.int')
  if(file.exists('NBU_END.pdb')) file.remove('NBU_END.pdb')
  if(file.exists('NBU_START.int')) file.remove('NBU_START.int')
  if(file.exists('NBU_START.pdb')) file.remove('NBU_START.pdb')
  if(file.exists('NBU_VIS.vmd')) file.remove('NBU_VIS.vmd')
  if(file.exists('NBU_traj.dcd')) file.remove('NBU_traj.dcd')
  if(file.exists('FRAMES_NBL.nc')) file.remove('FRAMES_NBL.nc')
  if(file.exists('ascii_based_analysis.tsv')) file.remove('ascii_based_analysis.tsv')
  if(file.exists('ascii_based_analysis.log')) file.remove('ascii_based_analysis.log')
  if(file.exists('ascii_based_analysis.key')) file.remove('ascii_based_analysis.key')
  if(file.exists('PROGIDX_000000000021.dat')) file.remove('PROGIDX_000000000021.dat')
  
  
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
  }
  
})