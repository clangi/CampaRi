context('run_campari')

test_that('Test run_campari from installation', {
  silent <- F
  plt_stff <- !silent
  if(!silent) {require(testthat); require(CampaRi)} 
  if(!dir.exists('to_delete_just_in_a_moment')) dir.create('to_delete_just_in_a_moment')
  expect_error(install_campari(installation_location = 'to_delete_just_in_a_moment', install_ncminer = T, install_threads = F, silent_built = F), NA)
  # expect_error(install_campari(install_ncminer = T, silent_built = silent), NA) # this way is broken if the wd is not writable
  # bin_dir <- system.file('extdata/for_campari/bin/', package = "CampaRi")
  # ca_exe <- paste0(bin_dir, dir(bin_dir)[1], '/', list.files(paste0(bin_dir, dir(bin_dir)[1]))[2])
  ca_exe <- system('ls to_delete_just_in_a_moment/for_campari/bin/*/campari', intern = T)[1]
  dir <- system(paste0('pwd ', ca_exe), intern = T)[1]
  ca_exe <- paste0(dir, '/', ca_exe)
  expect_true(is.character(ca_exe))
  expect_true(file.exists(ca_exe))
  # system('printf "NBU\nEND" &> nbu.in')
  data.table::fwrite(list('NBU'), file = 'nbu.in', row.names = F, col.names = F, verbose = !silent)
  data.table::fwrite(list('END'), file = 'nbu.in', append = T, row.names = F, col.names = F, verbose = !silent)
  expect_error(run_campari(FMCSC_SEQFILE="nbu.in", campari_exe = ca_exe, # you must have it defined according to CAMPARI's rules
                          # FMCSC_BASENAME="NBU", # lets try the base_name option
                          base_name = "NBU", print_status = !silent, # it will take 55 s in background ~
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
                          FMCSC_RSTOUT=20000000, silent = silent
  ), NA)
  
  trj <- data.table::fread("FYC.dat", header = F, skip = 1, data.table = FALSE, verbose = !silent)[,-1]
  trj <- sapply(trj, as.numeric) # always be sure that it is numeric!
  trj <- matrix(trj, nrow = 1000, ncol =3) # always be sure that it is numeric!
  expect_error(run_campari(trj = trj, base_name = "ascii_based_analysis", campari_exe = ca_exe,
                           FMCSC_CPROGINDMODE=1, #mst
                           FMCSC_CCOLLECT=1, print_status = !silent,
                           FMCSC_CMODE=4,
                           FMCSC_CDISTANCE=7, #rmsd without alignment 7 - dihedral distances need a complete analysis (pdb_format dcd pdb etc...) 
                           FMCSC_CPROGINDSTART=21, #starting snapshot 
                           # FMCSC_CPROGINDRMAX=1000, #search att
                           # FMCSC_BIRCHHEIGHT=2, #birch height
                           FMCSC_CMAXRAD=10880, #clustering
                           FMCSC_CRADIUS=10880,
                           FMCSC_CCUTOFF=10880, 
                           FMCSC_CPROGINDWIDTH=1000, silent = silent
                           ), NA) #local cut is automatically adjusted to 1/10 if it is too big (as here)
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
                           FMCSC_CPROGINDWIDTH=1000, silent = silent), NA)
  
  
  
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
  if(dir.exists('to_delete_just_in_a_moment')) unlink('to_delete_just_in_a_moment', recursive = T)
})