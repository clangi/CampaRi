#' @title plotting
#' @description
#'      plotting from file(X)
#'
#' @param fl file location
#' @param n_fold Number of links contractions (folds)
#' Default: \code{0}
#'
#'
#' @return tree: degree list, connectivity matrix and weights
#'
#' @export campari
campari<-function(nsnaps, wd, data_file, base_name, camp_home, 
                  cprogindstart = 1, pdb_format = 4, distance_met = 5,
                  birch_height = 5, cmaxrad = 2147483647, cradius = 2147483647,
                  cprogindwidth = 1,search_attempts = NULL, methodst = 1){
  # dirs/file definitions
  kfile <- paste0(base_name,".key")
  klog <-paste0(base_name,".log")
  seq_in <- paste0(base_name,".in")
  
  kchar<-readChar(paste0(wd,kfile), file.info(paste0(wd,kfile))$size)
  
  FMCSC_PDBANALYZE <- 1
  PARAMETERS <- paste0(camp_home,"params/abs3.2_opls.prm")  # file defining system energies. Irrelevant fuer blosse Analyse.
  FMCSC_BBSEGFILE <- paste0(camp_home,"data/bbseg2.dat")    # lookup table for secondary structure measures
  FMCSC_SEQFILE  <- paste0(wd,seq_in)                      # input file that defines the sequence of the molecule(s)
  
  FMCSC_PDB_FORMAT <- pdb_format 
  if(FMCSC_PDB_FORMAT==1) 
    FMCSC_PDBFILE <- paste0(wd,data_file)
  if(FMCSC_PDB_FORMAT==3) 
    FMCSC_XTCFILE <- paste0(wd,data_file)
  if(FMCSC_PDB_FORMAT==4) 
    FMCSC_DCDFILE <- paste0(wd,data_file)
  
  
  #1. CAMPARI expects a single trajectory file in pdb-format using the MODEL /ENDMDL syntax to denote the individual snapshots.
  #2. CAMPARI expects to find multiple pdb files with one snapshot each that are systematically numbered starting from the file provided via PDBFILE.
  #3. CAMPARI expects to find a single trajectory in binary xtc-format (GROMACS style).
  #4. CAMPARI expects to find a single trajectory in binary dcd-format (CHARMM/NAMD style).
  #5. CAMPARI expects to find a single trajectory in binary NetCDF-format (AMBER style). (reference)
  #Note that .xtc, .nc, and .dcd trajectory files are not annotated and that the order of atoms between
  #the file and CAMPARI's inner workings must be consistent. Since this is almost never true for binary
  #trajectory files obtained with other software, CAMPARI offers the user to provide a pdb template
  #which contains the order of atoms in the binary file in annotated form (see PDB_TEMPLATE).
  
  if(is.null(search_attempts)&&methodst==2) search_attempts <- nsnaps/10
  else search_attempts <- nsnaps / 10
  
  FMCSC_PDB_R_CONV <- 1 # different conventions for the formatting of PDB files.
  #1. CAMPARI
  #2. GROMOS
  #3. CHARMM
  #4. AMBER
  FMCSC_NRSTEPS <- nsnaps
  FMCSC_CDISTANCE <- distance_met # seq(0,10,1) #1-5 crashed because torsional
  
  FMCSC_BIRCHHEIGHT <- birch_height #12 #nothing is changing 1,51,2
  FMCSC_CMAXRAD <- cmaxrad  #nothing is changing 0,5,0.5
  FMCSC_CRADIUS <- cradius #nothing is changing seq(0.1,5.1,0.2)
  FMCSC_CCUTOFF <- 2147483647  #nothing seq(0,10) chardcut is different from cmaxrad which is used only in the case of sst (root lev)
  
  FMCSC_CPROGINDMODE <- methodst #exact MST = 1
  FMCSC_CPROGINDRMAX <- search_attempts #ONLY SST
  FMCSC_CPROGINDSTART <- cprogindstart #nothing is changing c(-1,seq(1,1801,200))
  FMCSC_CPROGINDWIDTH <-  cprogindwidth # 540
  
  FMCSC_CPROGMSTFOLD <- 0 # tree contraction
  FMCSC_PCAMODE      <-1   #no PCA is performed
  FMCSC_CREDUCEDIM   <-2 #ONLY PCA
  # this keyword allows the user to elect to run the clustering algorithm (→ CMODE) on a dataset of r
  # educed dimensionality that corresponds to the first NV data vectors in the transformed space, w
  # here NV is set by the choice for this keyword.
  FMCSC_CALIGN <- 0
  # this keyword can be used to specifically disable the alignment step that occurs before 
  # the actual RMSD of the two coordinate sets is computed. To achieve this, provide any 
  # value other than 1 (the default) for this on/off-type keyword.
  
  
  
  kchar<-sub(pattern = "PARAMETERS\\s[A-Za-z0-9_/.~]+", replacement = paste0("PARAMETERS ", PARAMETERS), x = kchar)
  kchar<-sub(pattern = "FMCSC_PDBANALYZE\\s[0-9]+", replacement = paste0("FMCSC_PDBANALYZE ", FMCSC_PDBANALYZE), x = kchar)
  kchar<-sub(pattern = "FMCSC_BBSEGFILE\\s[A-Za-z0-9_/.~]+", replacement = paste0("FMCSC_BBSEGFILE ", FMCSC_BBSEGFILE), x = kchar)
  kchar<-sub(pattern = "FMCSC_SEQFILE\\s[A-Za-z0-9_/.~]+", replacement = paste0("FMCSC_SEQFILE ", FMCSC_SEQFILE), x = kchar)
  if(FMCSC_PDB_FORMAT==1) 
    #if there is bla bla if not bla bla
    kchar<-sub(pattern = "FMCSC_PDBFILE\\s[A-Za-z0-9_/.~]+", replacement = paste0("FMCSC_PDBFILE ", FMCSC_PDBFILE), x = kchar)
  else
    kchar<-sub(pattern = "FMCSC_PDBFILE\\s[A-Za-z0-9_/.~]+", replacement = "", x = kchar)
  if(FMCSC_PDB_FORMAT==3) 
    kchar<-sub(pattern = "FMCSC_XTCFILE\\s[A-Za-z0-9_/.~]+", replacement = paste0("FMCSC_XTCFILE ", FMCSC_XTCFILE), x = kchar)
  else
    kchar<-sub(pattern = "FMCSC_XTCFILE\\s[A-Za-z0-9_/.~]+", replacement = "", x = kchar)
  if(FMCSC_PDB_FORMAT==4) 
    kchar<-sub(pattern = "FMCSC_DCDFILE\\s[A-Za-z0-9_/.~]+", replacement = paste0("FMCSC_DCDFILE ", FMCSC_DCDFILE), x = kchar)
  else
    kchar<-sub(pattern = "FMCSC_DCDFILE\\s[A-Za-z0-9_/.~]+", replacement = "", x = kchar)
  kchar<-sub(pattern = "FMCSC_PDB_FORMAT\\s[0-9]+", replacement = paste0("FMCSC_PDB_FORMAT ", FMCSC_PDB_FORMAT), x = kchar)
  kchar<-sub(pattern = "FMCSC_PDB_R_CONV\\s[0-9]+", replacement = paste0("FMCSC_PDB_R_CONV ", FMCSC_PDB_R_CONV), x = kchar)
  kchar<-sub(pattern = "FMCSC_NRSTEPS\\s[0-9]+", replacement = paste0("FMCSC_NRSTEPS ", FMCSC_NRSTEPS), x = kchar)
  kchar<-sub(pattern = "FMCSC_CDISTANCE\\s[0-9]+", replacement = paste0("FMCSC_CDISTANCE ",FMCSC_CDISTANCE),x = kchar)
  kchar<-sub(pattern = "FMCSC_BIRCHHEIGHT\\s[0-9]+", replacement = paste0("FMCSC_BIRCHHEIGHT ",FMCSC_BIRCHHEIGHT),x = kchar)
  kchar<-sub(pattern = "FMCSC_CMAXRAD\\s([0-9]+.[0-9]+|[0-9]+)",replacement = paste0("FMCSC_CMAXRAD ",FMCSC_CMAXRAD),x = kchar)
  kchar<-sub(pattern = "FMCSC_CRADIUS\\s([0-9]+.[0-9]+|[0-9]+)",replacement = paste0("FMCSC_CRADIUS ",FMCSC_CRADIUS),x = kchar)
  kchar<-sub(pattern = "FMCSC_CCUTOFF\\s([0-9]+.[0-9]+|[0-9]+)",replacement = paste0("FMCSC_CCUTOFF ",FMCSC_CCUTOFF),x = kchar)
  kchar<-sub(pattern = "FMCSC_CPROGINDMODE\\s[0-9]+",replacement = paste0("FMCSC_CPROGINDMODE ",FMCSC_CPROGINDMODE),x = kchar)
  kchar<-sub(pattern = "FMCSC_CPROGINDRMAX\\s[0-9]+",replacement = paste0("FMCSC_CPROGINDRMAX ",FMCSC_CPROGINDRMAX),x = kchar)
  kchar<-sub(pattern = "FMCSC_CPROGINDSTART\\s([0-9]+|-[0-9]+)",replacement = paste0("FMCSC_CPROGINDSTART ",FMCSC_CPROGINDSTART),x = kchar)
  kchar<-sub(pattern = "FMCSC_CPROGINDWIDTH\\s[0-9]+",replacement = paste0("FMCSC_CPROGINDWIDTH ",FMCSC_CPROGINDWIDTH),x = kchar)
  kchar<-sub(pattern = "FMCSC_CPROGMSTFOLD\\s[0-9]+",replacement = paste0("FMCSC_CPROGMSTFOLD ",FMCSC_CPROGMSTFOLD),x = kchar)
  kchar<-sub(pattern = "FMCSC_PCAMODE\\s[0-9]+",replacement = paste0("FMCSC_PCAMODE ",FMCSC_PCAMODE),x = kchar)
  kchar<-sub(pattern = "FMCSC_CREDUCEDIM\\s[0-9]+",replacement = paste0("FMCSC_CREDUCEDIM ",FMCSC_CREDUCEDIM),x = kchar)
  kchar<-sub(pattern = "FMCSC_CALIGN\\s[0-9]+",replacement = paste0("FMCSC_CALIGN ",FMCSC_CALIGN),x = kchar)
  fileConn<-file(kfile)
  writeLines(kchar, fileConn)
  close(fileConn)
  
  system(paste0(camp_home,"bin/x86_64/campari -k ",kfile, paste0(">& ",klog)))
  system(paste0("tail ",klog))
}
