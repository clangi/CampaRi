#' @title Original campari run
#' @description
#'      \code{campari} will call the original campari software, specified in the directory camp_home. 
#'      It is possible to use this function following the instructions in \url{http://campari.sourceforge.net/tutorial11.html}.
#'      For that purpose a keyfile and sequence file are present inside this package as an example (in \code{CampaRi/inst/extdata}).
#'
#' @param nsnaps Number of snapshots in the trajectory file.
#' @param working_dir Working directory.
#' @param data_file Input file (e.g. \code{trajectory.dcd}) location.
#' @param base_name This name should be used for the input files (\code{base_name.key} and \code{base_name.in}) in order to run campari.
#' @param camp_home Location of the installed campari software.
#' @param ... Analysis variables (similarly to \code{\link{mst_from_trj}}).
#' @details For details, please refer to the main documentation of the original campari software \url{http://campari.sourceforge.net/documentation.html}.
#' @seealso
#' \code{\link{mst_from_trj}}, \code{\link{gen_progindex}}, \code{\link{gen_annotation}}.
#' @examples 
#' \dontrun{
#' campari(nsnaps = 100, working_dir = "./", data_file = "file_dcd", camp_home = "/campari_home/",
#'         base_name = "nbu", pdb_format = 4,
#'         cprogindstart = 2,distance_met = 5, birch_height = 5, cmaxrad = 1, cradius = 0.1,
#'         cprogindwidth = 10,search_attempts = 10, methodst = 2)
#' }
#'
#' @export campari
campari<-function(nsnaps, working_dir, data_file, base_name, camp_home, ...){
  input_list <- list(...)
  if(!"cprogindstart" %in% names(input_list)) cprogindstart = 1
  else cprogindstart = as.numeric(input_list["cprogindstart"])
  if(!"pdb_format" %in% names(input_list)) pdb_format = 4
  else pdb_format = as.numeric(input_list["pdb_format"])
  if(!"distance_met" %in% names(input_list)) distance_met = 5
  else distance_met = as.numeric(input_list["distance_met"])
  if(!"birch_height" %in% names(input_list)) birch_height = 5
  else birch_height = as.numeric(input_list["birch_height"])
  if(!"cmaxrad" %in% names(input_list)) cmaxrad = 214748364
  else cmaxrad = as.numeric(input_list["cmaxrad"])
  if(!"cradius" %in% names(input_list)) cradius = 214748364
  else cradius = as.numeric(input_list["cradius"])
  if(!"cprogindwidth" %in% names(input_list)) cprogindwidth = as.integer(floor(nsnaps-1)/2)
  else cprogindwidth = as.numeric(input_list["cprogindwidth"])
  if(!"search_attempts" %in% names(input_list)) search_attempts = NULL
  else search_attempts = as.numeric(input_list["search_attempts"])
  if(!"methodst" %in% names(input_list)) methodst = 1
  else methodst = as.numeric(input_list["methodst"])
  # dirs/file definitions
  kfile <- paste0(working_dir,"/",base_name,".key")
  klog <-paste0(working_dir,"/",base_name,".log")
  if(!file.exists(kfile)) stop(kfile, " not in the working directory. ")

  
  kchar<-readChar(kfile, file.info(kfile)$size)
  
  FMCSC_PDBANALYZE <- 1
  PARAMETERS <- paste0(camp_home,"/params/abs3.2_opls.prm")  # file defining system energies. Irrelevant fuer blosse Analyse.
  FMCSC_BBSEGFILE <- paste0(camp_home,"/data/bbseg2.dat")    # lookup table for secondary structure measures
  if(!file.exists(PARAMETERS)) stop(PARAMETERS, " not present.")
  if(!file.exists(FMCSC_BBSEGFILE)) stop(FMCSC_BBSEGFILE, " not present.")
  
  FMCSC_PDB_FORMAT <- pdb_format 
  # pdb file
  if(FMCSC_PDB_FORMAT==1){
    FMCSC_PDBFILE <- paste0(working_dir, "/", data_file)
    seq_in <- paste0(working_dir,"/",strsplit(data_file, ".pdb")[[1]],".in")
  }
  # xtc file
  if(FMCSC_PDB_FORMAT==3) { 
    FMCSC_XTCFILE <- paste0(working_dir, "/", data_file)
    seq_in <- paste0(working_dir,"/",strsplit(data_file, ".xtc")[[1]],".in")
  }
  # dcd file
  if(FMCSC_PDB_FORMAT==4){
    FMCSC_DCDFILE <- paste0(working_dir, "/", data_file)
    seq_in <- paste0(working_dir,"/",strsplit(data_file, ".dcd")[[1]],".in")
  } 
  
  if(!file.exists(seq_in)) stop(seq_in, " not in the working directory. ")
  FMCSC_SEQFILE  <- seq_in                                  # input file that defines the sequence of the molecule(s)
  
  
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
  FMCSC_CCUTOFF <- 214748364  #nothing seq(0,10) chardcut is different from cmaxrad which is used only in the case of sst (root lev)
  
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
  
  first <- TRUE
  
  kchar<-sub(pattern = "PARAMETERS[ ]{1,}[A-Za-z0-9_/.~]+", replacement = paste0("PARAMETERS ", PARAMETERS), x = kchar)
  kchar<-sub(pattern = "FMCSC_PDBANALYZE[ ]{1,}[0-9]+", replacement = paste0("FMCSC_PDBANALYZE ", FMCSC_PDBANALYZE), x = kchar)
  kchar<-sub(pattern = "FMCSC_BBSEGFILE[ ]{1,}[A-Za-z0-9_/.~]+", replacement = paste0("FMCSC_BBSEGFILE ", FMCSC_BBSEGFILE), x = kchar)
  kchar<-sub(pattern = "FMCSC_SEQFILE[ ]{1,}[A-Za-z0-9_/.~]+", replacement = paste0("FMCSC_SEQFILE ", FMCSC_SEQFILE), x = kchar)
  
  # putting flagger to use the correct key word
  if(grepl(pattern = "FMCSC_PDBFILE[ ]{1,}[A-Za-z0-9_/.~]+", x = kchar)){
    if(first) kchar<-sub(pattern = "FMCSC_PDBFILE[ ]{1,}[A-Za-z0-9_/.~]+", replacement = "FMCSC_theFILE", x = kchar)
    else kchar<-sub(pattern = "FMCSC_PDBFILE[ ]{1,}[A-Za-z0-9_/.~]+", replacement = "", x = kchar)
    first <- FALSE 
  }
  if(grepl(pattern = "FMCSC_XTCFILE[ ]{1,}[A-Za-z0-9_/.~]+", x = kchar)){
    if(first) kchar<-sub(pattern = "FMCSC_XTCFILE[ ]{1,}[A-Za-z0-9_/.~]+", replacement = "FMCSC_theFILE", x = kchar)
    else kchar<-sub(pattern = "FMCSC_XTCFILE[ ]{1,}[A-Za-z0-9_/.~]+", replacement = "", x = kchar)
    first <- FALSE 
  }
  if(grepl(pattern = "FMCSC_DCDFILE[ ]{1,}[A-Za-z0-9_/.~]+", x = kchar)){
    if(first) kchar<-sub(pattern = "FMCSC_DCDFILE[ ]{1,}[A-Za-z0-9_/.~]+", replacement = "FMCSC_theFILE", x = kchar)
    else kchar<-sub(pattern = "FMCSC_DCDFILE[ ]{1,}[A-Za-z0-9_/.~]+", replacement = "", x = kchar)
    first <- FALSE 
  }
  
  if(first) kchar<-sub(pattern = "\n", replacement = "\nFMCSC_theFILE\n", x = kchar)
  
  
  if(FMCSC_PDB_FORMAT==1) 
    kchar<-sub(pattern = "FMCSC_theFILE", replacement = paste0("FMCSC_PDBFILE ", FMCSC_PDBFILE), x = kchar)
  else if(FMCSC_PDB_FORMAT==3) 
    kchar<-sub(pattern = "FMCSC_theFILE", replacement = paste0("FMCSC_XTCFILE ", FMCSC_XTCFILE), x = kchar)
  else if(FMCSC_PDB_FORMAT==4) 
    kchar<-sub(pattern = "FMCSC_theFILE", replacement = paste0("FMCSC_DCDFILE ", FMCSC_DCDFILE), x = kchar)
  else
    stop("File input not supported. We have 1/3/4")
  
  
  kchar<-sub(pattern = "FMCSC_PDB_FORMAT[ ]{1,}[0-9]+", replacement = paste0("FMCSC_PDB_FORMAT ", FMCSC_PDB_FORMAT), x = kchar)
  kchar<-sub(pattern = "FMCSC_PDB_R_CONV[ ]{1,}[0-9]+", replacement = paste0("FMCSC_PDB_R_CONV ", FMCSC_PDB_R_CONV), x = kchar)
  kchar<-sub(pattern = "FMCSC_NRSTEPS[ ]{1,}[0-9]+", replacement = paste0("FMCSC_NRSTEPS ", FMCSC_NRSTEPS), x = kchar)
  kchar<-sub(pattern = "FMCSC_CDISTANCE[ ]{1,}[0-9]+", replacement = paste0("FMCSC_CDISTANCE ",FMCSC_CDISTANCE),x = kchar)
  kchar<-sub(pattern = "FMCSC_BIRCHHEIGHT[ ]{1,}[0-9]+", replacement = paste0("FMCSC_BIRCHHEIGHT ",FMCSC_BIRCHHEIGHT),x = kchar)
  kchar<-sub(pattern = "FMCSC_CMAXRAD[ ]{1,}(-?[0-9]+.[0-9]+|[0-9]+)",replacement = paste0("FMCSC_CMAXRAD ",FMCSC_CMAXRAD),x = kchar)
  kchar<-sub(pattern = "FMCSC_CRADIUS[ ]{1,}(-?[0-9]+.[0-9]+|[0-9]+)",replacement = paste0("FMCSC_CRADIUS ",FMCSC_CRADIUS),x = kchar)
  kchar<-sub(pattern = "FMCSC_CCUTOFF[ ]{1,}([-?0-9]+.[0-9]+|[0-9]+)",replacement = paste0("FMCSC_CCUTOFF ",FMCSC_CCUTOFF),x = kchar)
  kchar<-sub(pattern = "FMCSC_CPROGINDMODE[ ]{1,}[0-9]+",replacement = paste0("FMCSC_CPROGINDMODE ",FMCSC_CPROGINDMODE),x = kchar)
  kchar<-sub(pattern = "FMCSC_CPROGINDRMAX[ ]{1,}[0-9]+",replacement = paste0("FMCSC_CPROGINDRMAX ",FMCSC_CPROGINDRMAX),x = kchar)
  kchar<-sub(pattern = "FMCSC_CPROGINDSTART[ ]{1,}([0-9]+|-[0-9]+)",replacement = paste0("FMCSC_CPROGINDSTART ",FMCSC_CPROGINDSTART),x = kchar)
  kchar<-sub(pattern = "FMCSC_CPROGINDWIDTH[ ]{1,}[0-9]+",replacement = paste0("FMCSC_CPROGINDWIDTH ",FMCSC_CPROGINDWIDTH),x = kchar)
  kchar<-sub(pattern = "FMCSC_CPROGMSTFOLD[ ]{1,}[0-9]+",replacement = paste0("FMCSC_CPROGMSTFOLD ",FMCSC_CPROGMSTFOLD),x = kchar)
  kchar<-sub(pattern = "FMCSC_PCAMODE[ ]{1,}[0-9]+",replacement = paste0("FMCSC_PCAMODE ",FMCSC_PCAMODE),x = kchar)
  kchar<-sub(pattern = "FMCSC_CREDUCEDIM[ ]{1,}[0-9]+",replacement = paste0("FMCSC_CREDUCEDIM ",FMCSC_CREDUCEDIM),x = kchar)
  kchar<-sub(pattern = "FMCSC_CALIGN[ ]{1,}[0-9]+",replacement = paste0("FMCSC_CALIGN ",FMCSC_CALIGN),x = kchar)
  fileConn<-file(kfile)
  writeLines(kchar, fileConn)
  close(fileConn)
  
  system(paste0(camp_home, "bin/x86_64/campari -k ", kfile, ">& ", klog))
  system(paste0("tail ", klog))
}
