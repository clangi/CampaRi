!--------------------------------------------------------------------------!
! LICENSE INFO:                                                            !
!--------------------------------------------------------------------------!
!    This file is part of CAMPARI.                                         !
!                                                                          !
!    Version 3.0                                                           !
!                                                                          !
!    Copyright (C) 2017, The CAMPARI development team (current and former  !
!                        contributors)                                     !
!                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang !
!                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     !
!                        Nicholas Lyle, Nicolas Bloechliger, Marco Bacci,  !
!                        Davide Garolini, Jiri Vymetal                     !
!                                                                          !
!    Website: http://sourceforge.net/projects/campari/                     !
!                                                                          !
!    CAMPARI is free software: you can redistribute it and/or modify       !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!    CAMPARI is distributed in the hope that it will be useful,            !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
!    GNU General Public License for more details.                          !
!                                                                          !
!    You should have received a copy of the GNU General Public License     !
!    along with CAMPARI.  If not, see <http://www.gnu.org/licenses/>.      !
!--------------------------------------------------------------------------!
! AUTHORSHIP INFO:                                                         !
!--------------------------------------------------------------------------!
!                                                                          !
! MAIN AUTHOR:   Andreas Vitalis                                           !
! CONTRIBUTIONS: Nicolas Bloechliger                                       !
!                                                                          !
!--------------------------------------------------------------------------!
!
#include "macros.i"
!
! cludata            : data extracted for clustering (all read into memory -> large)
! cstored            : number of snapshots stored for clustering (-> cludata)
! cmode              : choice of algorithm for clustering or similar tasks (1-5)
! cdis_crit          : the distance criterion underlying the clustering (1-8)
! cstorecalc         : frequency for collecting data for clustering
! cmaxsnaps          : allocation size for cludata in numbers of snapshots dimension
! calcsz,clstsz      : dimensionality size variables for stored data (-> cludata)
! cfilen             : input file for defining clustering subset
! cdofset/cl_imvec   : auxiliary variables for clustering data
! cradius            : general size threshold for clustering
! cprepmode          : options for signal (data) preprocessing
! cdistransform      : if cdis_crit 7-9, what functional transform to apply to raw features (before anything else)
! cdistrans_params   : parameters for this transform
! cchangeweights     : for eligible cdis_crit, replace weights for distance evaluations
! cl_mwvec           : keep track of mean weight per dimension
! cwcombination      : for locally adaptive weights, how to combine (L_? norm)
! cwwindowsz         : for locally adaptive weights, how big a window to use
! cdynwbuf           : for LAWs dependent on transition counts, buffer parameter for taking inverse
! cwacftau           : for ACF-based weights, fixed lag time
! csmoothie          : order for smoothing fxn (sets window size)
! cfeature_dump      : whether to dump features to a NetCDF file
! refine_clustering  : general flag to trigger refinement in clustering
! birchtree,scluster : arrays to hold clustering result
! sumssz             : global size variable for (t_scluster)%sums(:,:)
! cleadermode        : if cmode = 1/2/5(4): processing direction flags
! cmergecrit         : if cmode = 5(4): required cluster-merge criterion for threaded execution only
! clinkage           : if cmode = 3: linkage criterion (1-3)
! chardcut           : if cmode = 3(4): truncation cut for snapshot NB-list
! nblfilen,read_nbl_from_nc,cnblst,cnc_ids : if cmode = 3(4): snapshot NB-list variables
! clufilen           : pre-existing CG trajectory input file for modes 6,7
! clufcol            : column to select from clufilen
! resort_clustering  : instructs CAMPARI to resort ties in cluster size by centroid (not for file clustering)
! cequil             : mode for reweighting of statistics based on network
! caddlinks          : quantile for adding links based on distance (0 disables all augmentation)
! cequilbuf          : buffer parameter
! caddlkmd           : whether and how to augment network (add missing links)
! caddlinkwt         : what type of FP weight to use for added links
! csccwts            : original sampling weights per strongly connected component 
! cprogindex         : if cmode = 4: exact or approximate scheme (1-2)
! cprogindrmax       : if cmode = 4: number of nearest neighbor guesses for approximate scheme
! cprogindstart      : if cmode = 4: snapshot to start from
! cprogpwidth        : if cmode = 4: the (maximum) size of A/B partitions in localized cut
! cprogfold          : if cmode = 4: how often to collapse terminal vertices of the MST
! cprogrdepth        : if cmode = 4 and cprogindex = 2: auxiliary search depth
! cprogbatchsz       : if cmode = 4 and cprogindex = 2: batch size for random stretches 
! csivmax,csivmin    : if cmode = 4: used to identify starting points automatically
! c_nhier            : if cmode = 5(4): tree height
! cmaxrad            : if cmode = 5(4): root level size threshold
! c_multires         : if cmode = 5(4): how many levels to recluster in multi-pass and print
! ccfepmode          : if any clustering is performed, type of cFEP to produce
! align_for_clustering: if cdis_crit = 5/6 whether to do alignment
! cdofsbnds          : if cdis_crit = 6, selectors for separating sets for alignment and distance compu.
! pcamode            : whether to do PCA and what output to produce (1-3)
! reduced_dim_clustering: whether to proceed with clustering in reduced dimensional space after PCA
! align              : structure used for trajectory alignment (ALIGNFILE)
! csnap2tree         : temporary array for cmode = 4 and cprogindex = 2
! csnap2clus         : temporary array for cmode = 4 and cprogindex = 2 
! tmptree            : temporary array for cmode = 4 and cprogindex = 2
! sconnect           : snapshot connectivity map
! sconnectsz         : its size in first dimension
! synmode            : whether to do synth. trajs:  1 (hit clus), 2 (nstraj trajs. of nssnap), 3 as 2 and random ini clus
! nstraj             : number of synthetic trajectories to be generated
! refscc             : reference component for network analysis
! inissnap_usrslct   : snapshot from which all synthetic trajectories start (also to get minimal U set for pfold)
! endssnap_usrslct   : target snapshot for generation of synthetic trajectories (also to get minimal F set for pfold)
! nssnap             : max number of snapshots in a synthetic trajectory before the task of hitting endssnap_usrslct is abandoned
! nsskip             : output frequency for synthetic trajectory files
! ntbrks,ntbrks2,itbrklst : management of user-requested trajectory breaks (counts and counter)
! ntlnks             : management of user-requested link additions (count)
! remap_trajlnks     : whether to try to remap user-requested links if reference snapshots are missing
! cbirchbatch        : control parameter for parallel clustering
! hsl_verbosity      : verbosity setting for HSL routines
! eigvalmd           : how to rank the eigenvalues (0 is disabled, 1 is largest absolute val, 2 is rightmost, 3 imaginay part
! numeig             : number of eigevals (and potentially eigenvects) to be computed from spectral anal. of trans. mat.
! actnumeig          : real number of converged eigenvalues from spectral decomposition
! arnblks            : number of blocks in the arnoldi's method. The default choice is to use the arnblks = numeig + 2
! arnstps            : number of steps in the arnoldi's method. The default choice is arnstps = ceiling((numeig*8)/(numeig + 2))
! arnmaxitr          : max number of restarts in arnoldi's method
! arntol             : tolerance for the arnoldi's method
! tmat_report        : whether or not to write transition matrix(ces) to file(s).
! which_timeflow     : whether to study the forward time transition matrix (1), the backward time (2) or both (3)
! dopfold            : whether or not to compute pfold distribution on clusters
! dopfoldminus       : whether or not to compute also backward (minus (-)) committor probabilities
! pfold_report       : whether or not to write the matrix and right-hand side used to solve pfold and the list of intermed. clus
! clufold            : clusters that belong to the folded set
! cluunfold          : clusters that belong to the unfolded set
! clagt_msm          : lag time for sliding window in MSM. Currently campari does not crash only if a trace file is provided 
! brklnk_report      : report for breaks, links and PIGS-read-trace to standard output
! maxtime_eqmsm      : seconds to limit time execution of equilibrate_MSM
! maxtime_itmfpt     : seconds to limit time execution of iterative_mfpt
! pfolds             : store fwr +, bwr +, fwr -, bwr - committors in its 4 columns and nstruccls (nbasins) rows
!
module clusters
!
  type t_scluster
    integer nmbrs,center,geni,genidx,alsz,nb,ldis,active,rflw,nbalsz,chalsz,nchildren,parent,inscc,ix(2) ! misc.
    RTYPE radius,diam,nodewt(5),quality,sqsum ! properties
    RTYPE, ALLOCATABLE:: sums(:,:),lensnb(:,:),fewtsnb(:,:) ! centroid data and mean geometric length associated with edges
    integer, ALLOCATABLE:: snaps(:),tmpsnaps(:) ! snapshot lists
    integer, ALLOCATABLE:: map(:),children(:),wghtsnb(:,:),lstnb(:),flwnb(:,:) ! for graph
  end type t_scluster
!
  type t_ctree
    type(t_scluster), ALLOCATABLE:: cls(:) ! nothing but an array of clusters
    integer ncls,nclsalsz ! real and allocation sizes
    integer history       ! a helper for keeping track of paths through the tree
  end type t_ctree
!
  type t_cnblst
    integer nbs,alsz ! actual number of neighbors and current allocation size
    RTYPE, ALLOCATABLE:: dis(:) ! list of distances
    integer, ALLOCATABLE:: idx(:) ! and associated indices
    logical, ALLOCATABLE:: tagged(:) ! helper flag
  end type t_cnblst
!
  type t_progindextree
    integer nsnaps ! number of snapshots in tree
    integer nsiblings ! number of trees to merge with
    integer nsibalsz ! alloc size for that
    integer, ALLOCATABLE:: snaps(:) ! indices of snapshots in tree
    integer, ALLOCATABLE:: siblings(:) ! tree indices of tree to merge with
    integer mine(2) ! shortest edge leaving the tree
    RTYPE mind ! distance to nearest snapshot of tree, equals length(mine)
    integer ptr ! simple pointer
  end type t_progindextree
!
  type t_progindexcomponent
    integer ntrees ! number of trees in connected component
    integer, ALLOCATABLE:: trees(:) ! indices of trees in connected component
  end type t_progindexcomponent
!
  type t_adjlist
    integer deg ! degree of vertex
    integer alsz ! allocation size
    integer, ALLOCATABLE:: adj(:) ! list of adjacent vertices
    real(KIND=4), ALLOCATABLE:: dist(:) ! distance to the adjacent vertices
  end type t_adjlist
!
  type(t_cnblst), ALLOCATABLE:: cnblst(:)
  type(t_scluster), ALLOCATABLE:: scluster(:)
  type(t_ctree), ALLOCATABLE:: birchtree(:),thr_btree(:,:)
  type(t_adjlist), ALLOCATABLE:: approxmst(:)
  type(t_progindextree), ALLOCATABLE:: tmptree(:)
  RTYPE, ALLOCATABLE:: cludata(:,:),cl_imvec(:),cl_mwvec(:),csccwts(:)
  integer, ALLOCATABLE:: trbrkslst(:) ! a user-populated list of transitions to delete in graph-related analyses
  integer, ALLOCATABLE:: trlnkslst(:,:) ! a user-populated list of transitions to add in graph-related analyses
  integer cstored,cdis_crit,cstorecalc,calcsz,clstsz,sumssz,cmode,nstruccls,pcamode,cprogindex,cprogindrmax,cdofsbnds(4),cmaxsnaps
  integer csivmin,csivmax,cnc_ids(10),c_nhier,clinkage,cleadermode,reduced_dim_clustering,cprogindstart,cprogpwidth,ccfepmode
  integer ntbrks,ntbrks2,itbrklst,cequil,cchangeweights,cwacftau,cwcombination,cwwindowsz,csmoothie,cprepmode,cprogfold,caddlkmd
  integer(KIND=8) cdevalcnt,maxtime_eqmsm
  integer cprogrdepth,c_multires,cprogbatchsz,ntlnks,cbirchbatch,sconnectsz,clufcol,clagt_msm,eigvalmd,numeig,actnumeig,arnblks,&
 &brklnk_report,cdistransform
  integer arnmaxitr,arnstps,which_timeflow,refscc,synmode,nstraj,inissnap_usrslct,endssnap_usrslct,nssnap,nsskip,hsl_verbosity
  integer, ALLOCATABLE:: cdofset(:,:),sconnect(:,:)
  integer, ALLOCATABLE:: csnap2tree(:),csnap2clus(:,:),clufold(:),cluunfold(:)
  RTYPE cradius,cmaxrad,chardcut,cdynwbuf,cequilbuf,caddlinks,caddlinkwt,cmergecrit,arntol,cdistrans_params(2)
  character(MAXSTRLEN) cfilen,nblfilen,tbrkfilen,tlnkfilen,clufilen
  logical read_nbl_from_nc,align_for_clustering,refine_clustering,resort_clustering,remap_trajlnks
  logical tmat_report,doeigvec_msm,dopfold,dopfoldminus,pfold_report,cfeature_dump
  character(MAXSTRLEN) clufoldfile,cluunfoldfile
  RTYPE, ALLOCATABLE:: erp_msm_fwr(:),erp_msm_bwr(:),eip_msm_fwr(:),eip_msm_bwr(:) !real+imaginary part of tmat std state
  RTYPE, ALLOCATABLE:: eigvec_msm_fwr(:,:),eigvec_msm_bwr(:,:) !eigenvectors of transition matrix
  RTYPE, ALLOCATABLE:: pfolds(:,:) 
!
! other
  type t_align
    integer, ALLOCATABLE:: set(:),diffset(:) ! atom set
    RTYPE, ALLOCATABLE:: refxyz(:),curxyz(:),hlpxyz(:) ! coordinate helpers
    integer nr,diffnr,calc,mmol ! sizes, frequency, ref. mol.
    logical yes,refset,haveset,instrmsd
    character(MAXSTRLEN) filen ! for providing custom alignment set
  end type t_align
!
  type(t_align) align
!
end module clusters
!
