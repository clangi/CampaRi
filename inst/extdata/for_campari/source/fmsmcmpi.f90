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
!                                                                          !
!--------------------------------------------------------------------------!
#ifdef ENABLE_MPI
!
#include "macros.i"
!
!
subroutine MPI_StartMC()
!
  use iounit
  use keys
  use mpi
  use mpistuff
  use interfaces
  use, INTRINSIC:: ISO_C_BINDING
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer masterrank,ktag,initag,mstatus(MPI_STATUS_SIZE)
  integer i,j,dt(8),ierr,dummy(1),nktag,aone,kstag,atwo,athree,t1,t2
  integer, ALLOCATABLE:: keysz(:)
  character (len=12) rc(3)
  character(len=1,kind=c_char), dimension(MAXSTRLEN+1):: honam(0:MAXSTRLEN)
  character(MAXSTRLEN) hnam
  integer(c_int) hnerr
  integer(c_size_t) hnml
#ifdef ENABLE_THREADS
  integer tdtag,thrbu
  logical mainy
#endif
!
! some initialization
  dummy(1) = 0
  aone = 1
  atwo = 2
  athree = 3
  masterrank = 0
  kstag = 444
  nktag = 445
  initag = 555
  ktag = 556
!
! let's start up the universe
#ifdef ENABLE_THREADS
  tdtag = 443
  call find_nthreads() ! gets the desired number of threads into thrdat%maxn
  thrbu = thrdat%maxn
  call OMP_SET_NUM_THREADS(thrdat%maxn)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,ierr)
!$OMP MASTER
  call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,j,ierr)
  if (j.ne.MPI_THREAD_FUNNELED) then
    write(*,*) 'Fatal. MPI library does not allow required multithreaded execution model.'
    call fexit()
  end if
!$OMP END MASTER
!$OMP BARRIER
! test thread identiy 
!$OMP MASTER
  call MPI_IS_THREAD_MAIN(mainy,ierr)
  if (mainy.EQV..false.) then
    write(*,*) 'Fatal. When calling MPI functions from an OpenMP MASTER construct, the calling thread &
 &is not correctly identified or mapped.'
    call fexit()
  end if
!$OMP END MASTER
!$OMP BARRIER
!$OMP END PARALLEL
  call MPI_IS_THREAD_MAIN(mainy,ierr)
  if (mainy.EQV..false.) then
    write(*,*) 'Fatal. When calling MPI functions from outside of an OpenMP PARALLEL region, the calling thread &
 &is not correctly identified or mapped.'
    call fexit()
  end if
#else
  call MPI_INIT(ierr)
#endif
!
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_nodes,ierr)
  hnml = MAXSTRLEN
  honam(:) = ' '
  hnerr = handwrapped_gethostname(honam(0:(MAXSTRLEN-1)),hnml)
  do i=1,MAXSTRLEN
    hnam(i:i) = honam(i-1)
  end do
  call strlims(hnam,t1,t2)
  if (t2.gt.t1) then
    if (hnam(t2:t2).eq.C_NULL_CHAR) t2 = t2 - 1
  end if
  re_aux(10) = max(3,int(log10(1.0*mpi_nodes) + 1)) 
! 
!
! the master will distribute the keyword information and wait for each node
! successively to complete any parsing of the key-file
! there are still file-dependent setup procedures, but those don't seem to
! create conflicts
  if (myrank.eq.masterrank) then
!    write(*,*) 'M:Entered master.'
    call initial()
#ifdef ENABLE_THREADS
    thrdat%maxn = thrbu ! otherwise overwritten by initialization
#endif
    call makelogio(aone)
    call Date_and_Time(rc(1),rc(2),rc(3),dt)
    write(ilog,*)
    write(ilog,*) '----------- LOGFILE FOR MASTERNODE ------------'
    call write_license()
    write(ilog,*) 'Execution on host ',hnam(t1:t2),' started at ',rc(2)(1:2),':',rc(2)(3:4),' on ',&
 &                 rc(1)(5:6),'/',rc(1)(7:8),'/',rc(1)(1:4),'.'
    write(ilog,*)
    call getkey()
    call init_campprng(aone)
!   we parse the first two levels of the keyfile to determine scope
!   and sanity of the calculation
    write(ilog,*)
    write(ilog,*) '---   Now parsing keywords ...            ---'
    write(ilog,*)
    call parsekey(aone)
#ifdef ENABLE_THREADS
    thrdat%maxn = thrbu ! otherwise still overwritten if user request is unusable but system request was defined
#endif
    call parsekey(atwo)
!   a sanity check whether we have a sane MPI run at all
    if (((use_MPIAVG.EQV..false.).AND.(use_REMC.EQV..false.)).OR.&
 & ((use_REMC.EQV..true.).AND.(use_MPIAVG.EQV..true.))) then
      write(ilog,*) 'MPI runs must have either FMCSC_MPIAVG or FMCSC&
 &_REMC but not both set to true (1). Fatal exit (specify both keywords with the correct choices).'
      call fexit()
    end if
!   Sanity check against mismatch of nodes and replicas in MPI-RE
    if (use_REMC.EQV..true.) then
      if (mpi_nodes.ne.re_conditions) then
        write(ilog,*) 'Fatal: Number of nodes granted by mpirun (',mpi_nodes,') does not match number of replicas requested (',&
 &re_conditions,'). This is fatal in calculations using the replica exchange setup.'
        call fexit()
      end if
    else if (use_MPIAVG.EQV..true.) then
      if (mpi_nodes.ne.re_conditions) then
        write(ilog,*) 'Warning. Number of nodes granted by mpirun (',mpi_nodes,') does not match default or chosen value for &
 &FMCSC_REPLICAS (',re_conditions,') in the MPI averaging framework. Adjusting setting for FMCSC_REPLICAS.'
        re_conditions = mpi_nodes
      end if
    end if
    call parsekey(athree)
    write(ilog,*)
    write(ilog,*) '---   ... finished parsing keywords.      ---'
    write(ilog,*)
!   now we can allocate memory ...
    call allocate_mpistuff(aone)
!   and evtl. read the RE input file
    if (use_REMC.EQV..true.) then
      call read_refile()
    end if
    allocate(keysz(nkey))
    do j=1,nkey
      keysz(j) = size(key(j)%line)
    end do
!    write(*,*) 'M:Processed input.'
    do i=1,mpi_nodes-1
!      write(*,*) 'M:Sending to ',i
#ifdef ENABLE_THREADS
      call MPI_SEND(thrbu,1,MPI_INTEGER,i,tdtag,MPI_COMM_WORLD,ierr)
#endif
      call MPI_SEND(nkey,1,MPI_INTEGER,i,nktag,MPI_COMM_WORLD,ierr)
      call MPI_SEND(keysz,nkey,MPI_INTEGER,i,kstag,MPI_COMM_WORLD,ierr)
      do j=1,nkey
!        write(*,*) 'sending ',key(j)%line
        call MPI_SEND(key(j)%line,size(key(j)%line),MPI_CHARACTER,i,ktag+j,MPI_COMM_WORLD,ierr)
      end do
!      write(*,*) 'M:Waiting for completion from ',i
      call MPI_RECV(dummy,1,MPI_INTEGER,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.initag) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from slave ',i,'. Expected ',initag,'.'
        call fexit()
      end if
!      write(*,*) 'M:Done with ',i+1
    end do
    deallocate(keysz)
  else
!    write(*,*) 'S:Entered slave ',myrank,': Running initial.'
    call initial()
!    write(*,*) 'S:Slave ',myrank,': Waiting for master.'
#ifdef ENABLE_THREADS
    call MPI_RECV(thrbu,1,MPI_INTEGER,masterrank,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.tdtag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from master. Expected ',tdtag,'.'
      call fexit()
    end if
    thrdat%maxn = thrbu
#endif
    call MPI_RECV(dummy,1,MPI_INTEGER,masterrank,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.nktag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from master. Expected ',nktag,'.'
      call fexit()
    end if
    nkey = dummy(1)
    allocate(keysz(nkey))
    allocate(key(nkey))
    call MPI_RECV(keysz,nkey,MPI_INTEGER,masterrank,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.kstag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from master. Expected ',kstag,'.'
      call fexit()
    end if
    do j=1,nkey
      allocate(key(j)%line(keysz(j)))
      call MPI_RECV(key(j)%line,keysz(j),MPI_CHARACTER,masterrank,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.(ktag+j)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from master. Expected ',ktag+j,'.'
        call fexit()
      end if
    end do
    deallocate(keysz)
!    do j=1,nkey
!      write(*,*) key(j)%line
!    end do
!   log-file handles, initialization and input file processing
!    write(*,*) 'S:Slave ',myrank,' is processing input.'
    call makelogio(aone)
    call Date_and_Time(rc(1),rc(2),rc(3),dt)
    write(ilog,*)
    write(ilog,*) '--------- LOGFILE FOR NODE ',myrank+1,' ----------'
    write(ilog,*) 'Execution started on host ',hnam(t1:t2),' at ',rc(2)(1:2),':',rc(2)(3:4),' on ',&
 &rc(1)(5:6),'/',rc(1)(7:8),'/',rc(1)(1:4),'.'
    write(ilog,*)
    call init_campprng(aone)
    write(ilog,*)
    write(ilog,*) '---   Now parsing keywords ...            ---'
    write(ilog,*)
    call parsekey(aone)
#ifdef ENABLE_THREADS
    thrdat%maxn = thrbu ! otherwise still overwritten if user request is unusable but system request was defined
#endif
    call parsekey(atwo)
    call parsekey(athree)
    write(ilog,*)
    write(ilog,*) '---   ... finished parsing keywords.      ---'
    write(ilog,*)
    call allocate_mpistuff(aone)
    if (use_REMC.EQV..true.) then
      call read_refile()
    end if
    call MPI_SEND(dummy,1,MPI_INTEGER,masterrank,initag,&
 &MPI_COMM_WORLD,ierr)
!    write(*,*) 'S:Slave ',myrank,' done with setup.'
!   from here, the execution runs through the standard chainsaw code, MPI
!   will become relevant again when the first swap move is encountered
  end if
!
! now with mpi_granularity set, we can populate the communication flow structures
  call MPI_assigncommflow()
!   
end
!
!-----------------------------------------------------------------------------
!
subroutine MPI_StopMC()
!
  use iounit
  use mcsums
  use mpi
  use mpistuff
!
  implicit none
!
  integer i,j,masterrank,ierr,endtag
  integer mstatus(MPI_STATUS_SIZE)
!
  masterrank = 0
  endtag = 777
!
  if (myrank.eq.masterrank) then
!
    do i=1,mpi_nodes-1
!      write(*,*) 'sending end to ',i+1
      call MPI_Send(j,1,MPI_INTEGER,i,endtag,MPI_COMM_WORLD,ierr)
    end do
!
  else
!
    call MPI_RECV(j,1,MPI_INTEGER,masterrank,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.endtag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &' from master. Expected end signal (',endtag,').'
      call fexit()
    end if
!
  end if
!
end
!
!-----------------------------------------------------------------------
!
subroutine MPI_Synchronize()
!
  use iounit
  use mcsums
  use mpi
  use mpistuff
  use system
!
  implicit none
!
  integer i,j,masterrank,ierr,stptag,stps(mpi_nodes)
  integer mstatus(MPI_STATUS_SIZE),lowstp,aone,stptag2
!
  masterrank = 0
  aone = 1
  stptag = 888
  stptag2 = 889
!
  if (myrank.eq.masterrank) then
!
    stps(1) = nstep
    do i=1,mpi_nodes-1
      call MPI_RECV(j,aone,MPI_INTEGER,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.stptag) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &' from slave ',i,'. Expected step-signal (',stptag,').'
        call fexit()
      end if
      stps(i+1) = j
    end do
    lowstp = nsim+1
    do i=1,mpi_nodes
      if (stps(i).lt.lowstp) then
        lowstp = stps(i)
      end if
    end do
    nstep = lowstp
    do i=1,mpi_nodes-1
      call MPI_Send(lowstp,aone,MPI_INTEGER,i,stptag2,&
 &                                  MPI_COMM_WORLD,ierr)
    end do
!
  else
!
    call MPI_Send(nstep,aone,MPI_INTEGER,masterrank,stptag,&
 &                                 MPI_COMM_WORLD,ierr)
    call MPI_RECV(j,aone,MPI_INTEGER,masterrank,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.stptag2) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &' from master. Expected step2-signal (',stptag2,').'
      call fexit()
    end if
    nstep = j
!
  end if
!
  end
!
!-----------------------------------------------------------------------
!
subroutine MPI_SyncHybridCycle(istep)
!
  use mpi
  use mpistuff
  use movesets
!
  implicit none
!
  integer, INTENT(IN):: istep
!
  integer i,cycdims(4),azero,ierr,masterrank,cyctag(4)
  logical atrue
!
  atrue = .true.
  azero = 0
  masterrank = 0
  cyctag(1) = 887
  cyctag(2) = 886
  cyctag(3) = 885
  cyctag(4) = 884 
!
  if (myrank.eq.masterrank) then
    cycdims(1) = curcyc_start
    cycdims(2) = curcyc_end
  end if
!
  i = 4
  call MPI_ALLINTS1D(cyctag(1:4),azero,ierr,atrue,i,cycdims(1:4),use_MPIcolls)
!
  if (myrank.ne.masterrank) then
    curcyc_start = cycdims(3)
    curcyc_end = cycdims(4)
  end if
!
  end
!
!-----------------------------------------------------------------------
!
subroutine MPI_SyncInt(theint,themode)
!
  use mpi
  use mpistuff
  use system, ONLY: nsim
!
  implicit none
!
  integer, INTENT(INOUT):: theint
  integer, INTENT(IN):: themode
!
  integer i,nidims(2),ierr,masterrank,nstag(4)
  logical atrue
!
  atrue = .true.
  masterrank = 0
  nstag(1) = 893
  nstag(2) = 892
  nstag(3) = 891
  nstag(4) = 890 
!
  nidims(1) = theint
!
  i = 2
  call MPI_ALLINTS1D(nstag(1:4),themode,ierr,atrue,i,nidims(1:2),use_MPIcolls)
!
  theint = nidims(2)
!
  end

!
!---------------------------------------------------------------------------
! 
! this one collects the energy vs. parameter of hamiltonian information
! for every replica (i.e., a vector (single structure vs. all conditions)
! per replica)
! the size of the system (arrays) must be conserved
!
subroutine MPI_REMaster(istep,ndump,tpi)
!
  use iounit
  use atoms
  use accept
  use molecule
  use mpi
  use mpistuff
  use mcsums
  use energies
  use system
  use movesets
  use torsn
  use sequen
  use cutoffs
  use ems
  use zmatrix
!
  implicit none
!
  integer, INTENT(IN):: istep,ndump,tpi
!
  integer azero,atwo,athree
  integer i,j,k,ierr,reptag,masterrank,msgtag(4)
  integer mstatus(MPI_STATUS_SIZE)
  RTYPE dumes(3),boxovol,force3,frac1,frac2
  logical atrue,reolallbu,bcast,dynni,afalse
  integer(KIND=8) ee1,ee2
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
!
  if (tpi.le.0) then
    if (OMP_IN_PARALLEL().EQV..true.) then
      write(ilog,*) 'Fatal. When using hybrid MPI/multi-threaded code, MPI_REMaster(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
      call fexit()
    end if
  end if
#endif
!
  atrue = .true.
  afalse = .false.
  dynni = .true.
  if ((dyn_mode.eq.1).OR.((in_dyncyc.EQV..false.).AND.((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)))) dynni = .false.
  reptag = 10
  masterrank = 0
  azero = 0
  atwo = 2
  athree = 3
  boxovol = ens%insV
  dumes(:) = 0.0
!
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
!
  if (ens%flag.ne.1) then
    write(ilog,*) 'Fatal. Replica-exchange calculations are not supported for anything&
 & but the canonical (NVT) ensemble. This is most certainly a bug.'
    call fexit()
  end if
!
! increase the step counter
  nstep = nstep + 1
  mvcnt%nre = mvcnt%nre + 1
  reolallbu = reol_all
  if ((re_nbmode.eq.1).AND.(reol_all.EQV..false.)) then
    reol_all = .true.
  end if
  call System_Clock(ee1)
  mpi_evec(:) = 0.0
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
! every replica must recompute foreign energies as desired
  if (dynni.EQV..false.) then
    call lamenergy(mpi_evec,mpi_fvec,tpi)
  else
    call lamforce(mpi_evec,mpi_fvec,tpi)
  end if
!
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
 333 format(i10,100(g18.10,1x))
  reol_all = reolallbu
  call System_Clock(ee2)
  time_energy = time_energy + 1.0*(ee2 - ee1)
!
! now let's combine energy vectors -> note that this will not scale well with large numbers of replicas
  call System_Clock(ee1)
  if (myrank.eq.masterrank) then
    mpi_emat(myrank+1,:) = mpi_evec(:)
    do i=1,mpi_nodes-1
      if (re_nbmode.eq.2) then
        if (re_conditions.eq.2) then
          call MPI_RECV(mpi_evec(1:2),atwo,MPI_RTYPE,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
        else
          call MPI_RECV(mpi_evec(1:3),athree,MPI_RTYPE,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
        end if
      else
        call MPI_RECV(mpi_evec,re_conditions,MPI_RTYPE,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      end if
      if (mstatus(MPI_TAG).ne.reptag) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from slave ',i,'. Expected ',reptag,'.'
        call fexit()
      end if
      if (re_nbmode.eq.2) then
        if ((i.lt.(re_conditions-1)).AND.(i.gt.0)) then
          mpi_emat(i+1,i:(i+2)) = mpi_evec(1:3)
        else if (i.eq.0) then
          mpi_emat(i+1,(i+1):(i+2)) = mpi_evec(1:2)
        else if (re_conditions.eq.2) then
          mpi_emat(i+1,i:(i+1)) = mpi_evec(1:2)
        else
          mpi_emat(i+1,i:(i+1)) = mpi_evec(2:3)
        end if
      else
        mpi_emat(i+1,:) = mpi_evec(:)
      end if
    end do
  else
    if (re_nbmode.eq.2) then
      if (re_conditions.eq.2) then
        if (myrank.eq.0) then
          dumes(1:2) = mpi_evec((myrank+1):(myrank+2))
        else
          dumes(1:2) = mpi_evec((myrank):(myrank+1))
        end if
        call MPI_SEND(dumes(1:2),atwo,MPI_RTYPE,masterrank,reptag,MPI_COMM_WORLD,ierr)
      else
        if ((myrank.lt.(re_conditions-1)).AND.(myrank.gt.0)) then
          dumes(1:3) = mpi_evec(myrank:(myrank+2))
        else if (myrank.eq.0) then
          dumes(1:2) = mpi_evec((myrank+1):(myrank+2))
        else
          dumes(2:3) = mpi_evec((myrank):(myrank+1))
        end if
        call MPI_SEND(dumes,athree,MPI_RTYPE,masterrank,reptag,MPI_COMM_WORLD,ierr)
      end if
    else
      call MPI_SEND(mpi_evec,re_conditions,MPI_RTYPE,masterrank,reptag,MPI_COMM_WORLD,ierr)
    end if
  end if
  call System_Clock(ee2)
  time_comm = time_comm + 1.0*(ee2 - ee1)
!
! determine new mapping
  if (myrank.eq.masterrank) call RE_swap(mpi_emat,mpi_lmap(1:mpi_nodes,1))
!
! Broadcast mapping
  call System_Clock(ee1)
  msgtag(1) = 3367
  msgtag(2) = msgtag(1) + 1
  msgtag(3) = msgtag(1) + 2
  msgtag(4) = msgtag(1) + 3
  bcast = .true.
!
  k = 2*mpi_nodes
  call MPI_ALLINTS1D(msgtag(1:4),azero,ierr,bcast,k,mpi_lmap(:,1),use_MPIcolls) ! always uses MPI master as BC root
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
  j = mpi_lmap(mpi_nodes+myrank+1,1) - 1
  if (re_nbmode.eq.1) then
    do k=mpi_nodes+1,2*mpi_nodes
      if (mpi_lmap(k,1).eq.(myrank+1)) exit
    end do
    k = k - 1 - mpi_nodes
  else
    k = j
  end if
  if (j.gt.myrank) then
    call MPI_SendCoordinates(j,afalse,dynni,tpi)
    call MPI_ReceiveCoordinates(k,afalse,dynni,tpi)
    if (tpi.le.1) then
      call System_Clock(ee2)
      time_comm = time_comm + 1.0*(ee2 - ee1)
    end if
  else if (j.lt.myrank) then
    call MPI_ReceiveCoordinates(k,atrue,dynni,tpi)
    call MPI_SendCoordinates(j,atrue,dynni,tpi)
    if (tpi.le.1) then
      call System_Clock(ee2)
      time_comm = time_comm + 1.0*(ee2 - ee1)
    end if
  else ! no swap means we may have to restore energy
    if (tpi.le.1) then
      call System_Clock(ee2)
      time_comm = time_comm + 1.0*(ee2 - ee1)
      call System_Clock(ee1)
      if (use_EMICRO.EQV..true.) curmassm = .false.
    end if
    if (dynni.EQV..true.) then
#ifdef ENABLE_THREADS
      if (tpi.gt.0) then
        call force3_threads(esterms,esterms_tr,esterms_lr,atrue,esave)
      else
        esave = force3(esterms,esterms_tr,esterms_lr,atrue)
      end if
#else
      esave = force3(esterms,esterms_tr,esterms_lr,atrue)
#endif
    else
      call energy3(esterms,atrue,esave,tpi)
    end if
    if (tpi.le.1) then
      ens%insU = esave
      call System_Clock(ee2)
      time_energy = time_energy + 1.0*(ee2 - ee1)
    end if
  end if
!
! trace maintenance
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
  if (myrank.eq.masterrank) then
    do i=1,mpi_nodes
      mpi_lmap(i,1) = mpi_lmap(mpi_lmap(i+mpi_nodes,1),2)
    end do
    mpi_lmap(1:mpi_nodes,2) = mpi_lmap(1:mpi_nodes,1)
  end if
!
  if (mpi_lmap(myrank+mpi_nodes+1,1).ne.(myrank+1)) then
    acc%nre = acc%nre + 1
  end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
  if (dynni.EQV..true.) then
    if (fycxyz.eq.1) then
      call get_ensv(boxovol,tpi)
    else
      call get_cart_ensv(boxovol,tpi)
    end if
    call drift_removal(frac1,frac2,tpi)
  end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
  if (dynni.EQV..true.) call prt_curens(istep,boxovol,frac1,frac2)
  if ((myrank.eq.masterrank).AND.(inst_retr.EQV..true.)) call MPI_REPrtTrace(istep,iretr,mpi_lmap(1:mpi_nodes,2))
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
  if (tpi.le.1) call System_Clock(ee1)
  call mcstat(istep,ndump,tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
  if (tpi.le.1) then
    call System_Clock(ee2)
    time_analysis = time_analysis + 1.0*(ee2 - ee1)
  end if
!
end
!
!----------------------------------------------------------------------------------------------
!
! a universal fxn to send complete configuration to replica tnode (0:(N-1) numbering)
! logical frombuf selects sending from buffer [(:,2) vecs.] or from actual configuration arrays
! logical dynni restricts transfer to pure positional coordinates as for MC
!
subroutine MPI_SendCoordinates(tnode,frombuf,dynni,tpi)
!
  use mpi
  use mpistuff
  use mcsums
  use system
  use atoms
  use molecule
  use torsn
  use iounit
  use forces
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi,tnode
  logical, INTENT(IN):: frombuf,dynni
!
  integer aone,atwo,tortag,rbctag,xyztag,cmdtag,cldtag,imdtag,ildtag,ierr,bufix
  integer(KIND=8) ee1,ee2
#ifdef ENABLE_THREADS
  integer sta,sto
!
  if (tpi.gt.0) then
    sta = thr_limits(1,tpi)
    sto = thr_limits(2,tpi)
  else
    sta = 1
    sto = n
  end if
#endif
!
  if (tpi.le.1) then
    call System_Clock(ee1)
  end if
!
  tortag = 11
  rbctag = 12
  xyztag = 13
  aone = 1
  atwo = 2
  cmdtag = 14
  cldtag = 15
  imdtag = 16
  ildtag = 17
  bufix = 1
  if (frombuf.EQV..true.) bufix = 2
!
! send over structure under tags 11,12 or 13
  if ((fycxyz.eq.1).AND.(force_rexyz.EQV..false.)) then
    if (frombuf.EQV..false.) call transfer_torsn(aone,mpi_tvec(:,1),tpi)
    if (frombuf.EQV..false.) call transfer_rbc(aone,mpi_rbvec(:,1),tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    call MPI_SEND(mpi_tvec(:,bufix),ntorpuck,MPI_RTYPE,tnode,tortag,MPI_COMM_WORLD,ierr)
    call MPI_SEND(mpi_rbvec(:,bufix),9*nmol,MPI_RTYPE,tnode,rbctag,MPI_COMM_WORLD,ierr)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
  else if ((fycxyz.eq.2).OR.(force_rexyz.EQV..true.)) then
#ifdef ENABLE_THREADS
    if (frombuf.EQV..false.) mpi_cdsvec(sta:sto,1,1) = x(sta:sto)
    if (frombuf.EQV..false.) mpi_cdsvec(sta:sto,2,1) = y(sta:sto)
    if (frombuf.EQV..false.) mpi_cdsvec(sta:sto,3,1) = z(sta:sto)
#else
    if (frombuf.EQV..false.) mpi_cdsvec(1:n,1,1) = x(1:n)
    if (frombuf.EQV..false.) mpi_cdsvec(1:n,2,1) = y(1:n)
    if (frombuf.EQV..false.) mpi_cdsvec(1:n,3,1) = z(1:n)
#endif
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    call MPI_SEND(mpi_cdsvec(:,:,bufix),3*n,MPI_RTYPE,tnode,xyztag,MPI_COMM_WORLD,ierr)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
  else
    write(ilog,*) 'Fatal. Encountered unsupported choice of degrees of freedom (',fycxyz,') in&
 & MPI_SendCoordinates(...). This is most certainly a bug.'
    call fexit()
  end if
!
  if (dynni.EQV..true.) then
    if ((dyn_mode.eq.2).OR.(dyn_mode.eq.5)) then
      if (fycxyz.eq.1) then
        if (frombuf.EQV..false.) call manage_imds(mpi_dvec(:,:,1),aone,tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
        call MPI_SEND(mpi_dvec(:,:,bufix),(ndyntorsn+totrbd)*3,MPI_RTYPE,tnode,imdtag,MPI_COMM_WORLD,ierr)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
      else if (fycxyz.eq.2) then
        if (frombuf.EQV..false.) call manage_cmds(mpi_dvec(:,:,1),aone,tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
        call MPI_SEND(mpi_dvec(:,:,bufix),3*n,MPI_RTYPE,tnode,cmdtag,MPI_COMM_WORLD,ierr)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
      else
        write(ilog,*) 'Fatal. Encountered unsupported choice of degrees of freedom (',fycxyz,') in&
   & MPI_SendCoordinates(...). This is most certainly a bug.'
        call fexit()
      end if
    else if ((dyn_mode.eq.3).OR.(dyn_mode.eq.7)) then
      if (fycxyz.eq.1) then
        if (frombuf.EQV..false.) call manage_ilds(mpi_dvec(:,:,1),aone,tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
        call MPI_SEND(mpi_dvec(:,:,bufix),(ndyntorsn+totrbd)*4,MPI_RTYPE,tnode,ildtag,MPI_COMM_WORLD,ierr)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
      else if (fycxyz.eq.2) then
        if (frombuf.EQV..false.) call manage_clds(mpi_dvec(:,:,1),aone,tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
        call MPI_SEND(mpi_dvec(:,:,bufix),9*n,MPI_RTYPE,tnode,cldtag,MPI_COMM_WORLD,ierr)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
      else
        write(ilog,*) 'Fatal. Encountered unsupported choice of degrees of freedom (',fycxyz,') in&
   & MPI_SendCoordinates(...). This is most certainly a bug.'
        call fexit()
      end if
    else if (dyn_mode.eq.1) then
!     do nothing
    else
      write(ilog,*) 'Fatal. Encountered unsupported dynamics method (',dyn_mode,') in&
   & MPI_SendCoordinates(...). This is most certainly a bug.'
      call fexit()
    end if
  end if
!
  if (tpi.le.1) then
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
  end if
!
end
!
!----------------------------------------------------------------------------------------------
!
! a universal fxn to receive complete configuration from replica tnode (0:(N-1) numbering)
! recomputes energy/forces as necessary, but does not call any analysis 
! logical buf selects whether to back up configuration first (into (:,2) vecs.)
! logical dynni indicates whether or not to assume positional transfer only as in MC (use energy3 instead
! force3 if dynni is false)
!
subroutine MPI_ReceiveCoordinates(tnode,buf,dynni,tpi)
!
  use mpi
  use mpistuff
  use mcsums
  use system
  use atoms
  use molecule
  use torsn
  use iounit
  use cutoffs
  use sequen
  use energies
  use ems
  use forces
  use zmatrix
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi,tnode
  logical, INTENT(IN):: buf,dynni
!
  integer azero,aone,atwo,tortag,rbctag,xyztag,cmdtag,cldtag,imdtag,ildtag,ierr,i,dummy
  integer mstatus(MPI_STATUS_SIZE)
  RTYPE force3,otherT
  logical atrue
  integer(KIND=8) ee1,ee2
#ifdef ENABLE_THREADS
  integer sta,sto
  logical moljfl1(4),moljfl2(4)
!
  if (tpi.gt.0) then
    sta = thr_limits(1,tpi)
    sto = thr_limits(2,tpi)
  else
    sta = 1
    sto = n
  end if
#endif
!
  if (tpi.le.1) call System_Clock(ee1)
!
  tortag = 11
  rbctag = 12
  xyztag = 13
  azero = 0 
  aone = 1
  atwo = 2
  cmdtag = 14
  cldtag = 15
  imdtag = 16
  ildtag = 17
  atrue = .true.
!
  if ((fycxyz.eq.1).AND.(force_rexyz.EQV..false.)) then
    if (buf.EQV..true.) then
      call transfer_torsn(aone,mpi_tvec(:,2),tpi)
      call transfer_rbc(aone,mpi_rbvec(:,2),tpi)
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    call MPI_RECV(mpi_tvec(:,1),ntorpuck,MPI_RTYPE,tnode,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.tortag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from node ',tnode+1,'. Expected ',tortag,'.'
      call fexit()
    end if
    call MPI_RECV(mpi_rbvec(:,1),9*nmol,MPI_RTYPE,tnode,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.rbctag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from node ',tnode+1,'. Expected ',rbctag,'.'
      call fexit()
    end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!     let's update the configuration
!      write(*,*) 'Slave ',myrank,': Finalizing myself'
    call transfer_torsn(atwo,mpi_tvec(:,1),tpi)
    call transfer_rbc(atwo,mpi_rbvec(:,1),tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
    if (tpi.gt.0) then
      call fyczmat_threads(tpi)
!$OMP BARRIER
      do i=thr_limits(35,tpi),thr_limits(36,tpi)
        if (thr_mlgix(i).gt.0) cycle
!       we need to recompute internals
        call makexyz_formol(i)
        call update_rigid(i)
        call update_rigidm(i)
      end do
      moljfl1(:) = .false.
      moljfl1(1) = .true. ! only update_rigidm
      moljfl2(:) = moljfl1(:)
      moljfl2(3) = .true. ! only update_rigid
      do i=1,nmlgs
        call makexyz_threads(i,tpi) ! ends with a barrier
        call molops_threads(i,tpi,moljfl1)
        call molops_threads_geo(i,tpi,moljfl2)
      end do
! update pointer arrays such that analysis routines can work properly
!$OMP BARRIER
    else
#endif
    call fyczmat()
    do i=1,nmol
      call makexyz_formol(i)
      call update_rigid(i)
      call update_rigidm(i)
    end do
#ifdef ENABLE_THREADS
    end if
#endif
  else if ((fycxyz.eq.2).OR.(force_rexyz.EQV..true.)) then
    if (buf.EQV..true.) then
#ifdef ENABLE_THREADS
      mpi_cdsvec(sta:sto,1,2) = x(sta:sto)
      mpi_cdsvec(sta:sto,2,2) = y(sta:sto)
      mpi_cdsvec(sta:sto,3,2) = z(sta:sto)
#else
      mpi_cdsvec(1:n,1,2) = x(1:n)
      mpi_cdsvec(1:n,2,2) = y(1:n)
      mpi_cdsvec(1:n,3,2) = z(1:n)
#endif
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    call MPI_RECV(mpi_cdsvec(:,:,1),3*n,MPI_RTYPE,tnode,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.xyztag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from node ',tnode+1,'. Expected ',xyztag,'.'
      call fexit()
    end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
    x(sta:sto) = mpi_cdsvec(sta:sto,1,1)
    y(sta:sto) = mpi_cdsvec(sta:sto,2,1)
    z(sta:sto) = mpi_cdsvec(sta:sto,3,1)
!$OMP BARRIER
#else
    x(1:n) = mpi_cdsvec(1:n,1,1)
    y(1:n) = mpi_cdsvec(1:n,2,1)
    z(1:n) = mpi_cdsvec(1:n,3,1)
#endif
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      do i=thr_limits(35,tpi),thr_limits(36,tpi)
        if (thr_mlgix(i).gt.0) cycle
!       we need to recompute internals
        call genzmat(i)
        call update_rigid(i)
        call update_rigidm(i)
      end do
      moljfl1(:) = .true.
      moljfl1(3) = .false. ! all except c.o.m. velocity
      moljfl2(:) = .true.
      moljfl2(2) = .false.
      moljfl2(4) = .false. ! only update_rigid
      do i=1,nmlgs
!       covers update_rigidm, update_comv, genzmat, update_image
        call molops_threads(i,tpi,moljfl1)
        call molops_threads_geo(i,tpi,moljfl2)
      end do
! update pointer arrays such that analysis routines can work properly
!$OMP BARRIER
      call zmatfyc_threads(tpi,atwo)
    else
#endif
    do i=1,nmol
      call genzmat(i)
      call update_rigid(i)
      call update_rigidm(i)
    end do
    call zmatfyc2()
#ifdef ENABLE_THREADS
    end if
#endif
  else
    write(ilog,*) 'Fatal. Encountered unsupported choice of degrees of freedom (',fycxyz,') in&
 & MPI_ReceiveCoordinates(...). This is most certainly a bug.'
    call fexit()
  end if
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
  if (dynni.EQV..true.) then
    if ((dyn_mode.eq.2).OR.(dyn_mode.eq.5)) then
      if (fycxyz.eq.1) then
        if (buf.EQV..true.) then
          call manage_imds(mpi_dvec(:,:,2),aone,tpi)
        end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
        call MPI_RECV(mpi_dvec(:,:,1),3*(ndyntorsn+totrbd),MPI_RTYPE,tnode,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.imdtag) then
          write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from node ',tnode+1,'. Expected ',imdtag,'.'
          call fexit()
        end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
        call manage_imds(mpi_dvec(:,:,1),atwo,tpi)
      else if (fycxyz.eq.2) then
        if (buf.EQV..true.) then
          call manage_cmds(mpi_dvec(:,:,2),aone,tpi)
        end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
        call MPI_RECV(mpi_dvec(:,:,1),3*n,MPI_RTYPE,tnode,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.cmdtag) then
          write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from node ',tnode+1,'. Expected ',cmdtag,'.'
          call fexit()
        end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
        call manage_cmds(mpi_dvec(:,:,1),atwo,tpi)
      else
        write(ilog,*) 'Fatal. Encountered unsupported choice of degrees of freedom (',fycxyz,') in&
   & MPI_ReceiveCoordinates(...). This is most certainly a bug.'
        call fexit()
      end if
    else if ((dyn_mode.eq.3).OR.(dyn_mode.eq.7)) then
      if (fycxyz.eq.1) then
        if (buf.EQV..true.) then
          call manage_ilds(mpi_dvec(:,:,2),aone,tpi)
        end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
        call MPI_RECV(mpi_dvec(:,:,1),4*(ndyntorsn+totrbd),MPI_RTYPE,tnode,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.ildtag) then
          write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from node ',tnode+1,'. Expected ',ildtag,'.'
          call fexit()
        end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
        call manage_ilds(mpi_dvec(:,:,1),atwo,tpi)
      else if (fycxyz.eq.2) then
        if (buf.EQV..true.) then
          call manage_clds(mpi_dvec(:,:,2),aone,tpi)
        end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
        call MPI_RECV(mpi_dvec(:,:,1),9*n,MPI_RTYPE,tnode,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.cldtag) then
          write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),' from node ',tnode+1,'. Expected ',cldtag,'.'
          call fexit()
        end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
        call manage_clds(mpi_dvec(:,:,1),atwo,tpi)
      else
        write(ilog,*) 'Fatal. Encountered unsupported choice of degrees of freedom (',fycxyz,') in&
   & MPI_ReceiveCoordinates(...). This is most certainly a bug.'
        call fexit()
      end if
    else if (dyn_mode.eq.1) then
!     do nothing
    else
      write(ilog,*) 'Fatal. Encountered unsupported dynamics method (',dyn_mode,') in&
   & MPI_ReceiveCoordinates(...). This is most certainly a bug.'
      call fexit()
    end if
  end if
!
  if (tpi.le.1) then
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
  end if
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
  if (dynni.EQV..true.) then
    if (fycxyz.eq.1) then
      if (re_velmode.eq.1) then
        call randomize_velocities(tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
      else if ((use_REMC.EQV..true.).AND.(re_velmode.eq.2).AND.(Tthere.EQV..true.)) then
        dummy = tnode + 1
        otherT = re_mat(dummy,Tdim)
        call rescale_velocities(otherT,tpi)
      end if
    else if (fycxyz.eq.2) then
      if (re_velmode.eq.1) then
        call randomize_cart_velocities(tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
      else if ((use_REMC.EQV..true.).AND.(re_velmode.eq.2).AND.(Tthere.EQV..true.)) then
        dummy = tnode + 1
        otherT = re_mat(dummy,Tdim)
        call rescale_cart_velocities(otherT,tpi)
      end if
    end if
  end if
  if (tpi.le.1) then
    call System_Clock(ee1)
    if (use_EMICRO.EQV..true.) curmassm = .false.
  end if
  if (use_mcgrid.EQV..true.) then
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      sta = thr_limits(3,tpi)
      sto = thr_limits(4,tpi)
    else
      sta = 1
      sto = nseq
    end if
    do i=sta,sto
#else
    do i=1,nseq
#endif
      call updateresgp(i)
    end do
  end if
  if (dynni.EQV..true.) then
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      call force3_threads(esterms,esterms_tr,esterms_lr,atrue,esave)
    else
      esave = force3(esterms,esterms_tr,esterms_lr,atrue)
    end if
#else
    esave = force3(esterms,esterms_tr,esterms_lr,atrue)
#endif
  else
    call energy3(esterms,atrue,esave,tpi)
  end if
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
  ens%insU = esave
  if (dynni.EQV..true.) call ensv_helper()
 333 format(100(g18.10,1x))
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
  if (tpi.le.1) then
    call System_Clock(ee2)
    time_energy = time_energy + 1.0*(ee2 - ee1)
  end if

!
end
!
!---------------------------------------------------------------------------
! 
!
subroutine MPI_ASMaster(istep,tpi)
!
  use iounit
  use atoms
  use accept
  use molecule
  use mpi
  use mpistuff
  use mcsums
  use energies
  use torsn
  use system
  use movesets
  use forces
  use sequen
  use ems
  use cutoffs
  use zmatrix
  use clusters
  use fyoc
  use grandensembles
  use interfaces
#ifdef ENABLE_THREADS
  use threads
#endif
!
  implicit none
!
  integer, INTENT(IN):: tpi,istep
!
  integer tpn,azero,ixl,cpibu,medi(mpi_nodes,8)
  integer i,j,ierr,masterrank,k,modstep
  integer nmap(2*mpi_nodes),msgtag(4)
  RTYPE random,force3,eee,edum(MAXENERGYTERMS),edum2(MAXENERGYTERMS),edum3(MAXENERGYTERMS)
  logical bcast,dynni,afalse,sayyes
  integer, ALLOCATABLE:: progind(:),invvec(:),iv2(:)
  integer(KIND=8) ee1,ee2
  real(KIND=4), ALLOCATABLE:: distv(:)
  real(KIND=4) ddv(mpi_nodes,2)
#ifdef ENABLE_THREADS
  integer OMP_GET_NUM_THREADS
  logical OMP_IN_PARALLEL
!
  if (tpi.le.0) then
    if (OMP_IN_PARALLEL().EQV..true.) then
      write(ilog,*) 'Fatal. When using hybrid MPI/multi-threaded code, MPI_ASMaster(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
      call fexit()
    end if
  end if
#endif
!
  masterrank = 0
  azero = 0
  if (pdb_analyze.EQV..false.) then
    dynni = .true.
    if ((dyn_mode.eq.1).OR.((in_dyncyc.EQV..false.).AND.((dyn_mode.eq.5).OR.(dyn_mode.eq.7).OR.(dyn_mode.eq.8)))) dynni = .false.
  end if
  afalse = .false.
  sayyes = .true.
!
  if (myrank.eq.masterrank) then
#ifdef ENABLE_THREADS
    if (tpi.gt.0) then
      tpn = OMP_GET_NUM_THREADS()
    else
      tpn = 1
    end if
#else
    tpn = 1
#endif
!   analyze data
    if (pdb_analyze.EQV..false.) then
      ixl = cstored
    else
      ixl = cstored/mpi_nodes
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
    if (pdb_analyze.EQV..false.) cstored = cstored*mpi_nodes
!   fix some parameters if necessary
    if (csmoothie.ge.(cstored/2-1)) then
      write(ilog,*) 'Warning. Selected smoothing window size is too large (FMCSC_CSMOOTHORDER ',csmoothie,').'
      csmoothie = cstored/2-1
    end if
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
!
!   cluster new data
    if (cstored.gt.0) then
      if (tpi.le.1) call System_Clock(ee1)
      k = -2
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
      if (cprepmode.gt.0) call preprocess_cludata(azero)
      if (cstored.gt.4) call repopulate_weights()
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
#ifdef ENABLE_THREADS
        call birch_clustering_threads(k,nstruccls,tpi)
!$OMP BARRIER
#else
        call birch_clustering(k,nstruccls)
#endif
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!     get approximate MST
      call gen_MST_from_treeclustering(tpi)
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
      cpibu = cprogindstart
      allocate(distv(cstored))
      allocate(progind(cstored))
      allocate(invvec(cstored+2))
      allocate(iv2(cstored))
      cprogindstart = birchtree(c_nhier+1)%cls(1)%center ! center of largest cluster in hierarchical tree
!     generate progress index
      call gen_progind_from_adjlst(approxmst,cprogindstart,progind,distv,invvec,iv2)
      cprogindstart = cpibu
!     sort replicas' final snapshots by PI position, by PI distance from nearest other final ss, and by associated SST edge length
!     all sorts in decreasing order; maps are yielded in medi(:,5:7) with associated ranks per replica
!     from the ranks, get a flipped composite rank and recover map in medi(:,8)
      do i=ixl,cstored,ixl
        medi(i/ixl,1) = invvec(i+1)
        medi(i/ixl,2) = i/ixl
        write(*,*) 'Final pos ',i/ixl,medi(i/ixl,1) 
      end do
      i = 1
      j = mpi_nodes
      call merge_sort(ldim=mpi_nodes,up=afalse,list=medi(:,1),olist=medi(:,3),ilo=i,ihi=j,idxmap=medi(:,2),olist2=medi(:,5))
      do i=1,mpi_nodes
        medi(medi(i,5),2) = i
      end do
      medi(:,5) = medi(:,2)
      do i=1,mpi_nodes
        if (medi(i,5).eq.1) then
          medi(i,1) = invvec(i*ixl+1) - medi(medi(i,5)+1,3)
        else if (medi(i,5).eq.mpi_nodes) then
          medi(i,1) = medi(medi(i,5)-1,3) - invvec(i*ixl+1)
        else
         medi(i,1) = min(invvec(i*ixl+1)-medi(medi(i,5)+1,3),medi(medi(i,5)-1,3)-invvec(i*ixl+1))
        end if
        medi(i,2) = i
        write(*,*) 'Diff from ',i,medi(i,1)
      end do
      i = 1
      j = mpi_nodes
      call merge_sort(ldim=mpi_nodes,up=afalse,list=medi(:,1),olist=medi(:,3),ilo=i,ihi=j,idxmap=medi(:,2),olist2=medi(:,6))
      do i=1,mpi_nodes
        medi(medi(i,6),2) = i
      end do
      medi(:,6) = medi(:,2)
      do i=ixl,cstored,ixl
        ddv(i/ixl,1) = distv(invvec(i+1))
        medi(i/ixl,2) = i/ixl
        write(*,*) 'Length ',i/ixl,ddv(i/ixl,1)
      end do
      i = 1
      j = mpi_nodes
      call merge_sort(ldim=mpi_nodes,up=afalse,list=ddv(:,1),olist=ddv(:,2),ilo=i,ihi=j,idxmap=medi(:,2),olist2=medi(:,7))
      do i=1,mpi_nodes
        medi(medi(i,7),2) = i
      end do
      medi(:,7) = medi(:,2)
      do i=1,mpi_nodes
        medi(i,1) = sum(mpi_nodes+1-medi(i,5:7))
        medi(i,2) = i
      end do
      i = 1
      j = mpi_nodes
      call merge_sort(ldim=mpi_nodes,up=afalse,list=medi(:,1),olist=medi(:,3),ilo=i,ihi=j,idxmap=medi(:,2),olist2=medi(:,8))
!
!     now use the computed parameters to determine whether to reseed each replica
!     note that a replica being reseeded means that its own conformation is lost (no swaps allowed)
      do i=1,mpi_nodes
        nmap(i) = i-1
      end do
      do j=re_aux(7)+1,mpi_nodes
        i = medi(j,8)
        k = min(re_aux(7),ceiling(random()*re_aux(7)))
        write(*,*) 'proposing ',i,' to ',k
        if (random().lt.(1.0*(medi(k,3)-medi(j,3))/(1.0*(medi(1,3)-medi(mpi_nodes,3)+1)))) then
          nmap(i) = medi(k,8)-1
          write(*,*) 'conditional accept'
        end if
      end do
!
!     generate distribution parameters for positions of replica trajectories in PI
      medi(:,1:4) = 0
      do i=1,cstored
        j = (progind(i)+ixl-1)/ixl
        medi(j,1) = medi(j,1) + 1
        if (medi(j,1).le.(ixl/4))   medi(j,3) = i
        if (medi(j,1).le.(ixl/2))   medi(j,2) = i
        if (medi(j,1).le.3*(ixl/4)) medi(j,4) = i
      end do
!      write(*,*) '1/4IAN: ',medi(:,3)
!      write(*,*) 'MEDIAN: ',medi(:,2)
!      write(*,*) '3/4IAN: ',medi(:,4)
!      write(*,*) medi(:,4)-medi(:,3)
!
!     uniqueness overwrite
      do i=1,mpi_nodes
        if ((medi(i,4)-medi(i,3)).lt.ixl) nmap(i) = i-1
      end do
!
      deallocate(iv2)
      deallocate(invvec)
      deallocate(progind)
      deallocate(distv)
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
      if (pdb_analyze.EQV..false.) then
        if (allocated(approxmst).EQV..true.) then
          do i=max(tpi,1),size(approxmst),tpn
            if (allocated(approxmst(i)%adj).EQV..true.) deallocate(approxmst(i)%adj)
            if (allocated(approxmst(i)%dist).EQV..true.) deallocate(approxmst(i)%dist)
          end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
        end if
        if (allocated(birchtree).EQV..true.) then
          do i=1,c_nhier+1
            do j=max(tpi,1),birchtree(i)%ncls,tpn
              if (allocated(birchtree(i)%cls(j)%snaps).EQV..true.) deallocate(birchtree(i)%cls(j)%snaps)
              if (allocated(birchtree(i)%cls(j)%tmpsnaps).EQV..true.) deallocate(birchtree(i)%cls(j)%tmpsnaps)
              if (allocated(birchtree(i)%cls(j)%sums).EQV..true.) deallocate(birchtree(i)%cls(j)%sums)
              if (allocated(birchtree(i)%cls(j)%children).EQV..true.) deallocate(birchtree(i)%cls(j)%children)
            end do
          end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
          do i=max(tpi,1),c_nhier+1,tpn
            deallocate(birchtree(i)%cls)
          end do
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
        end if
        if (tpi.le.1) then
          deallocate(approxmst)
          deallocate(birchtree)
          call System_Clock(ee2)
          time_analysis = time_analysis + 1.0*(ee2 - ee1)
        end if
      end if
    else
      do i=1,mpi_nodes
        nmap(i) = i-1
      end do
    end if
!
!   print trace so breaks can be removed
    if (inst_retr.EQV..true.) then
#ifdef ENABLE_THREADS
!$OMP MASTER
#endif
      nmap(1:mpi_nodes) = nmap(1:mpi_nodes) + 1
      if (pdb_analyze.EQV..false.) then
        call MPI_REPrtTrace(istep,iretr,nmap(:))
      else
        i = re_aux(4)*re_aux(8)
        call MPI_REPrtTrace(i,iretr,nmap(:))
      end if
      nmap(1:mpi_nodes) = nmap(1:mpi_nodes) - 1
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
    end if
!
  end if
!
  if (pdb_analyze.EQV..false.) then
! to avoid confusion, reset cstored for all nodes
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    cstored = 0
!
!   Broadcast mapping
    call System_Clock(ee1)
    msgtag(1) = 3361
    msgtag(2) = msgtag(1) + 1
    msgtag(3) = msgtag(1) + 2
    msgtag(4) = msgtag(1) + 3
    bcast = .true.
!
    k = 2*mpi_nodes
    call MPI_ALLINTS1D(msgtag(1:4),azero,ierr,bcast,k,nmap,use_MPIcolls) ! always uses MPI master as BC root
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
#ifdef ENABLE_THREADS
    allocate(thr_hlper(mpi_nodes,1))
    thr_hlper(1:mpi_nodes,1) = nmap((mpi_nodes+1):(2*mpi_nodes))
!$OMP END MASTER
!$OMP BARRIER
    nmap((mpi_nodes+1):(2*mpi_nodes)) = thr_hlper(1:mpi_nodes,1)
#endif
!
    if (nmap(mpi_nodes+myrank+1).eq.myrank) then
      do i=1,mpi_nodes
        if (i.eq.(myrank+1)) cycle
        if (nmap(mpi_nodes+i).eq.myrank) then
          j = i-1
          call MPI_SendCoordinates(j,afalse,dynni,tpi)
        end if
      end do
    else
      j = nmap(mpi_nodes+myrank+1)
      call MPI_ReceiveCoordinates(j,afalse,dynni,tpi)
    end if
!
!   handle restart files separately from normal flow
    if (pdb_analyze.EQV..false.) then ! automatic
      modstep = mod(istep,rstout)
      if (modstep.eq.0) then

        if (tpi.le.1) then
          eee = esave
          edum(:) = esterms(:)
          edum2(:) = esterms_tr(:)
          edum3(:) = esterms_lr(:)
        end if
        if (dyn_mode.eq.1) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
          call energy(esterms,esave,tpi)
        else
#ifdef ENABLE_THREADS
!$OMP BARRIER
          call force3_threads(esterms,esterms_tr,esterms_lr,sayyes,esave)
#else
          esave = force3(esterms,esterms_tr,esterms_lr,sayyes)
#endif
        end if
#ifdef ENABLE_THREADS
!$OMP SINGLE
#endif
        call prt_restart(esave)
#ifdef ENABLE_THREADS
!$OMP END SINGLE
#endif
        if (tpi.le.1) then
          esave = eee
          esterms(:) = edum(:)
          esterms_tr(:) = edum2(:)
          esterms_lr(:) = edum3(:)
        end if
      end if
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
    deallocate(thr_hlper)
!$OMP END SINGLE
#endif
  end if
!
end
!
!---------------------------------------------------------------------------
! 
subroutine MPI_RETrajMode(istep)
!
  use iounit
  use atoms
  use accept
  use molecule
  use mpi
  use mpistuff
  use mcsums
!
  implicit none
!
  integer istep,azero,i,j,k,ierr,masterrank,msgtag(4),msgsz,xyztag
  integer mstatus(MPI_STATUS_SIZE)
  logical bcast
  integer(KIND=8) ee1,ee2
!
  if (re_aux(3).ne.1) return
!
  masterrank = 0
  xyztag = 13
  msgtag(1) = 61
  msgtag(2) = msgtag(1) + 1
  msgtag(3) = msgtag(1) + 2
  msgtag(4) = msgtag(1) + 3
  azero = 0
  bcast = .true.
!
  call System_Clock(ee2)
  mpi_lmap(1:mpi_nodes,1) = mpi_lmap(1:mpi_nodes,2)
  if (myrank.eq.masterrank) then
    call MPI_REReadTrace(istep,re_aux(6))
  end if
  msgsz = 2*mpi_nodes 
  call MPI_ALLINTS1D(msgtag(1:4),azero,ierr,bcast,msgsz,mpi_lmap(:,2),use_MPIcolls)
 !
  if (myrank.ne.masterrank) then
    mpi_lmap(1:mpi_nodes,2) = mpi_lmap((mpi_nodes+1):(2*mpi_nodes),2)
  end if
!
  if (mpi_lmap(myrank+1,2).ne.(myrank+1)) then
    do i=1,mpi_nodes
      if (mpi_lmap(i,2).eq.(myrank+1)) then
        j = i-1
        exit
      end if
    end do
    k = mpi_lmap(myrank+1,2) - 1
    if (mpi_lmap(myrank+1,2).lt.(myrank+1)) then
      mpi_cdsvec(1:n,1,2) = x(1:n)
      mpi_cdsvec(1:n,2,2) = y(1:n)
      mpi_cdsvec(1:n,3,2) = z(1:n)
      call MPI_RECV(mpi_cdsvec(:,:,1),3*n,MPI_RTYPE,j,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.xyztag) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from node ',j,'. Expected ',xyztag,'.'
        call fexit()
      end if
      x(1:n) = mpi_cdsvec(1:n,1,1)
      y(1:n) = mpi_cdsvec(1:n,2,1)
      z(1:n) = mpi_cdsvec(1:n,3,1)
      call MPI_SEND(mpi_cdsvec(:,:,2),3*n,MPI_RTYPE,k,xyztag,MPI_COMM_WORLD,ierr)
    else
      mpi_cdsvec(1:n,1,1) = x(1:n)
      mpi_cdsvec(1:n,2,1) = y(1:n)
      mpi_cdsvec(1:n,3,1) = z(1:n)
      call MPI_SEND(mpi_cdsvec(:,:,1),3*n,MPI_RTYPE,k,xyztag,MPI_COMM_WORLD,ierr)
      call MPI_RECV(mpi_cdsvec(:,:,2),3*n,MPI_RTYPE,j,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.xyztag) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from node ',j,'. Expected ',xyztag,'.'
        call fexit()
      end if
      x(1:n) = mpi_cdsvec(1:n,1,2)
      y(1:n) = mpi_cdsvec(1:n,2,2)
      z(1:n) = mpi_cdsvec(1:n,3,2)
    end if
  end if
  call System_Clock(ee1)
  time_comm = time_comm + 1.0*(ee2 - ee1)
  do i=1,nmol
    call genzmat(i)
    call update_rigid(i)
    call update_rigidm(i)
  end do
  call zmatfyc2()
!
end
!
!---------------------------------------------------------------------
!
subroutine MPI_REWritePB()
!
  use mpistuff
  use system
!
  implicit none
!
  integer freeunit,irepb,masterrank,i,j
  character(100) repbfile
  character(re_aux(10)) nod
  logical exists
!
  if (pdb_analyze.EQV..true.) return
  masterrank = 0
!
  if (myrank.ne.0) then
    write(*,*) 'Fatal. Called MPI_Write_REPB with slave node.'
    call fexit()
  end if
!
! setup filehandle
  call int2str(myrank,nod,re_aux(10))
  repbfile = 'N_'//nod(1:re_aux(10))//'_PROBABILITIES.dat'
  inquire(file=repbfile,exist=exists)
  if(exists) then
  irepb = freeunit()
    open (unit=irepb,file=repbfile,status='old')
    close(unit=irepb,status='delete')
  end if
  irepb = freeunit()
  open (unit=irepb,file=repbfile,status='new')
!
  do i=1,re_conditions
    do j=1,re_conditions
      if (i.eq.j) then
        re_probs(i,j) = 1.0
      else if (re_trans(i,j).gt.0) then
        re_probs(i,j) = re_probs(i,j)/(1.0*re_trans(i,j))
      else
        re_probs(i,j) = 0.0
      end if
    end do
  end do
!
 65   format(100000(g16.9,1x))
  do i=1,re_conditions
    write(irepb,65) (re_probs(i,j),j=1,re_conditions)
  end do
!
  close(unit=irepb)
!
end
!
!---------------------------------------------------------------------
!
subroutine MPI_REPrtTrace(istep,ihandle,prtmap)
!
  use mcsums
  use mpistuff
!
  implicit none
!
  integer bk,bk2,azero,istep,ihandle,prtmap(mpi_nodes)
  character(max(1,ceiling(log10(1.0*mpi_nodes+0.5)))) bnrs
  character(max(1,ceiling(log10(1.0*ceiling(log10(1.0*mpi_nodes+0.5)))))) bnrs2
  character(MAXSTRLEN) fmtstr
!
  azero = 0
!
  bk = max(1,ceiling(log10(1.0*mpi_nodes+0.5)))
  call int2str(mpi_nodes,bnrs,azero)
  bk2 = max(1,ceiling(log10(1.0*ceiling(log10(1.0*mpi_nodes+0.5)))))
  azero = 0
  call int2str(bk,bnrs2,azero)
  fmtstr='(i12,1x,'//bnrs(1:bk)//'(i'//bnrs2(1:bk2)//',1x))'
  write(ihandle,fmtstr) istep,prtmap
!
end
!
!---------------------------------------------------------------------
!
subroutine MPI_REReadTrace(istep,kstep)
!
  use mcsums
  use mpistuff
  use iounit
  use ISO_FORTRAN_ENV, ONLY: IOSTAT_END
!
  implicit none
!
  integer istep,kstep,lstep,iomess
!
  lstep = re_aux(5)/re_aux(4) ! integer division
  if (kstep.eq.(istep+lstep)*re_aux(4)) then
    mpi_lmap(1:mpi_nodes,2) = mpi_lmap((mpi_nodes+1):(2*mpi_nodes),1) ! we have read info for this exact step already: switch now
  end if
  do while (kstep.lt.(istep+lstep)*re_aux(4)) 
    read(iretr,*,iostat=iomess) kstep,mpi_lmap((mpi_nodes+1):(2*mpi_nodes),1)
    if (kstep.gt.(istep+lstep)*re_aux(4)) exit
    if (iomess.eq.IOSTAT_END) then
      exit
    else if (iomess.eq.2) then
      write(ilog,*) 'Fatal. File I/O error while reading replica exchange trace file (...).'
      call fexit()
    end if
    mpi_lmap(1:mpi_nodes,2) = mpi_lmap((mpi_nodes+1):(2*mpi_nodes),1)
  end do
!
end
!
!-------------------------------------------------------------------------------
!
subroutine MPI_REOverlap(tpi)
!
  use mpistuff
  use system
  use units
  use molecule
  use mcsums
  use energies
  use iounit
!
  implicit none
!
  integer, INTENT(IN):: tpi
!
  integer i,j,nioc,i_start,i_end
  RTYPE evec(mpi_nodes),fvec(mpi_nodes)
  logical prtout
#ifdef ENABLE_THREADS
  logical OMP_IN_PARALLEL
!
  if ((tpi.le.0).AND.(OMP_IN_PARALLEL().EQV..true.)) then
    write(ilog,*) 'Fatal. When using multi-threaded code, MPI_REOverlap(...) must always be called &
 &by all threads of the binding region (instead got a thread identifier of zero). This is a bug.'
    call fexit()
  end if
#endif
!
  if (reol_all.EQV..true.) then
    i_start = 1
    i_end = mpi_nodes
  else
    i_start = max((myrank+1)-1,1)
    i_end = min((myrank+1)+1,mpi_nodes)
  end if
  evec(:) = mpi_evec(:)
  fvec(:) = mpi_fvec(:)
#ifdef ENABLE_THREADS
!$OMP BARRIER
#endif
!
  if (Tonly.EQV..true.) then
!   we're just storing average energies
    if (tpi.le.1) re_olap(1,1) = re_olap(1,1) + esave
  else if (Tthere.EQV..true.) then
    if (dyn_mode.eq.1) then
      call lamenergy(mpi_evec,mpi_fvec,tpi)
    else
      call lamforce(mpi_evec,mpi_fvec,tpi)
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    do i=i_start,i_end !1,mpi_nodes
      re_olap(i,1) = re_olap(i,1) + mpi_evec(i)*re_mat(i,Tdim)*gasconst - mpi_evec(myrank+1)/invtemp
      re_olap(i,2) = re_olap(i,2) + exp(-invtemp*(mpi_evec(i)*re_mat(i,Tdim)*gasconst - mpi_evec(myrank+1)/invtemp))
    end do
    re_olap(1,3) = re_olap(1,3) + esave
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
  else if (noTI.EQV..true.) then
    if (dyn_mode.eq.1) then
      call lamenergy(mpi_evec,mpi_fvec,tpi)
    else
      call lamforce(mpi_evec,mpi_fvec,tpi)
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    do i=i_start,i_end !1,mpi_nodes
      re_olap(i,1) = re_olap(i,1) + (mpi_evec(i)-mpi_evec(myrank+1))/invtemp
      re_olap(i,2) = re_olap(i,2) + exp(-(mpi_evec(i)-mpi_evec(myrank+1)))
    end do
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
  else
    if (dyn_mode.eq.1) then
      call lamenergy(mpi_evec,mpi_fvec,tpi)
    else
      call lamforce(mpi_evec,mpi_fvec,tpi)
    end if
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    do i=i_start,i_end !1,mpi_nodes
      re_olap(i,1) = re_olap(i,1) + (mpi_evec(i)-mpi_evec(myrank+1))/invtemp
      re_olap(i,2) = re_olap(i,2) + exp(-(mpi_evec(i)-mpi_evec(myrank+1)))
      re_olap(i,3) = re_olap(i,3) + mpi_fvec(i)
    end do
#ifdef ENABLE_THREADS
!$OMP END MASTER
#endif
  end if
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP SINGLE
#endif
 62   format(1000000(g18.11,1x))
  nreolap = nreolap + 1
  prtout = .false.
  if (inst_reol.gt.0) then
    nioc = nstep/re_olcalc
    if (mod(nioc,inst_reol).eq.0) prtout = .true.
  end if
  if (prtout.EQV..true.) then
    if (Tthere.EQV..true.) then
      write(iremc,62)(re_mat(j,Tdim)*gasconst*mpi_evec(j),j=i_start,i_end)
    else
      write(iremc,62) (mpi_evec(j)/invtemp,j=i_start,i_end)
    end if
  end if
#ifdef ENABLE_THREADS
  mpi_evec(:) = evec(:)
  mpi_fvec(:) = fvec(:)
!$OMP END SINGLE
#endif
!
end
!
!------------------------------------------------------------------------------
!
subroutine MPI_WriteREOlap()
!
  use mpistuff
  use system
  use units
!
  implicit none
!
  integer OLNR
  parameter (OLNR=4)
  integer freeunit,ireol,i,j,i_start,i_end,ii,jj
  character(100) reolfile
  character(re_aux(10)) nod
  logical exists
!
 77   format(i6,1x,100(g14.7,1x))
!
  call int2str(myrank,nod,re_aux(10))
  reolfile = 'N_'//nod(1:re_aux(10))//'_OVERLAP.dat'
  call strlims(reolfile,ii,jj)
  inquire(file=reolfile(ii:jj),exist=exists)
  if (exists.EQV..true.) then
    ireol = freeunit()
    open(unit=ireol,file=reolfile(ii:jj),status='old')
    close(unit=ireol,status='delete')
  end if
!
  if (reol_all.EQV..true.) then
    i_start = 1
    i_end = mpi_nodes
  else
    i_start = max((myrank+1)-1,1)
    i_end = min((myrank+1)+1,mpi_nodes)
  end if
!
  if (nreolap.gt.0) then
    do i=1,mpi_nodes
      do j=1,OLNR
        re_olap(i,j) = re_olap(i,j)/(1.0*nreolap)
      end do
    end do
    ireol = freeunit()
    open(unit=ireol,file=reolfile(ii:jj),status='new')
    if (Tonly.EQV..true.) then
      write(ireol,77) myrank+1,kelvin,re_olap(1,1)
    else if (Tthere.EQV..true.) then
      do i=i_start,i_end !1,mpi_nodes
        if (i.eq.myrank+1) then
          write(ireol,77) i,re_mat(i,Tdim),re_olap(i,1),-log(re_olap(i,2))/invtemp,re_olap(1,3)
        else
          write(ireol,77) i,re_mat(i,Tdim),re_olap(i,1),-log(re_olap(i,2))/invtemp,0.0d0
        end if
      end do
    else if (noTI.EQV..true.) then
      do i=i_start,i_end !1,mpi_nodes
        write(ireol,77) i,re_olap(i,1),-log(re_olap(i,2))/invtemp
      end do
    else
      do i=i_start,i_end !1,mpi_nodes
        write(ireol,77) i,re_olap(i,1),-log(re_olap(i,2))/invtemp,re_olap(i,3)
      end do
    end if
    close(unit=ireol)
  end if
!
end
!
!------------------------------------------------------------------------------
!
! all the allocation statements must be robust and match the actual size
! i.e., values like nmoltyp, nmol, nseq, or even ampc%do_ss and the like
! MUST NOT CHANGE
!
subroutine MPI_AVGCollect()
!
  use iounit
  use system
  use sequen
  use energies
  use molecule
  use polyavg
  use paircorr
  use torsn
  use mpi
  use mpistuff
  use contacts
  use mcsums
  use atoms
  use dipolavg
  use diffrac
  use dssps
  use fos
  use ems
  use clusters
  use grandensembles
!
  implicit none
!
  integer MPI_POLSZ
  parameter (MPI_POLSZ=12)
!
  integer i,masterrank,tpi,rmgtag,rmstag,contag1,contag2,segtag1,savtag3
  integer mstatus(MPI_STATUS_SIZE),ierr,j,k,l,ii,idummr,aone,pertag1,pctag4
  integer pctag1,pctag2,pctag3,segtag2,pertag2,savtag,segtag3
  integer lcttag,rmmtag,pertag3,polatag,polatag2,mt,polatag3
  integer polatag4,polatag5,polatag6,polatag7,polatag8,clutag1
  integer clutag2,clutag3,pctag5,segtag4,polatag9,polatag10,diptag1
  integer diptag2,diptag3,diptag4,jctag,jctag2,diffrtag1
  integer diffrtag2,segtag5,inttag(4),dssptag1,dssptag2,dssptag3
  integer pctag6,pctag7,savtag2,contag3,segtag6,dssptag4,emtag1,emtag2,tbrktag,tbrktag2,nhtag1
  RTYPE, ALLOCATABLE:: map1(:),map2a(:,:),map2b(:,:),map2c(:,:)
  RTYPE, ALLOCATABLE:: map3a(:,:,:),map3b(:,:,:),map2d(:,:)
  RTYPE, ALLOCATABLE:: map3c(:,:,:),map3d(:,:,:)
  integer, ALLOCATABLE:: imap1a(:),imap2(:,:),imap1b(:)
  integer(KIND=8), ALLOCATABLE:: limap2(:,:),limap1(:)
  integer(KIND=8) ee1,ee2
  logical sayyes
#ifdef ENABLE_THREADS
  integer OMP_GET_THREAD_NUM
#endif
!
  sayyes = .true.
  aone = 1
  masterrank = 0
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  polatag = 92
  polatag2 = 93
  polatag3 = 95
  polatag4 = 101
  polatag5 = 97
  polatag6 = 98
  polatag7 = 99
  polatag8 = 100
  polatag9 = 107
  polatag10 = 108
  rmgtag = 80
  rmstag = 81
  rmmtag = 96
  segtag1 = 79
  segtag2 = 86
  contag1 = 78
  contag2 = 77
  pertag1 = 82
  pertag2 = 87
  pertag3 = 91
  pctag1 = 83
  pctag2 = 84
  pctag3 = 85
  savtag = 88
  segtag3 = 89
  pctag4 = 90
  lcttag = 94
  clutag1 = 102
  clutag2 = 103
  clutag3 = 104
  pctag5 = 105
  segtag4 = 106
  diptag1 = 109
  diptag2 = 110
  diptag3 = 111
  diptag4 = 112
  jctag = 113
  jctag2 = 114
  diffrtag1 = 115
  diffrtag2 = 116
  segtag5 = 117
  inttag(1) = 118
  inttag(2) = 119
  inttag(3) = 120
  inttag(4) = 121
  dssptag1 = 122
  dssptag2 = 123
  dssptag3 = 124
  pctag6 = 125
  pctag7 = 126
  savtag2 = 127
  contag3 = 128
  segtag6 = 129
  dssptag4 = 130
  savtag3 = 131
  emtag1 = 132
  emtag2 = 133
  tbrktag = 134
  tbrktag2 = 135
  nhtag1 = 136
!
  if (myrank.ne.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGCollect() with slave node.'
    call fexit()
  end if
!
! first average cheap (linear) polymer props
!
  if (polcalc.le.nsim) then
!   this re-scaling is really only necessary in the GC ensemble
    do mt=1,nangrps
      rgavg(mt) = npolavg(mt)*rgavg(mt)
      rg2avg(mt) = npolavg(mt)*rg2avg(mt)
      rgpcavg(mt,1) = npolavg(mt)*rgpcavg(mt,1)
      rgpcavg(mt,2) = npolavg(mt)*rgpcavg(mt,2)
      rgpcavg(mt,3) = npolavg(mt)*rgpcavg(mt,3)
      asphavg(mt) = npolavg(mt)*asphavg(mt)
      acylavg(mt) = npolavg(mt)*acylavg(mt)
      rgterm1avg(mt) = npolavg(mt)*rgterm1avg(mt)
      rgterm2avg(mt) = npolavg(mt)*rgterm2avg(mt)
      rgterm12avg(mt) = npolavg(mt)*rgterm12avg(mt)
!      dlavg(mt) = npolavg(mt)*dlavg(mt)
!      dlstavg(mt) = npolavg(mt)*dlstavg(mt)
      reteavg(mt) = npolavg(mt)*reteavg(mt)
      vtavg(mt) = npolavg(mt)*vtavg(mt)
    end do
    call System_Clock(ee1)
!   allocate
    allocate(map2a(nangrps,MPI_POLSZ))
    allocate(imap1a(nangrps))
    allocate(map3a(nangrps,TPBINS,DLSTBINS))
    allocate(map2b(nangrps,maxrgbins))
    do i=1,mpi_nodes-1
      call MPI_RECV(map2a,nangrps*MPI_POLSZ,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.polatag) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',polatag,'.'
        call fexit()
      end if
      call MPI_RECV(imap1a,nangrps,MPI_INTEGER,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.polatag2) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',polatag2,'.'
        call fexit()
      end if
      do mt=1,nangrps
        rgavg(mt) = rgavg(mt) + map2a(mt,1)
        rg2avg(mt) = rg2avg(mt) + map2a(mt,11)
        vtavg(mt) = vtavg(mt) + map2a(mt,12)
        rgpcavg(mt,1) = rgpcavg(mt,1) + map2a(mt,2)
        rgpcavg(mt,2) = rgpcavg(mt,2) + map2a(mt,3)
        rgpcavg(mt,3) = rgpcavg(mt,3) + map2a(mt,4)
        asphavg(mt) = asphavg(mt) + map2a(mt,5)
        acylavg(mt) = acylavg(mt) + map2a(mt,6)
        rgterm1avg(mt) = rgterm1avg(mt) + map2a(mt,7)
        rgterm2avg(mt) = rgterm2avg(mt) + map2a(mt,8)
        rgterm12avg(mt) = rgterm12avg(mt) + map2a(mt,9)
!        dlstavg(mt) = dlstavg(mt) + map2a(mt,7)
!        dlavg(mt) = dlavg(mt) + map2a(mt,8)
        reteavg(mt) = reteavg(mt) + map2a(mt,10)
        npolavg(mt) = npolavg(mt) + imap1a(mt)
      end do
      call MPI_RECV(map3a,nangrps*TPBINS*DLSTBINS,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.polatag3) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',polatag3,'.'
        call fexit()
      end if
      do mt=1,nangrps
        do j=1,TPBINS
          do k=1,DLSTBINS
            rdhist(mt,j,k) = rdhist(mt,j,k) + map3a(mt,j,k)
          end do
        end do
      end do
      call MPI_RECV(map2b,nangrps*maxrgbins,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.polatag4) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',polatag4,'.'
        call fexit()
      end if
      do mt=1,nangrps
        do j=1,maxrgbins
          retehist(mt,j) = retehist(mt,j) + map2b(mt,j)
        end do
      end do
      call MPI_RECV(map2b,nangrps*maxrgbins,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.polatag9) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',polatag9,'.'
        call fexit()
      end if
      do mt=1,nangrps
        do j=1,maxrgbins
          densproavg(mt,j) = densproavg(mt,j) + map2b(mt,j)
        end do
      end do
      call MPI_RECV(map2b,nangrps*maxrgbins,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.polatag10) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',polatag10,'.'
        call fexit()
      end if
      do mt=1,nangrps
        do j=1,maxrgbins
          rghist(mt,j) = rghist(mt,j) + map2b(mt,j)
        end do
      end do
    end do
    do mt=1,nangrps
      if (npolavg(mt).gt.0) then
        rgavg(mt) = rgavg(mt)/npolavg(mt)
        rg2avg(mt) = rg2avg(mt)/npolavg(mt)
        vtavg(mt) = vtavg(mt)/npolavg(mt)
        rgpcavg(mt,1) = rgpcavg(mt,1)/npolavg(mt)
        rgpcavg(mt,2) = rgpcavg(mt,2)/npolavg(mt)
        rgpcavg(mt,3) = rgpcavg(mt,3)/npolavg(mt)
        asphavg(mt) = asphavg(mt)/npolavg(mt)
        acylavg(mt) = acylavg(mt)/npolavg(mt)
        rgterm1avg(mt) = rgterm1avg(mt)/npolavg(mt)
        rgterm2avg(mt) = rgterm2avg(mt)/npolavg(mt)
        rgterm12avg(mt) = rgterm12avg(mt)/npolavg(mt)
        dlstavg(mt) = 1.0 - 3.0*rgterm12avg(mt)
        dlavg(mt) = 1.0 - 3.0*(rgterm1avg(mt)/rgterm2avg(mt))
!        dlstavg(mt) = dlstavg(mt)/npolavg(mt)
!        dlavg(mt) = dlavg(mt)/npolavg(mt)
        reteavg(mt) = reteavg(mt)/npolavg(mt)
      end if
    end do
!   de-allocate
    deallocate(map2a)
    deallocate(imap1a)
    deallocate(map3a)
    deallocate(map2b)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!   and finalize
    call prt_rdhist()
  end if
!
! expensive polymer stuff (quadratic: internal scaling, scattering)
!
  if (rhcalc.le.nsim) then
    call System_Clock(ee1)
!   allocate proper memory
    allocate(imap1a(nangrps))
    allocate(map2a(nangrps,nqv))
    allocate(map2b(nangrps,longestmol))
    do i=1,mpi_nodes-1
      call MPI_RECV(map2a,nangrps*nqv,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.polatag5) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',polatag5,'.'
        call fexit()
      end if
      call MPI_RECV(imap1a,nangrps,MPI_INTEGER,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.polatag6) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',polatag6,'.'
        call fexit()
      end if
      do mt=1,nangrps
        do j=1,nqv
          pscavg(mt,j) = pscavg(mt,j) + map2a(mt,j)
        end do
        nsctavg(mt) = nsctavg(mt) + imap1a(mt)
      end do
      call MPI_RECV(map2b,nangrps*longestmol,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.polatag7) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',polatag7,'.'
        call fexit()
      end if
      call MPI_RECV(imap1a,nangrps,MPI_INTEGER,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.polatag8) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',polatag8,'.'
        call fexit()
      end if
      do mt=1,nangrps
        do j=1,longestmol
          rhavg(mt,j) = rhavg(mt,j) + map2b(mt,j)
        end do
        nrhavg(mt) = nrhavg(mt) + imap1a(mt)
      end do
    end do
!   deallocate
    deallocate(imap1a)
    deallocate(map2a)
    deallocate(map2b)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!   and process
    call prt_hydrorad()
  end if
!
! backbone angular statistics maps
!
  if (angcalc.le.nsim) then
    call System_Clock(ee1)
!   first allocate proper memory
    allocate(map2a(torbins,torbins))
    allocate(imap1a(size(jccnt)))
    allocate(map1(nseq))
    if (nrspecrama.gt.0) then
      allocate(map3a(rama_alsz(1),torbins,torbins))
    end if
    if (nrmolrama.gt.0) then
      allocate(map3b(rama_alsz(2),torbins,torbins))
    end if
    do i=1,mpi_nodes-1
!     global
      call MPI_RECV(map2a,torbins*torbins,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.rmgtag) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',rmgtag,'.'
        call fexit()
      end if
!      write(*,*) 'SLRGRAM: ',rmap(1,3)
!      write(*,*) 'MAHGRAM: ',genfymap(1,3)
      do j=1,torbins
        do k=1,torbins
          genfymap(j,k) = genfymap(j,k) + map2a(j,k)
        end do
      end do
!     residue-specific maps
      if (nrspecrama.gt.0) then
        call MPI_RECV(map3a,rama_alsz(1)*torbins*torbins,&
 &MPI_RTYPE,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.rmstag) then
          write(*,*) 'Fatal. Received bad message (',&
 &mstatus(MPI_TAG),') from slave ',i,'. Expected ',rmstag,'.'
          call fexit()
        end if
!      write(*,*) 'SLRSRAM: ',smap(1,1,3),' ',smap(2,1,3)
!      write(*,*) 'MAHSRAM: ',specfymap(1,1,3),' ',specfymap(2,1,3)
        do ii=1,nrspecrama
          do j=1,torbins
            do k=1,torbins
              specfymap(ii,j,k) = specfymap(ii,j,k) + map3a(ii,j,k)
            end do
          end do
        end do
      end if
!     analysis group-specific maps
      if (nrmolrama.gt.0) then
        call MPI_RECV(map3b,rama_alsz(2)*torbins*torbins,&
 &MPI_RTYPE,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.rmmtag) then
          write(*,*) 'Fatal. Received bad message (',&
 &mstatus(MPI_TAG),') from slave ',i,'. Expected ',rmmtag,'.'
          call fexit()
        end if
!      write(*,*) 'SLRSRAM: ',smap(1,1,3),' ',smap(2,1,3)
!      write(*,*) 'MAHSRAM: ',specfymap(1,1,3),' ',specfymap(2,1,3)
        do ii=1,nrmolrama
          do j=1,torbins
            do k=1,torbins
              molfymap(ii,j,k) = molfymap(ii,j,k) + map3b(ii,j,k)
            end do
          end do
        end do
      end if
!     J-coupling data
      call MPI_RECV(map1,nseq,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.jctag) then
        write(*,*) 'Fatal. Received bad message (',&
 &mstatus(MPI_TAG),') from slave ',i,'. Expected ',jctag,'.'
        call fexit()
      end if
      call MPI_RECV(imap1a,size(jccnt),MPI_INTEGER,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.jctag2) then
        write(*,*) 'Fatal. Received bad message (',&
 &mstatus(MPI_TAG),') from slave ',i,'. Expected ',jctag2,'.'
        call fexit()
      end if
      jccnt(:) = jccnt(:) + imap1a(:)
      do ii=1,nseq
        jc_res(ii) = jc_res(ii) + map1(ii)
      end do
!
    end do
!   now deallocate memory
    deallocate(map2a)
    deallocate(map1)
    deallocate(imap1a)
    if (nrspecrama.gt.0) then
      deallocate(map3a)
    end if
    if (nrmolrama.gt.0) then
      deallocate(map3b)
    end if
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!   and process
    call printrama()
  end if
!
! next up: internal coordinate histograms
  if (intcalc.le.nsim) then
    do k=1,4
      if (do_ints(k).EQV..true.) then
        call System_Clock(ee1)
!       first allocate proper memory
        allocate(map2a(intszs(k,1),intszs(k,2)))
        do i=1,mpi_nodes-1
          call MPI_RECV(map2a,intszs(k,1)*intszs(k,2),MPI_RTYPE,&
 &i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
          if (mstatus(MPI_TAG).ne.inttag(k)) then
            write(*,*) 'Fatal. Received bad message (',&
 &mstatus(MPI_TAG),') from slave ',i,'. Expected ',inttag(k),'.'
            call fexit()
          end if
          if (k.eq.1) then
            bln_hist(:,:) = bln_hist(:,:) + map2a(:,:)
          else if (k.eq.2) then
            ang_hist(:,:) = ang_hist(:,:) + map2a(:,:)
          else if (k.eq.3) then
            imp_hist(:,:) = imp_hist(:,:) + map2a(:,:)
          else if (k.eq.4) then
            dih_hist(:,:) = dih_hist(:,:) + map2a(:,:)
          end if
        end do
        deallocate(map2a)
        call System_Clock(ee2)
        time_comm = time_comm + 1.0*(ee2 - ee1)
        time_analysis = time_analysis + 1.0*(ee1 - ee2)
      end if
    end do
    call prt_inthists()
  end if
!
! next we get backbone-segments
!
  if (segcalc.le.nsim) then
    call System_Clock(ee1)
!   first allocate proper memory
    allocate(map3d(NSEG,seg_alsz,maxseglen))
    allocate(imap1a(nangrps))
    allocate(map3a(nangrps,TPBINS,100))
    allocate(map3b(nangrps,2,100))
    allocate(map3c(nangrps,100,100))
    do i=1,mpi_nodes-1
      call MPI_RECV(map3d,NSEG*seg_alsz*maxseglen,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.segtag2) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',segtag2,'.'
        call fexit()
      end if
     call MPI_RECV(map3b,nangrps*200,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.segtag3) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',segtag3,'.'
        call fexit()
      end if
      call MPI_RECV(map3c,nangrps*100*100,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.segtag5) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',segtag5,'.'
        call fexit()
      end if
      call MPI_RECV(map3a,nangrps*100*100,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.segtag4) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',segtag4,'.'
        call fexit()
      end if
      call MPI_RECV(imap1a,nangrps,MPI_INTEGER,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.segtag6) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',segtag6,'.'
        call fexit()
      end if

      do j=1,nangrps
        sgprscnt(j) = sgprscnt(j) + imap1a(j)
      end do
      do j=1,NSEG
        do k=1,seg_alsz
          do l=1,maxseglen
            seg_perrs(j,k,l) = seg_perrs(j,k,l) + map3d(j,k,l)
          end do
        end do
      end do
      do j=1,nangrps
        do k=1,2
          do l=1,100
            z_hist(j,k,l) = z_hist(j,k,l) + map3b(j,k,l)
          end do
        end do
        do k=1,100
          do l=1,100
            z_hist2(j,k,l) = z_hist2(j,k,l) + map3c(j,k,l)
          end do
        end do
        do k=1,TPBINS
          do l=1,100
            z_2dhist(j,k,l) = z_2dhist(j,k,l) + map3a(j,k,l)
          end do
        end do
      end do
    end do
!   de-allocate
    deallocate(map3d)
    deallocate(map3a)
    deallocate(map3b)
    deallocate(map3c)
    deallocate(imap1a)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!   and process
    call prt_bb_segments()
  end if
!
! let's follow it up with DSSP-data
!
  if (dsspcalc.le.nsim) then
    call System_Clock(ee1)
!   first allocate proper memory
    allocate(map3b(pep_sz,11,pep_sz))
    allocate(imap1a(npepmol))
    allocate(map2a(4,100))
    allocate(map2b(100,100))
    do i=1,mpi_nodes-1
      call MPI_RECV(map3b,11*pep_sz*pep_sz,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.dssptag1) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',dssptag1,'.'
        call fexit()
      end if
      call MPI_RECV(imap1a,npepmol,MPI_INTEGER,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.dssptag4) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',dssptag4,'.'
        call fexit()
      end if
      call MPI_RECV(map2a,4*100,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.dssptag2) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',dssptag2,'.'
        call fexit()
      end if
      call MPI_RECV(map2b,100*100,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.dssptag3) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',dssptag3,'.'
        call fexit()
      end if
      do j=1,pep_sz
        do k=1,11
          do l=1,pep_sz
            dssp_perrs(j,k,l)=dssp_perrs(j,k,l) + map3b(j,k,l)
          end do
        end do
      end do
      do j=1,npepmol
        dssp_cnt(j) = dssp_cnt(j) + imap1a(j)
      end do
      do j=1,100
        do k=1,4
          dssp_hists(k,j) = dssp_hists(k,j) + map2a(k,j)
        end do
        do k=1,100
          dssp_2dhist(j,k) = dssp_2dhist(j,k) + map2b(j,k)
        end do
      end do
    end do
!   de-allocate
    deallocate(map3b)
    deallocate(imap1a)
    deallocate(map2a)
    deallocate(map2b)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!   and process
    call prt_dssp()
  end if
!
! and now contact histograms and maps as well as cluster statistics
!
  if (contactcalc.le.nsim) then
    call System_Clock(ee1)
!   first allocate proper memory for contact stuff
    allocate(map2b(nressolute,nressolute))
    allocate(imap2(nsolutes,nsolutes))
    allocate(map2a(2,20*nressolute))
    do i=1,mpi_nodes-1
      call MPI_RECV(map2b,nressolute*nressolute,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.contag1) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',contag1,'.'
        call fexit()
      end if
      call MPI_RECV(imap2,nsolutes*nsolutes,MPI_INTEGER,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.contag3) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',contag3,'.'
        call fexit()
      end if
!      write(*,*) 'SLRCOMAP: ',cmap(1,5,2),' ',cmap(2,5,2)
!      write(*,*) 'MAHCOMAP: ',contact_maps(1,5,2),' ',
! &contact_maps(2,5,2)
      do j=1,nressolute
        do k=1,nressolute
          contact_maps(j,k) = contact_maps(j,k) + map2b(j,k)
        end do
      end do
      do j=1,nsolutes
        ncontactavg(j,:) = ncontactavg(j,:) + imap2(j,:)
      end do
      call MPI_RECV(map2a,2*nressolute*20,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.contag2) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',contag2,'.'
        call fexit()
      end if
!      write(*,*) 'SLRCOMHI: ',chist(1,9),' ',chist(2,9)
!      write(*,*) 'MAHCOMHI: ',contact_hists(1,9),' ',
! &contact_hists(2,9)
      do j=1,nressolute*20
        contact_hists(1,j) = contact_hists(1,j) + map2a(1,j)
        contact_hists(2,j) = contact_hists(2,j) + map2a(2,j)
      end do
    end do
!   now de-allocate for contact-stuff
    deallocate(map2b)
    deallocate(imap2)
    deallocate(map2a)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!   and process
    call prt_rescontacts()
!
    if (clucalc*contactcalc.le.nsim) then
      call System_Clock(ee1)
!     first allocate proper memory for cluster stuff
      allocate(map1(nsolutes))
      allocate(map2a(nsolutes,nsolutes))
      allocate(map2b(nsolutes,100))
      do i=1,mpi_nodes-1
        call MPI_RECV(map1,nsolutes,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.clutag1) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',clutag1,'.'
          call fexit()
        end if
        do j=1,nsolutes
          clusys(j) = clusys(j) + map1(j)
        end do
        call MPI_RECV(map2a,nsolutes*nsolutes,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.clutag2) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',clutag2,'.'
          call fexit()
        end if
        do j=1,nsolutes
          do k=1,nsolutes
            clumol(j,k) = clumol(j,k) + map2a(j,k)
          end do
        end do
        call MPI_RECV(map2b,100*nsolutes,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.clutag3) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',clutag3,'.'
          call fexit()
        end if
        do j=1,nsolutes
          do k=1,100
            clucoo(j,k) = clucoo(j,k) + map2b(j,k)
          end do
        end do
      end do
!     de-allocate
      deallocate(map1)
      deallocate(map2a)
      deallocate(map2b)
      call System_Clock(ee2)
      time_comm = time_comm + 1.0*(ee2 - ee1)
      time_analysis = time_analysis + 1.0*(ee1 - ee2)
!     and process
      call prt_molclusters()
    end if
  end if
!
! next in the line is the PC analysis stuff
!
  if (pccalc.le.nsim) then
    call System_Clock(ee1)
!   first allocate memory
    if (gpc%nos.gt.0) then
      allocate(map2a(gpc%nos,MAXPCBINS))
      allocate(limap1(gpc%nlst))
    end if
    if (nsolutes.gt.1) then
      allocate(map3a(nangrps,nangrps,MAXPCBINS))
      allocate(limap2(nangrps,nangrps))
    end if
    if (ampc%do_bb.EQV..true.) then
      allocate(map2b(2,MAXPCBINS))
    end if
    if (ampc%do_ss.EQV..true.) then
      allocate(map2c(3,MAXPCBINS))
    end if
    if (ampc%do_bs.EQV..true.) then
      allocate(map2d(4,MAXPCBINS))
    end if
    do i=1,mpi_nodes-1
      if (ampc%do_bb.EQV..true.) then
        call MPI_RECV(map2b,2*MAXPCBINS,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.pctag1) then
          write(*,*) 'Fatal. Received bad message (',&
 &mstatus(MPI_TAG),') from slave ',i,'. Expected ',pctag1,'.'
          call fexit()
        end if
        do j=1,2
          do k=1,MAXPCBINS
            ampc%bb(j,k) = ampc%bb(j,k) + map2b(j,k)
          end do
        end do
      end if
      if (ampc%do_ss.EQV..true.) then
        call MPI_RECV(map2c,3*MAXPCBINS,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.pctag2) then
          write(*,*) 'Fatal. Received bad message (',&
 &mstatus(MPI_TAG),') from slave ',i,'. Expected ',pctag2,'.'
          call fexit()
        end if
        do j=1,3
          do k=1,MAXPCBINS
            ampc%ss(j,k) = ampc%ss(j,k) + map2c(j,k)
          end do
        end do
      end if
      if (ampc%do_bs.EQV..true.) then
        call MPI_RECV(map2d,4*MAXPCBINS,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.pctag3) then
          write(*,*) 'Fatal. Received bad message (',&
 &mstatus(MPI_TAG),') from slave ',i,'. Expected ',pctag3,'.'
          call fexit()
        end if
        do j=1,4
          do k=1,MAXPCBINS
            ampc%bs(j,k) = ampc%bs(j,k) + map2d(j,k)
          end do
        end do
      end if
      if (nsolutes.gt.1) then
        call MPI_RECV(map3a,nangrps*nangrps*MAXPCBINS,MPI_RTYPE&
 &,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.pctag5) then
          write(*,*) 'Fatal. Received bad message (',&
 &mstatus(MPI_TAG),') from slave ',i,'. Expected ',pctag5,'.'
          call fexit()
        end if
        call MPI_RECV(limap2,nangrps*nangrps,MPI_INTEGER8&
 &,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.pctag6) then
          write(*,*) 'Fatal. Received bad message (',&
 &mstatus(MPI_TAG),') from slave ',i,'. Expected ',pctag6,'.'
          call fexit()
        end if
        do j=1,nangrps
          do k=j,nangrps
            npcavg(j,k) = npcavg(j,k) + limap2(j,k)
            do l=1,MAXPCBINS
              rbc_pc(j,k,l) = rbc_pc(j,k,l) + map3a(j,k,l)
            end do
          end do
        end do
      end if
      if (gpc%nos.gt.0) then
        call MPI_RECV(map2a,gpc%nos*MAXPCBINS,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.pctag4) then
          write(*,*) 'Fatal. Received bad message (',&
 &mstatus(MPI_TAG),') from slave ',i,'. Expected ',pctag4,'.'
          call fexit()
        end if
        call MPI_RECV(limap1,gpc%nlst,MPI_INTEGER8,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.pctag7) then
          write(*,*) 'Fatal. Received bad message (',&
 &mstatus(MPI_TAG),') from slave ',i,'. Expected ',pctag7,'.'
          call fexit()
        end if
        do j=1,MAXPCBINS
          do k=1,gpc%nos
            gpc%pc(k,j) = gpc%pc(k,j) + map2a(k,j)
          end do
        end do
        do j=1,gpc%nlst
          gpc%navg(j) = gpc%navg(j) + limap1(j)
        end do
      end if
    end do
!   now de-allocate
    if (gpc%nos.gt.0) then
      deallocate(map2a)
      deallocate(limap1)
    end if
    if (nsolutes.gt.1) then
      deallocate(map3a)
      deallocate(limap2)
    end if
    if (ampc%do_bb.EQV..true.) deallocate(map2b)
    if (ampc%do_ss.EQV..true.) deallocate(map2c)
    if (ampc%do_bs.EQV..true.) deallocate(map2d)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!   and finalize
    if (do_amidepc.EQV..true.) call prt_amid_pc()
    call prt_rbc_pc()
    if (gpc%nos.gt.0) call prt_general_pc()
  end if
!
! now the persistence length (angular correlation fxn)
!
  if (polcalc.le.nsim) then
    call System_Clock(ee1)
!   first allocate memory
    allocate(imap1a(nangrps))
    allocate(map2b(nangrps,longestmol))
    allocate(map2a(nseq,9))
    do i=1,mpi_nodes-1
      call MPI_RECV(map2b,nangrps*longestmol,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.pertag1) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',pertag1,'.'
        call fexit()
      end if
      call MPI_RECV(map2a,nseq*9,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.pertag2) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',pertag2,'.'
        call fexit()
      end if
      call MPI_RECV(imap1a,nangrps,MPI_INTEGER,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.pertag3) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',pertag3,'.'
        call fexit()
      end if
!      write(*,*) 'SLRPERS: ',perv(1),' ',perv(3)
!      write(*,*) 'MAHPERS: ',persavg(1),' ',persavg(3)
      do j=1,nangrps
        npersavg(j) = npersavg(j) + imap1a(j)
        do k=1,longestmol
          persavg(j,k) = persavg(j,k) + map2b(j,k)
        end do
      end do
      do j=1,nseq
        do k=1,9
          turns_rs(j,k) = turns_rs(j,k) + map2a(j,k)
        end do
      end do
    end do
!   now de-allocate
    deallocate(imap1a)
    deallocate(map2b)
    deallocate(map2a)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!   and process
    call prt_pers()
  end if
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
! structural clustering is already preprocessed -> only communicate custom breaks
  if (use_MPIMultiSeed.EQV..false.) then
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    if (cstored.ge.0) then
      allocate(imap1a(cstored))
      allocate(imap1b(mpi_nodes*cstored))
      if (ntbrks2.gt.0) then
        imap1b(1:ntbrks2) = trbrkslst(1:ntbrks2)
      end if
      do i=1,mpi_nodes-1
        call MPI_RECV(idummr,aone,MPI_INTEGER,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.tbrktag) then
          write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from slave ',i,'. Expected ',tbrktag,'.'
          call fexit()
        end if
        if (idummr.gt.0) then
          call MPI_RECV(imap1a,cstored,MPI_INTEGER,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
          if (mstatus(MPI_TAG).ne.tbrktag2) then
            write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from slave ',i,'. Expected ',tbrktag2,'.'
            call fexit()
          end if
          imap1b((ntbrks2+1):(ntbrks2+idummr)) = imap1a(1:idummr) + i*cstored
          ntbrks2 = ntbrks2 + idummr
        end if
      end do
 
      if (allocated(trbrkslst).EQV..true.) deallocate(trbrkslst)
      if (ntbrks2.gt.0) then
        allocate(trbrkslst(ntbrks2))
        trbrkslst(1:ntbrks2) = imap1b(1:ntbrks2)
      end if
      deallocate(imap1b)
      deallocate(imap1a)
    end if
    cstored = cstored*mpi_nodes
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
    tpi = OMP_GET_THREAD_NUM() + 1
#else
    tpi = 0
#endif
    call do_clustering(tpi)
! in case of a PIGS analysis run, breaks are not supported, so we can save the time
  else
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
    cstored = cstored*mpi_nodes
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
    tpi = OMP_GET_THREAD_NUM() + 1
#else
    tpi = 0
#endif
    call do_clustering(tpi)
  end if
!
!
! here come atomic SAV's
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  if (savcalc.le.nsim) then
    call System_Clock(ee1)
!   allocate
    allocate(imap1a(nmol))
    allocate(map2a(2,n))
    if (savreq%nats.gt.0) allocate(map3a(savreq%nats,2,100))
!   remove normalization
    do j=1,nmol
      do k=atmol(j,1),atmol(j,2)
        atsavavg(:,k) = natsavavg(j)*atsavavg(:,k)
      end do
    end do
    do i=1,mpi_nodes-1
      call MPI_RECV(map2a,n*2,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.savtag) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',savtag,'.'
        call fexit()
      end if
      call MPI_RECV(imap1a,nmol,MPI_INTEGER,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.savtag2) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',savtag2,'.'
        call fexit()
      end if
      if (savreq%nats.gt.0) then
        call MPI_RECV(map3a,200*savreq%nats,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
        if (mstatus(MPI_TAG).ne.savtag3) then
          write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',savtag3,'.'
          call fexit()
        end if
      end if
!      write(*,*) 'SLRPERS: ',perv(1),' ',perv(3)
!      write(*,*) 'MAHPERS: ',persavg(1),' ',persavg(3)
      do j=1,nmol
        do k=atmol(j,1),atmol(j,2)
          atsavavg(:,k) = atsavavg(:,k) + imap1a(j)*map2a(:,k)
        end do
        natsavavg(j) = natsavavg(j) + imap1a(j)
      end do
      if (savreq%nats.gt.0) savreq%hists(:,:,:) = savreq%hists(:,:,:) + map3a(:,:,:)
    end do
!   de-allocate
    if (savreq%nats.gt.0) deallocate(map3a)
    deallocate(map2a)
    deallocate(imap1a)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
    do j=1,nmol
      do k=atmol(j,1),atmol(j,2)
        atsavavg(:,k) = atsavavg(:,k)/(1.0*natsavavg(j))
      end do
    end do
    call prt_sav()
  end if
!
! and dipole stuff ...
!
  if (dipcalc.le.nsim) then
    call System_Clock(ee1)
!   allocate
    allocate(map2a(nangrps,DIPDIM))
    allocate(map2b(nseq,DIPDIM))
    allocate(imap1a(nangrps))
    allocate(imap1b(nseq))
    do i=1,mpi_nodes-1
      call MPI_RECV(imap1a,nangrps,MPI_INTEGER,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.diptag1) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',diptag1,'.'
        call fexit()
      end if
      call MPI_RECV(map2a,nangrps*DIPDIM,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.diptag2) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',diptag2,'.'
        call fexit()
      end if
      do j=1,nangrps
        nmtdavg(j) = nmtdavg(j) + imap1a(j)
        do k=1,DIPDIM
          mtdavg(j,k) = mtdavg(j,k) + map2a(j,k)
        end do
      end do
!
      call MPI_RECV(imap1b,nseq,MPI_INTEGER,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.diptag3) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',diptag3,'.'
        call fexit()
      end if
      call MPI_RECV(map2b,nseq*DIPDIM,MPI_RTYPE,i,MPI_ANY_TAG,&
 &MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.diptag4) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',diptag4,'.'
        call fexit()
      end if
      do j=1,nseq
        nrsdavg(j) = nrsdavg(j) + imap1b(j)
        do k=1,DIPDIM
          rsdavg(j,k) = rsdavg(j,k) + map2b(j,k)
        end do
      end do
!
    end do
!   de-allocate
    deallocate(imap1a)
    deallocate(map2a)
    deallocate(imap1b)
    deallocate(map2b)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
    call prt_dipoles()
  end if
!
! diffraction maps
!
  if (diffrcalc.le.nsim) then
    call System_Clock(ee1)
!   first allocate proper memory
    allocate(map2a(rr_max,zz_max))
    do i=1,mpi_nodes-1
!     now get the map
      call MPI_RECV(map2a,rr_max*zz_max,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.diffrtag1) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',diffrtag1,'.'
        call fexit()
      end if
!     and the counter
      call MPI_RECV(idummr,aone,MPI_INTEGER,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.diffrtag2) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',diffrtag2,'.'
        call fexit()
      end if
      diffr_cnt = diffr_cnt + idummr
      do j=1,rr_max
        do k=1,zz_max
          am_diffr_map(j,k) = am_diffr_map(j,k) + map2a(j,k)
        end do
      end do
    end do
!   now deallocate memory
    deallocate(map2a)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!   and process
    call prt_diffraction()
  end if
!
! spatial density analysis
!
  if (emcalc.le.nsim) then
    call System_Clock(ee1)
!   first allocate proper memory
    allocate(map3a(emgrid%dim(1),emgrid%dim(2),emgrid%dim(3)))
    do i=1,mpi_nodes-1
!     now get the map
      call MPI_RECV(map3a,emgrid%dim(1)*emgrid%dim(2)*emgrid%dim(3),MPI_RTYPE,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.emtag1) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from slave ',i,'. Expected ',emtag1,'.'
        call fexit()
      end if
!     and the counter
      call MPI_RECV(idummr,aone,MPI_INTEGER,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.emtag2) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from slave ',i,'. Expected ',emtag2,'.'
        call fexit()
      end if
      emgrid%cnt = emgrid%cnt + idummr
      emgrid%avgmass(:,:,:) = emgrid%avgmass(:,:,:) + map3a(:,:,:)
    end do
!   now deallocate memory
    deallocate(map3a)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!   and process
    call prt_emmap(emgrid,1)
  end if
!
! particle number fluctuation analysis
!
  if (particlenumcalc.le.nsim) then
    call System_Clock(ee1)
!   first allocate proper memory
    allocate(imap2(size(numberhistogram,dim=1),size(numberhistogram,dim=2)))
    do i=1,mpi_nodes-1
!     now get the map
      call MPI_RECV(imap2,size(numberhistogram),MPI_INTEGER,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.nhtag1) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from slave ',i,'. Expected ',nhtag1,'.'
        call fexit()
      end if
      numberhistogram(:,:) = numberhistogram(:,:) + imap2(:,:)
    end do
!   now deallocate memory
    deallocate(imap2)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!   and process
    call prt_particlenumhistogram()
  end if
!
! (obsolete) LCT stuff
!
  if ((torlccalc.le.nsim).AND.(ntorlcs.gt.0)) then
    call System_Clock(ee1)
!   allocate
    allocate(map2a(ntorlcs,MAXTORLCBINS))
    do i=1,mpi_nodes-1
      call MPI_RECV(map2a,ntorlcs*MAXTORLCBINS,MPI_RTYPE,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
      if (mstatus(MPI_TAG).ne.lcttag) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),&
 &') from slave ',i,'. Expected ',lcttag,'.'
        call fexit()
      end if
!     now add those counts
      do j=1,ntorlcs
        do k=1,MAXTORLCBINS
          torlc_hist(j,k) = torlc_hist(j,k) + map2a(j,k)
        end do
      end do
    end do
!   de-allocate
    deallocate(map2a)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!   and finalize
    call prt_lc_tor()
  end if
!
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
end
!
!----------------------------------------------------------------------
!
! see above
!
subroutine MPI_AVGSend()
!
  use sequen
  use iounit
  use system
  use molecule
  use paircorr
  use polyavg
  use torsn
  use energies
  use mpi
  use mpistuff
  use contacts
  use mcsums
  use atoms
  use dipolavg
  use diffrac
  use dssps
  use fos
  use ems
  use clusters
  use grandensembles
!
  implicit none
!
  integer MPI_POLSZ
  parameter (MPI_POLSZ=12)
  integer aone,masterrank,rmgtag,rmstag,contag1,contag2,segtag1,polatag9
  integer ierr,pertag1,pctag1,pctag2,pctag3,segtag2,pertag2,savtag3
  integer savtag,segtag3,pctag4,lcttag,pertag3,polatag,polatag2
  integer rmmtag,mt,polatag3,polatag4,polatag5,polatag6,polatag10
  integer polatag7,polatag8,clutag1,clutag2,clutag3,pctag5,segtag4
  integer diptag1,diptag2,diptag3,diptag4,jctag,jctag2,diffrtag1
  integer diffrtag2,segtag5,inttag(4),dssptag1,dssptag2,dssptag3
  integer pctag6,pctag7,savtag2,contag3,segtag6,dssptag4,emtag1,emtag2,tbrktag,tbrktag2,nhtag1
  integer, ALLOCATABLE:: imap1a(:)
  RTYPE poly(nangrps,MPI_POLSZ)
  integer(KIND=8) ee1,ee2
!
#ifdef ENABLE_THREADS
!$OMP BARRIER
!$OMP MASTER
#endif
  masterrank = 0
  aone = 1
  polatag = 92
  polatag2 = 93
  polatag3 = 95
  polatag4 = 101
  polatag5 = 97
  polatag6 = 98
  polatag7 = 99
  polatag8 = 100
  polatag9 = 107
  polatag10 = 108
  rmgtag = 80
  rmstag = 81
  rmmtag = 96
  segtag1 = 79
  segtag2 = 86
  contag1 = 78
  contag2 = 77
  pertag1 = 82
  pertag2 = 87
  pertag3 = 91
  pctag1 = 83
  pctag2 = 84
  pctag3 = 85
  savtag = 88
  segtag3 = 89
  pctag4 = 90
  lcttag = 94
  clutag1 = 102
  clutag2 = 103
  clutag3 = 104
  pctag5 = 105
  segtag4 = 106
  diptag1 = 109
  diptag2 = 110
  diptag3 = 111
  diptag4 = 112
  jctag = 113
  jctag2 = 114
  diffrtag1 = 115
  diffrtag2 = 116
  segtag5 = 117
  inttag(1) = 118
  inttag(2) = 119
  inttag(3) = 120
  inttag(4) = 121
  dssptag1 = 122
  dssptag2 = 123
  dssptag3 = 124
  pctag6 = 125
  pctag7 = 126
  savtag2 = 127
  contag3 = 128
  segtag6 = 129
  dssptag4 = 130
  savtag3 = 131
  emtag1 = 132
  emtag2 = 133
  tbrktag = 134
  tbrktag2 = 135
  nhtag1 = 136
!
  if (myrank.eq.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGSend with master &
 &node.'
    call fexit()
  end if
!
! later, we will subtract the comm. time from the analysis time (accum. outside)
  call System_Clock(ee1)
!
! average polymer props
  if (polcalc.le.nsim) then
    do mt=1,nangrps
      poly(mt,1) = npolavg(mt)*rgavg(mt)
      poly(mt,2) = npolavg(mt)*rgpcavg(mt,1)
      poly(mt,3) = npolavg(mt)*rgpcavg(mt,2)
      poly(mt,4) = npolavg(mt)*rgpcavg(mt,3)
      poly(mt,5) = npolavg(mt)*asphavg(mt)
      poly(mt,6) = npolavg(mt)*acylavg(mt)
      poly(mt,7) = npolavg(mt)*rgterm1avg(mt)
      poly(mt,8) = npolavg(mt)*rgterm2avg(mt)
      poly(mt,9) = npolavg(mt)*rgterm12avg(mt)
!      poly(mt,7) = npolavg(mt)*dlstavg(mt)
!      poly(mt,8) = npolavg(mt)*dlavg(mt)
      poly(mt,10) = npolavg(mt)*reteavg(mt)
      poly(mt,11) = npolavg(mt)*rg2avg(mt)
      poly(mt,12) = npolavg(mt)*vtavg(mt)
    end do
    call MPI_SEND(poly,nangrps*MPI_POLSZ,MPI_RTYPE,&
 &masterrank,polatag,MPI_COMM_WORLD,ierr)
    call MPI_SEND(npolavg,nangrps,MPI_INTEGER,&
 &masterrank,polatag2,MPI_COMM_WORLD,ierr)
    call MPI_SEND(rdhist,nangrps*TPBINS*DLSTBINS,MPI_RTYPE,&
 &masterrank,polatag3,MPI_COMM_WORLD,ierr)
    call MPI_SEND(retehist,nangrps*maxrgbins,MPI_RTYPE,&
 &masterrank,polatag4,MPI_COMM_WORLD,ierr)
    call MPI_SEND(densproavg,nangrps*maxrgbins,MPI_RTYPE,&
 &masterrank,polatag9,MPI_COMM_WORLD,ierr)
    call MPI_SEND(rghist,nangrps*maxrgbins,MPI_RTYPE,&
 &masterrank,polatag10,MPI_COMM_WORLD,ierr)
  end if
! expensive average polymer props (internal scaling, scattering)
  if (rhcalc.le.nsim) then
    call MPI_SEND(pscavg,nangrps*nqv,MPI_RTYPE,&
 &masterrank,polatag5,MPI_COMM_WORLD,ierr)
    call MPI_SEND(nsctavg,nangrps,MPI_INTEGER,&
 &masterrank,polatag6,MPI_COMM_WORLD,ierr)
    call MPI_SEND(rhavg,nangrps*longestmol,MPI_RTYPE,&
 &masterrank,polatag7,MPI_COMM_WORLD,ierr)
    call MPI_SEND(nrhavg,nangrps,MPI_INTEGER,&
 &masterrank,polatag8,MPI_COMM_WORLD,ierr)
  end if
! global Ramachandran maps
  if (angcalc.le.nsim) then
    call MPI_SEND(genfymap,torbins*torbins,MPI_RTYPE,&
 &masterrank,rmgtag,MPI_COMM_WORLD,ierr)
!    write(*,*) 'SLSGRAM: ',genfymap(1,3)
!   specific Ramachandran maps
    if (nrspecrama.gt.0) then
      call MPI_SEND(specfymap,rama_alsz(1)*torbins*torbins,&
 &MPI_RTYPE,masterrank,rmstag,MPI_COMM_WORLD,ierr)
    end if
!    write(*,*) 'SLSSRAM: ',specfymap(1,1,3),' ',specfymap(2,1,3)
!   molecular Ramachandran maps
    if (nrmolrama.gt.0) then
      call MPI_SEND(molfymap,rama_alsz(2)*torbins*torbins,&
 &MPI_RTYPE,masterrank,rmmtag,MPI_COMM_WORLD,ierr)
    end if
    call MPI_SEND(jc_res,nseq,MPI_RTYPE,&
 &masterrank,jctag,MPI_COMM_WORLD,ierr)
    call MPI_SEND(jccnt,size(jccnt),MPI_INTEGER,&
 &masterrank,jctag2,MPI_COMM_WORLD,ierr)
  end if
! internal coordinate histograms
  if (intcalc.le.nsim) then
    if (do_ints(1).EQV..true.) then
      call MPI_SEND(bln_hist,intszs(1,1)*intszs(1,2),&
 &MPI_RTYPE,masterrank,inttag(1),MPI_COMM_WORLD,ierr)
    end if
    if (do_ints(2).EQV..true.) then
      call MPI_SEND(ang_hist,intszs(2,1)*intszs(2,2),&
 &MPI_RTYPE,masterrank,inttag(2),MPI_COMM_WORLD,ierr)
    end if
    if (do_ints(3).EQV..true.) then
      call MPI_SEND(imp_hist,intszs(3,1)*intszs(3,2),&
 &MPI_RTYPE,masterrank,inttag(3),MPI_COMM_WORLD,ierr)
    end if
    if (do_ints(4).EQV..true.) then
      call MPI_SEND(dih_hist,intszs(4,1)*intszs(4,2),&
 &MPI_RTYPE,masterrank,inttag(4),MPI_COMM_WORLD,ierr)
    end if
  end if
! backbone segments
  if (segcalc.le.nsim) then
    call MPI_SEND(seg_perrs,NSEG*seg_alsz*maxseglen,MPI_RTYPE,&
 &masterrank,segtag2,MPI_COMM_WORLD,ierr)
    call MPI_SEND(z_hist,nangrps*200,MPI_RTYPE,masterrank,&
 &segtag3,MPI_COMM_WORLD,ierr)
    call MPI_SEND(z_hist2,nangrps*100*100,MPI_RTYPE,masterrank,&
 &segtag5,MPI_COMM_WORLD,ierr)
    call MPI_SEND(z_2dhist,nangrps*TPBINS*100,MPI_RTYPE,&
 &masterrank,segtag4,MPI_COMM_WORLD,ierr)
    call MPI_SEND(sgprscnt,nangrps,MPI_INTEGER,&
 &masterrank,segtag6,MPI_COMM_WORLD,ierr)
  end if
! backbone segments
  if (dsspcalc.le.nsim) then
    call MPI_SEND(dssp_perrs,11*pep_sz*pep_sz,MPI_RTYPE,masterrank,dssptag1,MPI_COMM_WORLD,ierr)
    call MPI_SEND(dssp_cnt,npepmol,MPI_INTEGER,masterrank,dssptag4,MPI_COMM_WORLD,ierr)
    call MPI_SEND(dssp_hists,4*100,MPI_RTYPE,masterrank,dssptag2,MPI_COMM_WORLD,ierr)
    call MPI_SEND(dssp_2dhist,100*100,MPI_RTYPE,masterrank,dssptag3,MPI_COMM_WORLD,ierr)
  end if
! contact maps
  if (contactcalc.le.nsim) then
    call MPI_SEND(contact_maps,nressolute*nressolute,MPI_RTYPE,masterrank,&
 &contag1,MPI_COMM_WORLD,ierr)
    call MPI_SEND(ncontactavg,nsolutes*nsolutes,MPI_INTEGER,masterrank,&
 &contag3,MPI_COMM_WORLD,ierr)
!   contact histograms
    call MPI_SEND(contact_hists,40*nressolute,MPI_RTYPE,masterrank,&
 &contag2,MPI_COMM_WORLD,ierr)
!   cluster stats
    if (clucalc*contactcalc.le.nsim) then
      call MPI_SEND(clusys,nsolutes,MPI_RTYPE,masterrank,&
 &clutag1,MPI_COMM_WORLD,ierr)
      call MPI_SEND(clumol,nsolutes*nsolutes,MPI_RTYPE,masterrank,&
 &clutag2,MPI_COMM_WORLD,ierr)
      call MPI_SEND(clucoo,nsolutes*100,MPI_RTYPE,masterrank,&
 &clutag3,MPI_COMM_WORLD,ierr)
    end if
  end if
! pair correlation data 
  if (pccalc.le.nsim) then
    if (ampc%do_bb.EQV..true.) then
      call MPI_SEND(ampc%bb,2*MAXPCBINS,MPI_RTYPE,masterrank,&
 &pctag1,MPI_COMM_WORLD,ierr)
    end if
    if (ampc%do_ss.EQV..true.) then
      call MPI_SEND(ampc%ss,3*MAXPCBINS,MPI_RTYPE,masterrank,&
 &pctag2,MPI_COMM_WORLD,ierr)
    end if
    if (ampc%do_bs.EQV..true.) then
      call MPI_SEND(ampc%bs,4*MAXPCBINS,MPI_RTYPE,masterrank,&
 &pctag3,MPI_COMM_WORLD,ierr)
    end if
    if (nsolutes.gt.1) then
      call MPI_SEND(rbc_pc,nangrps*nangrps*MAXPCBINS,&
 &MPI_RTYPE,masterrank,pctag5,MPI_COMM_WORLD,ierr)
      call MPI_SEND(npcavg,nangrps*nangrps,&
 &MPI_INTEGER8,masterrank,pctag6,MPI_COMM_WORLD,ierr)
    end if
    if (gpc%nos.gt.0) then
      call MPI_SEND(gpc%pc,gpc%nos*MAXPCBINS,MPI_RTYPE,&
 &masterrank,pctag4,MPI_COMM_WORLD,ierr)
      call MPI_SEND(gpc%navg,gpc%nlst,MPI_INTEGER8,&
 &masterrank,pctag7,MPI_COMM_WORLD,ierr)
    end if
  end if
! bb angular correlation fxn
  if (polcalc.le.nsim) then
!    write(*,*) 'SLSPERS: ',persavg(1),' ',persavg(3)
    call MPI_SEND(persavg,nangrps*longestmol,MPI_RTYPE,masterrank,pertag1,&
 &MPI_COMM_WORLD,ierr)
    call MPI_SEND(turns_rs,nseq*9,MPI_RTYPE,masterrank,pertag2,&
 &MPI_COMM_WORLD,ierr)
    call MPI_SEND(npersavg,nangrps,MPI_INTEGER,masterrank,pertag3&
 &,MPI_COMM_WORLD,ierr)
  end if
! structural clustering is already preprocessed -> only communicate custom breaks
  if ((cstored.ge.0).AND.(use_MPIMultiSeed.EQV..false.)) then
    call MPI_SEND(ntbrks2,aone,MPI_INTEGER,masterrank,tbrktag,MPI_COMM_WORLD,ierr)
    if (ntbrks2.gt.0) then
      allocate(imap1a(cstored))
      imap1a(1:ntbrks2) = trbrkslst(1:ntbrks2)
      call MPI_SEND(imap1a,cstored,MPI_INTEGER,masterrank,tbrktag2,MPI_COMM_WORLD,ierr)
      deallocate(imap1a)
    end if
  end if
! atomic SAV's
  if (savcalc.le.nsim) then
    call MPI_SEND(atsavavg,2*n,MPI_RTYPE,masterrank,savtag,&
 &MPI_COMM_WORLD,ierr)
    call MPI_SEND(natsavavg,nmol,MPI_INTEGER,masterrank,savtag2,&
 &MPI_COMM_WORLD,ierr)
    if (savreq%nats.gt.0) then
      call MPI_SEND(savreq%hists,200*savreq%nats,MPI_RTYPE,masterrank,savtag3,&
 &MPI_COMM_WORLD,ierr)
    end if
  end if
! dipole stuff
  if (dipcalc.le.nsim) then
    call MPI_SEND(nmtdavg,nangrps,MPI_INTEGER,masterrank,diptag1&
 &,MPI_COMM_WORLD,ierr)
    call MPI_SEND(mtdavg,nangrps*DIPDIM,MPI_RTYPE,masterrank,diptag2&
 &,MPI_COMM_WORLD,ierr)
    call MPI_SEND(nrsdavg,nseq,MPI_INTEGER,masterrank,diptag3&
 &,MPI_COMM_WORLD,ierr)
    call MPI_SEND(rsdavg,nseq*DIPDIM,MPI_RTYPE,masterrank,diptag4&
 &,MPI_COMM_WORLD,ierr)
  end if
! diffraction stuff
  if (diffrcalc.le.nsim) then
    call MPI_SEND(am_diffr_map,rr_max*zz_max,MPI_RTYPE,masterrank&
 &,diffrtag1,MPI_COMM_WORLD,ierr)
    call MPI_SEND(diffr_cnt,1,MPI_INTEGER,masterrank,diffrtag2&
 &,MPI_COMM_WORLD,ierr)
  end if
! spatial density analysis
  if (emcalc.le.nsim) then
    call MPI_SEND(emgrid%avgmass,emgrid%dim(1)*emgrid%dim(2)*emgrid%dim(3),MPI_RTYPE,masterrank,emtag1,MPI_COMM_WORLD,ierr)
    call MPI_SEND(emgrid%cnt,1,MPI_INTEGER,masterrank,emtag2,MPI_COMM_WORLD,ierr)
  end if
! number histograms
  if (particlenumcalc.le.nsim) then
    call MPI_SEND(numberhistogram,size(numberhistogram),MPI_INTEGER,masterrank,nhtag1,MPI_COMM_WORLD,ierr)
  end if
! LCT analysis
  if ((torlccalc.le.nsim).AND.(ntorlcs.gt.0)) then
    call MPI_SEND(torlc_hist,ntorlcs*MAXTORLCBINS,MPI_RTYPE,&
 &masterrank,lcttag,MPI_COMM_WORLD,ierr)
  end if
!
  call System_Clock(ee2)
  time_comm = time_comm + 1.0*(ee2 - ee1)
  time_analysis = time_analysis + 1.0*(ee1 - ee2)
!
#ifdef ENABLE_THREADS
!$OMP END MASTER
!$OMP BARRIER
#endif
!
end
!
!---------------------------------------------------------------------------
!
subroutine MPI_AVGWriteEn()
!
  use iounit
  use mcsums
  use energies
  use mpi
  use mpistuff
!
  implicit none
!
  integer MPI_ENSZ
  parameter (MPI_ENSZ=MAXENERGYTERMS+1)
  integer masterrank,entag
  integer ierr,mstatus(MPI_STATUS_SIZE)
  RTYPE envals(MPI_ENSZ)
  integer(KIND=8) ee1,ee2
!
  masterrank = 0
  entag = 72
!
  if (myrank.ne.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGWriteEn with slave node (',myrank,').'
    call fexit()
  end if
!
 14   format(i12,1x,22(g14.7,1x))
  if (mpi_cnt_en.eq.0) then
    write(iene,14) nstep,esave,esterms(1),esterms(2),esterms(3),&
 &                 esterms(4),esterms(5),esterms(6),esterms(7),&
 &                 esterms(8),esterms(9),esterms(10),esterms(11),&
 &                 esterms(12),esterms(13),esterms(14),esterms(15),&
 &                 esterms(16),esterms(17),esterms(18),esterms(19),&
 &                 esterms(20),esterms(21)
  else
    call System_Clock(ee1)
    call MPI_RECV(envals,MPI_ENSZ,MPI_RTYPE,mpi_cnt_en,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.entag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from slave ',mpi_cnt_en,'. Expected ',entag,'.'
      call fexit()
    end if
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
    write(iene,14) nstep,envals(1),envals(2),envals(3),envals(4),&
 &                    envals(5),envals(6),envals(7),envals(8),&
 &                    envals(9),envals(10),envals(11),envals(12),&
 &                    envals(13),envals(14),envals(15),envals(16),&
 &                    envals(17),envals(18),envals(19),envals(20),&
 &                    envals(21),envals(22)
  end if
!
end
!
!---------------------------------------------------------------------------
!
subroutine MPI_AVGSendEn()
!
  use iounit
  use mcsums
  use energies
  use mpi
  use mpistuff
  use mcsums
!
  implicit none
!
  integer MPI_ENSZ
  parameter (MPI_ENSZ=MAXENERGYTERMS+1)
  integer masterrank,entag,i
  integer ierr
  RTYPE envals(MPI_ENSZ)
  integer(KIND=8) ee1,ee2
!
  masterrank = 0
  entag = 72
!
  if (myrank.eq.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGSendEn with master node.'
    call fexit()
  end if
!
  if (mpi_cnt_en.eq.myrank) then
!
    envals(1) = esave
    do i=1,MAXENERGYTERMS
      envals(i+1) = esterms(i)
    end do
!
    call System_Clock(ee1)
    call MPI_SEND(envals,MPI_ENSZ,MPI_RTYPE,masterrank,entag,MPI_COMM_WORLD,ierr)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!
  end if
!
end
!
!---------------------------------------------------------------------------
!
subroutine MPI_AVGWriteAcc()
!
  use iounit
  use mcsums
  use accept
  use mpi
  use mpistuff
  use mcsums
!
  implicit none
!
  integer MPI_ACCSZ
  parameter (MPI_ACCSZ=20)
  integer masterrank,acctag,i,j
  integer ierr,mstatus(MPI_STATUS_SIZE)
  integer accvals(MPI_ACCSZ),cumacv(MPI_ACCSZ)
  integer(KIND=8) ee1,ee2
!
  acctag = 73
  masterrank = 0
!
  if (myrank.ne.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGWriteAcc with slave node (',myrank,').'
    call fexit()
  end if
!
  cumacv(1)  = acc%naccept
  cumacv(2)  = acc%nfy
  cumacv(3)  = acc%nomega
  cumacv(4)  = acc%ndjcr
  cumacv(5)  = acc%ndocr
  cumacv(6)  = acc%nsjcr
  cumacv(7)  = acc%nujcr
  cumacv(8)  = acc%nchi
  cumacv(9)  = acc%nnuc
  cumacv(10) = acc%npucker
  cumacv(11) = acc%nnuccr
  cumacv(12) = acc%nrb
  cumacv(13) = acc%nrot
  cumacv(14) = acc%ntrans
  cumacv(15) = acc%nclurb
  cumacv(16) = acc%ninsert
  cumacv(17) = acc%ndelete
  cumacv(18) = acc%nidentitychange
  cumacv(19) = acc%nre
  cumacv(20) = acc%nother
!
  call System_Clock(ee1)
!
 24   format(i12,1x,20i11)
  do i=1,mpi_nodes-1
    call MPI_RECV(accvals,MPI_ACCSZ,MPI_INTEGER,i,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.acctag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from slave ',i,'. Expected ',acctag,'.'
      call fexit()
    end if
    do j=1,MPI_ACCSZ
      cumacv(j) = cumacv(j) + accvals(j)
    end do
  end do
!
  call System_Clock(ee2)
  time_comm = time_comm + 1.0*(ee2 - ee1)
  time_analysis = time_analysis + 1.0*(ee1 - ee2)
!
  write(iacc,24) nstep,cumacv(1),cumacv(2),cumacv(3),cumacv(4),&
 &cumacv(5),cumacv(6),cumacv(7),cumacv(8),cumacv(9),cumacv(10),&
 &cumacv(11),cumacv(12),cumacv(13),cumacv(14),cumacv(15),cumacv(16),&
 &cumacv(17),cumacv(18),cumacv(19),cumacv(20)
!
end
!
!---------------------------------------------------------------------------
!
subroutine MPI_AVGSendAcc()
!
  use iounit
  use mcsums
  use accept
  use mpi
  use mpistuff
  use mcsums
!
  implicit none
!
  integer MPI_ACCSZ
  parameter (MPI_ACCSZ=20)
  integer masterrank,acctag,ierr
  integer accvals(MPI_ACCSZ)
  integer(KIND=8) ee1,ee2
!
  masterrank = 0
  acctag = 73
!
  if (myrank.eq.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGSendAcc with master node.'
    call fexit()
  end if
!
  accvals(1)  = acc%naccept
  accvals(2)  = acc%nfy
  accvals(3)  = acc%nomega
  accvals(4)  = acc%ndjcr
  accvals(5)  = acc%ndocr
  accvals(6)  = acc%nsjcr
  accvals(7)  = acc%nujcr
  accvals(8)  = acc%nchi
  accvals(9)  = acc%nnuc
  accvals(10) = acc%npucker
  accvals(11) = acc%nnuccr
  accvals(12) = acc%nrb
  accvals(13) = acc%nrot
  accvals(14) = acc%ntrans
  accvals(15) = acc%nclurb
  accvals(16) = acc%ninsert
  accvals(17) = acc%ndelete
  accvals(18) = acc%nidentitychange
  accvals(19) = acc%nre
  accvals(20) = acc%nother
!
  call System_Clock(ee1)
  call MPI_SEND(accvals,MPI_ACCSZ,MPI_INTEGER,masterrank,acctag,MPI_COMM_WORLD,ierr)
  call System_Clock(ee2)
  time_comm = time_comm + 1.0*(ee2 - ee1)
  time_analysis = time_analysis + 1.0*(ee1 - ee2)
!
end
!
!---------------------------------------------------------------------------
!
subroutine MPI_AVGWriteEns()
!
  use iounit
  use mcsums
  use system
  use mpi
  use mpistuff
  use forces
!
  implicit none
!
  integer MPI_ENSSZ
  parameter (MPI_ENSSZ=11)
  integer masterrank,enstag
  integer ierr,mstatus(MPI_STATUS_SIZE)
  RTYPE envals(MPI_ENSSZ)
  integer(KIND=8) ee1,ee2
!
  masterrank = 0
  enstag = 76
!
  if (myrank.ne.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGWriteEn with slave node (',myrank,').'
    call fexit()
  end if
!
 14   format(i12,1x,15(g14.7,1x))
  if (mpi_cnt_ens.eq.0) then
    if ((ens%flag.eq.1).OR.(ens%flag.eq.2)) then
      if (fycxyz.eq.1) then
        write(iens,14) nstep,ens%insU,ens%insK,ens%insK+ens%insU,ens%insT,ens%insR(1),ens%insK2,ens%insK2+ens%insU
      else
        write(iens,14) nstep,ens%insU,ens%insK,ens%insK+ens%insU,ens%insT,ens%insR(1)
      end if
    else if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      if (fycxyz.eq.1) then
        write(iens,14) nstep,ens%insU,ens%insK,ens%insK+ens%insU,ens%insT,ens%insR(1),ens%insK2,ens%insK2+ens%insU,&
 &ens%insK+ens%insU+bnd_pV,ens%insP,ens%insV,bnd_pV
      else
        write(iens,14) nstep,ens%insU,ens%insK,ens%insK+ens%insU,ens%insT,ens%insR(1),&
 &ens%insK+ens%insU+bnd_pV,ens%insP,ens%insV,bnd_pV
      end if
    end if
  else
    call System_Clock(ee1)
    call MPI_RECV(envals,MPI_ENSSZ,MPI_RTYPE,mpi_cnt_ens,&
 &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.enstag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from slave ',mpi_cnt_ens,'. Expected ',enstag,'.'
      call fexit()
    end if
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
    if ((ens%flag.eq.1).OR.(ens%flag.eq.2)) then
      if (fycxyz.eq.1) then
        write(iens,14) nstep,envals(1),envals(2),envals(3),envals(5),envals(9),envals(10),envals(11)
      else
        write(iens,14) nstep,envals(1),envals(2),envals(3),envals(5),envals(9)
      end if
    else if ((ens%flag.eq.3).OR.(ens%flag.eq.4)) then
      if (fycxyz.eq.1) then
        write(iens,14) nstep,envals(1),envals(2),envals(3),envals(5),envals(9),envals(10),envals(11),&
 &envals(4),envals(6),envals(7),envals(8)
      else
        write(iens,14) nstep,envals(1),envals(2),envals(3),envals(5),envals(9),envals(4),envals(6),envals(7),envals(8)
      end if
    end if
  end if
!
end
!
!----------------------------------------------------------------------------
!
subroutine MPI_AVGSendEns()
!
  use iounit
  use mcsums
  use mpi
  use mpistuff
  use mcsums
  use system
  use forces
!
  implicit none
!
  integer MPI_ENSSZ
  parameter (MPI_ENSSZ=11)
  integer masterrank,enstag
  integer ierr
  RTYPE envals(MPI_ENSSZ)
  integer(KIND=8) ee1,ee2
!
  masterrank = 0
  enstag = 76
!
  if (myrank.eq.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGSendEns with master node.'
    call fexit()
  end if
!
  if (mpi_cnt_ens.eq.myrank) then
!
    envals(1) = ens%insU
    envals(2) = ens%insK
    envals(3) = ens%insK+ens%insU
    envals(4) = ens%insK+ens%insU+bnd_pV
    envals(5) = ens%insT
    envals(6) = ens%insP
    envals(7) = ens%insV
    envals(8) = bnd_pV
    envals(9) = ens%insR(1)
    envals(10) = ens%insK2
    envals(11) = ens%insK2+ens%insU
!
    call System_Clock(ee1)
    call MPI_SEND(envals,MPI_ENSSZ,MPI_RTYPE,masterrank,enstag,MPI_COMM_WORLD,ierr)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!
  end if
!
end
!
!---------------------------------------------------------------------------
!
subroutine MPI_AVGWriteTor()
!
  use sequen
  use mpi
  use mpistuff
  use mcsums
  use torsn
  use molecule
  use grandensembles
  use system
  use zmatrix
  use atoms
  use polypep
  use fyoc, ONLY: disulf
!
  implicit none
!
  integer masterrank,tortag,j,ttc,imol,shf,rs
  integer ierr,mstatus(MPI_STATUS_SIZE)
  RTYPE offset,getztor
  logical ismember
  integer(KIND=8) ee1,ee2
!
  tortag = 71
  masterrank = 0
  shf = 0
  if (ua_model.gt.0) shf = 1
!
  if (myrank.ne.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGWriteTor with slave node (',myrank,').'
    call fexit()
  end if
!
!
 10   format(i12,1x,20000000(f9.3,1x))
!
  if (mpi_cnt_tor.eq.0) then
    imol = 0
    do ttc=1,nzmlst
      if (molofrs(atmres(abs(torzmlst(ttc)))).ne.imol) then
        imol = molofrs(atmres(abs(torzmlst(ttc))))
        offset = 0.0
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
            if (ismember(ispresent,imol).EQV..false.) offset = 720.0
          end if
        end if
      end if
      if (torzmlst(ttc).gt.0) then
        curtvec(ttc) = ztor(torzmlst(ttc)) + offset
      else
        rs = atmres(abs(torzmlst(ttc)))
        if (disulf(rs).gt.0) then
          if (crosslink(crlk_idx(rs))%itstype.le.2) then
            if (abs(torzmlst(ttc)).eq.at(rs)%sc(2-shf)) then
              curtvec(ttc) = getztor(at(rs)%sc(2-shf),at(rs)%sc(3-shf),at(disulf(rs))%sc(3-shf),at(disulf(rs))%sc(2-shf))
            else if (abs(torzmlst(ttc)).eq.at(rs)%sc(3-shf)) then
              curtvec(ttc) = getztor(cai(disulf(rs)),at(disulf(rs))%sc(2-shf),at(disulf(rs))%sc(3-shf),at(rs)%sc(3-shf))
            else
              curtvec(ttc) = bang(abs(torzmlst(ttc))) + offset
            end if
          else
            curtvec(ttc) = bang(abs(torzmlst(ttc))) + offset
          end if
        else
          curtvec(ttc) = bang(abs(torzmlst(ttc))) + offset
        end if
      end if
    end do
!
    write(idihed,10) nstep,(curtvec(j),j=1,nzmlst)
  else
    call System_Clock(ee1)
    call MPI_RECV(curtvec(1:nzmlst),nzmlst,MPI_RTYPE,mpi_cnt_tor,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.tortag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from slave ',mpi_cnt_tor,'. Expected ',tortag,'.'
      call fexit()
    end if
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
    write(idihed,10) nstep,(curtvec(j),j=1,nzmlst)
  end if
!
end
!
!---------------------------------------------------------------------------
!
! MPI communication: ntorpuck needs to be a constant
!
subroutine MPI_AVGSendTor()
!
  use sequen
  use mpi
  use mpistuff
  use mcsums
  use torsn
  use molecule
  use grandensembles
  use system
  use zmatrix
  use atoms
  use polypep
  use fyoc, ONLY: disulf
!
  implicit none
!
  integer masterrank,tortag,ttc,imol,rs,shf
  integer ierr
  RTYPE offset,getztor
  logical ismember
  integer(KIND=8) ee1,ee2
!
  masterrank = 0
  tortag = 71
  shf = 0
  if (ua_model.gt.0) shf = 1
!
  if (myrank.eq.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGSendTor with master node.'
    call fexit()
  end if
!
  if (mpi_cnt_tor.eq.myrank) then
!
    imol = 0
    do ttc=1,nzmlst
      if (molofrs(atmres(abs(torzmlst(ttc)))).ne.imol) then
        imol = molofrs(atmres(abs(torzmlst(ttc))))
        offset = 0.0
        if ((ens%flag.eq.5).OR.(ens%flag.eq.6)) then
          if (ismember(fluctypes,moltypid(imol)).EQV..true.) then
            if (ismember(ispresent,imol).EQV..false.) offset = 720.0
          end if
        end if
      end if
      if (torzmlst(ttc).gt.0) then
        curtvec(ttc) = ztor(torzmlst(ttc)) + offset
      else
        rs = atmres(abs(torzmlst(ttc)))
        if (disulf(rs).gt.0) then
          if (crosslink(crlk_idx(rs))%itstype.le.2) then
            if (abs(torzmlst(ttc)).eq.at(rs)%sc(2-shf)) then
              curtvec(ttc) = getztor(at(rs)%sc(2-shf),at(rs)%sc(3-shf),at(disulf(rs))%sc(3-shf),at(disulf(rs))%sc(2-shf))
            else if (abs(torzmlst(ttc)).eq.at(rs)%sc(3-shf)) then
              curtvec(ttc) = getztor(cai(disulf(rs)),at(disulf(rs))%sc(2-shf),at(disulf(rs))%sc(3-shf),at(rs)%sc(3-shf))
            else
              curtvec(ttc) = bang(abs(torzmlst(ttc))) + offset
            end if
          else
            curtvec(ttc) = bang(abs(torzmlst(ttc))) + offset
          end if
        else
          curtvec(ttc) = bang(abs(torzmlst(ttc))) + offset
        end if
      end if
    end do
!
    call System_Clock(ee1)
    call MPI_SEND(curtvec(1:nzmlst),nzmlst,MPI_RTYPE,masterrank,tortag,MPI_COMM_WORLD,ierr)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!
  end if
!
end
!
!---------------------------------------------------------------------------
!
subroutine MPI_AVGWritePol()
!
  use iounit
  use mcsums
  use mpi
  use mpistuff
  use mcsums
!
  implicit none
!
  integer MPI_POLSZ
  parameter (MPI_POLSZ=8)
  integer masterrank,poltag
  integer ierr,mstatus(MPI_STATUS_SIZE)
  RTYPE polvals(MPI_POLSZ)
  RTYPE rgs,rgsten(3,3),rgsev(3),tps
  RTYPE asphs,acyls,dlts,dltsts
  integer(KIND=8) ee1,ee2
!
 10   format(i12,2x,g14.7,1x,f12.5,1x,f12.5,1x,f12.5,1x,f12.5,1x,g14.7,1x,g14.7,1x,g14.7)
!
  masterrank = 0
  poltag = 74
!
  if (myrank.ne.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGWritePol with slave node.'
    call fexit()
  end if
!
  if (mpi_cnt_pol.eq.0) then
!
    call gyrate2(rgs,tps,rgsten,rgsev,asphs,acyls,dlts,dltsts)
!
    write(ipolmr,10) nstep,rgs,tps,asphs/(rgs**2),acyls/(rgs**2),dltsts,rgsev(3),rgsev(2),rgsev(1)
!
  else
!
    call System_Clock(ee1)
    call MPI_RECV(polvals,MPI_POLSZ,MPI_RTYPE,mpi_cnt_pol,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.poltag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from slave ',mpi_cnt_pol,'. Expected ',poltag,'.'
      call fexit()
    end if
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
    write(ipolmr,10) nstep,polvals(1),polvals(2),polvals(3),polvals(4),polvals(5),polvals(6),polvals(7),polvals(8)
!
  end if
!
end
!
!---------------------------------------------------------------------------
!
subroutine MPI_AVGSendPol()
!
  use iounit
  use mpi
  use mpistuff
  use mcsums
!
  implicit none
!
  integer MPI_POLSZ
  parameter (MPI_POLSZ=8)
  integer masterrank,poltag
  integer ierr
  RTYPE polvals(MPI_POLSZ)
  RTYPE rgs,rgsten(3,3),rgsev(3),tps
  RTYPE asphs,acyls,dlts,dltsts
  integer(KIND=8) ee1,ee2
!
 10   format(i12,2x,g14.7,1x,f12.5,1x,f12.5,1x,f12.5,1x,f12.5,1x,g14.7,1x,g14.7,1x,g14.7)
!
  masterrank = 0
  poltag = 74
!
  if (myrank.eq.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGSendPol with master node.'
    call fexit()
  end if
!
  if (myrank.eq.mpi_cnt_pol) then
!
    call gyrate2(rgs,tps,rgsten,rgsev,asphs,acyls,dlts,dltsts)
!
    polvals(1) = rgs
    polvals(2) = tps
    polvals(3) = asphs/(rgs**2)
    polvals(4) = acyls/(rgs**2)
    polvals(5) = dltsts
    polvals(6) = rgsev(3)
    polvals(7) = rgsev(2)
    polvals(8) = rgsev(1)
!
    call System_Clock(ee1)
    call MPI_SEND(polvals,MPI_POLSZ,MPI_RTYPE,masterrank,poltag,MPI_COMM_WORLD,ierr)
    call System_Clock(ee2) 
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!
  end if
!
end
!
!-------------------------------------------------------------------------------------------
!
subroutine MPI_AVGWriteSAV()
!
  use iounit
  use mcsums
  use energies
  use mpi
  use mpistuff
  use atoms
  use fos
!
  implicit none
!
  integer masterrank,savtag,i
  integer ierr,mstatus(MPI_STATUS_SIZE)
  RTYPE savtmp(2+savreq%nats),cumsav,instsav(mpi_nodes),ainss(1+savreq%nats)
  RTYPE tsv,dummy
  logical sayno
  integer(KIND=8) ee1,ee2
!
  sayno = .false.
  savtag = 75
  masterrank = 0
!
  if (myrank.ne.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGWriteSAV with slave node (',myrank,').'
    call fexit()
  end if
!
  cumsav = savavg
!
  call System_Clock(ee1)
!
 24   format(i12,1x,5i11)
  do i=1,mpi_nodes-1
    call MPI_RECV(savtmp,2+savreq%nats,MPI_RTYPE,i,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.savtag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from slave ',i,'. Expected ',savtag,'.'
      call fexit()
    end if
    cumsav = cumsav + savtmp(1)
    instsav(i) = savtmp(2)
    if ((i.eq.mpi_cnt_sav).AND.(savreq%nats.gt.0)) then
      ainss(1:savreq%nats) = savtmp(3:(2+savreq%nats))
    end if
  end do
!
  call System_Clock(ee2)
  time_comm = time_comm + 1.0*(ee2 - ee1)
  time_analysis = time_analysis + 1.0*(ee1 - ee2)
!
 47   format(i10,10000000(1x,g15.7))
  if (mpi_cnt_sav.eq.0) then
!
    tsv = 0.0
    call get_molsav(tsv,dummy,sayno)
    if (savreq%nats.gt.0) then
      write(isasa,47) nstep,cumsav/mpi_nodes,tsv,atsav(savreq%idx(1:savreq%nats))/atbvol(savreq%idx(1:savreq%nats))
    else
      write(isasa,47) nstep,cumsav/mpi_nodes,tsv
    end if
!
  else
!
    if (savreq%nats.gt.0) then
      write(isasa,47) nstep,cumsav/mpi_nodes,instsav(mpi_cnt_sav),ainss(1:savreq%nats)
    else
      write(isasa,47) nstep,cumsav/mpi_nodes,instsav(mpi_cnt_sav)
    end if
!
  end if
!
end
!
!---------------------------------------------------------------------------
!
subroutine MPI_AVGSendSAV()
!
  use iounit
  use energies
  use mpi
  use mpistuff
  use mcsums
  use fos
  use atoms
!
  implicit none
!
  integer masterrank,savtag
  integer ierr
  RTYPE savtmp(2+savreq%nats),tsv,dummy
  logical sayno
  integer(KIND=8) ee1,ee2
!
  sayno = .false.
  masterrank = 0
  savtag = 75
!
  if (myrank.eq.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGSendSAV with master node.'
    call fexit()
  end if
!
  tsv = 0.0
  savtmp(1) = savavg
  call get_molsav(tsv,dummy,sayno)
  savtmp(2) = tsv
  if (savreq%nats.gt.0) then
    savtmp(3:(savreq%nats+2)) = atsav(savreq%idx(1:savreq%nats))/atbvol(savreq%idx(1:savreq%nats))
  end if
!
  call System_Clock(ee1)
  call MPI_SEND(savtmp,2+savreq%nats,MPI_RTYPE,masterrank,savtag,MPI_COMM_WORLD,ierr)
  call System_Clock(ee2)
  time_comm = time_comm + 1.0*(ee2 - ee1)
  time_analysis = time_analysis + 1.0*(ee1 - ee2)
!
end
!
!-------------------------------------------------------------------------------------------
!
subroutine MPI_AVGWaitStr()
!
  use mpi
  use mpistuff
  use mcsums
!
  implicit none
!
  integer masterrank,strtag
  integer ierr,dummy,aone
  integer mstatus(MPI_STATUS_SIZE)
!
  if (mpi_cnt_xyz.eq.myrank) return
!
  masterrank = 0
  aone = 1
  dummy = 0
  strtag = 3753753 + mpi_cnt_xyz
!
  call MPI_RECV(dummy,aone,MPI_INTEGER,mpi_cnt_xyz,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
  call MPI_SEND(dummy,aone,MPI_INTEGER,mpi_cnt_xyz,strtag,MPI_COMM_WORLD,ierr)
!
  if (mstatus(MPI_TAG).ne.strtag) then
    write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from slave ',mpi_cnt_xyz,'. Expected ',strtag,'.'
    call fexit()
  end if
!
end
!
!-----------------------------------------------------------------------------------------
!
subroutine MPI_AVGPingStr()
!
  use mpi
  use mpistuff
!
  implicit none
!
  integer masterrank,strtag
  integer ierr,dummy,aone
  integer mstatus(MPI_STATUS_SIZE)
!
  if (mpi_cnt_xyz.ne.myrank) return
!
  masterrank = 0
  aone = 1
  dummy = 0
  strtag = 3753753 + mpi_cnt_xyz
!
  call MPI_SEND(dummy,aone,MPI_INTEGER,masterrank,strtag,MPI_COMM_WORLD,ierr)
  call MPI_RECV(dummy,aone,MPI_INTEGER,masterrank,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
!
  if (mstatus(MPI_TAG).ne.strtag) then
    write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from master. Expected ',strtag,'.'
    call fexit()
  end if
!
end
!
!----------------------------------------------------------------------------------------
!
! MPI communication: ntorpuck needs to be a constant
!
subroutine MPI_AVGWriteTRCV(mt,dc,ttt)
!
  use iounit
  use torsn
  use math
  use mpi
  use mpistuff
  use mcsums
!
  implicit none
!
  integer masterrank,trcvtag,i,dc,mt
  integer ierr,mstatus(MPI_STATUS_SIZE)
  RTYPE ttt(ntorsn)
  integer(KIND=8) ee1,ee2
!
  trcvtag = 70
  masterrank = 0
!
  if (myrank.ne.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGWriteTRCV with slave node (',myrank,').'
    call fexit()
  end if
!
  if (mpi_cnt_trcv.eq.0) then
!   use the values passed through argument
  else
    call System_Clock(ee1)
    call MPI_RECV(ttt,dc,MPI_RTYPE,mpi_cnt_trcv,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
    if (mstatus(MPI_TAG).ne.trcvtag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from slave ',mpi_cnt_trcv,'. Expected ',trcvtag,'.'
      call fexit()
    end if
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
  end if
!
 86   format(2000000(g11.4,1x))
!
  if (covmode.eq.1) then
    write(itrcv(mt),86) (ttt(i)/RADIAN,i=1,dc)
  else if (covmode.eq.2) then
    write(itrcv(mt),86) (cos(ttt(i)/RADIAN),sin(ttt(i)/RADIAN),i=1,dc)
  else if (covmode.eq.3) then
    write(itrcv(mt),86) (ttt(i),i=1,dc)
  end if
!
end
!
!---------------------------------------------------------------------------
!
! MPI communication: ntorsn needs to be a constant
!
subroutine MPI_AVGSendTRCV(mt,dc,ttt)
!
  use iounit
  use torsn
  use mpi
  use mpistuff
  use mcsums
!
  implicit none
!
  integer masterrank,trcvtag,dc,mt
  integer ierr
  RTYPE ttt(ntorsn)
  integer(KIND=8) ee1,ee2
!
  masterrank = 0
  trcvtag = 70
!
  if (myrank.eq.masterrank) then
    write(*,*) 'Fatal. Called MPI_AVGSendTRCV with master node.'
    call fexit()
  end if
!
  if (mpi_cnt_trcv.eq.myrank) then
!
    call System_Clock(ee1)
    call MPI_SEND(ttt,dc,MPI_RTYPE,masterrank,trcvtag,MPI_COMM_WORLD,ierr)
    call System_Clock(ee2)
    time_comm = time_comm + 1.0*(ee2 - ee1)
    time_analysis = time_analysis + 1.0*(ee1 - ee2)
!
  end if
!
end
!
!-----------------------------------------------------------------------------------------
!
#ifdef LINK_NETCDF
!
subroutine MPI_AVGSyncNetCDF()
!
  use mpi
  use mpistuff
  use pdb
!
  implicit none
!
  integer masterrank,cdftag,i
  integer ierr,MPI_CDFSZ
  integer mstatus(MPI_STATUS_SIZE)
!
  masterrank = 0
  MPI_CDFSZ = 40
  cdftag = 376
!
  if (myrank.eq.masterrank) then
!
!   broadcast variable ID's
    do i=1,mpi_nodes-1
      call MPI_SEND(netcdf_ids,MPI_CDFSZ,MPI_INTEGER,i,cdftag,MPI_COMM_WORLD,ierr)
    end do
!
  else
    call MPI_RECV(netcdf_ids,MPI_CDFSZ,MPI_INTEGER,masterrank,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr)
!
    if (mstatus(MPI_TAG).ne.cdftag) then
      write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from master. Expected ',cdftag,'.'
      call fexit()
    end if
!
  end if
!
  netcdf_ids(40) = -1
!
end
!
#endif
!
!---------------------------------------------------------------------------
!
! the following are some replacement global MPI communication routines
! bc allows choosing between functionality of MPI_ALLREDUCE and MPI_REDUCE
! msgtag are preferably 4 unique values
! opcode is currently 1 for summation, 2 for max, 3 for min, and 0 for simple broadcast
! ierr holds the last non-zero error code (assuming routine exits cleanly)
! MPIcolls allows switching to MPI_(ALL)REDUCE and MPI_BCAST, note, however, that the root
! for BCAST and REDUCE cannot be altered (always masterrank = 0)
!
subroutine MPI_ALLINTS1D(msgtag,opcode,ierr,bc,msgsz2,intmsg1d,MPIcolls)
!
  use mpi
  use mpistuff
!
  implicit none
!
  integer opcode,ierr,ierr2,mstatus(MPI_STATUS_SIZE),msgsz,msgsz2,msgtag(4),i,masterrank
  integer intmsg1d(msgsz2)
  logical bc,MPIcolls
!
  ierr = 0
  masterrank = 0
!
  if ((opcode.eq.0).AND.(bc.EQV..false.)) return
!
  if (mod(msgsz2,2).ne.0) then
    write(*,*) 'Fatal. In MPI_ALLINTS1D, messages must be of even rank in the last dimension.'
    call fexit()
  end if
  msgsz = msgsz2/2
  intmsg1d((msgsz+1):(2*msgsz)) = intmsg1d(1:msgsz)
! block for using standard MPI collective communication syntax
  if (MPIcolls.EQV..true.) then
    if (opcode.eq.0) then
      call MPI_BCAST(intmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_INTEGER,masterrank,MPI_COMM_WORLD,ierr2)
    else if (bc.EQV..true.) then
      if (opcode.eq.1) then
        call MPI_ALLREDUCE(intmsg1d(1:msgsz),intmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.2) then
        call MPI_ALLREDUCE(intmsg1d(1:msgsz),intmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.3) then
        call MPI_ALLREDUCE(intmsg1d(1:msgsz),intmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr2)
      end if
    else
      if (opcode.eq.1) then
        call MPI_REDUCE(intmsg1d(1:msgsz),intmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_INTEGER,MPI_SUM,masterrank,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.2) then
        call MPI_REDUCE(intmsg1d(1:msgsz),intmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_INTEGER,MPI_MAX,masterrank,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.3) then
        call MPI_REDUCE(intmsg1d(1:msgsz),intmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_INTEGER,MPI_MIN,masterrank,MPI_COMM_WORLD,ierr2)
      end if
    end if
    if (ierr2.ne.0) ierr = ierr2
    return
  end if
  if (opcode.gt.0) then
    do i=1,mpi_granullst%n_r
      call MPI_RECV(intmsg1d(1:msgsz),msgsz,MPI_INTEGER,mpi_granullst%recvs(i),&
  &                 MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(1)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from granule node ',mpi_granullst%recvs(i)+1,'. &
  &Expected ',msgtag(1),'.'
        call fexit()
      end if
      if (opcode.eq.1) then ! sum
        intmsg1d((msgsz+1):(2*msgsz)) = intmsg1d((msgsz+1):(2*msgsz)) + intmsg1d(1:msgsz)
      else if (opcode.eq.2) then ! max
        intmsg1d((msgsz+1):(2*msgsz)) = max(intmsg1d((msgsz+1):(2*msgsz)),intmsg1d(1:msgsz)) ! op.s on vector
      else if (opcode.eq.3) then ! min
        intmsg1d((msgsz+1):(2*msgsz)) = min(intmsg1d((msgsz+1):(2*msgsz)),intmsg1d(1:msgsz)) ! op.s on vector
      end if
!     the original message is now lost in favor of the current composite message (if granule%n_r > 0)
    end do
    do i=1,mpi_granullst%n_s 
      call MPI_SEND(intmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_INTEGER,mpi_granullst%sends(i),&
  &                  msgtag(1),MPI_COMM_WORLD,ierr2) ! this sends either a local composite, or the untouched individual values
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   at this point the node with the connection leaving the granule should have the correct info for the granule
    do i=1,mpi_globlst%n_r
      call MPI_RECV(intmsg1d(1:msgsz),msgsz,MPI_INTEGER,mpi_globlst%recvs(i),&
  &                 MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(2)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from external node ',mpi_globlst%recvs(i)+1,'. &
  &Expected ',msgtag(2),'.'
        call fexit()
      end if
      if (opcode.eq.1) then ! sum
        intmsg1d((msgsz+1):(2*msgsz)) = intmsg1d((msgsz+1):(2*msgsz)) + intmsg1d(1:msgsz)
      else if (opcode.eq.2) then ! max
        intmsg1d((msgsz+1):(2*msgsz)) = max(intmsg1d((msgsz+1):(2*msgsz)),intmsg1d(1:msgsz))
      else if (opcode.eq.3) then ! min
        intmsg1d((msgsz+1):(2*msgsz)) = min(intmsg1d((msgsz+1):(2*msgsz)),intmsg1d(1:msgsz))
      end if
!     the granule composite is now lost in favor of the current, global composite message (if glob%n_r > 0)
    end do
    do i=1,mpi_globlst%n_s
      call MPI_SEND(intmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_INTEGER,mpi_globlst%sends(i),&
  &                  msgtag(2),MPI_COMM_WORLD,ierr2)! this sends either an intermediate composite, or the untouched granule vals
      if (ierr2.ne.0) ierr = ierr2 
    end do
  end if
! at this point, there should be exactly one node with the complete correct info in the second half of the array
  if (bc.EQV..true.) then
    do i=1,mpi_globlst%n_r2
      call MPI_RECV(intmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_INTEGER,mpi_globlst%recvs2(i),&
 &                  MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(3)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from external node ',mpi_globlst%recvs2(i)+1,'. &
 &Expected ',msgtag(3),'.'
        call fexit()
      end if
    end do
    do i=1,mpi_globlst%n_s2
      call MPI_SEND(intmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_INTEGER,mpi_globlst%sends2(i),&
 &                   msgtag(3),MPI_COMM_WORLD,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   at this point, all those nodes with mpi_granullst%n_r2 == 0 should have complete correct info
    do i=1,mpi_granullst%n_r2
      call MPI_RECV(intmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_INTEGER,mpi_granullst%recvs2(i),&
 &                  MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(4)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from granule node ',mpi_granullst%recvs2(i)+1,'. &
 &Expected ',msgtag(4),'.'
        call fexit()
      end if
    end do
    do i=1,mpi_granullst%n_s2
      call MPI_SEND(intmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_INTEGER,mpi_granullst%sends2(i),&
 &                   msgtag(4),MPI_COMM_WORLD,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   finally, every node should have complete correct info
  end if
!
end
!
!---------------------------------------------------------------------------------------------------------------------
!
subroutine MPI_ALLINTS2D(msgtag,opcode,ierr,bc,msgsz2,intmsg2d,MPIcolls)
!
  use mpi
  use mpistuff
!
  implicit none
!
  integer opcode,ierr,ierr2,mstatus(MPI_STATUS_SIZE),msgsz(2),msgsz2(2),msgtag(4),i,masterrank
  integer intmsg2d(msgsz2(1),msgsz2(2))
  logical bc,MPIcolls
!
  ierr = 0
  masterrank = 0
!
  if (mod(msgsz2(2),2).ne.0) then
    write(*,*) 'Fatal. In MPI_ALLINTS2D, messages must be of even rank in the last dimension.'
    call fexit()
  end if
  msgsz(:) = msgsz2(:)
  msgsz(2) = msgsz(2)/2
  intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) = intmsg2d(:,1:msgsz(2))
! block for using standard MPI collective communication syntax
  if (MPIcolls.EQV..true.) then
    if (opcode.eq.0) then
      call MPI_BCAST(intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),MPI_INTEGER,masterrank,MPI_COMM_WORLD,ierr2)
    else if (bc.EQV..true.) then
      if (opcode.eq.1) then
        call MPI_ALLREDUCE(intmsg2d(:,1:msgsz(2)),intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),&
 &                         MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.2) then
        call MPI_ALLREDUCE(intmsg2d(:,1:msgsz(2)),intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),&
 &                         MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.3) then
        call MPI_ALLREDUCE(intmsg2d(:,1:msgsz(2)),intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),&
 &                         MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr2)
      end if
    else
      if (opcode.eq.1) then
        call MPI_REDUCE(intmsg2d(:,1:msgsz(2)),intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),&
 &                      MPI_INTEGER,MPI_SUM,masterrank,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.2) then
        call MPI_REDUCE(intmsg2d(:,1:msgsz(2)),intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),&
 &                      MPI_INTEGER,MPI_MAX,masterrank,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.3) then
        call MPI_REDUCE(intmsg2d(:,1:msgsz(2)),intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),&
 &                      MPI_INTEGER,MPI_MIN,masterrank,MPI_COMM_WORLD,ierr2)
      end if
    end if
    if (ierr2.ne.0) ierr = ierr2
    return
  end if
  if (opcode.gt.0) then
    do i=1,mpi_granullst%n_r
      call MPI_RECV(intmsg2d(:,1:msgsz(2)),msgsz(1)*msgsz(2),MPI_INTEGER,mpi_granullst%recvs(i),&
  &                 MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(1)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from granule node ',mpi_granullst%recvs(i)+1,'. &
  &Expected ',msgtag(1),'.'
        call fexit()
      end if
      if (opcode.eq.1) then ! sum
        intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) = intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) + intmsg2d(:,1:msgsz(2))
      else if (opcode.eq.2) then ! max
        intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) = max(intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),intmsg2d(:,1:msgsz(2))) ! op.s on vector
      else if (opcode.eq.3) then ! min
        intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) = min(intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),intmsg2d(:,1:msgsz(2))) ! op.s on vector
      end if
!     the original message is now lost in favor of the current composite message (if granule%n_r > 0)
    end do
    do i=1,mpi_granullst%n_s 
      call MPI_SEND(intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),MPI_INTEGER,mpi_granullst%sends(i),&
  &                  msgtag(1),MPI_COMM_WORLD,ierr2) ! this sends either a local composite, or the untouched individual values
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   at this point the node with the connection leaving the granule should have the correct info for the granule
    do i=1,mpi_globlst%n_r
      call MPI_RECV(intmsg2d(:,1:msgsz(2)),msgsz(1)*msgsz(2),MPI_INTEGER,mpi_globlst%recvs(i),&
  &                 MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(2)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from external node ',mpi_globlst%recvs(i)+1,'. &
  &Expected ',msgtag(2),'.'
        call fexit()
      end if
      if (opcode.eq.1) then ! sum
        intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) = intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) + intmsg2d(:,1:msgsz(2))
      else if (opcode.eq.2) then ! max
        intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) = max(intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),intmsg2d(:,1:msgsz(2)))
      else if (opcode.eq.3) then ! min
        intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) = min(intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),intmsg2d(:,1:msgsz(2)))
      end if
!     the granule composite is now lost in favor of the current, global composite message (if glob%n_r > 0)
    end do
    do i=1,mpi_globlst%n_s
      call MPI_SEND(intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),MPI_INTEGER,mpi_globlst%sends(i),&
  &                  msgtag(2),MPI_COMM_WORLD,ierr2)! this sends either an intermediate composite, or the untouched granule vals
      if (ierr2.ne.0) ierr = ierr2 
    end do
  end if
! at this point, there should be exactly one node with the complete correct info in the second half of the array
  if (bc.EQV..true.) then
    do i=1,mpi_globlst%n_r2
      call MPI_RECV(intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),MPI_INTEGER,mpi_globlst%recvs2(i),&
 &                  MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(3)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from external node ',mpi_globlst%recvs2(i)+1,'. &
 &Expected ',msgtag(3),'.'
        call fexit()
      end if
    end do
    do i=1,mpi_globlst%n_s2
      call MPI_SEND(intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),MPI_INTEGER,mpi_globlst%sends2(i),&
 &                   msgtag(3),MPI_COMM_WORLD,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   at this point, all those nodes with mpi_granullst%n_r2 == 0 should have complete correct info
    do i=1,mpi_granullst%n_r2
      call MPI_RECV(intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),MPI_INTEGER,mpi_granullst%recvs2(i),&
 &                  MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(4)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from granule node ',mpi_granullst%recvs2(i)+1,'. &
 &Expected ',msgtag(4),'.'
        call fexit()
      end if
    end do
    do i=1,mpi_granullst%n_s2
      call MPI_SEND(intmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),MPI_INTEGER,mpi_granullst%sends2(i),&
 &                   msgtag(4),MPI_COMM_WORLD,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   finally, every node should have complete correct info
  end if
end
!
!--------------------------------------------------------------------------------------------------------------------
!
subroutine MPI_ALLINTS3D(msgtag,opcode,ierr,bc,msgsz2,intmsg3d,MPIcolls)
!
  use mpi
  use mpistuff
!
  implicit none
!
  integer opcode,ierr,ierr2,mstatus(MPI_STATUS_SIZE),msgsz(3),msgsz2(3),msgtag(4),i,masterrank
  integer intmsg3d(msgsz2(1),msgsz2(2),msgsz2(3))
  logical bc,MPIcolls
!
  ierr = 0
  masterrank = 0
!
  if ((opcode.eq.0).AND.(bc.EQV..false.)) return
!
  if (mod(msgsz2(3),2).ne.0) then
    write(*,*) 'Fatal. In MPI_ALLINTS3D, messages must be of even rank in the last dimension.'
    call fexit()
  end if
  msgsz(:) = msgsz2(:)
  msgsz(3) = msgsz(3)/2
  intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) = intmsg3d(:,:,1:msgsz(3))
! block for using standard MPI collective communication syntax
  if (MPIcolls.EQV..true.) then
    if (opcode.eq.0) then
      call MPI_BCAST(intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),MPI_INTEGER,masterrank,MPI_COMM_WORLD,ierr2)
    else if (bc.EQV..true.) then
      if (opcode.eq.1) then
        call MPI_ALLREDUCE(intmsg3d(:,:,1:msgsz(3)),intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),&
 &                         MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.2) then
        call MPI_ALLREDUCE(intmsg3d(:,:,1:msgsz(3)),intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),&
 &                         MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.3) then
        call MPI_ALLREDUCE(intmsg3d(:,:,1:msgsz(3)),intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),&
 &                         MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr2)
      end if
    else
      if (opcode.eq.1) then
        call MPI_REDUCE(intmsg3d(:,:,1:msgsz(3)),intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),&
 &                      MPI_INTEGER,MPI_SUM,masterrank,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.2) then
        call MPI_REDUCE(intmsg3d(:,:,1:msgsz(3)),intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),&
 &                      MPI_INTEGER,MPI_MAX,masterrank,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.3) then
        call MPI_REDUCE(intmsg3d(:,:,1:msgsz(3)),intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),&
 &                      MPI_INTEGER,MPI_MIN,masterrank,MPI_COMM_WORLD,ierr2)
      end if
    end if
    if (ierr2.ne.0) ierr = ierr2
    return
  end if
  if (opcode.gt.0) then
    do i=1,mpi_granullst%n_r
      call MPI_RECV(intmsg3d(:,:,1:msgsz(3)),msgsz(1)*msgsz(2)*msgsz(3),MPI_INTEGER,mpi_granullst%recvs(i),&
  &                 MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(1)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from granule node ',mpi_granullst%recvs(i)+1,'. &
  &Expected ',msgtag(1),'.'
        call fexit()
      end if
      if (opcode.eq.1) then ! sum
        intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) = intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) + intmsg3d(:,:,1:msgsz(3))
      else if (opcode.eq.2) then ! max
        intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) = max(intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),intmsg3d(:,:,1:msgsz(3)))
      else if (opcode.eq.3) then ! min
        intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) = min(intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),intmsg3d(:,:,1:msgsz(3)))
      end if
!     the original message is now lost in favor of the current composite message (if granule%n_r > 0)
    end do
    do i=1,mpi_granullst%n_s 
      call MPI_SEND(intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),MPI_INTEGER,mpi_granullst%sends(i),&
  &                  msgtag(1),MPI_COMM_WORLD,ierr2) ! this sends either a local composite, or the untouched individual values
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   at this point the node with the connection leaving the granule should have the correct info for the granule
    do i=1,mpi_globlst%n_r
      call MPI_RECV(intmsg3d(:,:,1:msgsz(3)),msgsz(1)*msgsz(2)*msgsz(3),MPI_INTEGER,mpi_globlst%recvs(i),&
  &                 MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(2)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from external node ',mpi_globlst%recvs(i)+1,'. &
  &Expected ',msgtag(2),'.'
        call fexit()
      end if
      if (opcode.eq.1) then ! sum
        intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) = intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) + intmsg3d(:,:,1:msgsz(3))
      else if (opcode.eq.2) then ! max
        intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) = max(intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),intmsg3d(:,:,1:msgsz(3)))
      else if (opcode.eq.3) then ! min
        intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) = min(intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),intmsg3d(:,:,1:msgsz(3)))
      end if
!     the granule composite is now lost in favor of the current, global composite message (if glob%n_r > 0)
    end do
    do i=1,mpi_globlst%n_s
      call MPI_SEND(intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),MPI_INTEGER,mpi_globlst%sends(i),&
  &                  msgtag(2),MPI_COMM_WORLD,ierr2)! this sends either an intermediate composite, or the untouched granule vals
      if (ierr2.ne.0) ierr = ierr2 
    end do
  end if
! at this point, there should be exactly one node with the complete correct info in the second half of the array
  if (bc.EQV..true.) then
    do i=1,mpi_globlst%n_r2
      call MPI_RECV(intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),MPI_INTEGER,mpi_globlst%recvs2(i),&
 &                  MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(3)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from external node ',mpi_globlst%recvs2(i)+1,'. &
 &Expected ',msgtag(3),'.'
        call fexit()
      end if
    end do
    do i=1,mpi_globlst%n_s2
      call MPI_SEND(intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),MPI_INTEGER,mpi_globlst%sends2(i),&
 &                   msgtag(3),MPI_COMM_WORLD,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   at this point, all those nodes with mpi_granullst%n_r2 == 0 should have complete correct info
    do i=1,mpi_granullst%n_r2
      call MPI_RECV(intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),MPI_INTEGER,mpi_granullst%recvs2(i),&
 &                  MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(4)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from granule node ',mpi_granullst%recvs2(i)+1,'. &
 &Expected ',msgtag(4),'.'
        call fexit()
      end if
    end do
    do i=1,mpi_granullst%n_s2
      call MPI_SEND(intmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),MPI_INTEGER,mpi_granullst%sends2(i),&
 &                   msgtag(4),MPI_COMM_WORLD,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   finally, every node should have complete correct info
  end if
!
end
!
!---------------------------------------------------------------------------------------------------------------------
!
subroutine MPI_ALLREALS1D(msgtag,opcode,ierr,bc,msgsz2,fltmsg1d,MPIcolls)
!
  use mpi
  use mpistuff
!
  implicit none
!
  integer opcode,ierr,ierr2,mstatus(MPI_STATUS_SIZE),msgsz,msgsz2,msgtag(4),i,masterrank
  RTYPE fltmsg1d(msgsz2)
  logical bc,MPIcolls
!
  ierr = 0
  masterrank = 0
!
  if ((opcode.eq.0).AND.(bc.EQV..false.)) return
!
  if (mod(msgsz2,2).ne.0) then
    write(*,*) 'Fatal. In MPI_ALLREALS1D, messages must be of even rank in the last dimension.'
    call fexit()
  end if
  msgsz = msgsz2/2
  fltmsg1d((msgsz+1):(2*msgsz)) = fltmsg1d(1:msgsz)
! block for using standard MPI collective communication syntax
  if (MPIcolls.EQV..true.) then
    if (opcode.eq.0) then
      call MPI_BCAST(fltmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_RTYPE,masterrank,MPI_COMM_WORLD,ierr2)
    else if (bc.EQV..true.) then
      if (opcode.eq.1) then
        call MPI_ALLREDUCE(fltmsg1d(1:msgsz),fltmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_RTYPE,MPI_SUM,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.2) then
        call MPI_ALLREDUCE(fltmsg1d(1:msgsz),fltmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_RTYPE,MPI_MAX,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.3) then
        call MPI_ALLREDUCE(fltmsg1d(1:msgsz),fltmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_RTYPE,MPI_MIN,MPI_COMM_WORLD,ierr2)
      end if
    else
      if (opcode.eq.1) then
        call MPI_REDUCE(fltmsg1d(1:msgsz),fltmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_RTYPE,MPI_SUM,masterrank,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.2) then
        call MPI_REDUCE(fltmsg1d(1:msgsz),fltmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_RTYPE,MPI_MAX,masterrank,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.3) then
        call MPI_REDUCE(fltmsg1d(1:msgsz),fltmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_RTYPE,MPI_MIN,masterrank,MPI_COMM_WORLD,ierr2)
      end if
    end if
    if (ierr2.ne.0) ierr = ierr2
    return
  end if
  if (opcode.gt.0) then
    do i=1,mpi_granullst%n_r
      call MPI_RECV(fltmsg1d(1:msgsz),msgsz,MPI_RTYPE,mpi_granullst%recvs(i),&
  &                 MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(1)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from granule node ',mpi_granullst%recvs(i)+1,'. &
  &Expected ',msgtag(1),'.'
        call fexit()
      end if
      if (opcode.eq.1) then ! sum
        fltmsg1d((msgsz+1):(2*msgsz)) = fltmsg1d((msgsz+1):(2*msgsz)) + fltmsg1d(1:msgsz)
      else if (opcode.eq.2) then ! max
        fltmsg1d((msgsz+1):(2*msgsz)) = max(fltmsg1d((msgsz+1):(2*msgsz)),fltmsg1d(1:msgsz)) ! op.s on vector
      else if (opcode.eq.3) then ! min
        fltmsg1d((msgsz+1):(2*msgsz)) = min(fltmsg1d((msgsz+1):(2*msgsz)),fltmsg1d(1:msgsz)) ! op.s on vector
      end if
!     the original message is now lost in favor of the current composite message (if granule%n_r > 0)
    end do
    do i=1,mpi_granullst%n_s 
      call MPI_SEND(fltmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_RTYPE,mpi_granullst%sends(i),&
  &                  msgtag(1),MPI_COMM_WORLD,ierr2) ! this sends either a local composite, or the untouched individual values
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   at this point the node with the connection leaving the granule should have the correct info for the granule
    do i=1,mpi_globlst%n_r
      call MPI_RECV(fltmsg1d(1:msgsz),msgsz,MPI_RTYPE,mpi_globlst%recvs(i),&
  &                 MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(2)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from external node ',mpi_globlst%recvs(i)+1,'. &
  &Expected ',msgtag(2),'.'
        call fexit()
      end if
      if (opcode.eq.1) then ! sum
        fltmsg1d((msgsz+1):(2*msgsz)) = fltmsg1d((msgsz+1):(2*msgsz)) + fltmsg1d(1:msgsz)
      else if (opcode.eq.2) then ! max
        fltmsg1d((msgsz+1):(2*msgsz)) = max(fltmsg1d((msgsz+1):(2*msgsz)),fltmsg1d(1:msgsz))
      else if (opcode.eq.3) then ! min
        fltmsg1d((msgsz+1):(2*msgsz)) = min(fltmsg1d((msgsz+1):(2*msgsz)),fltmsg1d(1:msgsz))
      end if
!     the granule composite is now lost in favor of the current, global composite message (if glob%n_r > 0)
    end do
    do i=1,mpi_globlst%n_s
      call MPI_SEND(fltmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_RTYPE,mpi_globlst%sends(i),&
  &                  msgtag(2),MPI_COMM_WORLD,ierr2)! this sends either an intermediate composite, or the untouched granule vals
      if (ierr2.ne.0) ierr = ierr2 
    end do
  end if
! at this point, there should be exactly one node with the complete correct info in the second half of the array
  if (bc.EQV..true.) then
    do i=1,mpi_globlst%n_r2
      call MPI_RECV(fltmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_RTYPE,mpi_globlst%recvs2(i),&
 &                  MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(3)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from external node ',mpi_globlst%recvs2(i)+1,'. &
 &Expected ',msgtag(3),'.'
        call fexit()
      end if
    end do
    do i=1,mpi_globlst%n_s2
      call MPI_SEND(fltmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_RTYPE,mpi_globlst%sends2(i),&
 &                   msgtag(3),MPI_COMM_WORLD,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   at this point, all those nodes with mpi_granullst%n_r2 == 0 should have complete correct info
    do i=1,mpi_granullst%n_r2
      call MPI_RECV(fltmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_RTYPE,mpi_granullst%recvs2(i),&
 &                  MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(4)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from granule node ',mpi_granullst%recvs2(i)+1,'. &
 &Expected ',msgtag(4),'.'
        call fexit()
      end if
    end do
    do i=1,mpi_granullst%n_s2
      call MPI_SEND(fltmsg1d((msgsz+1):(2*msgsz)),msgsz,MPI_RTYPE,mpi_granullst%sends2(i),&
 &                   msgtag(4),MPI_COMM_WORLD,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   finally, every node should have complete correct info
  end if
!
end
!
!---------------------------------------------------------------------------------------------------------------------
!
subroutine MPI_ALLREALS2D(msgtag,opcode,ierr,bc,msgsz2,fltmsg2d,MPIcolls)
!
  use mpi
  use mpistuff
!
  implicit none
!
  integer opcode,ierr,ierr2,mstatus(MPI_STATUS_SIZE),msgsz(2),msgsz2(2),msgtag(4),i,masterrank
  RTYPE fltmsg2d(msgsz2(1),msgsz2(2))
  logical bc,MPIcolls
!
  ierr = 0
  masterrank = 0
!
  if (mod(msgsz2(2),2).ne.0) then
    write(*,*) 'Fatal. In MPI_ALLREALS2D, messages must be of even rank in the last dimension.'
    call fexit()
  end if
  msgsz(:) = msgsz2(:)
  msgsz(2) = msgsz(2)/2
  fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) = fltmsg2d(:,1:msgsz(2))
  if (MPIcolls.EQV..true.) then
    if (opcode.eq.0) then
      call MPI_BCAST(fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),MPI_RTYPE,masterrank,MPI_COMM_WORLD,ierr2)
    else if (bc.EQV..true.) then
      if (opcode.eq.1) then
        call MPI_ALLREDUCE(fltmsg2d(:,1:msgsz(2)),fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),&
 &                         MPI_RTYPE,MPI_SUM,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.2) then
        call MPI_ALLREDUCE(fltmsg2d(:,1:msgsz(2)),fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),&
 &                         MPI_RTYPE,MPI_MAX,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.3) then
        call MPI_ALLREDUCE(fltmsg2d(:,1:msgsz(2)),fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),&
 &                         MPI_RTYPE,MPI_MIN,MPI_COMM_WORLD,ierr2)
      end if
    else
      if (opcode.eq.1) then
        call MPI_REDUCE(fltmsg2d(:,1:msgsz(2)),fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),&
 &                      MPI_RTYPE,MPI_SUM,masterrank,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.2) then
        call MPI_REDUCE(fltmsg2d(:,1:msgsz(2)),fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),&
 &                      MPI_RTYPE,MPI_MAX,masterrank,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.3) then
        call MPI_REDUCE(fltmsg2d(:,1:msgsz(2)),fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),&
 &                      MPI_RTYPE,MPI_MIN,masterrank,MPI_COMM_WORLD,ierr2)
      end if
    end if
    if (ierr2.ne.0) ierr = ierr2
    return
  end if
  if (opcode.gt.0) then
    do i=1,mpi_granullst%n_r
      call MPI_RECV(fltmsg2d(:,1:msgsz(2)),msgsz(1)*msgsz(2),MPI_RTYPE,mpi_granullst%recvs(i),&
  &                 MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(1)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from granule node ',mpi_granullst%recvs(i)+1,'. &
  &Expected ',msgtag(1),'.'
        call fexit()
      end if
      if (opcode.eq.1) then ! sum
        fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) = fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) + fltmsg2d(:,1:msgsz(2))
      else if (opcode.eq.2) then ! max
        fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) = max(fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),fltmsg2d(:,1:msgsz(2))) ! op.s on vector
      else if (opcode.eq.3) then ! min
        fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) = min(fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),fltmsg2d(:,1:msgsz(2))) ! op.s on vector
      end if
!     the original message is now lost in favor of the current composite message (if granule%n_r > 0)
    end do
    do i=1,mpi_granullst%n_s 
      call MPI_SEND(fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),MPI_RTYPE,mpi_granullst%sends(i),&
  &                  msgtag(1),MPI_COMM_WORLD,ierr2) ! this sends either a local composite, or the untouched individual values
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   at this point the node with the connection leaving the granule should have the correct info for the granule
    do i=1,mpi_globlst%n_r
      call MPI_RECV(fltmsg2d(:,1:msgsz(2)),msgsz(1)*msgsz(2),MPI_RTYPE,mpi_globlst%recvs(i),&
  &                 MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(2)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from external node ',mpi_globlst%recvs(i)+1,'. &
  &Expected ',msgtag(2),'.'
        call fexit()
      end if
      if (opcode.eq.1) then ! sum
        fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) = fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) + fltmsg2d(:,1:msgsz(2))
      else if (opcode.eq.2) then ! max
        fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) = max(fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),fltmsg2d(:,1:msgsz(2)))
      else if (opcode.eq.3) then ! min
        fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))) = min(fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),fltmsg2d(:,1:msgsz(2)))
      end if
!     the granule composite is now lost in favor of the current, global composite message (if glob%n_r > 0)
    end do
    do i=1,mpi_globlst%n_s
      call MPI_SEND(fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),MPI_RTYPE,mpi_globlst%sends(i),&
  &                  msgtag(2),MPI_COMM_WORLD,ierr2)! this sends either an intermediate composite, or the untouched granule vals
      if (ierr2.ne.0) ierr = ierr2 
    end do
  end if
! at this point, there should be exactly one node with the complete correct info in the second half of the array
  if (bc.EQV..true.) then
    do i=1,mpi_globlst%n_r2
      call MPI_RECV(fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),MPI_RTYPE,mpi_globlst%recvs2(i),&
 &                  MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(3)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from external node ',mpi_globlst%recvs2(i)+1,'. &
 &Expected ',msgtag(3),'.'
        call fexit()
      end if
    end do
    do i=1,mpi_globlst%n_s2
      call MPI_SEND(fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),MPI_RTYPE,mpi_globlst%sends2(i),&
 &                   msgtag(3),MPI_COMM_WORLD,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   at this point, all those nodes with mpi_granullst%n_r2 == 0 should have complete correct info
    do i=1,mpi_granullst%n_r2
      call MPI_RECV(fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),MPI_RTYPE,mpi_granullst%recvs2(i),&
 &                  MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(4)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from granule node ',mpi_granullst%recvs2(i)+1,'. &
 &Expected ',msgtag(4),'.'
        call fexit()
      end if
    end do
    do i=1,mpi_granullst%n_s2
      call MPI_SEND(fltmsg2d(:,(msgsz(2)+1):(2*msgsz(2))),msgsz(1)*msgsz(2),MPI_RTYPE,mpi_granullst%sends2(i),&
 &                   msgtag(4),MPI_COMM_WORLD,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   finally, every node should have complete correct info
  end if
end
!
!--------------------------------------------------------------------------------------------------------------------
!
subroutine MPI_ALLREALS3D(msgtag,opcode,ierr,bc,msgsz2,fltmsg3d,MPIcolls)
!
  use mpi
  use mpistuff
!
  implicit none
!
  integer opcode,ierr,ierr2,mstatus(MPI_STATUS_SIZE),msgsz(3),msgsz2(3),msgtag(4),i,masterrank
  RTYPE fltmsg3d(msgsz2(1),msgsz2(2),msgsz2(3))
  logical bc,MPIcolls
!
  ierr = 0
  masterrank = 0
!
  if ((opcode.eq.0).AND.(bc.EQV..false.)) return
!
  if (mod(msgsz2(3),2).ne.0) then
    write(*,*) 'Fatal. In MPI_ALLREALS3D, messages must be of even rank in the last dimension.'
    call fexit()
  end if
  msgsz(:) = msgsz2(:)
  msgsz(3) = msgsz(3)/2
  fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) = fltmsg3d(:,:,1:msgsz(3))
! block for using standard MPI collective communication syntax
  if (MPIcolls.EQV..true.) then
    if (opcode.eq.0) then
      call MPI_BCAST(fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),MPI_RTYPE,masterrank,MPI_COMM_WORLD,ierr2)
    else if (bc.EQV..true.) then
      if (opcode.eq.1) then
        call MPI_ALLREDUCE(fltmsg3d(:,:,1:msgsz(3)),fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),&
 &                         MPI_RTYPE,MPI_SUM,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.2) then
        call MPI_ALLREDUCE(fltmsg3d(:,:,1:msgsz(3)),fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),&
 &                         MPI_RTYPE,MPI_MAX,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.3) then
        call MPI_ALLREDUCE(fltmsg3d(:,:,1:msgsz(3)),fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),&
 &                         MPI_RTYPE,MPI_MIN,MPI_COMM_WORLD,ierr2)
      end if
    else
      if (opcode.eq.1) then
        call MPI_REDUCE(fltmsg3d(:,:,1:msgsz(3)),fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),&
 &                      MPI_RTYPE,MPI_SUM,masterrank,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.2) then
        call MPI_REDUCE(fltmsg3d(:,:,1:msgsz(3)),fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),&
 &                      MPI_RTYPE,MPI_MAX,masterrank,MPI_COMM_WORLD,ierr2)
      else if (opcode.eq.3) then
        call MPI_REDUCE(fltmsg3d(:,:,1:msgsz(3)),fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),&
 &                      MPI_RTYPE,MPI_MIN,masterrank,MPI_COMM_WORLD,ierr2)
      end if
    end if
    if (ierr2.ne.0) ierr = ierr2
    return
  end if
! now the custom version
  if (opcode.gt.0) then
    do i=1,mpi_granullst%n_r
      call MPI_RECV(fltmsg3d(:,:,1:msgsz(3)),msgsz(1)*msgsz(2)*msgsz(3),MPI_RTYPE,mpi_granullst%recvs(i),&
  &                 MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(1)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from granule node ',mpi_granullst%recvs(i)+1,'. &
  &Expected ',msgtag(1),'.'
        call fexit()
      end if
      if (opcode.eq.1) then ! sum
        fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) = fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) + fltmsg3d(:,:,1:msgsz(3))
      else if (opcode.eq.2) then ! max
        fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) = max(fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),fltmsg3d(:,:,1:msgsz(3))) 
      else if (opcode.eq.3) then ! min
        fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) = min(fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),fltmsg3d(:,:,1:msgsz(3)))
      end if
!     the original message is now lost in favor of the current composite message (if granule%n_r > 0)
    end do
    do i=1,mpi_granullst%n_s 
      call MPI_SEND(fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),MPI_RTYPE,mpi_granullst%sends(i),&
  &                  msgtag(1),MPI_COMM_WORLD,ierr2) ! this sends either a local composite, or the untouched individual values
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   at this point the node with the connection leaving the granule should have the correct info for the granule
    do i=1,mpi_globlst%n_r
      call MPI_RECV(fltmsg3d(:,:,1:msgsz(3)),msgsz(1)*msgsz(2)*msgsz(3),MPI_RTYPE,mpi_globlst%recvs(i),&
  &                 MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(2)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from external node ',mpi_globlst%recvs(i)+1,'. &
  &Expected ',msgtag(2),'.'
        call fexit()
      end if
      if (opcode.eq.1) then ! sum
        fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) = fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) + fltmsg3d(:,:,1:msgsz(3))
      else if (opcode.eq.2) then ! max
        fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) = max(fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),fltmsg3d(:,:,1:msgsz(3)))
      else if (opcode.eq.3) then ! min
        fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))) = min(fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),fltmsg3d(:,:,1:msgsz(3)))
      end if
!     the granule composite is now lost in favor of the current, global composite message (if glob%n_r > 0)
    end do
    do i=1,mpi_globlst%n_s
      call MPI_SEND(fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),MPI_RTYPE,mpi_globlst%sends(i),&
  &                  msgtag(2),MPI_COMM_WORLD,ierr2)! this sends either an intermediate composite, or the untouched granule vals
      if (ierr2.ne.0) ierr = ierr2 
    end do
  end if
! at this point, there should be exactly one node with the complete correct info in the second half of the array
  if (bc.EQV..true.) then
    do i=1,mpi_globlst%n_r2
      call MPI_RECV(fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),MPI_RTYPE,mpi_globlst%recvs2(i),&
 &                  MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(3)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from external node ',mpi_globlst%recvs2(i)+1,'. &
 &Expected ',msgtag(3),'.'
        call fexit()
      end if
    end do
    do i=1,mpi_globlst%n_s2
      call MPI_SEND(fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),MPI_RTYPE,mpi_globlst%sends2(i),&
 &                   msgtag(3),MPI_COMM_WORLD,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   at this point, all those nodes with mpi_granullst%n_r2 == 0 should have complete correct info
    do i=1,mpi_granullst%n_r2
      call MPI_RECV(fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),MPI_RTYPE,mpi_granullst%recvs2(i),&
 &                  MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
      if (mstatus(MPI_TAG).ne.msgtag(4)) then
        write(*,*) 'Fatal. Received bad message (',mstatus(MPI_TAG),') from granule node ',mpi_granullst%recvs2(i)+1,'. &
 &Expected ',msgtag(4),'.'
        call fexit()
      end if
    end do
    do i=1,mpi_granullst%n_s2
      call MPI_SEND(fltmsg3d(:,:,(msgsz(3)+1):(2*msgsz(3))),msgsz(1)*msgsz(2)*msgsz(3),MPI_RTYPE,mpi_granullst%sends2(i),&
 &                   msgtag(4),MPI_COMM_WORLD,ierr2)
      if (ierr2.ne.0) ierr = ierr2 
    end do
!   finally, every node should have complete correct info
  end if
!
end
!
!---------------------------------------------------------------------------------------------------------------------
!
subroutine MPI_assigncommflow()
!
  use iounit
  use mpistuff
!
  implicit none
!
  integer i,k,ki,MPI_bestk,grani,granc,granh
!
! sanity
  if (mpi_granularity.gt.mpi_nodes) then
    write(ilog,*) 'Warning. The specified granule size for the MPI communication layout exceeds the number of processes assigned &
 &by MPI. This means that a single granule is assumed.'
    mpi_granularity = mpi_nodes
  end if
!
! granule level communication (presumed fast)
  mpi_granullst%n_s = 0
  mpi_granullst%n_s2 = 0
  mpi_granullst%n_r = 0
  mpi_granullst%n_r2 = 0
  if (mpi_granularity.eq.1) then
    granc = mpi_nodes
  else
    if (mod(myrank+mpi_granularity,mpi_granularity).eq.0) then
      grani = ((myrank+mpi_granularity)/mpi_granularity)
    else
      grani = ceiling(dble(myrank)/dble(mpi_granularity))
    end if
    if (mod(mpi_nodes,mpi_granularity).eq.0) then
      granc = mpi_nodes/mpi_granularity
    else
      granc = ceiling(dble(mpi_nodes)/dble(mpi_granularity))
    end if
!   primitive for now
    granh = (grani-1)*mpi_granularity
    if ((granh.eq.grani).AND.(granh.eq.(mpi_nodes-1))) then
!     do nothing (single node overhanging as separate granule)
    else
      k = min(mpi_granularity-1,mpi_nodes-1-granh)
      if (granh.eq.myrank) then
        mpi_granullst%n_r = k
        mpi_granullst%n_s2 = k
        allocate(mpi_granullst%recvs(mpi_granullst%n_r))
        allocate(mpi_granullst%sends2(mpi_granullst%n_s2))
        do i=granh+1,granh+k
          mpi_granullst%recvs(i-granh) = i
          mpi_granullst%sends2(i-granh) = i
        end do
        mpi_granullst%n_s = 0
        mpi_granullst%n_r2 = 0
      else
        mpi_granullst%n_s = 1
        allocate(mpi_granullst%sends(mpi_granullst%n_s))
        mpi_granullst%sends(1) = granh
        mpi_granullst%n_r = 0
        mpi_granullst%n_s2 = 0
        mpi_granullst%n_r2 = 1
        allocate(mpi_granullst%recvs2(mpi_granullst%n_r2))
        mpi_granullst%recvs2(1) = granh
      end if
    end if
  end if
!
! external communication (presumed slow)
  mpi_globlst%n_s = 0
  mpi_globlst%n_s2 = 0
  mpi_globlst%n_r = 0
  mpi_globlst%n_r2 = 0
  if ((granc.gt.1).AND.(mod(myrank+mpi_granularity,mpi_granularity).eq.0)) then
!   dry run to get alloc sizes
    ki = granc-1
    k = MPI_bestk(ki)
    ki = 0
    do i=granc,1,-1
      grani = (i-1)*mpi_granularity
      if ((ki.lt.k).AND.(i.gt.1)) then 
        if (myrank.eq.grani) then
          mpi_globlst%n_r2 = mpi_globlst%n_r2 + 1
          mpi_globlst%n_s = mpi_globlst%n_s + 1
        else if (myrank.eq.(grani-mpi_granularity)) then
          mpi_globlst%n_s2 = mpi_globlst%n_s2 + 1
          mpi_globlst%n_r = mpi_globlst%n_r + 1
        end if
      else
        if (i.eq.1) cycle
        if (myrank.eq.grani) then
          mpi_globlst%n_r2 = mpi_globlst%n_r2 + 1
          mpi_globlst%n_s = mpi_globlst%n_s + 1
        else if (myrank.eq.0) then
          mpi_globlst%n_s2 = mpi_globlst%n_s2 + 1
          mpi_globlst%n_r = mpi_globlst%n_r + 1
        end if
      end if
      if (ki.eq.k) then
        k = k + 1
        ki = 0
      else
        ki = ki + 1
      end if
    end do
    if (mpi_globlst%n_s.gt.0) allocate(mpi_globlst%sends(mpi_globlst%n_s))
    if (mpi_globlst%n_s2.gt.0) allocate(mpi_globlst%sends2(mpi_globlst%n_s2))
    if (mpi_globlst%n_r.gt.0) allocate(mpi_globlst%recvs(mpi_globlst%n_r))
    if (mpi_globlst%n_r2.gt.0) allocate(mpi_globlst%recvs2(mpi_globlst%n_r2))
!   now the real run
    mpi_globlst%n_s = 0
    mpi_globlst%n_s2 = 0
    mpi_globlst%n_r = 0
    mpi_globlst%n_r2 = 0
    ki = granc-1
    k = MPI_bestk(ki)
    ki = 0
    do i=granc,1,-1
      grani = (i-1)*mpi_granularity
      if ((ki.lt.k).AND.(i.gt.1)) then 
        if (myrank.eq.grani) then
          mpi_globlst%n_r2 = mpi_globlst%n_r2 + 1
          mpi_globlst%n_s = mpi_globlst%n_s + 1
          mpi_globlst%recvs2(mpi_globlst%n_r2) = grani-mpi_granularity
          mpi_globlst%sends(mpi_globlst%n_s) = grani-mpi_granularity
        else if (myrank.eq.(grani-mpi_granularity)) then
          mpi_globlst%n_s2 = mpi_globlst%n_s2 + 1
          mpi_globlst%n_r = mpi_globlst%n_r + 1
          mpi_globlst%sends2(mpi_globlst%n_s2) = grani
          mpi_globlst%recvs(mpi_globlst%n_r) = grani
        end if
      else
        if (i.eq.1) cycle
        if (myrank.eq.grani) then
          mpi_globlst%n_r2 = mpi_globlst%n_r2 + 1
          mpi_globlst%n_s = mpi_globlst%n_s + 1
          mpi_globlst%recvs2(mpi_globlst%n_r2) = 0
          mpi_globlst%sends(mpi_globlst%n_s) = 0
        else if (myrank.eq.0) then
          mpi_globlst%n_s2 = mpi_globlst%n_s2 + 1
          mpi_globlst%n_r = mpi_globlst%n_r + 1
          mpi_globlst%sends2(mpi_globlst%n_s2) = grani
          mpi_globlst%recvs(mpi_globlst%n_r) = grani
        end if
      end if
      if (ki.eq.k) then
        k = k + 1
        ki = 0
      else
        ki = ki + 1
      end if
    end do
  end if
!
!  write(ilog,*) 'My rank: ',myrank
!  write(ilog,*) 'Granule: ',mpi_granullst%n_s,mpi_granullst%n_s2,mpi_granullst%n_r,mpi_granullst%n_r2
!  if (mpi_granullst%n_s.gt.0) write(ilog,*) mpi_granullst%sends(1:mpi_granullst%n_s)
!  if (mpi_granullst%n_s2.gt.0) write(ilog,*) mpi_granullst%sends2(1:mpi_granullst%n_s2)
!  if (mpi_granullst%n_r.gt.0) write(ilog,*) mpi_granullst%recvs(1:mpi_granullst%n_r)
!  if (mpi_granullst%n_r2.gt.0) write(ilog,*) mpi_granullst%recvs2(1:mpi_granullst%n_r2)
!  write(ilog,*) 'Global: ',mpi_globlst%n_s,mpi_globlst%n_s2,mpi_globlst%n_r,mpi_globlst%n_r2
!  if (mpi_globlst%n_s.gt.0) write(ilog,*) mpi_globlst%sends(1:mpi_globlst%n_s)
!  if (mpi_globlst%n_s2.gt.0) write(ilog,*) mpi_globlst%sends2(1:mpi_globlst%n_s2)
!  if (mpi_globlst%n_r.gt.0) write(ilog,*) mpi_globlst%recvs(1:mpi_globlst%n_r)
!  if (mpi_globlst%n_r2.gt.0) write(ilog,*) mpi_globlst%recvs2(1:mpi_globlst%n_r2)
!
end
!
!------------------------------------------------------------------------------------------------
!
function MPI_bestk(targetn)
!
  integer ki,k,targetn,MPI_bestk
!
  ki = 0
  k = 2
  do while (ki.lt.targetn)
    ki = ki + k
    k = k + 1
    if (ki.eq.(targetn)) then
      MPI_bestk = 1
      return
    else if (ki+1.eq.(targetn)) then
      MPI_bestk = 0
      return
    else if (ki.gt.(targetn)) then
      MPI_bestk = mod(targetn+1,2)
    end if
  end do
!
end
!
!---------------------------------------------------------------------------
!
#else
!
subroutine iamastupidmpisubstitute()
!
  implicit none
!
end
!
#endif
