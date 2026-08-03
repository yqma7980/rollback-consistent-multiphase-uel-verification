program wp10_i5_integrated_fd_harness
  use, intrinsic :: ieee_arithmetic
  use ifport
  implicit none

  integer, parameter :: ndofel=16, nrhs=1, nsvars=192
  integer, parameter :: nprops=27, mcrd=2, nnode=4
  integer, parameter :: nstate=13, ndir=3, nh=11, nblock=9
  integer, parameter :: iu(8)=(/1,2,5,6,9,10,13,14/)
  integer, parameter :: ip(4)=(/3,7,11,15/)
  integer, parameter :: isat(4)=(/4,8,12,16/)
  real(8), parameter :: sigma_global=-1.0d0
  real(8), parameter :: scale_u=1.0d-5
  real(8), parameter :: scale_p=1.0d7
  real(8), parameter :: scale_s=1.0d-1
  real(8), parameter :: kstab_num=1.0d-12
  real(8), parameter :: hvals(nh)=(/1.0d-2,3.0d-3,1.0d-3, &
       3.0d-4,1.0d-4,3.0d-5,1.0d-5,3.0d-6,1.0d-6, &
       3.0d-7,1.0d-7/)
  character(len=3), parameter :: block_names(nblock)=(/ &
       'KUU','KUP','KUS','KPU','KPP','KPS','KSU','KSP','KSS'/)

  real(8) :: props(nprops,nstate), coords(mcrd,nnode)
  real(8) :: ubase(ndofel,nstate), hcomm(nsvars,nstate)
  real(8) :: dir_u(8,ndir), dir_p(4,ndir), dir_s(4,ndir)
  real(8) :: rhs0(ndofel), jac0(ndofel,ndofel), hc0(nsvars)
  real(8) :: rhs_ref(ndofel,nstate),jac_ref(ndofel,ndofel,nstate)
  real(8) :: svars_ref(nsvars,nstate)
  real(8) :: jac_phys(ndofel,ndofel)
  real(8) :: rhsp(ndofel), jacp(ndofel,ndofel), hcp(nsvars)
  real(8) :: rhsm(ndofel), jacm(ndofel,ndofel), hcm(nsvars)
  real(8) :: uwork(ndofel), delta(ndofel), fd(8), kv(8), dvec(8)
  real(8) :: pnewdt, scale, norm_fd, norm_kv, err_l2, rel_l2
  real(8) :: err_inf, rel_inf, cosine
  integer :: rows(8), cols(8), nrows, ncols
  integer :: s,b,d,k,j,cp,cm,base_branch
  integer :: outunit,replayunit
  integer :: rhs_ok,jac_ok,svars_ok,branch_crossed,replay_ok(nstate)
  character(len=512) :: rawfile,dirfile,statefile,replayfile,run_id
  character(len=512) :: snapshotfile,subpatchfile,branchfile
  character(len=512) :: termfile,switchfile,fullfile,mixdirfile
  character(len=16) :: sid,did,bp,bm

  call get_command_argument(1,rawfile)
  call get_command_argument(2,dirfile)
  call get_command_argument(3,statefile)
  call get_command_argument(4,replayfile)
  call get_command_argument(5,run_id)
  call get_command_argument(6,snapshotfile)
  call get_command_argument(7,subpatchfile)
  call get_command_argument(8,branchfile)
  call get_command_argument(9,termfile)
  call get_command_argument(10,switchfile)
  call get_command_argument(11,fullfile)
  call get_command_argument(12,mixdirfile)
  if (len_trim(rawfile)==0 .or. len_trim(dirfile)==0 .or. &
      len_trim(statefile)==0 .or. len_trim(replayfile)==0 .or. &
      len_trim(snapshotfile)==0 .or. len_trim(subpatchfile)==0 .or. &
      len_trim(branchfile)==0 .or. len_trim(termfile)==0 .or. &
      len_trim(switchfile)==0 .or. len_trim(fullfile)==0 .or. &
      len_trim(mixdirfile)==0) stop 2
  if (len_trim(run_id)==0) run_id='run_unknown'

  call initialize_directions(dir_u,dir_p,dir_s)
  call initialize_states(props,coords,ubase,hcomm)
  call write_direction_manifest(trim(dirfile),dir_u,dir_p,dir_s)
  call write_mixed_direction_manifest(trim(mixdirfile),dir_u,dir_p,dir_s)

  open(newunit=replayunit,file=trim(replayfile),status='replace', &
       action='write')
  write(replayunit,'(a)') 'run_id,replay_stage,state_id,rhs_exact,'// &
       'amatrx_exact,'// &
       'physical_svars_exact,max_abs_rhs,max_abs_amatrx,'// &
       'max_abs_physical_svars,finite,baseline_branch,replay_status'

  do s=1,nstate
     call state_id(s,sid)
     call evaluate_uel(s,ubase(:,s),ubase(:,s),hcomm(:,s), &
          rhs0,jac0,hc0,pnewdt)
     rhs_ref(:,s)=rhs0
     jac_ref(:,:,s)=jac0
     svars_ref(:,s)=hc0
     call physical_jacobian(jac0,jac_phys)
     call fd_branch_signature(s,ubase(:,s),hcomm(:,s),props(:,s), &
          base_branch,bp)
     if (.not.used_finite(rhs0,jac0) .or. &
         .not.physical_svars_finite(hc0) .or. &
         base_branch/=expected_branch(s)) then
        write(replayunit,'(*(g0,:,","))') trim(run_id),'PRE_FD', &
             trim(sid),0,0,0, &
             0.0d0,0.0d0,0.0d0,0,base_branch,'BASELINE_INVALID'
        close(replayunit)
        stop 31
     end if

     uwork=ubase(:,s)
     uwork(iu(3))=uwork(iu(3))+2.0d-8
     uwork(ip(2))=uwork(ip(2))+2.0d4
     uwork(isat(4))=uwork(isat(4))+2.0d-4
     call evaluate_uel(s,uwork,ubase(:,s),hcomm(:,s), &
          rhsp,jacp,hcp,pnewdt)

     uwork=ubase(:,s)
     uwork(iu(6))=uwork(iu(6))-1.0d-8
     uwork(ip(4))=uwork(ip(4))-1.0d4
     uwork(isat(1))=uwork(isat(1))-1.0d-4
     call evaluate_uel(s,uwork,ubase(:,s),hcomm(:,s), &
          rhsm,jacm,hcm,pnewdt)

     call evaluate_uel(s,ubase(:,s),ubase(:,s),hcomm(:,s), &
          rhsp,jacp,hcp,pnewdt)
     replay_ok(s)=merge(1,0,all(rhs0==rhsp) .and. all(jac0==jacp) &
          .and. physical_svars_exact(hc0,hcp))
     write(replayunit,'(*(g0,:,","))') trim(run_id),'PRE_FD', &
          trim(sid), &
          merge(1,0,all(rhs0==rhsp)),merge(1,0,all(jac0==jacp)), &
          merge(1,0,physical_svars_exact(hc0,hcp)), &
          maxval(abs(rhs0-rhsp)),maxval(abs(jac0-jacp)), &
          max_physical_svars_difference(hc0,hcp), &
          merge(1,0,used_finite(rhsp,jacp) .and. &
          physical_svars_finite(hcp)),base_branch, &
          merge('PASS','FAIL',replay_ok(s)==1)
     if (replay_ok(s)/=1) then
        close(replayunit)
        stop 32
     end if
  end do
  close(replayunit)

  call write_state_manifest(trim(statefile),props,coords,ubase,hcomm)
  call write_snapshot(trim(snapshotfile),trim(run_id),props,ubase,hcomm)
  call write_branch_manifest(trim(branchfile),trim(run_id),props, &
       ubase,hcomm)
  call write_term_audit(trim(termfile),trim(run_id),props,ubase,hcomm)
  call write_switch_diagnostics(trim(switchfile),trim(run_id),props)

  open(newunit=outunit,file=trim(rawfile),status='replace',action='write')
  write(outunit,'(a)') 'run_id,state_id,block_id,direction_id,h,sigma,'// &
       'branch_id_plus,branch_id_minus,norm_fd,norm_kv,'// &
       'absolute_error_l2,relative_error_l2,absolute_error_inf,'// &
       'relative_error_inf,cosine_similarity,rhs_finite,'// &
       'amatrx_finite,candidate_svars_finite,replay_exact,branch_crossed'

  do s=1,nstate
     call state_id(s,sid)
     call evaluate_uel(s,ubase(:,s),ubase(:,s),hcomm(:,s), &
          rhs0,jac0,hc0,pnewdt)
     call physical_jacobian(jac0,jac_phys)
     call fd_branch_signature(s,ubase(:,s),hcomm(:,s),props(:,s), &
          base_branch,bp)
     do b=1,nblock
        call block_map(b,rows,cols,nrows,ncols,scale)
        do d=1,ndir
           write(did,'("D",i1)') d
           call get_direction(b,d,dir_u,dir_p,dir_s,dvec,ncols)
           do k=1,nh
              delta=0.0d0
              do j=1,ncols
                 delta(cols(j))=hvals(k)*scale*dvec(j)
              end do

              uwork=ubase(:,s)+delta
              call evaluate_uel(s,uwork,ubase(:,s),hcomm(:,s), &
                   rhsp,jacp,hcp,pnewdt)
              call fd_branch_signature(s,uwork,hcomm(:,s),props(:,s), &
                   cp,bp)

              uwork=ubase(:,s)-delta
              call evaluate_uel(s,uwork,ubase(:,s),hcomm(:,s), &
                   rhsm,jacm,hcm,pnewdt)
              call fd_branch_signature(s,uwork,hcomm(:,s),props(:,s), &
                   cm,bm)

               do j=1,nrows
                 fd(j)=(rhsp(rows(j))-rhsm(rows(j)))/(2.0d0*hvals(k))
                 kv(j)=0.0d0
                 kv(j)=sigma_global*sum(jac_phys(rows(j),cols(1:ncols)) &
                      * (scale*dvec(1:ncols)))
              end do

               norm_fd=l2norm(fd(1:nrows))
               norm_kv=l2norm(kv(1:nrows))
               err_l2=l2norm(fd(1:nrows)-kv(1:nrows))
              rel_l2=err_l2/max(norm_fd,norm_kv,1.0d-30)
               err_inf=maxval(abs(fd(1:nrows)-kv(1:nrows)))
               rel_inf=err_inf/max(maxval(abs(fd(1:nrows))), &
                    maxval(abs(kv(1:nrows))),1.0d-30)
              if (norm_fd>1.0d-300 .and. norm_kv>1.0d-300) then
                  cosine=sum(fd(1:nrows)*kv(1:nrows))/(norm_fd*norm_kv)
              else if (norm_fd<=1.0d-300 .and. norm_kv<=1.0d-300) then
                 cosine=1.0d0
              else
                 cosine=0.0d0
              end if

              rhs_ok=merge(1,0,vector_finite(rhsp) .and. &
                   vector_finite(rhsm))
              jac_ok=merge(1,0,matrix_finite(jac0) .and. &
                   matrix_finite(jacp) .and. matrix_finite(jacm))
              svars_ok=merge(1,0,physical_svars_finite(hcp) .and. &
                   physical_svars_finite(hcm))
              branch_crossed=merge(1,0,cp/=base_branch .or. &
                   cm/=base_branch .or. cp/=cm)

              write(outunit,'(*(g0,:,","))') trim(run_id),trim(sid), &
                   block_names(b),trim(did),hvals(k),-1,trim(bp),trim(bm), &
                   norm_fd,norm_kv,err_l2,rel_l2,err_inf,rel_inf,cosine, &
                   rhs_ok,jac_ok,svars_ok,replay_ok(s),branch_crossed
              if (rhs_ok/=1 .or. jac_ok/=1 .or. svars_ok/=1) then
                 close(outunit)
                 stop 33
              end if
           end do
        end do
     end do
  end do
  close(outunit)
  call write_saturation_subpatch(trim(subpatchfile),trim(run_id),props, &
       coords,ubase,hcomm)
  call write_full_matrix_fd(trim(fullfile),trim(run_id),props,ubase, &
       hcomm,dir_u,dir_p,dir_s,replay_ok)
  call append_post_full_replay(trim(replayfile),trim(run_id),props, &
       ubase,hcomm,rhs_ref,jac_ref,svars_ref)

contains

  subroutine initialize_states(p,c,u,h)
    real(8), intent(out) :: p(nprops,nstate),c(mcrd,nnode)
    real(8), intent(out) :: u(ndofel,nstate),h(nsvars,nstate)
    real(8) :: x,y,pwold,snold,sw,pc,peq,se,ds
    integer :: s,n,gp,off

    c(:,1)=(/0.0d0,0.0d0/)
    c(:,2)=(/1.0d0,0.0d0/)
    c(:,3)=(/1.0d0,1.0d0/)
    c(:,4)=(/0.0d0,1.0d0/)
    p=0.0d0
    do s=1,nstate
       p(1,s)=2.5d10; p(2,s)=0.3d0; p(3,s)=1.0d0
       p(4,s)=1.0d-6; p(5,s)=1.0d-3; p(6,s)=700.0d0
       p(7,s)=1.0d-3; p(8,s)=0.20d0; p(9,s)=0.05d0
       p(10,s)=1.0d4; p(11,s)=2.0d0; p(12,s)=2.0d0
       p(13,s)=2.0d0; p(14,s)=0.01d0; p(15,s)=0.10d0
       p(16,s)=6.0d-12; p(17,s)=0.0d0; p(18,s)=0.0d0
       p(19,s)=0.0d0; p(20,s)=0.0d0
       p(21,s)=1.0d-20; p(22,s)=1.0d-8
    end do
    p(18,3:nstate)=1.0d-8

    u=0.0d0
    do n=1,nnode
       x=c(1,n); y=c(2,n); off=4*(n-1)
       u(off+1,1)=1.0d-6*x
       u(off+2,1)=-5.0d-7*y
       u(off+3,1)=1.2d7+1.0d3*x+5.0d2*y
       u(off+4,1)=0.25d0+1.0d-3*x+5.0d-4*y

       u(off+1,2)=1.0d-5*x+2.0d-6*y
       u(off+2,2)=-8.0d-6*y
       u(off+3,2)=1.2d7+2.0d5*x-1.0d5*y
       u(off+4,2)=0.30d0+5.0d-3*x

       u(off+1,3)=5.0d-4*x+1.0d-4*y
       u(off+2,3)=-3.0d-4*y+5.0d-5*x
       u(off+3,3)=1.2d7+2.5d5*x-1.5d5*y
       u(off+4,3)=0.25d0+0.12d0*x+0.04d0*y
    end do
    do s=4,nstate
       u(:,s)=u(:,3)
    end do
    do n=1,nnode
       off=4*(n-1)
       u(off+4,4)=0.30d0
       u(off+4,5)=0.25d0
       u(off+4,6)=0.30d0
       u(off+4,7)=0.20d0
       u(off+4,8)=0.50d0
       u(off+4,9)=0.30d0
       u(off+4,10)=-0.10d0
       u(off+4,11)=1.10d0
       u(off+4,12)=0.95d0
       u(off+4,13)=0.01d0
    end do

    h=0.0d0
    do s=1,nstate
       select case(s)
       case(1); pwold=1.1999d7; snold=0.249d0
       case(2); pwold=1.1900d7; snold=0.285d0
       case(4); pwold=1.1800d7; snold=0.220d0
       case(5); pwold=1.1800d7; snold=0.245d0
       case(6); pwold=1.1800d7; snold=0.200d0
       case(7); pwold=1.1800d7; snold=0.300d0
       case(8); pwold=1.1800d7; snold=0.100d0
       case(9); pwold=1.1800d7; snold=0.700d0
       case(10); pwold=1.1800d7; snold=0.200d0
       case(11); pwold=1.1800d7; snold=0.800d0
       case(12); pwold=1.1800d7; snold=0.900d0
       case(13); pwold=1.1800d7; snold=0.020d0
       case default; pwold=1.1800d7; snold=0.220d0
       end select
       sw=1.0d0-snold
       ds=1.0d0-p(8,s)-p(9,s)
       se=max(1.0d-6,min(1.0d0,(sw-p(8,s))/ds))
       pc=p(10,s)*(se**(-1.0d0/p(11,s))-1.0d0)
       peq=pwold+(1.0d0-se)*pc
       do gp=1,4
          off=48*(gp-1)
          h(off+1,s)=pwold; h(off+2,s)=snold; h(off+3,s)=sw
          h(off+4,s)=pc; h(off+5,s)=peq; h(off+6,s)=0.0d0
          h(off+7,s)=0.15d0
          h(off+8,s)=max(p(21,s),min(p(16,s),p(22,s)))
          h(off+9,s)=1.0d-9; h(off+10:off+12,s)=0.0d0
          h(off+13,s)=peq-pwold; h(off+14,s)=1.0d0
       end do
    end do
  end subroutine initialize_states

  subroutine initialize_directions(du,dp,ds)
    real(8), intent(out) :: du(8,ndir),dp(4,ndir),ds(4,ndir)
    integer :: i,nseed
    integer, allocatable :: seed(:)
    real(8) :: r8(8),r4a(4),r4b(4)
    do i=1,8
       du(i,1)=1.0d0; du(i,2)=merge(1.0d0,-1.0d0,mod(i,2)==1)
    end do
    do i=1,4
       dp(i,1)=1.0d0; ds(i,1)=1.0d0
       dp(i,2)=merge(1.0d0,-1.0d0,mod(i,2)==1)
       ds(i,2)=merge(1.0d0,-1.0d0,mod(i,2)==1)
    end do
    call normalize(du(:,1)); call normalize(du(:,2))
    call normalize(dp(:,1)); call normalize(dp(:,2))
    call normalize(ds(:,1)); call normalize(ds(:,2))
    call random_seed(size=nseed); allocate(seed(nseed))
    do i=1,nseed
       seed(i)=20260723+104729*i
    end do
    call random_seed(put=seed)
    call random_number(r8); call random_number(r4a); call random_number(r4b)
    du(:,3)=2.0d0*r8-1.0d0
    dp(:,3)=2.0d0*r4a-1.0d0
    ds(:,3)=2.0d0*r4b-1.0d0
    call normalize(du(:,3)); call normalize(dp(:,3)); call normalize(ds(:,3))
    deallocate(seed)
  end subroutine initialize_directions

  subroutine normalize(vv)
    real(8), intent(inout) :: vv(:)
    real(8) :: nn
    nn=sqrt(sum(vv*vv))
    if (nn<=0.0d0) stop 34
    vv=vv/nn
  end subroutine normalize

  subroutine state_id(s,sname)
    integer, intent(in) :: s
    character(len=*), intent(out) :: sname
    select case(s)
    case(1:3); write(sname,'("I",i1)') s-1
    case(4); sname='LIMITER_OFF'
    case(5); sname='LIMITER_INACT'
    case(6); sname='RATIO_POS'
    case(7); sname='RATIO_NEG'
    case(8); sname='BETA_MIN_POS'
    case(9); sname='BETA_MIN_NEG'
    case(10); sname='FINAL_CLIP_LOW'
    case(11); sname='FINAL_CLIP_HIGH'
    case(12); sname='SE_CLIP_LOW'
    case default; sname='SE_CLIP_HIGH'
    end select
  end subroutine state_id

  integer function expected_branch(s)
    integer, intent(in) :: s
    select case(s)
    case(1:4); expected_branch=0
    case(5); expected_branch=1
    case(6); expected_branch=2
    case(7); expected_branch=3
    case(8); expected_branch=4
    case(9); expected_branch=5
    case(10); expected_branch=39
    case(11); expected_branch=24
    case(12); expected_branch=16
    case default; expected_branch=32
    end select
  end function expected_branch

  subroutine block_map(b,rr,cc,nr,nc,xscale)
    integer, intent(in) :: b
    integer, intent(out) :: rr(8),cc(8),nr,nc
    real(8), intent(out) :: xscale
    rr=0; cc=0
    select case(b)
    case(1); rr(1:8)=iu; cc(1:8)=iu; nr=8; nc=8; xscale=scale_u
    case(2); rr(1:8)=iu; cc(1:4)=ip; nr=8; nc=4; xscale=scale_p
    case(3); rr(1:8)=iu; cc(1:4)=isat; nr=8; nc=4; xscale=scale_s
    case(4); rr(1:4)=ip; cc(1:8)=iu; nr=4; nc=8; xscale=scale_u
    case(5); rr(1:4)=ip; cc(1:4)=ip; nr=4; nc=4; xscale=scale_p
    case(6); rr(1:4)=ip; cc(1:4)=isat; nr=4; nc=4; xscale=scale_s
    case(7); rr(1:4)=isat; cc(1:8)=iu; nr=4; nc=8; xscale=scale_u
    case(8); rr(1:4)=isat; cc(1:4)=ip; nr=4; nc=4; xscale=scale_p
    case default; rr(1:4)=isat; cc(1:4)=isat
       nr=4; nc=4; xscale=scale_s
    end select
  end subroutine block_map

  subroutine get_direction(b,d,du,dp,ds,vv,nc)
    integer, intent(in) :: b,d,nc
    real(8), intent(in) :: du(8,ndir),dp(4,ndir),ds(4,ndir)
    real(8), intent(out) :: vv(8)
    vv=0.0d0
    if (b==1 .or. b==4 .or. b==7) then
       vv(1:nc)=du(:,d)
    else if (b==2 .or. b==5 .or. b==8) then
       vv(1:nc)=dp(:,d)
    else
       vv(1:nc)=ds(:,d)
    end if
  end subroutine get_direction

  subroutine evaluate_uel(s,uin,uref,hin,rout,jout,hout,pnewdt_out)
    integer, intent(in) :: s
    real(8), intent(in) :: uin(ndofel),uref(ndofel),hin(nsvars)
    real(8), intent(out) :: rout(ndofel),jout(ndofel,ndofel)
    real(8), intent(out) :: hout(nsvars),pnewdt_out
    integer, parameter :: mlvarx=ndofel,mdload=1,npredf=1
    real(8) :: rhs(mlvarx,nrhs),amatrx(ndofel,ndofel)
    real(8) :: energy(8),ulocal(ndofel),duin(mlvarx,nrhs)
    real(8) :: v(ndofel),a(ndofel),time(2),params(3)
    real(8) :: adlmag(mdload,1),predef(2,npredf,nnode)
    real(8) :: ddlmag(mdload,1),pnewdt,period
    integer :: jdltyp(mdload,1),lflags(10),jprops(1)

    rhs=0.0d0; amatrx=0.0d0; hout=hin; energy=0.0d0
    ulocal=uin; duin=0.0d0; duin(:,1)=uin-uref
    v=0.0d0; a=0.0d0; time=(/600.0d0,1260.0d0/)
    params=0.0d0; adlmag=0.0d0; predef=0.0d0; ddlmag=0.0d0
    jdltyp=0; lflags=0; lflags(1)=1; lflags(3)=1; jprops=0
    pnewdt=1.0d6; period=600.0d0

    call configure_state_environment(s)
    call UEL(rhs,amatrx,hout,energy,ndofel,nrhs,nsvars, &
         props(:,s),nprops,coords,mcrd,nnode,ulocal,duin,v,a,1,time, &
         60.0d0,3,10,10287,params,0,jdltyp,adlmag,predef,npredf, &
         lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,1,period)
    rout=rhs(:,1); jout=amatrx; pnewdt_out=pnewdt
  end subroutine evaluate_uel

  subroutine configure_state_environment(s)
    integer, intent(in) :: s
    real(8) :: reg,tol,bmin
    call set_env_value('DBG_FIELD_CAPACITY_REG','0')
    call set_env_value('DBG_FIELD_ANCHOR_REG','0')
    call set_env_value('DBG_REG_S_DIAG','0')
    call set_env_value('DBG_PC_SLOPE_CAP_REG','0')
    call set_env_value('DBG_SN_BOUNDS','0')
    call raw_controls(s,reg,tol,bmin)
    if (reg>0.5d0) then
       call set_env_value('DBG_RAW_SN_JUMP_REG','1')
    else
       call set_env_value('DBG_RAW_SN_JUMP_REG','0')
    end if
    if (s==8 .or. s==9) then
       call set_env_value('DBG_RAW_SN_JUMP_TOL','1.0E-2')
       call set_env_value('DBG_RAW_SN_JUMP_BETA_MIN','1.0E-1')
    else
       call set_env_value('DBG_RAW_SN_JUMP_TOL','5.0E-2')
       call set_env_value('DBG_RAW_SN_JUMP_BETA_MIN','5.0E-2')
    end if
  end subroutine configure_state_environment

  subroutine raw_controls(s,reg,tol,bmin)
    integer, intent(in) :: s
    real(8), intent(out) :: reg,tol,bmin
    reg=0.0d0; tol=5.0d-2; bmin=5.0d-2
    if (s>=5 .and. s<=9) reg=1.0d0
    if (s==8 .or. s==9) then
       tol=1.0d-2
       bmin=1.0d-1
    end if
  end subroutine raw_controls

  subroutine set_env_value(name,value)
    character(len=*), intent(in) :: name,value
    logical(4) :: ok
    ok=setenvqq(trim(name)//'='//trim(value))
    if (.not.ok) stop 61
  end subroutine set_env_value

  subroutine branch_signature(hh,p,code,label)
    real(8), intent(in) :: hh(nsvars),p(nprops)
    integer, intent(out) :: code
    character(len=*), intent(out) :: label
    integer :: gp,off
    real(8) :: se,kabs,ds
    code=0; ds=1.0d0-p(8)-p(9)
    do gp=1,4
       off=48*(gp-1)
       se=(1.0d0-hh(off+2)-p(8))/ds
       kabs=hh(off+8)
       if (se<=1.0d-5) code=ior(code,1)
       if (se>=1.0d0-1.0d-5) code=ior(code,2)
       if (kabs<=p(21)*(1.0d0+1.0d-8)) code=ior(code,4)
       if (kabs>=p(22)*(1.0d0-1.0d-8)) code=ior(code,8)
    end do
    if (code==0) then
       label='SMOOTH'
    else
       write(label,'("BRANCH_",i0)') code
    end if
  end subroutine branch_signature

  subroutine fd_branch_signature(s,uin,hin,p,code,label)
    integer, intent(in) :: s
    real(8), intent(in) :: uin(ndofel),hin(nsvars),p(nprops)
    integer, intent(out) :: code
    character(len=*), intent(out) :: label
    real(8) :: sn,snold,delta,ad,reg,tol,bmin,beta,snuse
    real(8) :: ds,se,seps
    integer :: gp,off,rawbranch,pcbranch
    sn=sum(uin(isat))/4.0d0
    snold=0.0d0
    do gp=1,4
       off=48*(gp-1)
       snold=snold+hin(off+2)/4.0d0
    end do
    call raw_controls(s,reg,tol,bmin)
    delta=sn-snold; ad=abs(delta); beta=1.0d0; rawbranch=0
    if (reg>0.5d0) then
       rawbranch=1
       if (ad>tol) then
          beta=tol/(ad+1.0d-30)
          rawbranch=merge(2,3,delta>0.0d0)
          if (beta<bmin) then
             beta=bmin
             rawbranch=merge(4,5,delta>0.0d0)
          end if
          if (beta>1.0d0) then
             beta=1.0d0
             rawbranch=6
          end if
       end if
    end if
    snuse=snold+beta*delta
    if (snuse<0.0d0) rawbranch=7
    if (snuse>1.0d0) rawbranch=8
    ds=max(1.0d-12,1.0d0-p(8)-p(9)); seps=1.0d-6
    se=(1.0d0-sn-p(8))/ds; pcbranch=0
    if (se<=seps) pcbranch=1
    if (se>=1.0d0-seps) pcbranch=2
    code=rawbranch+16*pcbranch
    write(label,'("I3B_",i0)') code
  end subroutine fd_branch_signature

  subroutine write_direction_manifest(path,du,dp,ds)
    character(len=*), intent(in) :: path
    real(8), intent(in) :: du(8,ndir),dp(4,ndir),ds(4,ndir)
    integer :: unit,d,i
    character(len=2) :: dd
    open(newunit=unit,file=path,status='replace',action='write')
    write(unit,'(a)') 'field_id,direction_id,component,'// &
         'scaled_component,physical_basis_component,characteristic_scale,seed'
    do d=1,ndir
       write(dd,'("D",i1)') d
       do i=1,8
          write(unit,'(*(g0,:,","))') 'u',trim(dd),i,du(i,d), &
               scale_u*du(i,d),scale_u,20260723
       end do
       do i=1,4
          write(unit,'(*(g0,:,","))') 'p_w',trim(dd),i,dp(i,d), &
               scale_p*dp(i,d),scale_p,20260723
          write(unit,'(*(g0,:,","))') 'S_n',trim(dd),i,ds(i,d), &
               scale_s*ds(i,d),scale_s,20260723
       end do
    end do
    close(unit)
  end subroutine write_direction_manifest

  subroutine write_state_manifest(path,p,c,u,h)
    character(len=*), intent(in) :: path
    real(8), intent(in) :: p(nprops,nstate),c(mcrd,nnode)
    real(8), intent(in) :: u(ndofel,nstate),h(nsvars,nstate)
    integer :: unit,s,i,j,code
    character(len=16) :: ss,label
    open(newunit=unit,file=path,status='replace',action='write')
    write(unit,'(a)') 'state_id,record_type,index1,index2,name,value,units'
    do s=1,nstate
       call state_id(s,ss)
       do i=1,nprops
          write(unit,'(*(g0,:,","))') trim(ss),'PROPS',i,0,'PROP',p(i,s),'-'
       end do
       do j=1,nnode
          do i=1,mcrd
             write(unit,'(*(g0,:,","))') trim(ss),'COORDS',i,j,'COORD',c(i,j),'m'
          end do
       end do
       do i=1,ndofel
          write(unit,'(*(g0,:,","))') trim(ss),'U',i,0,'TRIAL_BASE',u(i,s),'-'
          write(unit,'(*(g0,:,","))') trim(ss),'DU',i,0,'BASE_INCREMENT',0.0d0,'-'
       end do
       do i=1,nsvars
          write(unit,'(*(g0,:,","))') trim(ss),'SVARS',i,0,'COMMITTED',h(i,s),'-'
       end do
       write(unit,'(*(g0,:,","))') trim(ss),'TIME',1,0,'STEP_TIME',600.0d0,'s'
       write(unit,'(*(g0,:,","))') trim(ss),'TIME',2,0,'TOTAL_TIME',1260.0d0,'s'
       write(unit,'(*(g0,:,","))') trim(ss),'DTIME',0,0,'DTIME',60.0d0,'s'
       write(unit,'(*(g0,:,","))') trim(ss),'SCALE',1,0,'X_U',scale_u,'m'
       write(unit,'(*(g0,:,","))') trim(ss),'SCALE',2,0,'X_P',scale_p,'Pa'
       write(unit,'(*(g0,:,","))') trim(ss),'SCALE',3,0,'X_S',scale_s,'-'
       call evaluate_uel(s,u(:,s),u(:,s),h(:,s),rhs0,jac0,hc0,pnewdt)
       call branch_signature(hc0,p(:,s),code,label)
       write(unit,'(*(g0,:,","))') trim(ss),'BRANCH',0,0,trim(label),dble(code),'-'
    end do
    close(unit)
  end subroutine write_state_manifest

  subroutine write_snapshot(path,run_label,p,u,h)
    character(len=*), intent(in) :: path,run_label
    real(8), intent(in) :: p(nprops,nstate),u(ndofel,nstate)
    real(8), intent(in) :: h(nsvars,nstate)
    integer :: unit,s,i,j,gp,off,code
    character(len=16) :: ss,label
    open(newunit=unit,file=path,status='replace',action='write')
    write(unit,'(a)') 'run_id,state_id,record_type,index1,index2,value'
    do s=1,nstate
       call state_id(s,ss)
       call evaluate_uel(s,u(:,s),u(:,s),h(:,s),rhs0,jac0,hc0,pnewdt)
       call branch_signature(hc0,p(:,s),code,label)
       do i=1,ndofel
          write(unit,'(*(g0,:,","))') trim(run_label),trim(ss), &
               'RHS',i,0,rhs0(i)
          do j=1,ndofel
             write(unit,'(*(g0,:,","))') trim(run_label),trim(ss), &
                  'AMATRX',i,j,jac0(i,j)
          end do
       end do
       do gp=1,4
          off=48*(gp-1)
          do i=1,14
             write(unit,'(*(g0,:,","))') trim(run_label),trim(ss), &
                  'PHYSICAL_SVARS',gp,i,hc0(off+i)
          end do
          write(unit,'(*(g0,:,","))') trim(run_label),trim(ss), &
               'CLOSURE_PW',gp,0,hc0(off+1)
          write(unit,'(*(g0,:,","))') trim(run_label),trim(ss), &
               'CLOSURE_SN',gp,0,hc0(off+2)
          write(unit,'(*(g0,:,","))') trim(run_label),trim(ss), &
               'CLOSURE_PC',gp,0,hc0(off+4)
          write(unit,'(*(g0,:,","))') trim(run_label),trim(ss), &
               'CLOSURE_PEQ',gp,0,hc0(off+5)
          write(unit,'(*(g0,:,","))') trim(run_label),trim(ss), &
               'KABS',gp,0,hc0(off+8)
       end do
       write(unit,'(*(g0,:,","))') trim(run_label),trim(ss), &
            'BRANCH',0,0,dble(code)
    end do
    close(unit)
  end subroutine write_snapshot

  subroutine write_branch_manifest(path,run_label,p,u,h)
    character(len=*), intent(in) :: path,run_label
    real(8), intent(in) :: p(nprops,nstate),u(ndofel,nstate)
    real(8), intent(in) :: h(nsvars,nstate)
    integer :: unit,s,gp,off,code,sn_branch,sn_switch
    integer :: pc_branch,pc_switch
    real(8) :: ds,se_raw,kabs,se_low_dist,se_high_dist
    real(8) :: k_low_dist,k_high_dist
    real(8) :: sn_gp_mean,sn_old_mean,sn_use,beta,dsnuse
    real(8) :: sn_dist,d2pc,pc_dist,reg,tol,bmin
    character(len=16) :: ss,label
    open(newunit=unit,file=path,status='replace',action='write')
    write(unit,'(a)') 'run_id,state_id,gp,branch_id,branch_label,'// &
         'se_raw,se_low_distance,se_high_distance,kabs,'// &
         'k_lower_distance,k_upper_distance,sn_use_branch,'// &
         'sn_branch_distance,sn_switch,pc2_branch,'// &
         'pc_branch_distance,pc_switch'
    do s=1,nstate
       call state_id(s,ss)
       call evaluate_uel(s,u(:,s),u(:,s),h(:,s),rhs0,jac0,hc0,pnewdt)
       call branch_signature(hc0,p(:,s),code,label)
       sn_gp_mean=sum(u(isat,s))/4.0d0
       sn_old_mean=0.0d0
       do gp=1,4
          off=48*(gp-1)
          sn_old_mean=sn_old_mean+h(off+2,s)/4.0d0
       end do
       call raw_controls(s,reg,tol,bmin)
#ifdef WP10_I3_BUILD
       call UPS_EVAL_SN_RAW_USE_JET(sn_gp_mean,sn_old_mean,reg,tol, &
            bmin,sn_use,beta,dsnuse,sn_branch,sn_dist,sn_switch)
       call UPS_PC_SECOND_DERIVATIVE(sn_gp_mean,p(:,s),d2pc, &
            pc_branch,pc_dist,pc_switch)
#else
       sn_use=sn_gp_mean; beta=1.0d0; dsnuse=1.0d0
       sn_branch=-999; sn_dist=-1.0d0; sn_switch=0
       d2pc=0.0d0; pc_branch=-999; pc_dist=-1.0d0; pc_switch=0
#endif
       ds=max(1.0d-12,1.0d0-p(8,s)-p(9,s))
       do gp=1,4
          off=48*(gp-1)
          se_raw=(1.0d0-hc0(off+2)-p(8,s))/ds
          kabs=hc0(off+8)
          se_low_dist=se_raw-1.0d-6
          se_high_dist=(1.0d0-1.0d-6)-se_raw
          k_low_dist=kabs-p(21,s)
          k_high_dist=p(22,s)-kabs
          write(unit,'(*(g0,:,","))') trim(run_label),trim(ss),gp, &
                code,trim(label),se_raw,se_low_dist,se_high_dist,kabs, &
                k_low_dist,k_high_dist,sn_branch,sn_dist,sn_switch, &
                pc_branch,pc_dist,pc_switch
       end do
    end do
    close(unit)
  end subroutine write_branch_manifest

  subroutine write_term_audit(path,run_label,p,u,h)
    character(len=*), intent(in) :: path,run_label
    real(8), intent(in) :: p(nprops,nstate),u(ndofel,nstate)
    real(8), intent(in) :: h(nsvars,nstate)
    real(8) :: kt(4,4,8),ksum(4,4),kret(4,4)
    real(8) :: frob,maxabs,signedsum,maxdiff
    integer :: unit,s,t,i,j,finite_ok
    character(len=16) :: ss
    character(len=20), parameter :: tname(8)=(/ &
         'KSS_STORAGE         ','KSS_STRAIN_SAT      ', &
         'KSS_MOBILITY        ','KSS_PC_FIRST        ', &
         'KSS_PC_SECOND       ','KSS_PAIRED_CAPACITY ', &
         'KSS_PAIRED_ANCHOR   ','KSS_PAIRED_DIAG     '/)
    open(newunit=unit,file=path,status='replace',action='write')
    write(unit,'(a)') 'run_id,state_id,term_id,frobenius_norm,'// &
         'max_abs,signed_sum,units,finite,sum_to_returned_max_abs'
    do s=1,nstate
       call state_id(s,ss)
       call evaluate_uel(s,u(:,s),u(:,s),h(:,s),rhs0,jac0,hc0,pnewdt)
       kt=0.0d0
#ifdef WP10_I3_BUILD
       call UPS_WP10I3_GET_KSS_AUDIT(kt(:,:,1),kt(:,:,2),kt(:,:,3), &
            kt(:,:,4),kt(:,:,5),kt(:,:,6),kt(:,:,7),kt(:,:,8))
#endif
       ksum=0.0d0
       do t=1,8
          ksum=ksum+kt(:,:,t)
          frob=sqrt(sum(kt(:,:,t)*kt(:,:,t)))
          maxabs=maxval(abs(kt(:,:,t)))
          signedsum=sum(kt(:,:,t))
          finite_ok=merge(1,0,matrix_finite(kt(:,:,t)))
          write(unit,'(*(g0,:,","))') trim(run_label),trim(ss), &
               trim(tname(t)),frob,maxabs,signedsum,'RS_per_S', &
               finite_ok,0.0d0
       end do
        do i=1,4
           do j=1,4
              kret(i,j)=jac0(isat(i),isat(j))
           end do
           kret(i,i)=kret(i,i)-kstab_num
        end do
       maxdiff=maxval(abs(ksum-kret))
       write(unit,'(*(g0,:,","))') trim(run_label),trim(ss), &
            'KSS_TOTAL',sqrt(sum(ksum*ksum)),maxval(abs(ksum)), &
            sum(ksum),'RS_per_S',merge(1,0,matrix_finite(ksum)),maxdiff
    end do
    close(unit)
  end subroutine write_term_audit

  subroutine write_switch_diagnostics(path,run_label,p)
    character(len=*), intent(in) :: path,run_label
    real(8), intent(in) :: p(nprops,nstate)
    real(8) :: sn,snold,reg,tol,bmin,snuse,beta,dsn,dist
    real(8) :: d2pc,pcdist,ds,seps
    integer :: unit,k,bid,sw,pcbid,pcsw
    character(len=20) :: name
    open(newunit=unit,file=path,status='replace',action='write')
    write(unit,'(a)') 'run_id,diagnostic_id,sn,sn_old,reg,tol,'// &
         'beta_min,sn_use,beta,dsn_use_dsn,sn_branch,sn_distance,'// &
         'sn_switch,d2pc_dsn2,pc_branch,pc_distance,pc_switch'
    ds=1.0d0-p(8,1)-p(9,1)
    seps=1.0d-6
    do k=1,6
       reg=0.0d0; tol=5.0d-2; bmin=1.0d-1
       select case(k)
       case(1)
          name='ABS_DELTA_EQ_TOL'; reg=1.0d0
          snold=0.20d0; sn=snold+tol
       case(2)
          name='BETA_EQ_BETA_MIN'; reg=1.0d0
          snold=0.20d0; sn=snold+tol/bmin
       case(3)
          name='SN_USE_EQ_ZERO'; snold=0.20d0; sn=0.0d0
       case(4)
          name='SN_USE_EQ_ONE'; snold=0.80d0; sn=1.0d0
       case(5)
          name='SE_EQ_LOW'
          snold=0.70d0; sn=1.0d0-p(8,1)-ds*seps
       case default
          name='SE_EQ_HIGH'
          snold=0.10d0; sn=1.0d0-p(8,1)-ds*(1.0d0-seps)
       end select
#ifdef WP10_I3_BUILD
       call UPS_EVAL_SN_RAW_USE_JET(sn,snold,reg,tol,bmin,snuse, &
            beta,dsn,bid,dist,sw)
       call UPS_PC_SECOND_DERIVATIVE(sn,p(:,1),d2pc,pcbid,pcdist,pcsw)
#else
       snuse=sn; beta=1.0d0; dsn=1.0d0
       bid=-999; dist=-1.0d0; sw=0
       d2pc=0.0d0; pcbid=-999; pcdist=-1.0d0; pcsw=0
#endif
       write(unit,'(*(g0,:,","))') trim(run_label),trim(name),sn, &
            snold,reg,tol,bmin,snuse,beta,dsn,bid,dist,sw,d2pc, &
            pcbid,pcdist,pcsw
    end do
    close(unit)
  end subroutine write_switch_diagnostics

  subroutine write_saturation_subpatch(path,run_label,p,c,u,h)
    character(len=*), intent(in) :: path,run_label
    real(8), intent(in) :: p(nprops,nstate),c(mcrd,nnode)
    real(8), intent(in) :: u(ndofel,nstate),h(nsvars,nstate)
    real(8), parameter :: hpatch=1.0d-4
    real(8) :: dphys(ndofel),up(ndofel),um(ndofel)
    real(8) :: rp(ndofel),rm(ndofel),jp(ndofel,ndofel)
    real(8) :: jphys(ndofel,ndofel)
    real(8) :: hp(nsvars),hm(nsvars),fdv(4),kvv(4),pnd
    real(8) :: nfd,nkv,el2,rl2,einf,rinf,cosv,tmp(4)
    integer :: unit,m,j,n,off,finite_ok
    character(len=20) :: mode_id
    call evaluate_uel(3,u(:,3),u(:,3),h(:,3),rhs0,jac0,hc0,pnewdt)
    call physical_jacobian(jac0,jphys)
    open(newunit=unit,file=path,status='replace',action='write')
    write(unit,'(a)') 'run_id,state_id,mode_id,h,norm_fd,norm_kv,'// &
         'absolute_error_l2,relative_error_l2,absolute_error_inf,'// &
         'relative_error_inf,cosine_similarity,finite'
    do m=1,2
       dphys=0.0d0
       if (m==1) then
          mode_id='S_UNIFORM'
          tmp=1.0d0
       else
          mode_id='S_AFFINE'
          do n=1,nnode
             tmp(n)=c(1,n)+c(2,n)-1.0d0
          end do
       end if
       tmp=scale_s*tmp/sqrt(sum(tmp*tmp))
       do n=1,nnode
          off=4*(n-1)
          dphys(off+4)=tmp(n)
       end do
       up=u(:,3)+hpatch*dphys
       um=u(:,3)-hpatch*dphys
       call evaluate_uel(3,up,u(:,3),h(:,3),rp,jp,hp,pnd)
       call evaluate_uel(3,um,u(:,3),h(:,3),rm,jp,hm,pnd)
       do j=1,4
          fdv(j)=(rp(isat(j))-rm(isat(j)))/(2.0d0*hpatch)
          kvv(j)=sigma_global*sum(jphys(isat(j),:)*dphys)
       end do
       nfd=l2norm(fdv); nkv=l2norm(kvv); el2=l2norm(fdv-kvv)
       rl2=el2/max(nfd,nkv,1.0d-30)
       einf=maxval(abs(fdv-kvv))
       rinf=einf/max(maxval(abs(fdv)),maxval(abs(kvv)),1.0d-30)
       if (nfd>1.0d-300 .and. nkv>1.0d-300) then
          cosv=sum(fdv*kvv)/(nfd*nkv)
       else if (nfd<=1.0d-300 .and. nkv<=1.0d-300) then
          cosv=1.0d0
       else
          cosv=0.0d0
       end if
       finite_ok=merge(1,0,vector_finite(rp) .and. vector_finite(rm) &
            .and. matrix_finite(jac0) .and. physical_svars_finite(hp) &
            .and. physical_svars_finite(hm))
       write(unit,'(*(g0,:,","))') trim(run_label),'I2',trim(mode_id), &
            hpatch,nfd,nkv,el2,rl2,einf,rinf,cosv,finite_ok
       if (finite_ok/=1) then
          close(unit)
          stop 62
       end if
    end do
    close(unit)
  end subroutine write_saturation_subpatch

  subroutine write_affine_patch(path,run_label,p,c,u,h)
    character(len=*), intent(in) :: path,run_label
    real(8), intent(in) :: p(nprops,nstate),c(mcrd,nnode)
    real(8), intent(in) :: u(ndofel,nstate),h(nsvars,nstate)
    integer, parameter :: nmode=8
    real(8), parameter :: hpatch=1.0d-4
    real(8) :: dphys(ndofel),up(ndofel),um(ndofel)
    real(8) :: rp(ndofel),rm(ndofel),jp(ndofel,ndofel)
    real(8) :: hp(nsvars),hm(nsvars),fdru(8),kvru(8)
    real(8) :: nfd,nkv,el2,rl2,einf,rinf,cosv,pnd
    integer :: unit,m,j,cp,cm,base_code,finite_ok,branch_crossed
    character(len=16) :: mode_id,field_id,bp,bm,base_label
    call evaluate_uel(2,u(:,2),u(:,2),h(:,2),rhs0,jac0,hc0,pnewdt)
    call branch_signature(hc0,p(:,2),base_code,base_label)
    open(newunit=unit,file=path,status='replace',action='write')
    write(unit,'(a)') 'run_id,state_id,mode_id,field_id,h,'// &
         'norm_fd,norm_kv,absolute_error_l2,relative_error_l2,'// &
         'absolute_error_inf,relative_error_inf,cosine_similarity,'// &
         'branch_id_plus,branch_id_minus,branch_crossed,finite'
    do m=1,nmode
       call affine_direction(m,c,dphys,mode_id,field_id)
       up=u(:,2)+hpatch*dphys
       um=u(:,2)-hpatch*dphys
       call evaluate_uel(2,up,u(:,2),h(:,2),rp,jp,hp,pnd)
       call branch_signature(hp,p(:,2),cp,bp)
       call evaluate_uel(2,um,u(:,2),h(:,2),rm,jp,hm,pnd)
       call branch_signature(hm,p(:,2),cm,bm)
       do j=1,8
          fdru(j)=(rp(iu(j))-rm(iu(j)))/(2.0d0*hpatch)
          kvru(j)=sigma_global*sum(jac0(iu(j),:)*dphys)
       end do
       nfd=l2norm(fdru); nkv=l2norm(kvru)
       el2=l2norm(fdru-kvru)
       rl2=el2/max(nfd,nkv,1.0d-30)
       einf=maxval(abs(fdru-kvru))
       rinf=einf/max(maxval(abs(fdru)),maxval(abs(kvru)),1.0d-30)
       if (nfd>1.0d-300 .and. nkv>1.0d-300) then
          cosv=sum(fdru*kvru)/(nfd*nkv)
       else if (nfd<=1.0d-300 .and. nkv<=1.0d-300) then
          cosv=1.0d0
       else
          cosv=0.0d0
       end if
       branch_crossed=merge(1,0,cp/=base_code .or. cm/=base_code .or. cp/=cm)
       finite_ok=merge(1,0,vector_finite(rp) .and. vector_finite(rm) &
            .and. matrix_finite(jac0) .and. physical_svars_finite(hp) &
            .and. physical_svars_finite(hm))
       write(unit,'(*(g0,:,","))') trim(run_label),'I1',trim(mode_id), &
            trim(field_id),hpatch,nfd,nkv,el2,rl2,einf,rinf,cosv, &
            trim(bp),trim(bm),branch_crossed,finite_ok
       if (finite_ok/=1) then
          close(unit)
          stop 41
       end if
    end do
    close(unit)
  end subroutine write_affine_patch

  subroutine affine_direction(mode,c,dv,mode_id,field_id)
    integer, intent(in) :: mode
    real(8), intent(in) :: c(mcrd,nnode)
    real(8), intent(out) :: dv(ndofel)
    character(len=*), intent(out) :: mode_id,field_id
    integer :: n,off
    real(8) :: x,y
    dv=0.0d0; field_id='u'
    select case(mode)
    case(1); mode_id='TX'
       do n=1,nnode; off=4*(n-1); dv(off+1)=scale_u; end do
    case(2); mode_id='TY'
       do n=1,nnode; off=4*(n-1); dv(off+2)=scale_u; end do
    case(3); mode_id='ROT'
       do n=1,nnode
          off=4*(n-1); x=c(1,n); y=c(2,n)
          dv(off+1)=-scale_u*y; dv(off+2)=scale_u*x
       end do
    case(4); mode_id='EXX'
       do n=1,nnode; off=4*(n-1); dv(off+1)=scale_u*c(1,n); end do
    case(5); mode_id='EYY'
       do n=1,nnode; off=4*(n-1); dv(off+2)=scale_u*c(2,n); end do
    case(6); mode_id='GXY'
       do n=1,nnode
          off=4*(n-1); x=c(1,n); y=c(2,n)
          dv(off+1)=0.5d0*scale_u*y
          dv(off+2)=0.5d0*scale_u*x
       end do
    case(7); mode_id='P_UNIFORM'; field_id='p_w'
       do n=1,nnode; off=4*(n-1); dv(off+3)=scale_p; end do
    case default; mode_id='S_UNIFORM'; field_id='S_n'
       do n=1,nnode; off=4*(n-1); dv(off+4)=scale_s; end do
    end select
  end subroutine affine_direction


  subroutine write_mixed_direction_manifest(path,du,dp,ds)
    character(len=*), intent(in) :: path
    real(8), intent(in) :: du(8,ndir),dp(4,ndir),ds(4,ndir)
    real(8) :: dscaled(ndofel),dphys(ndofel)
    integer :: unit,m,i,local_index
    character(len=2) :: mid,source_id
    character(len=3) :: field_id
    open(newunit=unit,file=path,status='replace',action='write')
    write(unit,'(a)') 'direction_id,component,field_id,field_component,'// &
         'scaled_component,physical_component,characteristic_scale,'// &
         'source_direction,seed'
    do m=1,ndir
       call mixed_direction(m,du,dp,ds,dscaled,dphys)
       write(mid,'("M",i1)') m
       do i=1,ndofel
          call dof_metadata(i,field_id,local_index)
          call mixed_source_id(m,field_id,source_id)
          write(unit,'(*(g0,:,","))') trim(mid),i,trim(field_id), &
               local_index,dscaled(i),dphys(i),dof_scale(i), &
               trim(source_id),20260723
       end do
    end do
    close(unit)
  end subroutine write_mixed_direction_manifest

  subroutine mixed_direction(m,du,dp,ds,dscaled,dphys)
    integer, intent(in) :: m
    real(8), intent(in) :: du(8,ndir),dp(4,ndir),ds(4,ndir)
    real(8), intent(out) :: dscaled(ndofel),dphys(ndofel)
    real(8) :: factor
    integer :: ku,kp,ks,j
    factor=1.0d0/sqrt(3.0d0)
    select case(m)
    case(1); ku=1; kp=2; ks=3
    case(2); ku=2; kp=3; ks=1
    case default; ku=3; kp=1; ks=2
    end select
    dscaled=0.0d0
    do j=1,8
       dscaled(iu(j))=factor*du(j,ku)
    end do
    do j=1,4
       dscaled(ip(j))=factor*dp(j,kp)
       dscaled(isat(j))=factor*ds(j,ks)
    end do
    do j=1,ndofel
       dphys(j)=dscaled(j)*dof_scale(j)
    end do
  end subroutine mixed_direction

  subroutine mixed_source_id(m,field_id,source_id)
    integer, intent(in) :: m
    character(len=*), intent(in) :: field_id
    character(len=*), intent(out) :: source_id
    integer :: idx
    if (field_id=='u') then
       select case(m); case(1); idx=1; case(2); idx=2
       case default; idx=3; end select
    else if (field_id=='p_w') then
       select case(m); case(1); idx=2; case(2); idx=3
       case default; idx=1; end select
    else
       select case(m); case(1); idx=3; case(2); idx=1
       case default; idx=2; end select
    end if
    write(source_id,'("D",i1)') idx
  end subroutine mixed_source_id

  subroutine dof_metadata(i,field_id,local_index)
    integer, intent(in) :: i
    character(len=*), intent(out) :: field_id
    integer, intent(out) :: local_index
    integer :: j
    field_id=''
    local_index=0
    do j=1,8
       if (i==iu(j)) then
          field_id='u'; local_index=j; return
       end if
    end do
    do j=1,4
       if (i==ip(j)) then
          field_id='p_w'; local_index=j; return
       end if
       if (i==isat(j)) then
          field_id='S_n'; local_index=j; return
       end if
    end do
    stop 71
  end subroutine dof_metadata

  real(8) function dof_scale(i)
    integer, intent(in) :: i
    if (any(i==iu)) then
       dof_scale=scale_u
    else if (any(i==ip)) then
       dof_scale=scale_p
    else if (any(i==isat)) then
       dof_scale=scale_s
    else
       stop 72
    end if
  end function dof_scale

  subroutine write_full_matrix_fd(path,run_label,p,u,h,du,dp,ds,replay)
    character(len=*), intent(in) :: path,run_label
    real(8), intent(in) :: p(nprops,nstate),u(ndofel,nstate)
    real(8), intent(in) :: h(nsvars,nstate)
    real(8), intent(in) :: du(8,ndir),dp(4,ndir),ds(4,ndir)
    integer, intent(in) :: replay(nstate)
    real(8) :: dscaled(ndofel),dphys(ndofel),up(ndofel),um(ndofel)
    real(8) :: rp(ndofel),rm(ndofel),jp(ndofel,ndofel),jm(ndofel,ndofel)
    real(8) :: hp(nsvars),hm(nsvars),jphys(ndofel,ndofel)
    real(8) :: fd(ndofel),kv(ndofel),sfd(ndofel),skv(ndofel)
    real(8) :: row_scale(ndofel),err(ndofel)
    real(8) :: nfd,nkv,errl2,rell2,errinf,relinf,cosine,pnd
    integer :: unit,ss,m,k,i,j,cp,cm,base_branch
    integer :: rhsok,jacok,svarsok,crossed
    character(len=16) :: state_name,bp,bm
    character(len=2) :: mid
    character(len=3) :: field_id
    integer :: local_index
    open(newunit=unit,file=path,status='replace',action='write')
    write(unit,'(a)') 'run_id,state_id,direction_id,h,sigma,component,'// &
         'row_field,row_field_component,row_scale,fd_component,'// &
         'kv_component,scaled_fd_component,scaled_kv_component,'// &
         'norm_scaled_fd,norm_scaled_kv,absolute_error_scaled_l2,'// &
         'relative_error_scaled_l2,absolute_error_scaled_inf,'// &
         'relative_error_scaled_inf,cosine_similarity,'// &
         'branch_id_plus,branch_id_minus,rhs_finite,amatrx_finite,'// &
         'candidate_svars_finite,replay_exact,branch_crossed'
    do ss=1,nstate
       call state_id(ss,state_name)
       call evaluate_uel(ss,u(:,ss),u(:,ss),h(:,ss), &
            rhs0,jac0,hc0,pnewdt)
       call physical_jacobian(jac0,jphys)
       call fd_branch_signature(ss,u(:,ss),h(:,ss),p(:,ss), &
            base_branch,bp)
       do i=1,ndofel
          row_scale(i)=0.0d0
          do j=1,ndofel
             row_scale(i)=row_scale(i)+(jphys(i,j)*dof_scale(j))**2
          end do
          row_scale(i)=max(sqrt(row_scale(i)),1.0d-30)
       end do
       do m=1,ndir
          write(mid,'("M",i1)') m
          call mixed_direction(m,du,dp,ds,dscaled,dphys)
          do k=1,nh
             up=u(:,ss)+hvals(k)*dphys
             um=u(:,ss)-hvals(k)*dphys
             call evaluate_uel(ss,up,u(:,ss),h(:,ss),rp,jp,hp,pnd)
             call fd_branch_signature(ss,up,h(:,ss),p(:,ss),cp,bp)
             call evaluate_uel(ss,um,u(:,ss),h(:,ss),rm,jm,hm,pnd)
             call fd_branch_signature(ss,um,h(:,ss),p(:,ss),cm,bm)
             fd=(rp-rm)/(2.0d0*hvals(k))
             kv=sigma_global*matmul(jphys,dphys)
             sfd=fd/row_scale
             skv=kv/row_scale
             err=sfd-skv
             nfd=l2norm(sfd); nkv=l2norm(skv)
             errl2=l2norm(err)
             rell2=errl2/max(nfd,nkv,1.0d-30)
             errinf=maxval(abs(err))
             relinf=errinf/max(maxval(abs(sfd)),maxval(abs(skv)),1.0d-30)
             if (nfd>1.0d-300 .and. nkv>1.0d-300) then
                cosine=sum(sfd*skv)/(nfd*nkv)
             else if (nfd<=1.0d-300 .and. nkv<=1.0d-300) then
                cosine=1.0d0
             else
                cosine=0.0d0
             end if
             rhsok=merge(1,0,vector_finite(rp) .and. vector_finite(rm))
             jacok=merge(1,0,matrix_finite(jac0) .and. &
                  matrix_finite(jp) .and. matrix_finite(jm) .and. &
                  matrix_finite(jphys))
             svarsok=merge(1,0,physical_svars_finite(hp) .and. &
                  physical_svars_finite(hm))
             crossed=merge(1,0,cp/=base_branch .or. cm/=base_branch &
                  .or. cp/=cm)
             do i=1,ndofel
                call dof_metadata(i,field_id,local_index)
                write(unit,'(*(g0,:,","))') trim(run_label), &
                     trim(state_name),trim(mid),hvals(k),-1,i, &
                     trim(field_id),local_index,row_scale(i),fd(i), &
                     kv(i),sfd(i),skv(i),nfd,nkv,errl2,rell2,errinf, &
                     relinf,cosine,trim(bp),trim(bm),rhsok,jacok, &
                     svarsok,replay(ss),crossed
             end do
             if (rhsok/=1 .or. jacok/=1 .or. svarsok/=1) stop 73
          end do
       end do
    end do
    close(unit)
  end subroutine write_full_matrix_fd

  subroutine append_post_full_replay(path,run_label,p,u,h,rref,jref,href)
    character(len=*), intent(in) :: path,run_label
    real(8), intent(in) :: p(nprops,nstate),u(ndofel,nstate)
    real(8), intent(in) :: h(nsvars,nstate)
    real(8), intent(in) :: rref(ndofel,nstate)
    real(8), intent(in) :: jref(ndofel,ndofel,nstate)
    real(8), intent(in) :: href(nsvars,nstate)
    integer :: unit,ss,code,ok,finite_flag
    character(len=16) :: state_name,label
    open(newunit=unit,file=path,status='old',position='append',action='write')
    do ss=1,nstate
       call state_id(ss,state_name)
       call evaluate_uel(ss,u(:,ss),u(:,ss),h(:,ss), &
            rhs0,jac0,hc0,pnewdt)
       call fd_branch_signature(ss,u(:,ss),h(:,ss),p(:,ss),code,label)
       ok=merge(1,0,all(rref(:,ss)==rhs0) .and. &
            all(jref(:,:,ss)==jac0) .and. &
            physical_svars_exact(href(:,ss),hc0))
       finite_flag=merge(1,0,used_finite(rhs0,jac0) .and. &
            physical_svars_finite(hc0))
       write(unit,'(*(g0,:,","))') trim(run_label),'POST_FULL', &
            trim(state_name),merge(1,0,all(rref(:,ss)==rhs0)), &
            merge(1,0,all(jref(:,:,ss)==jac0)), &
            merge(1,0,physical_svars_exact(href(:,ss),hc0)), &
            maxval(abs(rref(:,ss)-rhs0)), &
            maxval(abs(jref(:,:,ss)-jac0)), &
            max_physical_svars_difference(href(:,ss),hc0), &
            finite_flag,code,merge('PASS','FAIL',ok==1)
       if (ok/=1 .or. finite_flag/=1) then
          close(unit)
          stop 74
       end if
    end do
    close(unit)
  end subroutine append_post_full_replay

  subroutine physical_jacobian(jin,jout)
    real(8), intent(in) :: jin(ndofel,ndofel)
    real(8), intent(out) :: jout(ndofel,ndofel)
    integer :: ii
    jout=jin
    do ii=1,ndofel
       jout(ii,ii)=jout(ii,ii)-kstab_num
    end do
  end subroutine physical_jacobian

  real(8) function l2norm(vv)
    real(8), intent(in) :: vv(:)
    l2norm=sqrt(sum(vv*vv))
  end function l2norm

  logical function vector_finite(vv)
    real(8), intent(in) :: vv(:)
    vector_finite=all(ieee_is_finite(vv)) .and. all(abs(vv)<=1.0d100)
  end function vector_finite

  logical function matrix_finite(mm)
    real(8), intent(in) :: mm(:,:)
    matrix_finite=all(ieee_is_finite(mm)) .and. all(abs(mm)<=1.0d100)
  end function matrix_finite

  logical function used_finite(rr,jj)
    real(8), intent(in) :: rr(ndofel),jj(ndofel,ndofel)
    used_finite=vector_finite(rr) .and. matrix_finite(jj)
  end function used_finite

  logical function physical_svars_finite(hh)
    real(8), intent(in) :: hh(nsvars)
    integer :: gp,off
    physical_svars_finite=.true.
    do gp=1,4
       off=48*(gp-1)
       if (.not.vector_finite(hh(off+1:off+14))) &
            physical_svars_finite=.false.
    end do
  end function physical_svars_finite

  logical function physical_svars_exact(a,b)
    real(8), intent(in) :: a(nsvars),b(nsvars)
    integer :: gp,off
    physical_svars_exact=.true.
    do gp=1,4
       off=48*(gp-1)
       if (.not.all(a(off+1:off+14)==b(off+1:off+14))) &
            physical_svars_exact=.false.
    end do
  end function physical_svars_exact

  real(8) function max_physical_svars_difference(a,b)
    real(8), intent(in) :: a(nsvars),b(nsvars)
    integer :: gp,off
    max_physical_svars_difference=0.0d0
    do gp=1,4
       off=48*(gp-1)
       max_physical_svars_difference=max(max_physical_svars_difference, &
            maxval(abs(a(off+1:off+14)-b(off+1:off+14))))
    end do
  end function max_physical_svars_difference

end program wp10_i5_integrated_fd_harness
