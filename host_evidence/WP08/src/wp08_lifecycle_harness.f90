program wp08_lifecycle_harness
  use, intrinsic :: ieee_arithmetic
  implicit none

  integer, parameter :: ndofel=16, nrhs=1, nsvars=192
  integer, parameter :: nprops=27, mcrd=2, nnode=4
  real(8), parameter :: tol=1.0d-12
  real(8) :: props(nprops), coords(mcrd,nnode)
  real(8) :: u0(ndofel), u_trial(ndofel), h0(nsvars), hread(nsvars)
  real(8) :: r_ref(ndofel), j_ref(ndofel,ndofel), hc_ref(nsvars)
  real(8) :: r_test(ndofel), j_test(ndofel,ndofel), hc_test(nsvars)
  real(8) :: r_tmp(ndofel), j_tmp(ndofel,ndofel), hc_tmp(nsvars)
  real(8) :: pnewdt
  real(8) :: m0_y0, m0_ytrial, m0_yreplay
  real(8) :: m1_y0, m1_ytrial, m1_yreplay
  character(len=512) :: outfile, statefile
  integer :: outunit, stateunit
  logical :: all_pass

  call get_command_argument(1,outfile)
  if (len_trim(outfile)==0) outfile='wp08_lifecycle_harness.csv'
  call get_command_argument(2,statefile)
  if (len_trim(statefile)==0) statefile='wp08_svars.bin'

  call initialize_model(props,coords,u0,h0)
  u_trial=u0
  u_trial(5)=u_trial(5)+2.0d-5
  u_trial(7)=u_trial(7)+2.5d5
  u_trial(8)=min(0.45d0,u_trial(8)+0.08d0)
  u_trial(11)=u_trial(11)-1.5d5
  u_trial(12)=min(0.45d0,u_trial(12)+0.05d0)

  open(newunit=outunit,file=trim(outfile),status='replace',action='write')
  write(outunit,'(a)') 'test_id,model,rhs_rel,jac_rel,physical_svars_rel,expected,result'
  all_pass=.true.

  call evaluate_uel(u0,h0,600.0d0,1260.0d0,60.0d0,3,10, &
       r_ref,j_ref,hc_ref,pnewdt)

  ! RC-T1: an intervening trial at the same declared committed state
  ! must not alter the replayed baseline output.
  call evaluate_uel(u_trial,h0,600.0d0,1260.0d0,60.0d0,3,10, &
       r_tmp,j_tmp,hc_tmp,pnewdt)
  call evaluate_uel(u0,h0,600.0d0,1260.0d0,60.0d0,3,10, &
       r_test,j_test,hc_test,pnewdt)
  call record_comparison(outunit,'RC-T1_ORDER','M1',r_ref,j_ref,hc_ref, &
       r_test,j_test,hc_test,'PASS',tol,all_pass)

  ! RC-T2: emulate a rejected trial/cutback by changing trial state,
  ! increment metadata and DTIME, discarding all returned candidate state,
  ! then replaying the original declared state.
  call evaluate_uel(u_trial,h0,615.0d0,1275.0d0,15.0d0,3,11, &
       r_tmp,j_tmp,hc_tmp,pnewdt)
  call evaluate_uel(u0,h0,600.0d0,1260.0d0,60.0d0,3,10, &
       r_test,j_test,hc_test,pnewdt)
  call record_comparison(outunit,'RC-T2_REJECT','M1',r_ref,j_ref,hc_ref, &
       r_test,j_test,hc_test,'PASS',tol,all_pass)

  ! RC-T3: insert a KINC=0 terminal/output-style call. The host-owned
  ! committed state is deliberately restored before replay.
  call evaluate_uel(u_trial,h0,600.0d0,1260.0d0,1.0d-12,3,0, &
       r_tmp,j_tmp,hc_tmp,pnewdt)
  call evaluate_uel(u0,h0,600.0d0,1260.0d0,60.0d0,3,10, &
       r_test,j_test,hc_test,pnewdt)
  call record_comparison(outunit,'RC-T3_TERMINAL','M1',r_ref,j_ref,hc_ref, &
       r_test,j_test,hc_test,'PASS',tol,all_pass)

  ! RC-T4 micro serialization: this verifies the declared SVARS payload,
  ! not Abaqus restart handling of primary DOFs.
  open(newunit=stateunit,file=trim(statefile),status='replace', &
       access='stream',form='unformatted',action='write')
  write(stateunit) h0
  close(stateunit)
  hread=0.0d0
  open(newunit=stateunit,file=trim(statefile),status='old', &
       access='stream',form='unformatted',action='read')
  read(stateunit) hread
  close(stateunit)
  call evaluate_uel(u0,hread,600.0d0,1260.0d0,60.0d0,3,10, &
       r_test,j_test,hc_test,pnewdt)
  call record_comparison(outunit,'RC-T4_SVARS_SERIALIZE','M1', &
       r_ref,j_ref,hc_ref,r_test,j_test,hc_test,'PASS',tol,all_pass)

  call evaluate_uel(u0,h0,600.0d0,1260.0d0,60.0d0,3,10, &
       r_test,j_test,hc_test,pnewdt)
  call record_comparison(outunit,'RC-T5_REPEAT_PROCESS','M1', &
       r_ref,j_ref,hc_ref,r_test,j_test,hc_test,'PASS',tol,all_pass)

  ! Minimal, explicit M0 counterexample. This microkernel is not the
  ! production UEL; it reproduces the rejected-trial SAVE-cache mechanism.
  call legacy_capref_kernel(0,0.0d0,0.0d0,m0_y0)
  call legacy_capref_kernel(1,10.0d0,10.0d0,m0_y0)
  call legacy_capref_kernel(1,20.0d0,10.0d0,m0_ytrial)
  call legacy_capref_kernel(1,10.0d0,10.0d0,m0_yreplay)
  call record_scalar(outunit,'RC-T2_LEGACY_COUNTEREXAMPLE','M0', &
       abs(m0_yreplay-m0_y0),'EXPECTED_FAIL', &
       abs(m0_yreplay-m0_y0)>tol,all_pass)

  call safe_capref_kernel(10.0d0,10.0d0,m1_y0)
  call safe_capref_kernel(20.0d0,10.0d0,m1_ytrial)
  call safe_capref_kernel(10.0d0,10.0d0,m1_yreplay)
  call record_scalar(outunit,'RC-T2_SAFE_MICROKERNEL','M1', &
       abs(m1_yreplay-m1_y0),'PASS', &
       abs(m1_yreplay-m1_y0)<=tol,all_pass)

  close(outunit)
  if (.not.all_pass) stop 2

contains

  subroutine initialize_model(p,c,u,h)
    real(8), intent(out) :: p(nprops),c(mcrd,nnode),u(ndofel),h(nsvars)
    real(8) :: swr,snr,ds,se,pc,peq
    integer :: node,gp,off

    p=0.0d0
    p(1)=1.0d9
    p(2)=0.25d0
    p(3)=0.80d0
    p(4)=1.0d-9
    p(5)=1.0d-3
    p(6)=700.0d0
    p(7)=1.0d-3
    p(8)=0.20d0
    p(9)=0.05d0
    p(10)=1.0d4
    p(11)=2.0d0
    p(12)=2.0d0
    p(13)=2.0d0
    p(14)=5.0d-3
    p(15)=0.10d0
    p(16)=1.0d-14
    p(17)=0.0d0
    p(18)=0.0d0
    p(19)=0.0d0
    p(20)=0.0d0
    p(21)=1.0d-20
    p(22)=1.0d-8

    c(:,1)=(/0.0d0,0.0d0/)
    c(:,2)=(/1.0d0,0.0d0/)
    c(:,3)=(/1.0d0,1.0d0/)
    c(:,4)=(/0.0d0,1.0d0/)

    u=0.0d0
    do node=1,nnode
       off=4*(node-1)
       u(off+1)=1.0d-5*c(1,node)
       u(off+2)=-5.0d-6*c(2,node)
       u(off+3)=1.0d7+1.0d5*c(1,node)
       u(off+4)=0.10d0+0.01d0*c(1,node)
    end do

    swr=p(8)
    snr=p(9)
    ds=1.0d0-swr-snr
    se=(1.0d0-0.10d0-swr)/ds
    pc=p(10)*(se**(-1.0d0/p(11))-1.0d0)
    peq=1.0d7+(1.0d0-se)*pc
    h=0.0d0
    do gp=1,4
       off=48*(gp-1)
       h(off+1)=1.0d7
       h(off+2)=0.10d0
       h(off+3)=0.90d0
       h(off+4)=pc
       h(off+5)=peq
       h(off+6)=0.0d0
       h(off+7)=0.15d0
       h(off+8)=1.0d-14
       h(off+9)=0.0d0
       h(off+10:off+12)=0.0d0
       h(off+13)=peq-1.0d7
       h(off+14)=1.0d0
    end do
  end subroutine initialize_model

  subroutine evaluate_uel(uin,hin,tstep,ttotal,dtime,kstep,kinc, &
       rout,jout,hout,pnewdt_out)
    real(8), intent(in) :: uin(ndofel),hin(nsvars),tstep,ttotal,dtime
    integer, intent(in) :: kstep,kinc
    real(8), intent(out) :: rout(ndofel),jout(ndofel,ndofel)
    real(8), intent(out) :: hout(nsvars),pnewdt_out
    integer, parameter :: mlvarx=ndofel, mdload=1, npredf=1
    real(8) :: rhs(mlvarx,nrhs),amatrx(ndofel,ndofel)
    real(8) :: svars_local(nsvars),energy(8),ulocal(ndofel)
    real(8) :: du(mlvarx,nrhs),v(ndofel),a(ndofel),time(2)
    real(8) :: params(3),adlmag(mdload,1),predef(2,npredf,nnode)
    real(8) :: ddlmag(mdload,1),pnewdt,period
    integer :: jdltyp(mdload,1),lflags(10),jprops(1)
    integer :: jtype,jelem,ndload,njprop

    rhs=0.0d0
    amatrx=0.0d0
    svars_local=hin
    energy=0.0d0
    ulocal=uin
    du=0.0d0
    v=0.0d0
    a=0.0d0
    time=(/tstep,ttotal/)
    params=0.0d0
    adlmag=0.0d0
    predef=0.0d0
    ddlmag=0.0d0
    jdltyp=0
    lflags=0
    lflags(1)=1
    lflags(3)=1
    jprops=0
    jtype=1
    jelem=10287
    ndload=0
    njprop=1
    pnewdt=1.0d6
    period=600.0d0

    call UEL(rhs,amatrx,svars_local,energy,ndofel,nrhs,nsvars, &
         props,nprops,coords,mcrd,nnode,ulocal,du,v,a,jtype,time,dtime, &
         kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf, &
         lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop,period)

    rout=rhs(:,1)
    jout=amatrx
    hout=svars_local
    pnewdt_out=pnewdt
  end subroutine evaluate_uel

  subroutine record_comparison(unit,test_id,model,ra,ja,ha,rb,jb,hb, &
       expected,threshold,global_pass)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: test_id,model,expected
    real(8), intent(in) :: ra(ndofel),ja(ndofel,ndofel),ha(nsvars)
    real(8), intent(in) :: rb(ndofel),jb(ndofel,ndofel),hb(nsvars)
    real(8), intent(in) :: threshold
    logical, intent(inout) :: global_pass
    real(8) :: dr,dj,dh
    logical :: passed

    dr=relative_vector_difference(ra,rb)
    dj=relative_matrix_difference(ja,jb)
    dh=relative_physical_state_difference(ha,hb)
    passed=dr<=threshold .and. dj<=threshold .and. dh<=threshold
    if (.not.passed) global_pass=.false.
    write(unit,'(a,",",a,3(",",es24.16e3),",",a,",",a)') &
         trim(test_id),trim(model),dr,dj,dh,trim(expected), &
         merge('PASS','FAIL',passed)
  end subroutine record_comparison

  subroutine record_scalar(unit,test_id,model,diff,expected,passed, &
       global_pass)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: test_id,model,expected
    real(8), intent(in) :: diff
    logical, intent(in) :: passed
    logical, intent(inout) :: global_pass
    if (.not.passed) global_pass=.false.
    write(unit,'(a,",",a,3(",",es24.16e3),",",a,",",a)') &
         trim(test_id),trim(model),diff,0.0d0,0.0d0,trim(expected), &
         merge('PASS','FAIL',passed)
  end subroutine record_scalar

  real(8) function relative_vector_difference(a1,a2)
    real(8), intent(in) :: a1(:),a2(:)
    real(8) :: scale
    if (.not.all(ieee_is_finite(a1)) .or. &
        .not.all(ieee_is_finite(a2))) then
       relative_vector_difference=huge(1.0d0)
       return
    end if
    scale=max(1.0d0,maxval(abs(a1)),maxval(abs(a2)))
    relative_vector_difference=maxval(abs(a1-a2))/scale
  end function relative_vector_difference

  real(8) function relative_matrix_difference(a1,a2)
    real(8), intent(in) :: a1(:,:),a2(:,:)
    real(8) :: scale
    if (.not.all(ieee_is_finite(a1)) .or. &
        .not.all(ieee_is_finite(a2))) then
       relative_matrix_difference=huge(1.0d0)
       return
    end if
    scale=max(1.0d0,maxval(abs(a1)),maxval(abs(a2)))
    relative_matrix_difference=maxval(abs(a1-a2))/scale
  end function relative_matrix_difference

  real(8) function relative_physical_state_difference(a1,a2)
    real(8), intent(in) :: a1(nsvars),a2(nsvars)
    real(8) :: scale,diff
    integer :: gp,first,last
    scale=1.0d0
    diff=0.0d0
    do gp=1,4
       first=48*(gp-1)+1
       last=first+13
       if (.not.all(ieee_is_finite(a1(first:last))) .or. &
           .not.all(ieee_is_finite(a2(first:last)))) then
          relative_physical_state_difference=huge(1.0d0)
          return
       end if
       scale=max(scale,maxval(abs(a1(first:last))), &
            maxval(abs(a2(first:last))))
       diff=max(diff,maxval(abs(a1(first:last)-a2(first:last))))
    end do
    relative_physical_state_difference=diff/scale
  end function relative_physical_state_difference

  subroutine legacy_capref_kernel(action,x,committed,y)
    integer, intent(in) :: action
    real(8), intent(in) :: x,committed
    real(8), intent(out) :: y
    real(8), save :: cache=0.0d0
    logical, save :: initialized=.false.
    if (action==0) then
       cache=committed
       initialized=.false.
       y=0.0d0
       return
    end if
    if (.not.initialized) then
       cache=committed
       initialized=.true.
    end if
    y=x-cache
    cache=x
  end subroutine legacy_capref_kernel

  subroutine safe_capref_kernel(x,committed,y)
    real(8), intent(in) :: x,committed
    real(8), intent(out) :: y
    y=x-committed
  end subroutine safe_capref_kernel

end program wp08_lifecycle_harness
