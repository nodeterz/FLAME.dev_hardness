!*****************************************************************************************
subroutine cal_ann_cent3(parini,atoms,symfunc,ann_arr)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, convert_ann_epotd
    use mod_symfunc, only: typ_symfunc
    use mod_electrostatics, only: typ_poisson
    use mod_linked_lists, only: typ_pia_arr
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_poisson):: poisson
    !local variables
    type(typ_pia_arr):: pia_arr_tmp
    integer:: iat, i, j, ng, ng_2
    real(8):: epot_c, out_ann, out_ann_2
    real(8):: time1, time2, time3, time4, time5, time6, time7, time8
    real(8):: tt1, tt2, tt3, fx_es, fy_es, fz_es, hinv(3,3), vol
    real(8):: tt4, tt5, tt6
    call f_routine(id='cal_ann_cent3')
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%chi_i(1:atoms%nat))
        allocate(ann_arr%chi_o(1:atoms%nat))
        allocate(ann_arr%chi_d(1:atoms%nat))
        allocate(ann_arr%fat_chi(1:3,1:atoms%nat))
        allocate(ann_arr%hardness_i(1:atoms%nat))
        allocate(ann_arr%hardness_o(1:atoms%nat))
        allocate(ann_arr%hardness_d(1:atoms%nat))
        allocate(ann_arr%fat_hardness(1:3,1:atoms%nat))
        allocate(ann_arr%a(1:(atoms%nat+1)*(atoms%nat+1)))
        ann_arr%fat_chi=0.d0
        ann_arr%chi_i=0.d0
        ann_arr%chi_o=0.d0
        ann_arr%chi_d=0.d0
        ann_arr%fat_hardness=0.d0
        ann_arr%hardness_i=0.d0
        ann_arr%hardness_o=0.d0
        ann_arr%hardness_d=0.d0
        ann_arr%a=0.d0
    else
        ann_arr%fat_chi=0.d0
        ann_arr%chi_i=0.d0
        ann_arr%chi_o=0.d0
        ann_arr%chi_d=0.d0
        ann_arr%fat_hardness=0.d0
        ann_arr%hardness_i=0.d0
        ann_arr%hardness_o=0.d0
        ann_arr%hardness_d=0.d0
        ann_arr%a=0.d0
    endif
    if(parini%iverbose>=2) call cpu_time(time1)
    if(parini%iverbose>=2) call cpu_time(time2)
    if(ann_arr%compute_symfunc) then
        call symmetry_functions(parini,ann_arr,atoms,symfunc,.true.)
    else
        symfunc%linked_lists%rcut=ann_arr%rcut
        symfunc%linked_lists%triplex=.true.
        call call_linkedlist(parini,atoms,.true.,symfunc%linked_lists,pia_arr_tmp)
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%fatpq(1:3,1:symfunc%linked_lists%maxbound_rad))
        allocate(ann_arr%stresspq(1:3,1:3,1:symfunc%linked_lists%maxbound_rad))
        allocate(ann_arr%fatpq_hardness(1:3,1:symfunc%linked_lists%maxbound_rad))
        allocate(ann_arr%stresspq_hardness(1:3,1:3,1:symfunc%linked_lists%maxbound_rad))
    endif
    if(parini%iverbose>=2) call cpu_time(time3)
    over_iat: do iat=1,atoms%nat
        i=atoms%itypat(iat)
        j=parini%ntypat+atoms%itypat(iat)
        ng=ann_arr%ann(i)%nn(0)
        ng_2=ann_arr%ann(j)%nn(0)
        ann_arr%ann(i)%y(1:ng,0)=symfunc%y(1:ng,iat)
        ann_arr%ann(j)%y(1:ng_2,0)=symfunc%y(1:ng_2,iat)
        if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
            call cal_architecture(ann_arr%ann(i),out_ann)
            call cal_architecture(ann_arr%ann(j),out_ann_2)
            call cal_force_chi_part1(parini,symfunc,iat,atoms,out_ann,ann_arr)
            call cal_force_hardness_part1(parini,symfunc,iat,atoms,out_ann_2,ann_arr)
            ann_arr%hardness_i(iat)=out_ann_2
            tt4=tanh(ann_arr%ann(j)%prefactor_hardness*out_ann_2)
            ann_arr%hardness_o(iat)=ann_arr%ann(j)%ampl_hardness*tt4+ann_arr%ann(j)%hardness0
        elseif(trim(ann_arr%event)=='train') then
            call cal_architecture_der(ann_arr%ann(i),out_ann)
            ann_arr%chi_i(iat)=out_ann
            tt1=tanh(ann_arr%ann(i)%prefactor_chi*out_ann)
            ann_arr%chi_o(iat)=ann_arr%ann(i)%ampl_chi*tt1+ann_arr%ann(i)%chi0
            call convert_ann_epotd(ann_arr%ann(i),ann_arr%num(i),ann_arr%g_per_atom_chi(1,iat))
            ann_arr%g_per_atom_chi(1:ann_arr%num(1),iat)=ann_arr%g_per_atom_chi(1:ann_arr%num(1),iat)*ann_arr%ann(i)%ampl_chi*ann_arr%ann(i)%prefactor_chi*(1.d0-tt1**2)

            call cal_architecture_der(ann_arr%ann(j),out_ann_2)
            ann_arr%hardness_i(iat)=out_ann_2
            tt4=tanh(ann_arr%ann(j)%prefactor_hardness*out_ann_2)
            ann_arr%hardness_o(iat)=ann_arr%ann(j)%ampl_hardness*tt4+ann_arr%ann(j)%hardness0
            call convert_ann_epotd(ann_arr%ann(j),ann_arr%num(i),ann_arr%g_per_atom_hardness(1,iat))
            ann_arr%g_per_atom_hardness(1:ann_arr%num(1),iat)=parini%ampl_grad_hardness*ann_arr%g_per_atom_hardness(1:ann_arr%num(1),iat)*&
                                                              ann_arr%ann(j)%ampl_hardness*ann_arr%ann(j)%prefactor_hardness*(1.d0-tt4**2)
        else
            stop 'ERROR: undefined content for ann_arr%event'
        endif
    enddo over_iat
    write(98,'(a,3es14.6)')'J_max,J_min_J_var',maxval(ann_arr%hardness_o(1:parini%ntypat)),minval(ann_arr%hardness_o(1:parini%ntypat)),&
                                               maxval(ann_arr%hardness_o(1:parini%ntypat))-minval(ann_arr%hardness_o(1:parini%ntypat))
    call init_electrostatic_cent3(parini,atoms,ann_arr,ann_arr%a,poisson)
    !This must be here since contribution from coulomb
    !interaction is calculated during the process of charge optimization.
    if(parini%iverbose>=2) call cpu_time(time4)
    call get_qat_from_chi_cent1(parini,ann_arr,atoms,poisson,ann_arr%a)
    if(parini%iverbose>=2) call cpu_time(time5)
    atoms%stress(1:3,1:3)=0.d0
    atoms%fat(1:3,1:atoms%nat)=0.d0
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        call cal_force_chi_part2(parini,symfunc,atoms,ann_arr)
        call cal_force_hardness_part2(parini,symfunc,atoms,ann_arr)
    endif !end of if for potential
    if(parini%iverbose>=2) call cpu_time(time6)
         call get_electrostatic_cent3(parini,atoms,ann_arr,epot_c,ann_arr%a,poisson)
    if(parini%iverbose>=2) then
        call cpu_time(time7)
        call yaml_mapping_open('Timing of CENT1')
        call yaml_map('initialize matrix',time2-time1)
        call yaml_map('calculation of symfunc',time3-time2)
        call yaml_map('neural network process',time4-time3)
        call yaml_map('linear equations solver',time5-time4)
        call yaml_map('force (SR term)',time6-time5)
        call yaml_map('energy (SR+LR), force (LR)',time7-time6)
        call yaml_map('total time',time7-time1)
        call yaml_mapping_close()
        !write(*,'(a,f8.3)') 'Timing:cent3: initialize matrix          ',time2-time1
        !write(*,'(a,f8.3)') 'Timing:cent3: calculation of symfunc     ',time3-time2
        !write(*,'(a,f8.3)') 'Timing:cent3: neural network process     ',time4-time3
        !write(*,'(a,f8.3)') 'Timing:cent3: linear equations solver    ',time5-time4
        !write(*,'(a,f8.3)') 'Timing:cent3: force (SR term)            ',time6-time5
        !write(*,'(a,f8.3)') 'Timing:cent3: energy (SR+LR), force (LR) ',time7-time6
        !write(*,'(a,f8.3)') 'Timing:cent3: total time                 ',time7-time1
    endif !end of if for printing out timing.
    atoms%epot=epot_c
    if(trim(ann_arr%event)=='evalu') then
        tt1=0.d0
        tt2=0.d0
        tt3=0.d0
        tt4=0.d0
        tt5=0.d0
        tt6=0.d0
        do iat=1,atoms%nat
            fx_es=atoms%fat(1,iat)-ann_arr%fat_chi(1,iat)
            fy_es=atoms%fat(2,iat)-ann_arr%fat_chi(2,iat)
            fz_es=atoms%fat(3,iat)-ann_arr%fat_chi(3,iat)
            tt1=tt1+fx_es**2+fy_es**2+fz_es**2
            tt2=tt2+ann_arr%fat_chi(1,iat)**2+ann_arr%fat_chi(2,iat)**2+ann_arr%fat_chi(3,iat)**2
            tt3=tt3+fx_es*ann_arr%fat_chi(1,iat)+fy_es*ann_arr%fat_chi(2,iat)+fz_es*ann_arr%fat_chi(3,iat)

            fx_es=atoms%fat(1,iat)-ann_arr%fat_hardness(1,iat)
            fy_es=atoms%fat(2,iat)-ann_arr%fat_hardness(2,iat)
            fz_es=atoms%fat(3,iat)-ann_arr%fat_hardness(3,iat)
            tt4=tt4+fx_es**2+fy_es**2+fz_es**2
            tt5=tt5+ann_arr%fat_hardness(1,iat)**2+ann_arr%fat_hardness(2,iat)**2+ann_arr%fat_hardness(3,iat)**2
            tt6=tt6+fx_es*ann_arr%fat_hardness(1,iat)+fy_es*ann_arr%fat_hardness(2,iat)+fz_es*ann_arr%fat_hardness(3,iat)
        enddo
        tt1=sqrt(tt1)
        tt2=sqrt(tt2)
        ann_arr%fchi_angle=tt3/(tt1*tt2)
        ann_arr%fchi_norm=tt2/tt1
        tt4=sqrt(tt4)
        tt5=sqrt(tt5)
        ann_arr%fhardness_angle=tt6/(tt4*tt5)
        ann_arr%fhardness_norm=tt5/tt4
        write(*,'(a,4es14.6)') 'f',tt2,tt5
    endif
    call fini_electrostatic_cent1(parini,atoms,poisson)
    !call repulsive_potential_cent(parini,atoms,ann_arr)
    call getvol_alborz(atoms%cellvec,vol)
    call invertmat_alborz(atoms%cellvec,hinv)
    !The following line is inconsistent with the definition of stress tensor
    atoms%stress(1:3,1:3)=atoms%stress(1:3,1:3)*vol
    do i=1,3
    do j=1,3
        atoms%celldv(i,j)=vol*(atoms%stress(i,1)*hinv(j,1)+atoms%stress(i,2)*hinv(j,2)+atoms%stress(i,3)*hinv(j,3))
    enddo
    enddo

    deallocate(symfunc%linked_lists%prime_bound)
    deallocate(symfunc%linked_lists%bound_rad)
    deallocate(symfunc%linked_lists%bound_ang)
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        call f_free(symfunc%y)
        call f_free(symfunc%y0d)
        call f_free(symfunc%y0dr)
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        deallocate(ann_arr%chi_i)
        deallocate(ann_arr%chi_o)
        deallocate(ann_arr%chi_d)
        deallocate(ann_arr%hardness_i)
        deallocate(ann_arr%hardness_o)
        deallocate(ann_arr%hardness_d)
        deallocate(ann_arr%a)
        deallocate(ann_arr%fat_chi)
        deallocate(ann_arr%fat_hardness)
        deallocate(ann_arr%fatpq)
        deallocate(ann_arr%stresspq)
        deallocate(ann_arr%fatpq_hardness)
        deallocate(ann_arr%stresspq_hardness)
    endif
    call f_release_routine()
end subroutine cal_ann_cent3
!*****************************************************************************************
subroutine init_electrostatic_cent3(parini,atoms,ann_arr,a,poisson)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    type(typ_poisson), intent(inout):: poisson
    real(8),allocatable :: gausswidth(:)
    !local variables
    integer:: iat, jat
    real(8):: vol, c
    real(8):: dx, dy, dz, r, tt1, tt2, pi, beta_iat, beta_jat, gama, ttf
    associate(epot_es=>ann_arr%epot_es)
    pi=4.d0*atan(1.d0)
    ann_arr%ener_ref=0.d0
    do iat=1,atoms%nat
        ann_arr%ener_ref=ann_arr%ener_ref+ann_arr%ann(atoms%itypat(iat))%ener_ref
    enddo
    if (.not. parini%ewald) then 
        poisson%alpha = maxval(ann_arr%ann(:)%gausswidth)
    else 
        if (parini%alpha_ewald<0.d0) then
            call getvol_alborz(atoms%cellvec,vol)
            c=2.2d0
            poisson%alpha = 1.d0/(c*sqrt(pi)*(atoms%nat/vol**2)**(1.d0/6.d0))
            write(*,*)"optimized alpha = ", poisson%alpha
        else
            poisson%alpha=parini%alpha_ewald
        endif
    end if
    if(trim(parini%syslinsolver_ann)=='direct' .or. trim(parini%syslinsolver_ann)=='apply_matrix') then
        if(trim(atoms%boundcond)/='free') then
            write(*,*) 'ERROR: syslinsolver=direct can be used only for free BC.'
        endif
        call get_amat_cent3(atoms,ann_arr,a)
    elseif(trim(parini%syslinsolver_ann)=='operator') then
        if(trim(atoms%boundcond)=='bulk' .or. trim(atoms%boundcond)=='slab') then
            allocate(gausswidth(atoms%nat))
            gausswidth(:)=ann_arr%ann(atoms%itypat(:))%gausswidth
            poisson%task_finit="alloc_rho:set_ngp"
            call init_hartree(parini,atoms,poisson,gausswidth)
            deallocate(gausswidth)
        else
            write(*,*) 'ERROR: currently syslinsolver=operator only for BC=bulk/slab.'
            stop
        endif
    else
        write(*,*) 'ERROR: unknown value for syslinsolver',trim(parini%syslinsolver_ann)
        stop
    endif
    end associate
end subroutine init_electrostatic_cent3
!*****************************************************************************************
subroutine get_amat_cent3(atoms,ann_arr,a)
    use mod_interface
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    !local variables
    integer:: iat, jat
    real(8):: dx, dy, dz, r, pi, beta_iat, beta_jat, gama
    pi=4.d0*atan(1.d0)
    call update_ratp(atoms)
    do iat=1,atoms%nat
        a(iat,atoms%nat+1)=1.d0
        a(atoms%nat+1,iat)=1.d0
        beta_iat=ann_arr%ann(atoms%itypat(iat))%gausswidth
        gama=1.d0/sqrt(beta_iat**2+beta_iat**2)
        a(iat,iat)=gama*2.d0/sqrt(pi)+ann_arr%hardness_o(iat)
        do jat=iat+1,atoms%nat
            dx=atoms%ratp(1,jat)-atoms%ratp(1,iat)
            dy=atoms%ratp(2,jat)-atoms%ratp(2,iat)
            dz=atoms%ratp(3,jat)-atoms%ratp(3,iat)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            beta_jat=ann_arr%ann(atoms%itypat(jat))%gausswidth
            gama=1.d0/sqrt(beta_iat**2+beta_jat**2)
            a(iat,jat)=erf(gama*r)/r
            a(jat,iat)=a(iat,jat)
        enddo
    enddo
    a(atoms%nat+1,atoms%nat+1)=0.d0
end subroutine get_amat_cent3
!*****************************************************************************************
subroutine cal_force_hardness_part1(parini,symfunc,iat,atoms,out_ann,ann_arr)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_symfunc), intent(in):: symfunc
    integer, intent(in):: iat
    type(typ_atoms), intent(in):: atoms
    real(8), intent(in):: out_ann
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    real(8):: sxx, sxy, sxz, syx, syy, syz, szx, szy, szz
    real(8):: ttx, tty, ttz, tt1, tt2
    integer:: ib, i, j
    i=parini%ntypat+atoms%itypat(iat)
    ann_arr%hardness_i(iat)=out_ann
    tt1=tanh(ann_arr%ann(i)%prefactor_hardness*out_ann)
    ann_arr%hardness_o(iat)=ann_arr%ann(i)%ampl_hardness*tt1+ann_arr%ann(i)%hardness0
    if(trim(ann_arr%event)/='train') then
        tt2=ann_arr%ann(i)%ampl_hardness*ann_arr%ann(i)%prefactor_hardness*(1.d0-tt1**2)
        do ib=symfunc%linked_lists%prime_bound(iat),symfunc%linked_lists%prime_bound(iat+1)-1
            ttx=0.d0 ; tty=0.d0 ; ttz=0.d0
            do j=1,ann_arr%ann(i)%nn(0)
                ttx=ttx+ann_arr%ann(i)%d(j)*symfunc%y0d(j,1,ib)
                tty=tty+ann_arr%ann(i)%d(j)*symfunc%y0d(j,2,ib)
                ttz=ttz+ann_arr%ann(i)%d(j)*symfunc%y0d(j,3,ib)
            enddo
            ann_arr%fatpq_hardness(1,ib)=ttx*tt2
            ann_arr%fatpq_hardness(2,ib)=tty*tt2
            ann_arr%fatpq_hardness(3,ib)=ttz*tt2
        enddo
    endif
    if(trim(ann_arr%event)=='potential') then
        do ib=symfunc%linked_lists%prime_bound(iat),symfunc%linked_lists%prime_bound(iat+1)-1
            sxx=0.d0 ; sxy=0.d0 ; sxz=0.d0
            syx=0.d0 ; syy=0.d0 ; syz=0.d0
            szx=0.d0 ; szy=0.d0 ; szz=0.d0
            do j=1,ann_arr%ann(i)%nn(0)
                sxx=sxx+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,1,ib)
                sxy=sxy+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,2,ib)
                sxz=sxz+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,3,ib)
                syx=syx+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,4,ib)
                syy=syy+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,5,ib)
                syz=syz+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,6,ib)
                szx=szx+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,7,ib)
                szy=szy+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,8,ib)
                szz=szz+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,9,ib)
            enddo
            ann_arr%stresspq_hardness(1,1,ib)=-sxx*tt2
            ann_arr%stresspq_hardness(2,1,ib)=-syx*tt2
            ann_arr%stresspq_hardness(3,1,ib)=-szx*tt2
            ann_arr%stresspq_hardness(1,2,ib)=-sxy*tt2
            ann_arr%stresspq_hardness(2,2,ib)=-syy*tt2
            ann_arr%stresspq_hardness(3,2,ib)=-szy*tt2
            ann_arr%stresspq_hardness(1,3,ib)=-sxz*tt2
            ann_arr%stresspq_hardness(2,3,ib)=-syz*tt2
            ann_arr%stresspq_hardness(3,3,ib)=-szz*tt2
        enddo
    endif
end subroutine cal_force_hardness_part1
!*****************************************************************************************
subroutine cal_force_hardness_part2(parini,symfunc,atoms,ann_arr)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_symfunc), intent(in):: symfunc
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    real(8):: ttx, tty, ttz, qnet, hinv(3,3), vol
    integer:: ib, i, j, iat, jat
    if(trim(ann_arr%event)/='train') then
        do ib=1,symfunc%linked_lists%maxbound_rad
            iat=symfunc%linked_lists%bound_rad(1,ib)
            jat=symfunc%linked_lists%bound_rad(2,ib)
            if(trim(ann_arr%approach)=='eem1' .or. trim(ann_arr%approach)=='cent1' .or. trim(ann_arr%approach)=='cent3') then
                qnet=atoms%qat(iat)
            elseif(trim(ann_arr%approach)=='cent2') then
                qnet=atoms%zat(iat)+atoms%qat(iat)
            else
                write(*,'(2a)') 'ERROR: unknown approach in ANN, ',trim(ann_arr%approach)
                stop
            endif
            ttx=ann_arr%fatpq_hardness(1,ib)*qnet
            tty=ann_arr%fatpq_hardness(2,ib)*qnet
            ttz=ann_arr%fatpq_hardness(3,ib)*qnet
            ann_arr%fat_hardness(1,iat)=ann_arr%fat_hardness(1,iat)+ttx
            ann_arr%fat_hardness(2,iat)=ann_arr%fat_hardness(2,iat)+tty
            ann_arr%fat_hardness(3,iat)=ann_arr%fat_hardness(3,iat)+ttz
            ann_arr%fat_hardness(1,jat)=ann_arr%fat_hardness(1,jat)-ttx
            ann_arr%fat_hardness(2,jat)=ann_arr%fat_hardness(2,jat)-tty
            ann_arr%fat_hardness(3,jat)=ann_arr%fat_hardness(3,jat)-ttz
            if(trim(ann_arr%event)=='potential') then
            atoms%stress(1,1)=atoms%stress(1,1)+ann_arr%stresspq(1,1,ib)*qnet
            atoms%stress(2,1)=atoms%stress(2,1)+ann_arr%stresspq(2,1,ib)*qnet
            atoms%stress(3,1)=atoms%stress(3,1)+ann_arr%stresspq(3,1,ib)*qnet
            atoms%stress(1,2)=atoms%stress(1,2)+ann_arr%stresspq(1,2,ib)*qnet
            atoms%stress(2,2)=atoms%stress(2,2)+ann_arr%stresspq(2,2,ib)*qnet
            atoms%stress(3,2)=atoms%stress(3,2)+ann_arr%stresspq(3,2,ib)*qnet
            atoms%stress(1,3)=atoms%stress(1,3)+ann_arr%stresspq(1,3,ib)*qnet
            atoms%stress(2,3)=atoms%stress(2,3)+ann_arr%stresspq(2,3,ib)*qnet
            atoms%stress(3,3)=atoms%stress(3,3)+ann_arr%stresspq(3,3,ib)*qnet
            endif
        enddo
        do iat=1,atoms%nat
            atoms%fat(1,iat)=atoms%fat(1,iat)+ann_arr%fat_hardness(1,iat)
            atoms%fat(2,iat)=atoms%fat(2,iat)+ann_arr%fat_hardness(2,iat)
            atoms%fat(3,iat)=atoms%fat(3,iat)+ann_arr%fat_hardness(3,iat)
        enddo
    endif
    call getvol_alborz(atoms%cellvec,vol)
    if(trim(ann_arr%event)=='potential') then
    atoms%stress(1:3,1:3)=atoms%stress(1:3,1:3)/vol !*atoms%nat !not certain if this is needed!!!
    endif
end subroutine cal_force_hardness_part2
!*****************************************************************************************
subroutine get_electrostatic_cent3(parini,atoms,ann_arr,epot_c,a,poisson)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(out):: epot_c
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iat, jat
    real(8):: vol, c
    real(8):: dx, dy, dz, r, tt1, tt2, pi, beta_iat, beta_jat, gama, ttf
    associate(epot_es=>ann_arr%epot_es)
    tt1=0.d0
    tt2=0.d0
    do iat=1,atoms%nat
        tt1=tt1+ann_arr%chi_o(iat)*atoms%qat(iat)
        tt2=tt2+atoms%qat(iat)**2*0.5d0*ann_arr%hardness_o(iat)
    enddo
    call cal_electrostatic_ann(parini,atoms,ann_arr,a,poisson)
    epot_c=epot_es+tt1+tt2+ann_arr%ener_ref
    end associate
end subroutine get_electrostatic_cent3
