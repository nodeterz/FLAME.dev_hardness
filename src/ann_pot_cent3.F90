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
    integer:: iat, i, j, ng
    real(8):: epot_c, out_ann
    real(8):: time1, time2, time3, time4, time5, time6, time7, time8
    real(8):: tt1, tt2, tt3, fx_es, fy_es, fz_es, hinv(3,3), vol
    call f_routine(id='cal_ann_cent3')
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%fat_chi(1:3,1:atoms%nat))
        allocate(ann_arr%chi_i(1:atoms%nat))
        allocate(ann_arr%chi_o(1:atoms%nat))
        allocate(ann_arr%chi_d(1:atoms%nat))
        allocate(ann_arr%a(1:(atoms%nat+1)*(atoms%nat+1)))
        ann_arr%fat_chi=0.d0
        ann_arr%chi_i=0.d0
        ann_arr%chi_o=0.d0
        ann_arr%chi_d=0.d0
        ann_arr%a=0.d0
    else
        ann_arr%fat_chi=0.d0
        ann_arr%chi_i=0.d0
        ann_arr%chi_o=0.d0
        ann_arr%chi_d=0.d0
        ann_arr%a=0.d0
    endif
    if(parini%iverbose>=2) call cpu_time(time1)
    call init_electrostatic_cent1(parini,atoms,ann_arr,ann_arr%a,poisson)
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
    endif
    if(parini%iverbose>=2) call cpu_time(time3)
    over_iat: do iat=1,atoms%nat
        i=atoms%itypat(iat)
        ng=ann_arr%ann(i)%nn(0)
        ann_arr%ann(i)%y(1:ng,0)=symfunc%y(1:ng,iat)
        if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
            call cal_architecture(ann_arr%ann(i),out_ann)
            call cal_force_chi_part1(parini,symfunc,iat,atoms,out_ann,ann_arr)
        elseif(trim(ann_arr%event)=='train') then
            call cal_architecture_der(ann_arr%ann(i),out_ann)
            ann_arr%chi_i(iat)=out_ann
            tt1=tanh(ann_arr%ann(i)%prefactor_chi*out_ann)
            ann_arr%chi_o(iat)=ann_arr%ann(i)%ampl_chi*tt1+ann_arr%ann(i)%chi0
            call convert_ann_epotd(ann_arr%ann(i),ann_arr%num(i),ann_arr%g_per_atom(1,iat))
            ann_arr%g_per_atom(1:ann_arr%num(1),iat)=ann_arr%g_per_atom(1:ann_arr%num(1),iat)*ann_arr%ann(i)%ampl_chi*ann_arr%ann(i)%prefactor_chi*(1.d0-tt1**2)
        else
            stop 'ERROR: undefined content for ann_arr%event'
        endif
    enddo over_iat
    !This must be here since contribution from coulomb
    !interaction is calculated during the process of charge optimization.
    if(parini%iverbose>=2) call cpu_time(time4)
    call get_qat_from_chi_cent1(parini,ann_arr,atoms,poisson,ann_arr%a)
    if(parini%iverbose>=2) call cpu_time(time5)
    atoms%stress(1:3,1:3)=0.d0
    atoms%fat(1:3,1:atoms%nat)=0.d0
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        call cal_force_chi_part2(parini,symfunc,atoms,ann_arr)
    endif !end of if for potential
    if(parini%iverbose>=2) call cpu_time(time6)
    call get_electrostatic_cent1(parini,atoms,ann_arr,epot_c,ann_arr%a,poisson)
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
        do iat=1,atoms%nat
            fx_es=atoms%fat(1,iat)-ann_arr%fat_chi(1,iat)
            fy_es=atoms%fat(2,iat)-ann_arr%fat_chi(2,iat)
            fz_es=atoms%fat(3,iat)-ann_arr%fat_chi(3,iat)
            tt1=tt1+fx_es**2+fy_es**2+fz_es**2
            tt2=tt2+ann_arr%fat_chi(1,iat)**2+ann_arr%fat_chi(2,iat)**2+ann_arr%fat_chi(3,iat)**2
            tt3=tt3+fx_es*ann_arr%fat_chi(1,iat)+fy_es*ann_arr%fat_chi(2,iat)+fz_es*ann_arr%fat_chi(3,iat)
        enddo
        tt1=sqrt(tt1)
        tt2=sqrt(tt2)
        ann_arr%fchi_angle=tt3/(tt1*tt2)
        ann_arr%fchi_norm=tt2/tt1
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
        deallocate(ann_arr%a)
        deallocate(ann_arr%fat_chi)
        deallocate(ann_arr%fatpq)
        deallocate(ann_arr%stresspq)
    endif
    call f_release_routine()
end subroutine cal_ann_cent3
!*****************************************************************************************
