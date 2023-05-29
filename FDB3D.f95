!-FDB_beta optimization
!--Program written by Pau Besalu and Guillem Pey* in 2022-2023
Program FDBb
double precision :: E_0,E_p,E_f,E_r,E_fdb
double precision, dimension(3) :: mu,mu_p,mu_f,mu_r,F
double precision, dimension(3,3) :: alpha,alpha_p,alpha_f,alpha_r,alpha_tmp
double precision, dimension(3,3,3) :: beta,beta_p,beta_f,beta_r
character*80 hermenegildo
character*10 ET
character*3 axis_scan
character*6 approximation
double precision x0,y0,z0,R,c_i,target_barrier,potential
integer i,j,k,d,approx

        COMMON /values/ E_0,mu,alpha,beta,mu_f,alpha_f,beta_f
        COMMON /properties/ target_barrier,x0,y0,z0,R,c_i,d,potential
        COMMON /keywords/ approximation,axis_scan,ET !(Dipole/dipole or 1, Alpha/alpha or 2, Beta/Beta or 3)

call readvalues()

call minimize()

End
!-------------------------------------!
Function E_fdb(F,approx)
        implicit none
        double precision :: E_0,E_fdb,E_mu,E_alpha
        double precision, dimension(3) :: mu,F,mu_f
        double precision, dimension(3,3) :: alpha,alpha_f
        double precision, dimension(3,3,3) :: beta,beta_f
        integer i,j,k,d,approx
        double precision target_barrier,x0,y0,z0,R,c_i,potential
        character*10 ET
        character*3 axis_scan
        character*6 approximation
        COMMON /values/ E_0,mu,alpha,beta,mu_f,alpha_f,beta_f
        COMMON /properties/ target_barrier,x0,y0,z0,R,c_i,d,potential
        COMMON /keywords/ approximation,axis_scan,ET

        E_fdb=E_0-potential

        do i=1,3
                E_fdb=E_fdb-mu(i)*F(i)
        end do

        !if (mu_f(1) < 0) then
        !   write(*,*) "WARNING. LIGAND NOT PROPERLY ORIENTED. STOP"
        !   stop
        !end if

        E_fdb=E_fdb-mu_f(1)*sqrt(F(1)**2+F(2)**2+F(3)**2) !adds solvent correction
        E_mu=E_fdb

        do i=1,3
                do j=1,3
                        E_fdb=E_fdb-0.5d0*alpha(i,j)*F(i)*F(j)
                end do
        end do
        E_fdb=E_fdb-0.5d0*alpha_f(1,1)*sqrt(F(1)**2+F(2)**2+F(3)**2)**2 !adds solvent correction
        E_alpha=E_fdb

        do i=1,3
                do j=1,3
                        do k=1,3
                                E_fdb=E_fdb+1.0d0/6.0d0*beta(i,j,k)*F(i)*F(j)*F(k)
                        end do
                end do
        end do
        E_fdb=E_fdb+1.0d0/6.0d0*beta_f(1,1,1)*sqrt(F(1)**2+F(2)**2+F(3)**2)**3 !adds solvent correction

        if (approx==0) then
                write(*,*) "You must choose an approximation to optimise the system. Stop"
                stop
        else if(approx==1) then
                E_fdb=E_mu
        else if(approx==2) then
                E_fdb=E_alpha
        else if(approx==3) then
                continue
        else
                write(*,*) "I dont know which approximation to use in FDB. Stop"
                stop
        end if
        return
End
!-------------------------------------!
Subroutine readvalues()
        implicit none
        ! Intrinsic values of the system
        double precision :: E_0,E_p,E_f,E_r
        double precision, dimension(3) :: mu,mu_p,mu_f,mu_r
        double precision, dimension(3,3) :: alpha,alpha_p,alpha_f,alpha_r,alpha_tmp
        double precision, dimension(3,3,3) :: beta,beta_p,beta_f,beta_r

        character*80 hermenegildo
        ! Values declared by the user itself
        character*10 ET
        character*3 axis_scan
        character*6 approximation
        double precision x0,y0,z0,potential
        double precision R,c_i,target_barrier
        integer d,approx,e_ox,e_red,coef
        ! Dummy variables for loops
        integer i,j,k
        ! Declaring common variables for subroutines so they do not have to be declared
        COMMON /values/ E_0,mu,alpha,beta,mu_f,alpha_f,beta_f
        COMMON /properties/ target_barrier,x0,y0,z0,R,c_i,d,potential
        COMMON /keywords/ approximation,axis_scan,ET
        ! Reading input values from the file
        read(*,*) hermenegildo
        read(*,*) axis_scan,target_barrier
        read(*,*) approximation
        read(*,*) ET,potential
        write(*,*) "Prepotential", potential

                ! Consideration of oxidation of reductions
                coef=0
        e_ox=index(ET,"oxidation")
        e_ox=e_ox+index(ET,"Oxidation")
        e_ox=e_ox+index(ET,"Ox")
        e_ox=e_ox+index(ET,"ox")
                if(e_ox.gt.0) coef=-1
        e_red=index(ET,"reduction")
        e_red=e_red+index(ET,"Reduction")
        e_red=e_red+index(ET,"Red")
        e_red=e_red+index(ET,"red")
                if(e_red.gt.0) coef=1

        potential=(-potential*9.6485d4)/4184.0d0-98.6991 ! Conversion from volts to kcal/mol
                                                         ! 98.6991 ----> Conversion of -4.28eV to kcal/mol
        potential=potential*coef/627.51d0 ! Now "potential" is the energy correction in A.U to be added in the FDB
        write(*,*) "Postpotential", potential
        read(*,*) x0,y0,z0;x0=x0*1.0d-4;y0=y0*1.0d-4;z0=z0*1.0d-4
        read(*,*) R,c_i
        R=R*1.0d-4
        read(*,*) d
        ! Echoing the input
        write(*,*) ! Just so it can be better screenshot
        write(*,*) "The input for this run is:"
        write(*,'(" Central point",xF8.3,xF8.3,xF8.3)') x0,y0,z0
        write(*,'(" Radius ",f7.3," and confidence ",f4.1,"%")') R,c_i
        write(*,'(" Number of points ",i10," (1D) ",i10," (2D) ",i10," (3D)")') d,d**2,d**3
        write(*,'(" The computed potential is: ",F8.3)') potential
        if (d > 1E3) write(*,*) "WARNING: THE GRID MIGHT BE TOO DENSE FOR 3D"
        if (d < 1E2) write(*,*) "WARNING: THE GRID MIGHT BE TOO LOW FOR EACH AXIS"
        read(*,*) hermenegildo; read(*,*) hermenegildo; read(*,*) hermenegildo
!---------------------------------REACTANTS-------------------------------------!
        read(*,*) E_r
        read(*,*) hermenegildo
        read(*,*) mu_r(1),mu_r(2),mu_r(3)
        read(*,*) hermenegildo
        read(*,*) alpha_r(1,1),alpha_r(1,2),alpha_r(2,2),alpha_r(1,3),alpha_r(2,3)
        read(*,*) alpha_r(3,3)
                ! Conversion of alpha to its transposed values
                alpha_r(2,1)=alpha_r(1,2);alpha_r(3,1)=alpha_r(1,3);alpha_r(3,2)=alpha(2,3)

        read(*,*) hermenegildo
        read(*,*) beta_r(1,1,1),beta_r(1,1,2),beta_r(1,2,2),beta_r(2,2,2),beta_r(1,1,3)
        read(*,*) beta_r(1,2,3),beta_r(2,2,3),beta_r(1,3,3),beta_r(2,3,3),beta_r(3,3,3)
                ! Conversion of beta to its transposed values
                                        !xxy
                beta_r(2,1,1)=beta_r(1,1,2);beta_r(1,2,1)=beta_r(1,1,2)
                                        !xyy
                beta_r(2,1,2)=beta_r(1,2,2);beta_r(2,2,1)=beta_r(1,2,2)
                                        !xxz
                beta_r(3,1,1)=beta_r(1,1,3);beta_r(1,3,1)=beta_r(1,1,3)
                                        !yyz
                beta_r(3,2,2)=beta_r(2,2,3);beta_r(2,3,2)=beta_r(2,2,3)
                                        !xzz
                beta_r(3,1,3)=beta_r(1,3,3);beta_r(3,3,1)=beta_r(1,3,3)
                                        !yzz
                beta_r(3,3,2)=beta_r(2,3,3);beta_r(3,2,3)=beta_r(2,3,3)
                                        !xyz
                beta_r(1,3,2)=beta_r(1,2,3);beta_r(2,1,3)=beta_r(1,2,3)
                beta_r(2,3,1)=beta_r(1,2,3);beta_r(3,1,2)=beta_r(1,2,3)
                beta_r(3,2,1)=beta_r(1,2,3)
        read(*,*) hermenegildo
        read(*,*) alpha_tmp(1,1),alpha_tmp(1,2),alpha_tmp(1,3)
        read(*,*) alpha_tmp(2,1),alpha_tmp(2,2),alpha_tmp(2,3)
        read(*,*) alpha_tmp(3,1),alpha_tmp(3,2),alpha_tmp(3,3)
        do i=1,3
                do j=1,3
                        alpha_r(i,j)=alpha_r(i,j)+alpha_tmp(i,j)
                end do
        end do
       ! write(*,*) "=========================================="
       ! write(*,*) "Properties for reactants (R)"
       ! write(*,*) "------------------------------------------"
       ! write(*,*) "Dipole moment"
       ! write(*,*) (mu_r(i),i=1,3)
       ! write(*,*) "Polarizability matrix"
       ! do i=1,3
       !         write(*,*) (alpha_r(i,j),j=1,3)
       ! end do
       ! write(*,*) "First hyperpolarizability tensor"
       ! write(*,*) "X _ _"
       ! do j=1,3
       !         write(*,*) (beta_r(1,j,k),k=1,3)
       ! end do
       ! write(*,*) "Y _ _"
       ! do j=1,3
       !         write(*,*) (beta_r(2,j,k),k=1,3)
       ! end do
       ! write(*,*) "Z _ _"
       ! do j=1,3
       !         write(*,*) (beta_r(3,j,k),k=1,3)
       ! end do
       ! write(*,*) "=========================================="
!---------------------------------PRODUCTS------------------------------!
        read(*,*) hermenegildo; read(*,*) hermenegildo
        read(*,*) E_p
        read(*,*) hermenegildo
        read(*,*) mu_p(1),mu_p(2),mu_p(3)
        read(*,*) hermenegildo
        read(*,*) alpha_p(1,1),alpha_p(1,2),alpha_p(2,2),alpha_p(1,3),alpha_p(2,3)
        read(*,*) alpha_p(3,3)
                ! Conversion of alpha to its transposed values
                alpha_p(2,1)=alpha_p(1,2);alpha_p(3,1)=alpha_p(1,3);alpha_p(3,2)=alpha_p(2,3)
        read(*,*) hermenegildo
        read(*,*) beta_p(1,1,1),beta_p(1,1,2),beta_p(1,2,2),beta_p(2,2,2),beta_p(1,1,3)
        read(*,*) beta_p(1,2,3),beta_p(2,2,3),beta_p(1,3,3),beta_p(2,3,3),beta_p(3,3,3)
                ! Conversion of beta to its transposed values
                                        !xxy
                beta_p(2,1,1)=beta_p(1,1,2);beta_p(1,2,1)=beta_p(1,1,2)
                                        !xyy
                beta_p(2,1,2)=beta_p(1,2,2);beta_p(2,2,1)=beta_p(1,2,2)
                                        !xxz
                beta_p(3,1,1)=beta_p(1,1,3);beta_p(1,3,1)=beta_p(1,1,3)
                                        !yyz
                beta_p(3,2,2)=beta_p(2,2,3);beta_p(2,3,2)=beta_p(2,2,3)
                                        !xzz
                beta_p(3,1,3)=beta_p(1,3,3);beta_p(3,3,1)=beta_p(1,3,3)
                                        !yzz
                beta_p(3,3,2)=beta_p(2,3,3);beta_p(3,2,3)=beta_p(2,3,3)
                                        !xyz
                beta_p(1,3,2)=beta_p(1,2,3);beta_p(2,1,3)=beta_p(1,2,3)
                beta_p(2,3,1)=beta_p(1,2,3);beta_p(3,1,2)=beta_p(1,2,3)
                beta_p(3,2,1)=beta_p(1,2,3)
        read(*,*) hermenegildo
        read(*,*) alpha_tmp(1,1),alpha_tmp(1,2),alpha_tmp(1,3)
        read(*,*) alpha_tmp(2,1),alpha_tmp(2,2),alpha_tmp(2,3)
        read(*,*) alpha_tmp(3,1),alpha_tmp(3,2),alpha_tmp(3,3)
        do i=1,3
                do j=1,3
                        alpha_p(i,j)=alpha_p(i,j)+alpha_tmp(i,j)
                end do
        end do
       ! write(*,*)
       ! write(*,*) "Properties for products (P)"
       ! write(*,*) "------------------------------------------"
       ! write(*,*) "Dipole moment"
       ! write(*,*) (mu_p(i),i=1,3)
       ! write(*,*) "Polarizability matrix"
       ! do i=1,3
       !         write(*,*) (alpha_p(i,j),j=1,3)
       ! end do
       ! write(*,*) "First hyperpolarizability tensor"
       ! write(*,*) "X _ _"
       ! do j=1,3
       !         write(*,*) (beta_p(1,j,k),k=1,3)
       ! end do
       ! write(*,*) "Y _ _"
       ! do j=1,3
       !         write(*,*) (beta_p(2,j,k),k=1,3)
       ! end do
       ! write(*,*) "Z _ _"
       ! do j=1,3
       !         write(*,*) (beta_p(3,j,k),k=1,3)
       ! end do
       ! write(*,*) "=========================================="
!---------------------------LIGANDS----------------------------------------!
        read(*,*) hermenegildo; read(*,*) hermenegildo
        read(*,*) E_f
        read(*,*) hermenegildo
        read(*,*) mu_f(1),mu_f(2),mu_f(3)
                ! Unique conversion for ligands
                                mu_f(1)=mu_f(3);mu_f(2)=mu_f(3)
        read(*,*) hermenegildo
        read(*,*) alpha_f(1,1),alpha_f(1,2),alpha_f(2,2),alpha_f(1,3),alpha_f(2,3)
        read(*,*) alpha_f(3,3)
                ! Unique conversion for ligands
                alpha_f(1,1)=alpha_f(3,3);alpha_f(2,2)=alpha_f(3,3)
                ! Conversion of alpha to its transposed values
                alpha_f(2,1)=alpha_f(1,2);alpha_f(3,1)=alpha_f(1,3);alpha_f(3,2)=alpha_f(2,3)
        read(*,*) hermenegildo
        read(*,*) beta_f(1,1,1),beta_f(1,1,2),beta_f(1,2,2),beta_f(2,2,2),beta_f(1,1,3)
        read(*,*) beta_f(1,2,3),beta_f(2,2,3),beta_f(1,3,3),beta_f(2,3,3),beta_f(3,3,3)
                ! Unique conversion for ligands
                beta_f(1,1,1)=beta_f(3,3,3);beta_f(2,2,2)=beta_f(3,3,3)
                ! Conversion of beta to its transposed values
                                        !xxy
                beta_f(2,1,1)=beta_f(1,1,2);beta_f(1,2,1)=beta_f(1,1,2)
                                        !xyy
                beta_f(2,1,2)=beta_f(1,2,2);beta_f(2,2,1)=beta_f(1,2,2)
                                        !xxz
                beta_f(3,1,1)=beta_f(1,1,3);beta_f(1,3,1)=beta_f(1,1,3)
                                        !yyz
                beta_f(3,2,2)=beta_f(2,2,3);beta_f(2,3,2)=beta_f(2,2,3)
                                        !xzz
                beta_f(3,1,3)=beta_f(1,3,3);beta_f(3,3,1)=beta_f(1,3,3)
                                        !yzz
                beta_f(3,3,2)=beta_f(2,3,3);beta_f(3,2,3)=beta_f(2,3,3)
                                        !xyz
                beta_f(1,3,2)=beta_f(1,2,3);beta_f(2,1,3)=beta_f(1,2,3)
                beta_f(2,3,1)=beta_f(1,2,3);beta_f(3,1,2)=beta_f(1,2,3)
                beta_f(3,2,1)=beta_f(1,2,3)
        read(*,*) hermenegildo
        read(*,*) alpha_tmp(1,1),alpha_tmp(1,2),alpha_tmp(1,3)
        read(*,*) alpha_tmp(2,1),alpha_tmp(2,2),alpha_tmp(2,3)
        read(*,*) alpha_tmp(3,1),alpha_tmp(3,2),alpha_tmp(3,3)
        do i=1,3
                do j=1,3
                        alpha_f(i,j)=alpha_f(i,j)+alpha_tmp(i,j)
                end do
        end do
       ! write(*,*)
       ! write(*,*) "Properties for fakes (F)"
       ! write(*,*) "------------------------------------------"
       ! write(*,*) "Dipole moment"
       ! write(*,*) (mu_f(i),i=1,3)
       ! write(*,*) "Polarizability matrix"
       ! do i=1,3
       !         write(*,*) (alpha_f(i,j),j=1,3)
       ! end do
       ! write(*,*) "First hyperpolarizability tensor"
       ! write(*,*) "X _ _"
       ! do j=1,3
       !         write(*,*) (beta_f(1,j,k),k=1,3)
       ! end do
       ! write(*,*) "Y _ _"
       ! do j=1,3
       !         write(*,*) (beta_f(2,j,k),k=1,3)
       ! end do
       ! write(*,*) "Z _ _"
       ! do j=1,3
       !         write(*,*) (beta_f(3,j,k),k=1,3)
       ! end do
       ! write(*,*) "=========================================="
!-----------------------------SYSTEM PROPERTIES-------------------------!
        mu=0.0; alpha=0.0; beta=0.0
        ! This very next line must be modified by hand to consider the stoichiometry of the reaction
        E_0=E_p+E_f-E_r !in atomic units
        do i=1,3
                mu(i)=mu_p(i)-mu_r(i)
                do j=1,3
                        alpha(i,j)=alpha_p(i,j)-alpha_r(i,j)
                        do k=1,3
                                beta(i,j,k)=beta_p(i,j,k)-beta_r(i,j,k)
                        end do
                end do
        end do
        write(*,*)
        write(*,*) "============================================="
        write(*,*) "Global properties of the system"
        write(*,*) "------------------------------------------"
        write(*,'(x"Gibbs free energy (dG)",xF9.4,x"kcal/mol")') E_0*6.2751d2
        write(*,*) "Dipole moment vector"
        write(*,'(xE11.4,xxxxxE11.4,xxxxxE11.4)') (mu(i),i=1,3)
        write(*,*) "Polarizability matrix"
        do i=1,3
                write(*,'(xE11.4,xxxxxE11.4,xxxxxE11.4)') (alpha(i,j),j=1,3)
        end do
        write(*,*) "First hyperpolarizability tensor"
        write(*,*) "X _ _"
        do j=1,3
                write(*,'(xE11.4,xxxxxE11.4,xxxxxE11.4)') (beta(1,j,k),k=1,3)
        end do
        write(*,*) "Y _ _"
        do j=1,3
                write(*,'(xE11.4,xxxxxE11.4,xxxxxE11.4)') (beta(2,j,k),k=1,3)
        end do
        write(*,*) "Z _ _"
        do j=1,3
                write(*,'(xE11.4,xxxxxE11.4,xxxxxE11.4)') (beta(3,j,k),k=1,3)
        end do
        write(*,*) "============================================="
End
!--------------------------------------------------------------!
Subroutine minimize()
        implicit none
        double precision target_barrier
        double precision x0,y0,z0
        double precision R,c_i,E_0,potential
        integer d
        integer x_axis,y_axis,z_axis,dirscan
        integer K_dipole,K_alpha,K_beta,approx
        character*10 ET
        character*3 axis_scan
        character*6 approximation
        COMMON /properties/ target_barrier,x0,y0,z0,R,c_i,d,potential
        COMMON /keywords/ approximation,axis_scan,ET

                ! To know what approximation will be used (dipole moment, alpha or beta)
        K_dipole=index(approximation,"dipole")
        K_dipole=K_dipole+index(approximation,"Dipole")
                if (K_dipole.gt.0) approx=1
        K_alpha=index(approximation,"alpha")
        K_alpha=K_alpha+index(approximation,"Alpha")
                if (K_alpha.gt.0) approx=2
        K_beta=index(approximation,"beta")
        K_beta=K_beta+index(approximation,"Beta")
                if (K_beta.gt.0) approx=3

                ! To know how many scans will be performed
        x_axis=scan(axis_scan,"x")
        x_axis=x_axis+scan(axis_scan,"X")
        if (x_axis<0) then
                x_axis=1
        end if
        y_axis=scan(axis_scan,"y")
        y_axis=y_axis+scan(axis_scan,"Y")
        if (y_axis.gt.0) then
                y_axis=1
        end if
        z_axis=scan(axis_scan,"z")
        z_axis=z_axis+scan(axis_scan,"Z")
        if (z_axis.gt.0) then
                z_axis=1
        end if
        !write(*,*) x_axis,y_axis,z_axis
        dirscan=x_axis+y_axis+z_axis

                ! Scan selector
        if(dirscan==1) then
                call min_1D(approx,x_axis,y_axis,z_axis)
        else if(dirscan==2) then
                if (x_axis==0) then
                   z_axis=0
                   call min_1D(approx,x_axis,y_axis,z_axis)
                   y_axis=0;z_axis=1
                   call min_1D(approx,x_axis,y_axis,z_axis)
                   z_axis=1;y_axis=1
                   call min_2D(approx,x_axis,y_axis,z_axis)
                else if (y_axis==0) then
                   z_axis=0
                   call min_1D(approx,x_axis,y_axis,z_axis)
                   x_axis=0;z_axis=1
                   call min_1D(approx,x_axis,y_axis,z_axis)
                   x_axis=1;z_axis=1
                   call min_2D(approx,x_axis,y_axis,z_axis)
                else if (z_axis==0) then
                   y_axis=0
                   call min_1D(approx,x_axis,y_axis,z_axis)
                   y_axis=1;x_axis=0
                   call min_1D(approx,x_axis,y_axis,z_axis)
                   x_axis=1;y_axis=1
                   call min_2D(approx,x_axis,y_axis,z_axis)
                end if
        else if(dirscan==3) then ! in order
                y_axis=0;z_axis=0 ! X scan
                !call min_1D(approx,x_axis,y_axis,z_axis)
                y_axis=1;x_axis=0 ! Y scan
                !call min_1D(approx,x_axis,y_axis,z_axis)
                y_axis=0;z_axis=1 ! Z scan
                !call min_1D(approx,x_axis,y_axis,z_axis)
                x_axis=1 ! XZ scan
                !call min_2D(approx,x_axis,y_axis,z_axis)
                x_axis=0;y_axis=1 ! YZ scan
                !call min_2D(approx,x_axis,y_axis,z_axis)
                x_axis=1;z_axis=0  ! XY scan
                !call min_2D(approx,x_axis,y_axis,z_axis)
                z_axis=1 ! XYZ scan
                !call min_3D(approx,x_axis,y_axis,z_axis)
                ! Locate the minimum field for a given barrier
                call minfield_3D(approx,x_axis,y_axis,z_axis)
        else
                write(*,*) "WARNING! BAD INPUT!"
                stop
        end if
End
!--------------------------------------------------------------!
Subroutine min_1D(approx,x_axis,y_axis,z_axis)
        implicit none
        character*10 ET
        character*3 axis_scan
        character*6 approximation
        double precision distance,distance_save
        double precision, dimension(3) :: F
        double precision, dimension(0:10,4) :: E_min
        double precision E_fdb,E_0
        integer x_axis,y_axis,z_axis,approx
        double precision x0,y0,z0
        double precision c_i,target_barrier,R,R_amp
        double precision step_i,tolerance,compare,potential
        integer i,p,pp,d,min_zone
        COMMON /properties/ target_barrier,x0,y0,z0,R,c_i,d,potential
        COMMON /keywords/ approximation,axis_scan,ET

        compare=R*1.0d-1 ! compare= minimum separation to consider two basins as different minimum regions
        R_amp=R*(1+c_i*0.01);F=0.0;E_min=9.9e9;tolerance=0.1 !(in kcal/mol)

        ! Initialization of F
        F(1)=x0+x_axis*(x0-R_amp)
        F(2)=y0+y_axis*(y0-R_amp)
        F(3)=z0+z_axis*(z0-R_amp)

        ! Inicialization of the scan
        step_i=2*abs(x0-R_amp)/(d-1)
        p=0;pp=0 ! Points in the grid, positive values of E_fdb
        min_zone=0 ! Indexes for E_min matrix
        write(*,*)
        write(*,*) "                    1 D RESULTS",x_axis,y_axis,z_axis
        write(*,*) "========================================================================================================"
        do i=0,(d-1)
                if(distance(F).le.R) then
                        p=p+1
                        !                        "Visual" scan
                        !write(*,*) "E_fdb",E_fdb(F,approx)*6.2751d2,"F:",F(1)*1.0d4,F(2)*1.0d4,F(3)*1.0d4
                        !write(*,*) abs(E_fdb(F,approx)*6.2751d2 - target_barrier)
                        if(abs(E_fdb(F,approx)*6.2751d2 - target_barrier).lt.tolerance) then
                            if (E_fdb(F,approx).ge.0.0d0) then
                                pp=pp+1

                                !change the zone (basin) where you are locating
                                !the minima with threshold for considering two
                                !basins as different called "compare". (change
                                !it at the beggining if you want to).
                                if (abs(F(1)-E_min(min_zone,2)).gt. compare .or. abs(F(2)-E_min(min_zone,3)) &
                                &.gt. compare .or. abs(F(3)-E_min(min_zone,4)).gt. compare) min_zone=min_zone+1

                                if (abs(E_min(min_zone,1)-target_barrier) > abs(E_fdb(F,approx)-target_barrier)) then
                                      !within a zone, gets the minimum energy.
                                      E_min(min_zone,1)=E_fdb(F,approx)
                                      E_min(min_zone,2)=F(1)
                                      E_min(min_zone,3)=F(2)
                                      E_min(min_zone,4)=F(3)
                                      distance_save=distance(F)
                                end if
                            end if
                        end if
                end if
                !x_ y_ and z_ axis are pseudo-logical parameters to select the proper
                !direction so that only one direction will be updated
                F(1)=F(1)+x_axis*(step_i)
                F(2)=F(2)+y_axis*(step_i)
                F(3)=F(3)+z_axis*(step_i)
        end do
        write(*,'(I8,x"points have been evaluated (",I3,x"of them meet the criteria)")'),p,pp
        write(*,*) "========================================================================================================"
        if(pp.gt.0) then
        write(*,'(x"Ideal field for"x,F7.4,x"kcal/mol for a density of points"x,I5)') target_barrier,d
        do i=1,min_zone
        write(*,'(x"F_x",ES11.3,x"F_y",ES11.3,x"F_z",ES11.3,x"Strenght",F11.5,x"u.a.",xx"Energy (kcal/mol)",F14.8)') &
        & E_min(i,2)*1.0d4,E_min(i,3)*1.0d4,E_min(i,4)*1.0d4,sqrt((E_min(i,2)*1.0d4)**2+(E_min(i,3)*1.0d4)**2 + &
        & (E_min(i,4)*1.0d4)**2),E_min(i,1)*6.2751d2
        end do
        write(*,*) "========================================================================================================"
        write(*,*) ! Just so it can be better screenshot
        end if
End
!--------------------------------------------------------------!
Subroutine min_2D(approx,x_axis,y_axis,z_axis)
        implicit none
        integer,parameter :: nbasins=10
        character*10 ET
        character*3 axis_scan
        character*6 approximation
        double precision, dimension(3) :: F
        double precision, dimension(0:nbasins,4) :: E_min
        double precision distance,distance_save,E_fdb,E_0,currentE
        double precision step_i,step_j,tolerance,compare
        double precision x0,y0,z0,c_i,target_barrier,R,R_amp,potential
        integer i,j,p,pp,d,min_zone
        integer x_axis,y_axis,z_axis,approx
        COMMON /properties/ target_barrier,x0,y0,z0,R,c_i,d,potential
        COMMON /keywords/ approximation,axis_scan,ET

        compare=R*1.0d-1        !compare= minimum separation to consider two basins as different minimum regions
        R_amp=R*(1+c_i*0.01);F=0.0;E_min=9.9e9;tolerance=0.1 !(in kcal/mol)

        write(*,*) "                    2 D ANALYSIS"
        write(*,*) "========================================================================================================"
        99      continue
        ! Inicialization of F
        F(1)=x0+x_axis*(x0-R_amp)
        F(2)=y0+y_axis*(y0-R_amp)
        F(3)=z0+z_axis*(z0-R_amp)
        ! Inicialization of the scan. Calculation of the stepsize
        if (z_axis==0) then !x_axis & y_axis == 1
                step_i=2*abs(x0-R_amp)/(d-1)
                step_j=2*abs(y0-R_amp)/(d-1)
        else if (y_axis==0) then !x_axis & z_axis == 1
                step_i=2*abs(x0-R_amp)/(d-1)
                step_j=2*abs(z0-R_amp)/(d-1)
        else if (x_axis==0) then !y_axis & z_axis == 1
                step_i=2*abs(y0-R_amp)/(d-1)
                step_j=2*abs(z0-R_amp)/(d-1)
        else
                write(*,*) "Unable to detect the two axes on 2D minimization.Stop"
                stop
        end if
        p=0;pp=0 !points in the grid, positive values of Efdb
        min_zone=0 !indexes for E_min matrix the i index
        do i=0,d-1
                do j=0,d-1
                        currentE=E_fdb(F,approx)*6.2751d2
                        if(distance(F).le.R) then
                                p=p+1
                        !                        "Visual" scan
                        !"E_fdb",currentE,"F:",F(1)*1.0d4,F(2)*1.0d4,F(3)*1.0d4
                        !write(*,*) abs(currentE - target_barrier)
                                if(abs(currentE - target_barrier).lt.tolerance) then
                                        if (currentE.ge.0.0d0) then
                                                pp=pp+1
                                                !change the zone (basin) where you are locating
                                                !the minima with threshold for considering two
                                                !basins as different called "compare". (change
                                                !it at the beggining if you want to).
                                                if (abs(F(1)-E_min(min_zone,2)).gt.compare .or. abs(F(2)-E_min(min_zone,3)) &
                                                & .gt.compare .or. abs(F(3)-E_min(min_zone,4)).gt.compare) then ! Change basin
                                                        min_zone=min_zone+1
                                                        if (min_zone.gt.nbasins)then ! Must restart and increase the "compare"
                                                                compare=compare*2.0d0
                                                                write(*,*) "Too many basins"
                                                                write(*,*) "The &
& separation within basins is increased to repeat the analysis"
                                                                write(*,'("New separation is",xES11.4," a.u &
&of field strenght")') compare
                                                                write(*,*)
                                                                goto 99
                                                        end if

                                                        if (abs(E_min(min_zone,1)-target_barrier).gt. &
                                                        & abs(currentE-target_barrier)) then  ! Within a zone, get the minimum energy.
                                                                E_min(min_zone,1)=currentE
                                                                E_min(min_zone,2)=F(1)
                                                                E_min(min_zone,3)=F(2)
                                                                E_min(min_zone,4)=F(3)
                                                                distance_save=distance(F)
                                                        end if
                                                end if
                                        end if
                                end if
                        end if
                        !x_ y_ and z_ axis are pseudo-logical parameters to select the proper
                        !direction so that only one direction will be updated. Update
                        !the one associated with j dummy index (second coordinate)
                        if (z_axis==0) then
                                F(2)=F(2)+y_axis*(step_j)
                        else
                                F(3)=F(3)+z_axis*(step_j)
                        end if
                end do
                !x_ y_ and z_ axis are pseudo-logical parameters to select the proper
                !direction so that only one direction will be updated. Update the one
                !associated with i dummy index (first coordinate)
                if (x_axis==0) then !YZ scan
                     F(2)=F(2)+step_i
                     !restart the other field (Z)
                     F(3)=z0+z_axis*(z0-R_amp)
                else !XY or XZ scan
                     F(1)=F(1)+step_i
                     !restart the other field (Y or Z)
                     if (y_axis==1) then
                       F(2)=y0+y_axis*(y0-R_amp)
                     else !restart Z
                       F(3)=z0+z_axis*(z0-R_amp)
                     end if
                end if

        end do
        write(*,*) "                    2 D RESULTS",x_axis,y_axis,z_axis
        write(*,*) "========================================================================================================"
        write(*,'(I8,x"points have been evaluated (",I8,x"of them meet the criteria)")'),p,pp
        write(*,*) "========================================================================================================"
        if(pp.gt.0) then
        write(*,'(x"Ideal field for"x,F7.4,x"kcal/mol for a density of points",xI4)') target_barrier,d
        do i=1,min_zone
        write(*,'(x"F_x",ES11.3,x"F_y",ES11.3,x"F_z",ES11.3,x"Strenght",F11.5,x"u.a.",xx"Energy (kcal/mol)",F14.8)') &
        & E_min(i,2)*1.0d4,E_min(i,3)*1.0d4,E_min(i,4)*1.0d4,sqrt((E_min(i,2)*1.0d4)**2+(E_min(i,3)*1.0d4)**2 + &
        & (E_min(i,4)*1.0d4)**2),E_min(i,1)
        end do
        write(*,*) "========================================================================================================"
        write(*,*) ! Just so it can be better screenshot
        end if
End

!----------------------------------------------------------------------!
Subroutine min_3D(approx,x_axis,y_axis,z_axis)
        implicit none
        integer,parameter :: nbasins=10
        character*10 ET
        character*3 axis_scan
        character*6 approximation
        double precision, dimension(3) :: F
        double precision, dimension(0:nbasins,4) :: E_min
        double precision distance,distance_save,E_fdb,E_0,currentE
        double precision step_i,step_j,step_k,tolerance,compare
        double precision x0,y0,z0,c_i,target_barrier,R,R_amp,potential
        integer i,j,k,p,pp,d,min_zone
        integer x_axis,y_axis,z_axis,approx
        COMMON /properties/ target_barrier,x0,y0,z0,R,c_i,d,potential
        COMMON /keywords/ approximation,axis_scan,ET

        compare=R*1.0d-1        !compare= minimum separation to consider two basins as different minimum regions
        R_amp=R*(1+c_i*0.01);F=0.0;E_min=9.9e9;tolerance=0.1 !(in kcal/mol)

        write(*,*) "                    3 D ANALYSIS"
        write(*,*) "========================================================================================================"
        99      continue
        ! Inicialization of F
        F(1)=x0+x_axis*(x0-R_amp)
        F(2)=y0+y_axis*(y0-R_amp)
        F(3)=z0+z_axis*(z0-R_amp)
        ! Inicialization of the scan. Calculation of the stepsize
        step_i=2*abs(x0-R_amp)/(d-1)
        step_j=2*abs(y0-R_amp)/(d-1)
        step_k=2*abs(z0-R_amp)/(d-1)

        p=0;pp=0 !points in the grid, positive values of Efdb
        min_zone=0 !indexes for E_min matrix the i index
        do i=0,d-1
                do j=0,d-1
                        do k=0,d-1
                                currentE=E_fdb(F,approx)*6.2751d2
                                if(distance(F).le.R) then
                                        p=p+1
                                        !                        "Visual" scan
                                        !"E_fdb",currentE,"F:",F(1)*1.0d4,F(2)*1.0d4,F(3)*1.0d4
                                        if(abs(currentE - target_barrier).lt.tolerance) then
                                                if (currentE.ge.0.0d0) then
                                                        pp=pp+1
                                                        !change the zone (basin) where you are locating
                                                        !the minima with threshold for considering two
                                                        !basins as different called "compare"
                                                        !(change it at the beggining if you want to)
                                                        if (abs(F(1)-E_min(min_zone,2)).gt.compare &
                                                        & .or. abs(F(2)-E_min(min_zone,3)).gt.compare &
                                                        & .or.abs(F(3)-E_min(min_zone,4)).gt.compare) then ! Change basin
                                                                min_zone=min_zone+1
                                                                if(min_zone.gt.nbasins)then ! Must restart and increase the "compare"
                                                                        compare=compare*2.0d0
                                                                        write(*,*) "Too many basins"
                                                                write(*,*)  "The separation within basins &
                                                                & is increased to repeat the analysis"
                                                                        write(*,'(" New separation is",xES11.4," a.u &
                                                                        & of field strenght")') compare
                                                                        write(*,*)
                                                                        goto 99
                                                                end if
                                                                if(abs(E_min(min_zone,1)-target_barrier).gt. &
                                                                & abs(currentE-target_barrier)) then  ! Within a zone, get the minimum energy.
                                                                        E_min(min_zone,1)=currentE
                                                                        E_min(min_zone,2)=F(1)
                                                                        E_min(min_zone,3)=F(2)
                                                                        E_min(min_zone,4)=F(3)
                                                                        distance_save=distance(F)
                                                                end if
                                                        end if
                                                end if
                                        end if
                                end if
                                F(3)=F(3)+z_axis*(step_k)
                        end do
                        F(2)=F(2)+y_axis*(step_j)
                        F(3)=z0+z_axis*(z0-R_amp)
                end do
                F(1)=F(1)+x_axis*(step_i)
                F(2)=y0+y_axis*(y0-R_amp)
                F(3)=z0+z_axis*(z0-R_amp)
        end do
        write(*,*) "                    3 D RESULTS"
        write(*,*) "========================================================================================================"
        write(*,'(I10,x"points have been evaluated (",I8,x"of them meet the criteria)")'),p,pp
        write(*,*) "========================================================================================================"
        if(pp.gt.0) then
        write(*,'(x"Ideal field for"x,F7.4,x"kcal/mol for a density of points",xI4)') target_barrier,d
        do i=1,min_zone
        write(*,'(x"F_x",ES11.3,x"F_y",ES11.3,x"F_z",ES11.3,x"Strenght",F11.5,x"u.a.",xx"Energy (kcal/mol)",F14.8)') &
        & E_min(i,2)*1.0d4,E_min(i,3)*1.0d4,E_min(i,4)*1.0d4,sqrt((E_min(i,2))**2+(E_min(i,3))**2 + &
        & (E_min(i,4))**2),E_min(i,1)
        end do
        write(*,*) "========================================================================================================"
        write(*,*) ! Just so it can be better screenshot
        end if
End
!---------------------------------------------
Subroutine minfield_3D(approx,x_axis,y_axis,z_axis)
!for the targeted barrier. Search the lowest possible F-strength.
        implicit none
        integer,parameter :: nbasins=10
        character*10 ET
        character*3 axis_scan
        character*6 approximation
        double precision, dimension(3) :: F
        double precision, dimension(4) :: E_min
        double precision distance,distance_save,E_fdb,E_0,currentE
        double precision step_i,step_j,step_k,tolerance,compare
        double precisionx0,y0,z0,c_i,target_barrier,R,R_amp,min_module
        double precision potential
        integer i,j,k,p,pp,d
        integer x_axis,y_axis,z_axis,approx
        COMMON /properties/ target_barrier,x0,y0,z0,R,c_i,d,potential
        COMMON /keywords/ approximation,axis_scan,ET

        compare=R*1.0d-1        !compare= minimum separation to consider two basins as different minimum regions
        R_amp=R*(1+c_i*0.01);F=0.0;E_min=9.9e9;tolerance=0.1 !(in kcal/mol)

        write(*,*) "                    3 D ANALYSIS MINIMUM FIELD"
        write(*,*) "========================================================================================================"
        ! Inicialization of F
        F(1)=x0+x_axis*(x0-R_amp)
        F(2)=y0+y_axis*(y0-R_amp)
        F(3)=z0+z_axis*(z0-R_amp)
        ! Inicialization of the scan. Calculation of the stepsize
        step_i=2*abs(x0-R_amp)/(d-1)
        step_j=2*abs(y0-R_amp)/(d-1)
        step_k=2*abs(z0-R_amp)/(d-1)

        p=0;pp=0 !points in the grid, positive values of Efdb
        min_module=1E9
        do i=0,d-1
                do j=0,d-1
                        do k=0,d-1
                                currentE=E_fdb(F,approx)*6.2751d2
                                if(distance(F).le.R) then
                                        p=p+1
                                        !                        "Visual" scan
                                        !"E_fdb",currentE,"F:",F(1)*1.0d4,F(2)*1.0d4,F(3)*1.0d4
                                        if(abs(currentE - target_barrier).lt.tolerance) then
                                                if (currentE.ge.0.0d0) then
                                                     if(sqrt(F(1)**2+F(2)**2+F(3)**2)<min_module) then
                                                          E_min(1)=currentE
                                                          E_min(2)=F(1)
                                                          E_min(3)=F(2)
                                                          E_min(4)=F(3)
                                                          !write(*,*) (E_min(pp),pp=2,4),E_min(1)
                                                          min_module=sqrt(F(1)**2+F(2)**2+F(3)**2)
                                                     end if
                                                end if
                                        end if
                                end if
                                F(3)=F(3)+z_axis*(step_k)
                        end do
                        F(2)=F(2)+y_axis*(step_j)
                        F(3)=z0+z_axis*(z0-R_amp)
                end do
                F(1)=F(1)+x_axis*(step_i)
                F(2)=y0+y_axis*(y0-R_amp)
                F(3)=z0+z_axis*(z0-R_amp)
        end do
        write(*,*) "                    3 D RESULTS"
        write(*,*) "========================================================================================================"
        write(*,'(I10,x"points have been evaluated ")') p
        write(*,*) "========================================================================================================"
        write(*,'(x"Ideal field for"x,F7.4,x"kcal/mol for a density of points",xI4)') target_barrier,d
        write(*,'(x"F_x",ES11.3,x"F_y",ES11.3,x"F_z",ES11.3,x"Strenght",F11.5,x"u.a. (",i3.3,")",xx"Energy (kcal/mol)",F14.8)') &
        & E_min(2)*1.0d4,E_min(3)*1.0d4,E_min(4)*1.0d4,sqrt((E_min(2))**2+(E_min(3))**2+ (E_min(4))**2), &
        & int(sqrt((E_min(2))**2+(E_min(3))**2+ (E_min(4))**2)*10000),E_min(1)
        write(*,*) "========================================================================================================"
        write(*,*) ! Just so it can be better screenshot
End
!----------------------------------------------------------------------!
Function distance(F)
        implicit none
        double precision distance
        double precision potential
        double precision,dimension(3) :: F
                ! COMMON variables
        double precision target_barrier,x0,y0,z0,R,c_i
        integer d,approx
        COMMON /properties/ target_barrier,x0,y0,z0,R,c_i,d,potential

        distance=sqrt((x0-F(1))**2+(y0-F(2))**2+(z0-F(3))**2)

        return
End
