!-FDB_b optimization
!--Program written by Pau Besalu and Guillem Pey in 2022-2023
Program FDBb
double precision :: E_0,E_p,E_f,E_r,E_fdb
double precision, dimension(3) :: mu,mu_p,mu_f,mu_r,F
double precision, dimension(3,3) :: alpha,alpha_p,alpha_f,alpha_r,alpha_tmp
double precision, dimension(3,3,3) :: beta,beta_p,beta_f,beta_r
character*80 hermenegildo
character*3 key1
double precision x0,y0,z0
double precision R,c_i,d
double precision target_barrier
integer i,j,k

COMMON /values/ E_0,mu,alpha,beta
COMMON /properties/target_barrier,x0,y0,z0,R,c_i,d,key1

! If one has to modify the default stoichiometry: "/stoichiometry"
! Falta fer la tria dels scans segons la keyword key1
! Falta programar els tipus de scan

call readvalues()

End
!-------------------------------------!
Function E_fdb(F)
        implicit none
        double precision :: E_0,E_fdb
        double precision, dimension(3) :: mu,F
        double precision, dimension(3,3) :: alpha
        double precision, dimension(3,3,3) :: beta
        integer i,j,k
        COMMON /values/ E_0,mu,alpha,beta

        E_fdb=E_0
        do i=1,3
                E_fdb=E_fdb-mu(i)*F(i)
                do j=1,3
                        E_fdb=E_fdb-alpha(i,j)*F(i)*F(j)
                        do k=1,3
                                E_fdb=E_fdb-beta(i,j,k)*F(i)*F(j)*F(k)
                        end do
                end do
        end do
End
!-------------------------------------!
Subroutine readvalues()
        implicit none
! Intrinsic values of the system
        double precision :: E_0,E_p,E_f,E_r
        double precision, dimension(3) :: mu,mu_p,mu_f,mu_r
        double precision, dimension(3,3) :: alpha,alpha_p,alpha_f,alpha_r,alpha_tmp
        double precision, dimension(3,3,3) :: beta,beta_p,beta_f,beta_r
        ! Defining them as zero so they can be operated over

        character*80 hermenegildo
! Values declared by the user itself
        character*3 key1
        double precision x0,y0,z0
        double precision R,c_i,d
        double precision target_barrier
! Other
        integer i,j,k
! Declaring common variables for subroutines so they do not have to be declared
        COMMON /values/ E_0,mu,alpha,beta
        COMMON /properties/target_barrier,x0,y0,z0,R,c_i,d,key1

! Reading input values from the file
        read(*,*) hermenegildo
        read(*,*) key1,target_barrier
        read(*,*) x0,y0,z0
        read(*,*) R,c_i
        read(*,*) d
        ! Echoing the input
        write(*,*) "The input for this run is:"
        write(*,*) "Central point",x0,y0,z0
        write(*,'(" Radius ",f7.3," and confidence ",f4.1,"%")') R,c_i
        write(*,*) "Density of grid", d
        read(*,*) hermenegildo; read(*,*) hermenegildo; read(*,*) hermenegildo
        ! Reactant values
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
        ! Product values
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
        read(*,*) hermenegildo; read(*,*) hermenegildo
        read(*,*) E_f
        read(*,*) hermenegildo
        read(*,*) mu_f(1),mu_f(2),mu_f(3)
        read(*,*) hermenegildo
        read(*,*) alpha_f(1,1),alpha_f(1,2),alpha_f(2,2),alpha_f(1,3),alpha_f(2,3)
        read(*,*) alpha_f(3,3)
                ! Conversion of alpha to its transposed values
                alpha_f(2,1)=alpha_f(1,2);alpha_f(3,1)=alpha_f(1,3);alpha_f(3,2)=alpha_f(3,2)
        read(*,*) hermenegildo
        read(*,*) beta_f(1,1,1),beta_f(1,1,2),beta_f(1,2,2),beta_f(2,2,2),beta_f(1,1,3)
        read(*,*) beta_f(1,2,3),beta_f(2,2,3),beta_f(1,3,3),beta_f(2,3,3),beta_f(3,3,3)
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
        mu=0.0; alpha=0.0; beta=0.0
        ! This very next line must be modified by hand to consider the stoichiometry of the reaction
        E_0=E_p+E_f-E_r
        do i=1,3
                mu(i)=mu_p(i)+mu_f(i)-mu_r(i)
                do j=1,3
                        alpha(i,j)=alpha_p(i,j)+alpha_f(i,j)-alpha_r(i,j)
                        do k=1,3
                                beta(i,j,k)=beta_p(i,j,k)+beta_f(i,j,k)-beta_r(i,j,k)
                        end do
                end do
        end do
        write(*,*)
        write(*,*)
        write(*,*) "=========================================="
        write(*,*) "Global properties of the system"
        write(*,*) "------------------------------------------"
        write(*,*) "Gibbs free energy (dG)",E_0*627.51,"kcal/mol"
        write(*,*) "Dipole moment vector"
        write(*,*) (mu(i),i=1,3)
        write(*,*) "Polarizability matrix"
        do i=1,3
                write(*,*) (alpha(i,j),j=1,3)
        end do
        write(*,*) "First hyperpolarizability tensor"
        write(*,*) "X _ _"
        do j=1,3
                write(*,*) (beta(1,j,k),k=1,3)
        end do
        write(*,*) "Y _ _"
        do j=1,3
                write(*,*) (beta(2,j,k),k=1,3)
        end do
        write(*,*) "Z _ _"
        do j=1,3
                write(*,*) (beta(3,j,k),k=1,3)
        end do
        write(*,*) "=========================================="
End
