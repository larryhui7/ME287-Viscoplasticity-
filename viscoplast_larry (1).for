      SUBROUTINE VUMAT (
!      Read only (unmodifiable) variables :-
     +                    NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV,
     +                    NPROPS, LANNEAL, STEP_TIME, TOTAL_TIME,
     +                    DT, CMNAME, COORD_MP, CHAR_LENGTH, PROPS,
     +                    DENSITY, STRAIN_INC, REL_SPIN_INC,
     +                    TEMP_OLD, STRETCH_OLD, DEFGRAD_OLD,
     +                    FIELD_OLD, STRESS_OLD, STATE_OLD, 
     +                    ENER_INTERN_OLD, ENER_INELAS_OLD, TEMP_NEW,
     +                    STRETCH_NEW, DEFGRAD_NEW, FIELD_NEW, 
!      Read and ! write (modifiable) variables :-
     +                    STRESS_NEW, STATE_NEW, ENER_INTERN_NEW,
     +                    ENER_INELAS_NEW)
     
      INCLUDE 'VABA_PARAM.INC'
      
      DIMENSION COORD_MP(NBLOCK,*),CHAR_LENGTH(NBLOCK), PROPS(NPROPS),
     +       DENSITY(NBLOCK), STRAIN_INC(NBLOCK,NDIR+NSHR), 
     +       REL_SPIN_INC(NBLOCK,NSHR), TEMP_OLD(NBLOCK), 
     +       STRETCH_OLD(NBLOCK,NDIR+NSHR),
     +       DEFGRAD_OLD(NBLOCK,NDIR+NSHR+NSHR), 
     +       FIELD_OLD(NBLOCK,NFIELDV), STRESS_OLD(NBLOCK,NDIR+NSHR),
     +       STATE_OLD(NBLOCK,NSTATEV), ENER_INTERN_OLD(NBLOCK),
     +       ENER_INELAS_OLD(NBLOCK), TEMP_NEW(NBLOCK), 
     +       STRETCH_NEW(NBLOCK,NDIR+NSHR), 
     +       DEFGRAD_NEW(NBLOCK,NDIR+NSHR+NSHR), 
     +       FIELD_NEW(NBLOCK,NFIELDV), STRESS_NEW(NBLOCK,NDIR+NSHR),
     +       STATE_NEW(NBLOCK,NSTATEV), ENER_INTERN_NEW(NBLOCK),
     +       ENER_INELAS_NEW(NBLOCK)

      character*8 CMNAME

      integer  km




 !!!!  Don't change anything above this line.  Start modifying below.



      real*8 F_t(3,3),F_tau(3,3),Fdot(3,3),U_tau(3,3),U_inv(3,3), det_U,R_tau(3,3)
      real*8 T_tau(3,3), Ee(3,3), det_F, E_tau(3,3), Re_tau(3,3)
      real*8 Fp_t(3,3), Fp_tau(3,3), Fe(3,3), eigvecs(3,3), eigvals(3) 
      real*8 nu_p, tau,  I_1(3,3), Dp(3,3), Te_0(3,3), stretches(3,3)
      real*8 Fp_inv(3,3), Ee0(3,3), tr_Ee, dtmp, D_dev_tau(3,3)
      real*8 L_tau(3,3), D_tau(3,3), tr_D_tau, F_tau_inv(3,3)
      real*8 det_Fe, zero_m(3,3),  snake, z, Ue_inv(3,3),Tvis_tau(3,3)
      real*8 mu, kappa, dens_p,  excess, Te_sph, youngs, nu, tens_yield
      real*8 F_inv(3,3), shear_yield, visc_ext, Ep(3,3), visc
      real*8 Finc(3,3), Linc(3,3), Dinc(3,3), tr_Dinc, D_dev_inc(3,3)
      real*8 Delta_Tvis(3,3), Tvis_old(3,3), F_avg(3,3), F_avg_inv(3,3)
      parameter(ONE=1., ONE_HALF=0.5, TWO=2., ZERO=0.,
     +                ONE_THIRD=1./3., TWO_THIRD=2./3.,
     +                THREE_HALF=1.5, THREE=3.)

     
        I_1 = reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1 /), (/3,3/))
	Zero_m = reshape((/0, 0, 0, 0, 0, 0, 0, 0, 0 /), (/3,3/))
	z=0.d0

      ! Get properties defined in input file
      ! For VISCOELASTICITY (no plasticity):
      !
        youngs             		= props(01)
        nu				= props(02)
	ka			= props(03)
	ma			= props(04)
	
	mu=youngs/(2.d0*(1.d0+nu))
	kappa=youngs/(3.d0*(1.d0-2.d0*nu))
	

      do km = 1,NBLOCK

        !Copy the old and new Deformation gradients into F_T and F_tau 
        !respectively
        !
        F_T(1,1) = DEFGRAD_OLD(KM,1)
        F_T(2,2) = DEFGRAD_OLD(KM,2)
        F_T(3,3) = DEFGRAD_OLD(KM,3)
        F_T(1,2) = DEFGRAD_OLD(KM,4)

        F_tau(1,1) = DEFGRAD_NEW(KM,1)
        F_tau(2,2) = DEFGRAD_NEW(KM,2)
        F_tau(3,3) = DEFGRAD_NEW(KM,3)
        F_tau(1,2) = DEFGRAD_NEW(KM,4)

        U_tau(1,1) = STRETCH_NEW(KM,1)
        U_tau(2,2) = STRETCH_NEW(KM,2)
        U_tau(3,3) = STRETCH_NEW(KM,3)
        U_tau(1,2) = STRETCH_NEW(KM,4)

        if (NSHR .eq. 1) then
          F_T(2,3) = ZERO
          F_T(3,1) = ZERO
          F_T(2,1) = DEFGRAD_OLD(KM,5)
          F_T(3,2) = ZERO
          F_T(1,3) = ZERO

          F_tau(2,3) = ZERO
          F_tau(3,1) = ZERO
          F_tau(2,1) = DEFGRAD_NEW(KM,5)
          F_tau(3,2) = ZERO
          F_tau(1,3) = ZERO

          U_tau(2,3) = ZERO
          U_tau(3,1) = ZERO
          U_tau(2,1) = U_tau(1,2)
          U_tau(3,2) = ZERO
          U_tau(1,3) = ZERO
        else
          F_T(2,3) = DEFGRAD_OLD(KM,5)
          F_T(3,1) = DEFGRAD_OLD(KM,6)
          F_T(2,1) = DEFGRAD_OLD(KM,7)
          F_T(3,2) = DEFGRAD_OLD(KM,8)
          F_T(1,3) = DEFGRAD_OLD(KM,9)

          F_tau(2,3) = DEFGRAD_NEW(KM,5)
          F_tau(3,1) = DEFGRAD_NEW(KM,6)
          F_tau(2,1) = DEFGRAD_NEW(KM,7)
          F_tau(3,2) = DEFGRAD_NEW(KM,8)
          F_tau(1,3) = DEFGRAD_NEW(KM,9)

          U_tau(2,3) = STRETCH_NEW(KM,5)
          U_tau(3,1) = STRETCH_NEW(KM,6)
          U_tau(2,1) = U_tau(1,2)
          U_tau(3,2) = U_tau(2,3)
          U_tau(1,3) = U_tau(3,1)
        endif
	 
		 
	 if(total_time.eq.0 .and. step_time.eq.0) then
          !
          ! If first dummy step: Initialize state vars
          !
          State_old(km,1:9)   = (/1, 0, 0, 0, 1, 0, 0, 0, 1 /) ! Fp
          State_old(km,10:18) = (/0, 0, 0, 0, 0, 0, 0, 0, 0 /)  ! Dev part of Mandel
          State_old(km,21:29) = (/0,0,0,0,0,0,0,0,0/)




	call SPECTRAL(matmul(transpose(F_tau), F_tau), eigvals, eigvecs)

           
	 E_tau = matmul(matmul(eigvecs, 
     + 0.5d0*reshape((/ dlog(eigvals(1)),z,z,z,dlog(eigvals(2)),z,z,z, 
     + dlog(eigvals(3)) /), (/3, 3/))), transpose(eigvecs))

          tr_Ee = E_tau(1,1)+E_tau(2,2)+E_tau(3,3) 
          Ee0 = E_tau - (tr_Ee/3.d0) * I_1
	
	    ! Griffin suggestion of midpoint
	    do i = 1, 3
        do j = 1, 3
          Finc(i,j) = F_tau(i,j) - F_T(i,j)
          F_avg(i,j) = 0.5d0 * (F_tau(i,j) + F_T(i,j))
        end do
      end do
      !was call matinv(F_tau, F_tau_inv, dtmp)
          call matinv(F_avg, F_avg_inv, dtmp)
    	do i = 1, 3
         do j = 1, 3
            Linc(i,j) = 0.0d0
            do k = 1, 3
               !Linc(i,j) = Linc(i,j) + Finc(i,k) * F_tau_inv(k,j)
               Linc(i,j) = Linc(i,j) + Finc(i,k) * F_avg_inv(k,j)
            end do
         end do
      end do
      
    	Dinc=0.5d0*(Linc+transpose(Linc))
    	tr_Dinc=(Dinc(1,1))+(Dinc(2,2))+(Dinc(3,3))
    	D_dev_inc=Dinc-(tr_Dinc)/3.d0 *I_1

      Tvis_old = zero_m
      Delta_Tvis = 2.d0*ma*D_dev_inc + ka*tr_Dinc*I_1
      Tvis_tau = Tvis_old + Delta_Tvis

      T_tau = 2.d0 * mu * Ee0 + kappa * tr_Ee * I_1+Tvis_tau

      state_new(km,21:29) = reshape(Tvis_tau, (/9/))

      else

          !
          ! Get state vars from last step
          !

          ! Get viscous stress from last step
          Tvis_old = reshape(state_old(km,21:29), (/3,3/))   	   					          
          
      ! Compute midpoint
      do i = 1, 3
        do j = 1, 3
          Finc(i,j) = F_tau(i,j) - F_T(i,j)
        end do
      end do
      
      do i = 1, 3
        do j = 1, 3
          F_avg(i,j) = 0.5d0 * (F_tau(i,j) + F_T(i,j))
        end do
      end do
      
      call matinv(F_avg, F_avg_inv, dtmp)
      
    	do i = 1, 3
         do j = 1, 3
            Linc(i,j) = 0.d0
            do k = 1, 3
               Linc(i,j) = Linc(i,j) + Finc(i,k) * F_avg_inv(k,j)
            end do
         end do
      end do
      
    	Dinc=0.5d0*(Linc+transpose(Linc))
    	tr_Dinc=(Dinc(1,1))+(Dinc(2,2))+(Dinc(3,3))
    	D_dev_inc=Dinc-(tr_Dinc)/3.d0 *I_1
    	
    	Delta_Tvis=2.d0*ma*D_dev_inc+ka*tr_Dinc*I_1
    	Tvis_tau=Tvis_old+Delta_Tvis

      ! No plastic deformation for viscoelasticity
      Fp_tau = I_1
      nu_p = 0.d0
      state_new(km,19)=0.d0
      state_new(km,20)=0.d0
		
      ! Since Fp = I, we have Fe = F
      Fe = F_tau
 	
       call SPECTRAL(matmul(transpose(Fe), Fe), eigvals, eigvecs)
           
	Ee = matmul(matmul(eigvecs, 
     + 0.5d0*reshape((/dlog(eigvals(1)),z,z,z,dlog(eigvals(2)),z,z,z, 
     + dlog(eigvals(3)) /), (/3, 3/))), transpose(eigvecs))
	   
	   tr_Ee = Ee(1,1)+Ee(2,2)+Ee(3,3)
	   Ee0 = Ee - (tr_Ee/3.d0) * I_1

	det_Fe=dexp(tr_Ee)

	! Neo-Hookean model
	Te_sph = kappa * tr_Ee
	Te_0 = 2.d0 * mu * Ee0

	Ue_inv = matmul(matmul(eigvecs, 
     + reshape((/ eigvals(1)**(-.5d0),z,z,z,eigvals(2)**(-.5d0),z,z,z, 
     + eigvals(3)**(-.5d0) /), (/3, 3/))), transpose(eigvecs))
          
       Re_tau = matmul(F_tau, Ue_inv)

	! Total Cauchy stress = Hyperelastic + Viscous
	T_tau = matmul(Re_tau, matmul(Te_0+I_1*Te_sph, 
     +   transpose(Re_tau)))/det_Fe + Tvis_tau

	state_new(km,1:9) = reshape(I_1, (/9/))         
	state_new(km,10:18) =  reshape(Te_0, (/9/))    
        state_new(km,21:29) = reshape(Tvis_tau, (/9/)) 

      end if ! if(total_time.eq.0 .and. step_time.eq.0)
	
!   Do not alter anything below this line.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



         ! Rotate Cauchy stress to ABAQUS stress
         !
        call MATINV(U_tau, U_inv, det_U)
        R_tau = matmul(F_tau, U_inv)
        T_tau = matmul(transpose(R_tau),matmul(T_tau,R_tau))
	 
	   
	  
        
	  ! Update ABAQUS stresses 
        !
        do i = 1,NDIR
            STRESS_NEW(KM,i) = T_tau(i,i)
        end do
        if (NSHR .ne. 0) then
          STRESS_NEW(KM,NDIR+1) = T_tau(1,2)
          if (NSHR .ne. 1) then
            STRESS_NEW(KM,NDIR+2) = T_tau(2,3)
            if (NSHR .ne. 2) then
              STRESS_NEW(KM,NDIR+3) = T_tau(1,3)
            endif
          endif
        endif


        !
        ! Specific internal energy
        !
        stress_power = ONE_HALF * (
     +          ( stress_old(km,1)+stress_new(km,1) )*strain_inc(km,1) +
     +          ( stress_old(km,2)+stress_new(km,2) )*strain_inc(km,2) +
     +          ( stress_old(km,3)+stress_new(km,3) )*strain_inc(km,3))
          if(nshr .eq. 1) then
             stress_power = stress_power + one_half*(
     +       TWO*(stress_old(km,4)+stress_new(km,4))*strain_inc(km,4) )
          else
             stress_power = stress_power + one_half*(
     +       TWO*(stress_old(km,4)+stress_new(km,4))*strain_inc(km,4) +
     +       TWO*(stress_old(km,5)+stress_new(km,5))*strain_inc(km,5) +
     +       TWO*(stress_old(km,6)+stress_new(km,6))*strain_inc(km,6) )
          endif
          ener_intern_new(km) = ener_intern_old(km)
     +                          + stress_power/density(km)

      end do  ! loop over km

      return
      end

C*****************************************************************
      SUBROUTINE MATINV(A,A_INV,DET_A_INV)
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION A(3,3), A_INV(3,3)

      PARAMETER(ZERO=0., ONE=1.)

      DET_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))

      IF (DET_A .LE. ZERO) THEN
        ! write(10,*) 'WARNING: SUBROUTINE MATINV:' 
        ! write(10,*) 'WARNING: DET of MAT is zero/negative!!'
      ENDIF

      DET_A_INV = ONE/DET_A

      A_INV(1,1) = DET_A_INV*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_INV(1,2) = DET_A_INV*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_INV(1,3) = DET_A_INV*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_INV(2,1) = DET_A_INV*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_INV(2,2) = DET_A_INV*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_INV(2,3) = DET_A_INV*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_INV(3,1) = DET_A_INV*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_INV(3,2) = DET_A_INV*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_INV(3,3) = DET_A_INV*(A(1,1)*A(2,2)-A(2,1)*A(1,2))

      RETURN
      END
	
C****************************************************************************
C**********************************************************************
C	THE FOLLOWING SUBROUTINES CALCULATE THE SPECTRAL
C	DECOMPOSITION OF A SYMMETRIC THREE BY THREE MATRIX
C**********************************************************************
	SUBROUTINE SPECTRAL(A,D,V)
C
C	THIS SUBROUTINE CALCULATES THE EIGENVALUES AND EIGENVECTORS OF
C	A SYMMETRIC 3 BY 3 MATRIX [A]. 
C
C	THE OUTPUT CONSISTS OF A VECTOR D CONTAINING THE THREE
C	EIGENVALUES IN ASCENDING ORDER, AND
C	A MATRIX [V] WHOSE COLUMNS CONTAIN THE CORRESPONDING
C	EIGENVECTORS.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER(NP=3)
	DIMENSION D(NP),V(NP,NP)
	DIMENSION A(3,3),E(NP,NP)

	DO 2 I = 1,3
          DO 1 J= 1,3
            E(I,J) = A(I,J)
1	  CONTINUE
2	CONTINUE

	CALL JACOBI(E,3,NP,D,V,NROT)
	CALL EIGSRT(D,V,3,NP)

	RETURN
	END

C**********************************************************************
	SUBROUTINE JACOBI(A,N,NP,D,V,NROT)

C	COMPUTES ALL EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C	MATRIX [A], WHICH IS OF SIZE N BY N, STORED IN A PHYSICAL 
C	NP BY BP ARRAY. ON OUTPUT, ELEMENTS OF [A] ABOVE THE DIAGONAL 
C	ARE DESTROYED, BUT THE DIAGONAL AND SUB-DIAGONAL ARE UNCHANGED
C	AND GIVE FULL INFORMATION ABOUT THE ORIGINAL SYMMETRIC MATRIX.
C	VECTOR D RETURNS THE EIGENVALUES OF [A] IN ITS FIRST N ELEMENTS.
C	[V] IS A MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS AS
C	[A] WHOSE COLUMNS CONTAIN, ON OUTPUT, THE NORMALIZED
C	EIGENVECTORSOF [A]. NROT RETURNS THE NUMBER OF JACOBI ROTATIONS
C	WHICH WERE REQUIRED.

C	THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", PAGE 346.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER (NMAX =100)
	DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)

C	INITIALIZE [V] TO THE IDENTITY MATRIX

	DO 12 IP = 1,N	
	  DO 11 IQ = 1,N
	    V(IP,IQ) = 0.D0
11        CONTINUE
          V(IP,IP) = 1.D0
12	CONTINUE

C	INITIALIZE [B] AND [D] TO THE DIAGONAL OF [A], AND Z TO ZERO.
C	THE VECTOR Z WILL ACCUMULATE TERMS OF THE FORM T*A_PQ AS
C	IN EQUATION (11.1.14)

	DO 13 IP = 1,N
	  B(IP) = A(IP,IP)
	  D(IP) = B(IP)
	  Z(IP) = 0.D0
13	CONTINUE
C
	NROT = 0
	DO 24 I = 1,50

C	SUM OFF-DIAGONAL ELEMENTS

          SM = 0.D0
          DO 15 IP = 1, N-1
            DO 14 IQ = IP + 1, N
	      SM = SM + DABS ( A(IP,IQ ))
14          CONTINUE
15        CONTINUE

C	IF SUM = 0., THEN RETURN. THIS IS THE NORMAL RETURN
C	WHICH RELIES ON QUADRATIC CONVERGENCE TO MACHINE 
C	UNDERFLOW.

          IF ( SM .EQ. 0.D0) RETURN
C
C	  IF ( SM .LT. 1.0D-15) RETURN

C	IN THE FIRST THREE SWEEPS CARRY OUT THE PQ ROTATION ONLY IF
C	|A_PQ| > TRESH, WHERE TRESH IS SOME THRESHOLD VALUE, 
C	SEE EQUATION (11.1.25). THEREAFTER TRESH = 0.

          IF ( I .LT. 4) THEN
            TRESH = 0.2D0*SM/N**2
          ELSE
            TRESH = 0.D0
          ENDIF
C
          DO 22 IP = 1, N-1
            DO 21 IQ = IP+1,N
              G = 100.D0*DABS(A(IP,IQ))

C	AFTER FOUR SWEEPS, SKIP THE ROTATION IF THE
C	OFF-DIAGONAL ELEMENT IS SMALL.

	      IF ((I .GT. 4) .AND. (DABS(D(IP))+G .EQ. DABS(D(IP)))
     +            .AND. ( DABS(D(IQ))+G .EQ. DABS(D(IQ)))) THEN
                A(IP,IQ) = 0.D0
              ELSE IF ( DABS(A(IP,IQ)) .GT. TRESH) THEN
                H = D(IQ) - D(IP)
                IF (DABS(H)+G .EQ. DABS(H)) THEN

C	T = 1./(2.*THETA), EQUATION(11.1.10)

	          T =A(IP,IQ)/H
	        ELSE
	          THETA = 0.5D0*H/A(IP,IQ)
	          T =1.D0/(DABS(THETA)+DSQRT(1.D0+THETA**2))
	          IF (THETA .LT. 0.D0) T = -T
	        ENDIF
	        C = 1.D0/DSQRT(1.D0 + T**2)
	        S = T*C
	        TAU = S/(1.D0 + C)
	        H = T*A(IP,IQ)
	        Z(IP) = Z(IP) - H
	        Z(IQ) = Z(IQ) + H
	        D(IP) = D(IP) - H
	        D(IQ) = D(IQ) + H
	        A(IP,IQ) = 0.D0

C	CASE OF ROTATIONS 1 <= J < P
				
	        DO 16 J = 1, IP-1
	          G = A(J,IP)
	          H = A(J,IQ)
	          A(J,IP) = G - S*(H + G*TAU)
	          A(J,IQ) = H + S*(G - H*TAU)
16	        CONTINUE

C	CASE OF ROTATIONS P < J < Q

	        DO 17 J = IP+1, IQ-1
	          G = A(IP,J)
	          H = A(J,IQ)
	          A(IP,J) = G - S*(H + G*TAU)
	          A(J,IQ) = H + S*(G - H*TAU)
17	        CONTINUE

C	CASE OF ROTATIONS Q < J <= N

	        DO 18 J = IQ+1, N
                  G = A(IP,J)
	          H = A(IQ,J)
	          A(IP,J) = G - S*(H + G*TAU)
	          A(IQ,J) = H + S*(G - H*TAU)
18	        CONTINUE
	        DO 19 J = 1,N
	          G = V(J,IP)
	          H = V(J,IQ)
	          V(J,IP) = G - S*(H + G*TAU)
	          V(J,IQ) = H + S*(G - H*TAU)
19	        CONTINUE
	        NROT = NROT + 1
              ENDIF
21	    CONTINUE
22	  CONTINUE

C	UPDATE D WITH THE SUM OF T*A_PQ, AND REINITIALIZE Z

	  DO 23 IP = 1, N
	    B(IP) = B(IP) + Z(IP)
	    D(IP) = B(IP)
	    Z(IP) = 0.D0
23	  CONTINUE
24	CONTINUE

C	IF THE ALGORITHM HAS REACHED THIS STAGE, THEN
C	THERE ARE TOO MANY SWEEPS, PRINT A DIAGNOSTIC
C	AND STOP.

	! write (80,'(/1X,A/)') '50 ITERS IN JACOBI SHOULD NEVER HAPPEN'

	RETURN
	END

C**********************************************************************
	SUBROUTINE EIGSRT(D,V,N,NP)

C	GIVEN THE EIGENVALUES [D] AND EIGENVECTORS [V] AS OUTPUT FROM
C	JACOBI, THIS ROUTINE SORTS THE EIGENVALUES INTO ASCENDING ORDER, 
C	AND REARRANGES THE COLUMNS OF [V] ACCORDINGLY.

C	THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", P. 348.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION D(NP),V(NP,NP)

	DO 13 I = 1,N-1
	  K = I
	  P = D(I)
	  DO 11 J = I+1,N
	    IF (D(J) .GE. P) THEN
	      K = J
	      P = D(J)
	    END IF
11	  CONTINUE
	  IF (K .NE. I) THEN
	    D(K) = D(I)
	    D(I) = P
	    DO 12 J = 1,N
	      P = V(J,I)
	      V(J,I) = V(J,K)
	      V(J,K) = P
12	    CONTINUE
  	  ENDIF
13	CONTINUE

	RETURN
	END

	
