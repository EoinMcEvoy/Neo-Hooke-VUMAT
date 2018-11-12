      SUBROUTINE vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      INCLUDE 'vaba_param.inc'
C
c      DEFAULT
      DIMENSION props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
	 
c     ADDITIONAL
      DIMENSION dgrad9(ndir+nshr+nshr),dgrad(ndir,ndir),bmat(ndir,ndir),
     1 b2mat(ndir,ndir),xkirch(ndir,ndir),xkirchV(ndir+nshr)
C
      CHARACTER*80 cmname
C
C Neo-hooke hyperelasticity for 3D elements
C Props(1) = C10
C Props(2) = D1
C
C The state variables are stored as:
C      STATE(*,1) = invariant I1
C
      PARAMETER( zero = 0., one = 1., two = 2., three = 3.,
     1  four=4., five=5., six=6.)
C
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  
      c10 = props(1)
      d1 = props(2)
	  
C    START NBLOCK LOOP
      DO km = 1, nblock
	  
C    CONVERT DEFORMATION GRADIENT TO 3x3 MATRIX
      DO i = 1,ndir+2*nshr
         dgrad9(i) = defgradNew(km,i)
      END DO	  
      CALL k9vector2matrix(dgrad9, dgrad, nshr)
	  
C    CALCULATE DETERMINANT OF DEFORMATION GRADIENT

      xJ = dgrad(1,1) * dgrad(2,2) * dgrad(3,3) -
     1 dgrad(1,2) * dgrad(2,1) * dgrad(3,3)
      xJ = xJ + dgrad(1,2) * dgrad(2,3) * dgrad(3,1) +
     1 dgrad(1,3) * dgrad(3,2) * dgrad(2,1) -
     2 dgrad(1,3) * dgrad(3,1) * dgrad(2,2) -
     3 dgrad(2,3) * dgrad(3,2) * dgrad(1,1)  
		 
C      ZERO MATRICES
        DO I=1,3
           DO J=1,3
              bmat(I,J) = zero  	
           END DO
        END DO 
		
C      CALCULATE LEFT CAUCHY-GREEN DEFORMATION TENSOR       
        DO I=1,3
            DO J=1,3
                DO K = 1,3				
                    bmat(I,J) = bmat(I,J) + dgrad(I,K)*dgrad(J,K)				
                END DO
            END DO
        END DO
		
C-----------------------------------------------------------------------
C     CALCULCATE THE INVARIANTS                
C-----------------------------------------------------------------------
C        
        xI1 = bmat(1,1)+bmat(2,2)+bmat(3,3)
C        
C       kmtms (M, N, L, A, KA, B, KB, C, KC)
        CALL kmtms(3, 3, 3, bmat, 3, bmat, 3, b2mat, 3)
C        
        trb2 = b2mat(1,1)+b2mat(2,2)+b2mat(3,3)
C        
        xI2 = ((xI1**two) - trb2)/two
C        
        xI3=xJ**two
	 
C--------------------------------------------------------------------
C     CALCULATE KIRCH STRESS
C--------------------------------------------------------------------
C
C     PART 1
      coeff1 = two*c10/(xJ**(two/three))
	  
C     KIRCHHOFF STRESS PART 2
      trbmat= coeff1*(bmat(1,1)+bmat(2,2)+bmat(3,3))/three
	  
C     KIRCHHOFF STRESS PART 1   
      DO I=1,3
        DO J=1,3		
            xkirch(I,J) = bmat(I,J) * coeff1			
        END DO
      END DO
	  
C     SUBTRACT THE PART 2   
      DO I = 1,3	  
        xkirch(I,I) = xkirch(I,I) - trbmat		
      END DO
	  
C     FORM PART 3      
      coeff3 = two*(xJ-one)*xJ/d1
	  
C     ADD TO THE PREVIOUS PARTS      
      DO I = 1,3  
        xkirch(I,I) = xkirch(I,I) + coeff3	
      END DO
	     
C    EXPLICIT REQUIRES STRESS AS (11,22,33,12,23,31)
      CALL kmatrix2vector_explicit(xkirch, xkirchV, nshr)
	  
C    CONVERT TO CAUCHY STRESS
      DO I = 1,ndir+nshr  
         stressNew(km,I) = (one/xJ)*xkirchV(I)
      END DO
	  
      stateNew(km,1) = xI1
	  
      xener=c10*(xI1-three)+(one/d1)*((xJ-one)**two)
      enerInternNew(km) = xener/density(km)
		
C	  END NBLOCK LOOP
      END DO
	  

      RETURN
      CONTAINS
	  
	  
	  
	  
C------------------------------------------------------------------
C           SUBROUTINES 
C------------------------------------------------------------------
C
C * kmatrix2vector_explicit  -  Convert a 3x3 matrix to a 6x1 vector
C
C * k9vector2matrix  -  Convert a 9x1 vector to a 3x3 matrix
C 
C * kmtms  -  Multiply two 2nd order tensors
c-----------------------------------------------------------------------------------------------

      SUBROUTINE kmatrix2vector_explicit(xmat, vec, nshr)
      
      INCLUDE 'vaba_param.inc'
      
      INTENT(IN) :: xmat, nshr
      INTENT(OUT):: vec
	  
C    Explicit requires stress as (11,22,33,12,23,31)
      
      DIMENSION xmat(3,3), vec(6)
  
        DO i=1,3
            vec(i) = xmat(i,i)
        END DO
               
        vec(4) = xmat(1,2)
        
        IF (nshr==3) THEN
            vec(5) = xmat(2,3)
            vec(6) = xmat(3,1)
        END IF
      
      END SUBROUTINE kmatrix2vector_explicit	
	  
c-----------------------------------------------------------------------------------------------	  
c-----------------------------------------------------------------------------------------------

      SUBROUTINE k9vector2matrix(vec, xmat, nshr)
      
      INCLUDE 'vaba_param.inc'
      
      INTENT(IN) :: vec, nshr
      INTENT(OUT):: xmat
      
      DIMENSION xmat(3,3), vec(9)
  
        DO i=1,3
            xmat(i,i) = vec(i)
        END DO
               
        xmat(1,2) = vec(4)
		
        IF (nshr==1) THEN
           xmat(2,1) = vec(5)
        ELSE
           xmat(2,3) = vec(5)
           xmat(3,1) = vec(6)
           xmat(2,1) = vec(7)
           xmat(3,2) = vec(8)
           xmat(1,3) = vec(9)
        END IF
      
      END SUBROUTINE k9vector2matrix	
	  
c-----------------------------------------------------------------------------------------------	  
c-----------------------------------------------------------------------------------------------
	       
      SUBROUTINE kmtms (M, N, L, A, KA, B, KB, C, KC)
      
      INCLUDE 'vaba_param.inc'
C      
      INTENT(IN) :: M, N, L, A, KA, B, KB, KC
      INTENT(OUT):: C      
C      
C    PRODUCT OF REAL MATRICES
C
      DIMENSION A(KA,N), B(KB,L), C(KC,L)
      DOUBLE PRECISION W
C       
C
      DO 30 J = 1,L
         DO 20 I = 1,M
            W = 0.D0
            DO 10 K = 1,N
               W = W + A(I,K) * B(K,J)
   10       CONTINUE
            C(I,J) = W
   20    CONTINUE
   30 CONTINUE
   
      RETURN
	  
      END SUBROUTINE kmtms
	  
c-----------------------------------------------------------------------------------------------

 
      END