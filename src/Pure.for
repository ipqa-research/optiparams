      SUBROUTINE paramsfromdelta1(del1,Tc,Pc,OM,ac,b,rk,Dc)
      implicit double precision (A-H,O-Z)
      PARAMETER (RGAS=0.08314472d0)
C	Critical constants must be given in K and bar
C	b will be in L/mol and ac in bar*(L/mol)**2
      PARAMETER (A0=0.0017,B0=1.9681,C0=-2.7238)
      PARAMETER (A1=-2.4407,B1=7.4513,C1=12.504)
	COMMON /ABd1/ a,bb,d	

	d1=(1+del1**2)/(1+del1)
	y=1+(2*(1+del1))**(1.0d0/3)+(4/(1+del1))**(1.0d0/3)
	OMa=(3*y*y+3*y*d1+d1**2+d1-1.0d0)/(3*y+d1-1.0d0)**2
	OMb=1/(3*y+d1-1.0d0)
	Zc=y/(3*y+d1-1.0d0)
	RT=RGAS*Tc
	ac=OMa*RT**2/Pc
	b=OMb*RT/Pc
	Dc=Pc/(Zc*RT)
c	IF(nmodel.EQ.3)THEN
		rk=(A1*Zc+A0)*OM**2+(B1*Zc+B0)*OM+(C1*Zc+C0) ! initial guess for k parameter
		Tr=0.7D0	! Change here to use another Pv than the one at Tr 0.7
		Pvdat=Pc*10**-(1.0D0+OM)
c            if(IVAP==1)then ! added 29/06/2013 in order to allow for better reproductions of Pv curves
c			    READ(NIN,*)T,Pvdat
c               Tr=T/Tc
c            end if
		a=ac*(3/(2+Tr))**rk
	 	bb=b
	 	d=del1
	 	T=Tr*Tc
	 	CALL VaporPressure(T,Pvdat,Pv,Dc,RHOL,RHOV,phiL)
		if(Pv>Pvdat)then
			dk = 0.1
		else
			dk = -0.1
		end if
		err = 1.0
		do while (err > 0.005)
			Pold = Pv
			oldk = rk
			rk = rk + dk
  			a=ac*(3/(2+Tr))**rk
			CALL VaporPressure(T,Pvdat,Pv,Dc,RHOL,RHOV,phiL)
			dk = -(Pv-Pvdat)*(rk-oldk)/(Pv-Pold)
			err = abs(Pv-Pvdat)/Pvdat
		end do
		rk = rk + dk
c		ELSE
c			Zc=y/(3*y+d1-1.0d0)
c			Vceos=Zc*RGAS*Tc/Pc
c			IF(nmodel.eq.1)THEN
c				rm=0.48+1.574*OM-0.175*OM**2  ! m from SRK
c			ELSE
c				rm=0.37464+1.54226*OM-0.26992*OM**2  ! m from PR
c			END IF
c		END IF

      RETURN
      END
c
	SUBROUTINE VaporPressure(T,PVini,Pv,Dc,RHOL,RHOV,phiL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ERRMAX=1.D-8)
	dphi = 0.0D0
	P = PVini
	n=1
 30	call VCALCp(1,T,P,V)
	RHOL = 1/V
	call VCALCp(-1,T,P,V) ! SOLVE for vapor density
	RHOV = 1/V
	if(RHOL.LT.0.9*dc) then
		P=1.01*P
		go to 30
	else if(RHOV.GT.dc) then
		P=0.99*P
		go to 30
	end if
	call FUG_CALC(T,P,1/RHOL,phi)
	phiL = phi
	call FUG_CALC(T,P,V,phi)
	phiV = phi
	dphiold = dphi
	dphi = phiV - phiL
	IF (ABS(dphi).gt.ERRMAX) THEN
		Pold = Plast
		Plast = P
	if(dphiold.eq.0.0D0) then  ! .or.Tr.gt.0.975
		P = P * (phiL/phiV)
	else
		P = Plast - dphi*(Plast-Pold)/(dphi-dphiold)
	end if
c		n=n+1
		GO TO 30
	END IF
	PV = P
	return
	END

      SUBROUTINE VCALCp(ITYP,T,P,V)
C
C     ROUTINE FOR CALCULATION OF VOLUME, GIVEN PRESSURE
C
C     INPUT:
C
C     ITYP:        TYPE OF ROOT DESIRED
C     T:           TEMPERATURE
C     P:           PRESSURE
C
C     OUTPUT:
C
C     V:           VOLUME
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (RGAS=0.08314472d0)
      LOGICAL FIRST_RUN
	COMMON /ABd1/ a,b,d1	
      FIRST_RUN = .TRUE.
	VCP = b
      S3R = 1.D0/VCP
      ITER = 0
C
      ZETMIN = 0.D0
      ZETMAX = .99D0
      IF (ITYP .GT. 0) THEN
         ZETA = .5D0
         ELSE
C..............IDEAL GAS ESTIMATE
         ZETA = MIN (.5D0,VCP*P/(RGAS*T))
      ENDIF
  100 CONTINUE
C	WRITE(*,*)'ZETA',ZETA
      V = VCP/ZETA
      ITER = ITER + 1
	CALL vdWg_Derivs(1,T,V,F,F_V,F_2V,F_N)
      PCALC = RGAS*T*(1/V - F_V)
C	WRITE(*,*)'PCALC',PCALC
      IF (PCALC .GT. P) THEN
         ZETMAX = ZETA
         ELSE
         ZETMIN = ZETA
      ENDIF
c	write(*,*)'VCALC V=',V
      AT  = F - LOG(V) + V*P/(T*RGAS)
      DER = RGAS*T*(F_2V+1.D0)*S3R
      DEL = -(PCALC-P)/DER
      ZETA = ZETA + MAX (MIN(DEL,0.1D0),-.1D0)
      IF (ZETA .GT. ZETMAX .OR. ZETA .LT. ZETMIN)
     &      ZETA = .5D0*(ZETMAX+ZETMIN)
      IF (ABS(DEL) .GT. 1D-10) GOTO 100
      IF (ITYP .EQ. 0 ) THEN
C
C FIRST RUN WAS VAPOUR; RERUN FOR LIQUID
C
         IF (FIRST_RUN) THEN
            VVAP = V
            AVAP = AT
            FIRST_RUN = .FALSE.
            ZETA = 0.5D0
            GOTO 100
            ELSE
            IF (AT .GT. AVAP) V = VVAP
          ENDIF
      ENDIF
	return
      END
C
	SUBROUTINE FUG_CALC(T,P,V,phi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (RGAS=0.08314472d0)
	RT = RGAS*T
	Z = P*V/RT
	CALL vdWg_Derivs(2,T,V,F,F_V,F_2V,F_N)
      PHI=EXP(F_N)/Z
	return
	END
C
C
	subroutine vdWg_Derivs(NDER,T,V,F,F_V,F_2V,F_N)
c	
C     THE SUBROUTINE CALCULATES THE CONTRIBUTION TO THE RESIDUAL,
C     REDUCED HELMHOLZ ENERGY (F) AND
C     ITS FIRST AND SECOND DERIVATIVE WRT V
C
C     INPUT:
C	NDER:		 indicates which derivatives are required.
C				 1 is used for density calculation and 2 for fugacity
c		NDER = 1: CALCULATES F, F_V AND F_2V      
c 		NDER = 2: CALCULATES F AND F_N 
C     T:           TEMPERATURE
C     V:           VOLUME (ML/MOL) or (ML) for checking n-derivatives
C
C     OUTPUT:	   NDER
C     F:			5	A^RES/RT CONTRIBUTION (DIMENSIONLESS) or (MOLES)
C     F_V:		5	1ST V-DERIVATIVE OF F
C     F_2V:			1ST V-DERIVATIVE OF F_V  (*V**2)
C     F_N:			1ST N-DERIVATIVE OF F
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (RGAS=0.08314472d0)
	COMMON /ABd1/ a,b,d
	C = (1-d)/(1+d)
	aRT = a / (RGAS*T)
      ETA = 0.25 * b / V
	SUMC = c*b+V
	SUMD = d*b+V
      REP = -log(1-4*ETA)
	ATT = aRT*LOG(SUMD/SUMC)/(b*(C-D))
	ATTV = aRT/SUMC/SUMD
      REPV = 1/(1-4*ETA)-1
      REP2V = 1/(1-4*ETA)**2-1
      ATT2V = aRT*V**2*(1/SUMD**2-1/SUMC**2)/(b*(C-D))
	F = REP+ATT
	F_V = (-REPV/V+ATTV)
	IF (NDER.EQ.1) THEN
	F_2V = REP2V-ATT2V
	ELSE
		F_N = REP + ATT - V*F_V
	END IF
	return
	end