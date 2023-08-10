	subroutine readcomp(nin,nout)
      implicit DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2,RGAS=0.08314472d0)
C	Critical constants must be given in K and bar
C	b will be in L/mol and ac in bar*(L/mol)**2
	DOUBLE PRECISION Kij(nco,nco),lij(nco,nco),mij(nco,nco)
	DOUBLE PRECISION Kinf,Kinf1,Kinf2,K01,K02,lijk(nco,nco,nco)
c     PARAMETER (A0=0.0017,B0=1.9681,C0=-2.7238)
c     PARAMETER (A1=-2.4407,B1=7.4513,C1=12.504)
	dimension ac(nco),b(nco),del1(nco),rk(nco),diam(nco),vc(nco)
	dimension D(6),OM(nco),Vceos(nco)
      CHARACTER*10 fluid(nco)
	COMMON/names/fluid
      COMMON/CRIT/TC(nco),PC(nco),DCeos(nco)
	COMMON /COMPONENTS/ ac,b,del1,rk,Kij,NTDEP
	COMMON /K12/ K12
	COMMON/COVOL/bb1(2)
	COMMON /rule/ncomb,iRuleDel
	COMMON /SpecRep/arepfr,arep,dardT,dardT2        ! August 2016
	COMMON /bcross/bij(nco,nco)
	COMMON /D1cross/D1ij(nco,nco)
	COMMON /Tdep/ Kinf,Tstar
	COMMON /Kcubic/Kinf1,Kinf2,K01,K02,Tstar1,Tstar2,C1,C2
	COMMON /bcrosscub/bijk(nco,nco,nco)
c	D=[0.428363,18.496215,0.338426,0.660,789.723105,2.512392]
	NC=2
	read(NIN,*) ncomb,NTDEP   !,iRuleDel
	iRuleDel = 1
	do i=1,nc
	READ(NIN,'(A)')fluid(i)
	READ(NIN,*)Tc(i),Pc(i),OM(i),Vceos(i),Zrat
	RT=RGAS*Tc(i)
	Zc=Pc(i)*Vceos(i)/RT
	Zcin=Zc/Zrat
	Vc(i)=Vceos(i)/Zrat
	dceos(i)=1/Vceos(i)
	READ(NIN,*)ac(i),b(i),del1(i),rk(i)
 4	bb1(i)=b(i)
	write(nout,'(A)')fluid(i)
	write(nout,1)Tc(i),Pc(i),Vc(i),OM(i)
	write(nout,3)Zcin,Zrat,Zc,Vceos(i)
	write(nout,2)ac(i),b(i),del1(i),rk(i)
	Kij(i,i)=0.0D0
	Lij(i,i)=0.0D0
	mij(i,i)=0.0D0
	IF(i.gt.1)then
		if(ncomb.lt.2)then
			READ(NIN,*) (Kij(j,i),j=1,i-1)
			if(NTDEP.ge.1)READ(NIN,*) Kinf
			if(NTDEP.eq.1)READ(NIN,*)Tstar
			READ(NIN,*) (lij(j,i),j=1,i-1)
		    if(ncomb == -1)then ! Added August 2016 for CO2 new approach
		        READ(NIN,*) arepfr
		    end if
		else
			READ(NIN,*) K01,K02
			if(NTDEP.ge.1)READ(NIN,*) Kinf1,Kinf2
			if(NTDEP.eq.1)READ(NIN,*)Tstar1,Tstar2
			if(NTDEP.eq.2)READ(NIN,*)C1,C2
			READ(NIN,*) Lijk(1,1,2),Lijk(1,2,2)
		end if
	if(iRuleDel==2) READ(NIN,*) (mij(j,i),j=1,i-1)  ! added 24/08/2017
	ENDIF
	end do
c	if(NTDEP.eq.1)Tstar=258.0
c	if(NTDEP.eq.2)READ(NIN,*)Tstar
	write(NOUT,*)
	write(nout,*)'Tc, Pc and Vc are given in K, bar and L/mol respectively'
c	Lij(1,2)=0.10d0		!  delete or comment!!!!!!!!!!!!!!!!!
c	Kij(1,2)=0.10d0		!  delete or comment!!!!!!!!!!!!!!!!!
	K12=Kij(1,2)
 1	FORMAT('Tc=',F9.4,'   Pc =',F9.4,'   Vc =',F8.4,'   OM =',F7.4)
 3	FORMAT('Zc=',F9.4,' Zcrat=',F9.4,' Zceos=',F8.4,' Vceos=',F7.4)
 2	FORMAT('ac=',F9.4,'    b =',F9.4,'  del1=',F8.4,'    k =',F7.4)
	write(NOUT,*)
	if(ncomb.lt.2)then
		if(NTDEP.EQ.0)then
			write(NOUT,*)' K12 = ',Kij(1,2)
			write(NOUT,*)
		else
			write(NOUT,*)' K012 = ',Kij(1,2)
			write(NOUT,*)
			write(NOUT,*)' Kinf = ',Kinf
			write(NOUT,*)
			write(NOUT,*)'Tstar = ',Tstar
			write(NOUT,*)
		end if
c		write(NOUT,*)'  KIJ MATRIX'
c		DO I=1,NC
c		write(NOUT,6)FLUID(I),(Kij(j,i),j=1,i-1)
c		END DO
c		write(NOUT,*)
		write(NOUT,*)'  LIJ MATRIX'
		DO I=1,NC
		write(NOUT,6)FLUID(I),(Lij(j,i),j=1,i-1)
		END DO
		if(iRuleDel==2)then
		    write(NOUT,*)'  MIJ MATRIX'
		    DO I=1,NC
		    write(NOUT,6)FLUID(I),(mij(j,i),j=1,i-1)
		    END DO
		end if
		if(ncomb == -1)then ! Added August 2016 for CO2 new approach
		    write(NOUT,*) 'Repulsive fraction of ac in CO2: ',arepfr
		end if
	else
		if(NTDEP.EQ.0)then
			write(NOUT,*)' Kijk:     112      122'
			write(NOUT,7)K01,K02
			write(NOUT,*)
		else
			write(NOUT,*)' K0ijk:    112      122'
			write(NOUT,7)K01,K02
			write(NOUT,*)
			write(NOUT,*)'Kinfijk:   112      122'
			write(NOUT,7)Kinf1,Kinf2
			write(NOUT,*)
			write(NOUT,*)'Tstar  :   112      122'
			write(NOUT,8)Tstar1,Tstar2
			write(NOUT,*)
		end if
		if(NTDEP.EQ.2)then
			write(NOUT,*)' Cijk:     112      122'
			write(NOUT,7)C1,C2
			write(NOUT,*)
		end if
		write(NOUT,*)' Lijk:     112      122'
		write(NOUT,7)Lijk(1,1,2),Lijk(1,2,2)
		write(NOUT,*)
	end if
	write(NOUT,*)
	write(NOUT,*)' Combining rules:'
	if(ncomb<=0)then
	write(NOUT,*)' 0: Classical or van der Waals '
	  if(ncomb==-1)write(NOUT,*)' -1: Specific Repulsion for CO2 '
		do i=1,nc
		do j=i,nc
		bij(i,j)=(1-lij(i,j))*(b(i)+b(j))/2
		bij(j,i)=bij(i,j)
		end do
		end do
	else if(ncomb.eq.3)then
		do i=1,nc
		bijk(i,i,i)=b(i)
		do j=i+1,nc
			bijk(i,i,j)=(1-lijk(i,i,j))*(2*b(i)+b(j))/3
			bijk(i,j,i)=bijk(i,i,j)
			bijk(j,i,i)=bijk(i,i,j)
			bijk(i,j,j)=(1-lijk(i,j,j))*(b(i)+2*b(j))/3
			bijk(j,i,j)=bijk(i,j,j)
			bijk(j,j,i)=bijk(i,j,j)
			do k=j+1,nc	! only possible with three or more components
				bijk(i,j,k)=(1-lijk(i,j,k))*(b(i)+b(j)+b(k))/3
				bijk(j,i,k)=bijk(i,j,k)
				bijk(i,k,j)=bijk(i,j,k)
				bijk(j,k,i)=bijk(i,j,k)
				bijk(k,i,j)=bijk(i,j,k)
				bijk(k,j,i)=bijk(i,j,k)
			end do
		end do
		end do
	else
	write(NOUT,*)' 1: Lorentz-Berthelot'
		third=1.0d0/3
		do i=1,nc
		diam(i)=b(i)**third
		end do
		do i=1,nc
		do j=i,nc
		bij(i,j)=((1-lij(i,j))*(diam(i)+diam(j))/2)**3
		bij(j,i)=bij(i,j)
		end do
		end do
	end if
	if(iRuleDel==2)then
		do i=1,nc   ! added 24/08/2017
		do j=i,nc
		D1ij(i,j)=(1-mij(i,j))*SQRT(del1(i)*del1(j))
		D1ij(j,i)=D1ij(i,j)
		end do
		end do
      end if
 6	FORMAT(A10,4F8.5)
 7	FORMAT(9x,F7.4,2x,F7.4)
 8	FORMAT(9x,F7.2,2x,F7.2)
	end
c
c	Then Kij values will be called indicating the lower index first, e.g. Kij(1,3)
c
c
	subroutine aTder(ac,Tc,rk,T,a,dadT,dadT2)
      implicit DOUBLE PRECISION (A-H,O-Z)
c	Given ac,Tc and the k parameter of the RKPR correlation, as well as the actual T,
c	this subroutine calculates a(T) and its first and second derivatives with T.
	COMMON /MODEL/ NMODEL
	Tr=T/Tc
	IF(NMODEL.LE.2)THEN
		rm=rk
		a=ac*(1+rm*(1-sqrt(Tr)))**2
		dadT=ac*rm*(rm-(rm+1)/sqrt(Tr))/Tc
		dadT2=ac*rm*(rm+1)/(2*Tc**2*Tr**1.5D0)
	ELSE
		a=ac*(3/(2+Tr))**rk
		dadT=-rk*a/Tc/(2+Tr)
		dadT2=-(rk+1)*dadT/Tc/(2+Tr)
	END IF
	end
C
	subroutine aijTder(NTD,nc,T,aij,daijdT,daijdT2)
      implicit DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2)
	DOUBLE PRECISION Kinf,Kij0(nco,nco),Kij(nco,nco)
	dimension ai(nco),daidT(nco),daidT2(nco)
	dimension aij(nco,nco),daijdT(nco,nco),daijdT2(nco,nco)
      COMMON/CRIT/TC(nco),PC(nco),DC(nco)
	COMMON /COMPONENTS/ ac(nco),b(nco),d1(nco),rk(nco),Kij0,NTDEP
	COMMON /bcross/bij(nco,nco)
	COMMON /rule/ncomb,iRuleDel
	COMMON /SpecRep/arepfr,arep,dardT,dardT2        ! August 2016
	COMMON /Tdep/ Kinf,Tstar
	IF(NTDEP.GE.1)THEN
		Kij=0.0D0
		Kij(1,2)=Kinf+Kij0(1,2)*exp(-T/Tstar)
		Kij(2,1)=Kij(1,2)
C	ELSE IF(NTDEP.EQ.2)THEN
C		Kij=0.0D0
C		Kij(1,2)=Kij0(1,2)*exp(-T/Tstar)
C		Kij(2,1)=Kij(1,2)
	ELSE
		Kij=Kij0
	END IF
	DO i=1,nc
	call aTder(ac(i),Tc(i),rk(i),T,ai(i),daidT(i),daidT2(i))
        if(i==1.and.ncomb==-1)then !  August 2016 for CO2 new approach
            arep = arepfr*ai(1)
            ai(i)=ai(i)+arep
            dardT = arepfr*daidT(1)
            daidT(i)=daidT(i)+dardT
            dardT2 = arepfr*daidT2(1)
            daidT2(i)=daidT2(i)+dardT2
        end if
	aij(i,i)=ai(i)
	daijdT(i,i)=daidT(i)
	daijdT2(i,i)=daidT2(i)
	IF (i.gt.1) THEN
	do j=1,i-1
	aij(j,i)=sqrt(ai(i)*ai(j))*(1-Kij(j,i))
	aij(i,j)=aij(j,i)
	if(NTD.EQ.1)then
		daijdT(j,i)=(1-Kij(j,i))*(sqrt(ai(i)/ai(j))*daidT(j)+
     &						  sqrt(ai(j)/ai(i))*daidT(i))/2
		daijdT(i,j)=daijdT(j,i)
		daijdT2(j,i)=(1-Kij(j,i))*(daidT(j)*daidT(i)/sqrt(ai(i)*ai(j))
     &		+sqrt(ai(i)/ai(j))*(daidT2(j)-daidT(j)**2/(2*ai(j)))
     &		+sqrt(ai(j)/ai(i))*(daidT2(i)-daidT(i)**2/(2*ai(i))))/2
		daijdT2(i,j)=daijdT2(j,i)
	end if
	end do
	END IF
	END DO
	if (ncomb.eq.1) then
		DO i=1,nc-1
		DO j=i+1,nc
			barrgij=bij(i,j)/sqrt(b(i)*b(j))
			aij(i,j)=barrgij*aij(i,j)
			aij(j,i)=aij(i,j)
			daijdT(i,j)=barrgij*daijdT(i,j)
			daijdT(j,i)=daijdT(i,j)
			daijdT2(i,j)=barrgij*daijdT2(i,j)
			daijdT2(j,i)=daijdT2(i,j)
		END DO
		END DO
	end if
	IF(NTDEP.ge.1.and.NTD.EQ.1)THEN
		aux=daijdT(1,2)
		ratK=Kij(1,2)/(1-Kij(1,2))/Tstar
		daijdT(1,2)=aux+aij(1,2)*ratK
		daijdT(2,1)=daijdT(1,2)
		daijdT2(1,2)=daijdT2(1,2)+(2*aux-aij(1,2)/Tstar)*ratK  ! 2* was missing (before aux)
c												since implementation in 2005	(06-03-2008)
		daijdT2(2,1)=daijdT2(1,2)
	END IF
	end
C
	subroutine aijkTder(NTD,nc,T,aijk,daijkdT,daijkdT2)
      implicit DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2)
	DOUBLE PRECISION Kinf1,Kinf2,K01,K02,Kijk(nco,nco,nco),Kij(nco,nco)
	dimension ai(nco),daidT(nco),daidT2(nco)
	dimension aijk(nco,nco,nco),daijkdT(nco,nco,nco),daijkdT2(nco,nco,nco)
      COMMON/CRIT/TC(nco),PC(nco),DC(nco)
	COMMON /COMPONENTS/ ac(nco),b(nco),d1(nco),rk(nco),Kij,NTDEP
c	COMMON /Tdep/ Tstar
	COMMON /Kcubic/Kinf1,Kinf2,K01,K02,Tstar1,Tstar2,C1,C2
	COMMON /bcross/bij(nco,nco)
	Kijk=0.0D0
	IF(NTDEP.EQ.1)THEN
		Kijk(1,1,2)=Kinf1+K01*exp(-T/Tstar1)
		Kijk(1,2,1)=Kijk(1,1,2)
		Kijk(2,1,1)=Kijk(1,1,2)
		Kijk(2,2,1)=Kinf2+K02*exp(-T/Tstar2)
		Kijk(2,1,2)=Kijk(2,2,1)
		Kijk(1,2,2)=Kijk(2,2,1)
	ELSE IF(NTDEP.EQ.2)THEN
		A1=K01
		B1=Kinf1
		SUM1=(C1+T)
		rT1=T/SUM1
		Kijk(1,1,2)=A1+B1*rT1
		Kijk(1,2,1)=Kijk(1,1,2)
		Kijk(2,1,1)=Kijk(1,1,2)
		A2=K02
		B2=Kinf2
		SUM2=(C2+T)
		rT2=T/SUM2
		Kijk(2,2,1)=A2+B2*rT2
		Kijk(2,1,2)=Kijk(2,2,1)
		Kijk(1,2,2)=Kijk(2,2,1)
	ELSE
		Kijk(1,1,2)=K01
		Kijk(1,2,1)=K01
		Kijk(2,1,1)=K01
		Kijk(2,2,1)=K02
		Kijk(2,1,2)=K02
		Kijk(1,2,2)=K02
	END IF
	third=1.0d0/3
	DO i=1,nc
	call aTder(ac(i),Tc(i),rk(i),T,ai(i),daidT(i),daidT2(i))
	aijk(i,i,i)=ai(i)
	daijkdT(i,i,i)=daidT(i)
	daijkdT2(i,i,i)=daidT2(i)
	IF (i.gt.1) THEN
		do j=1,i-1
			aijk(j,j,i)=(1-Kijk(j,j,i))*(ai(i)*ai(j)*ai(j))**third
			aijk(j,i,j)=aijk(j,j,i)
			aijk(i,j,j)=aijk(j,j,i)
			aijk(j,i,i)=(1-Kijk(j,i,i))*(ai(i)*ai(i)*ai(j))**third
			aijk(i,i,j)=aijk(j,i,i)
			aijk(i,j,i)=aijk(j,i,i)
		  if(NTD.EQ.1)then
			daijkdT(j,j,i)=(2*daidT(j)/ai(j)+daidT(i)/ai(i))*aijk(j,j,i)/3
			daijkdT(j,i,j)=daijkdT(j,j,i)
			daijkdT(i,j,j)=daijkdT(j,j,i)
			daijkdT(j,i,i)=(daidT(j)/ai(j)+2*daidT(i)/ai(i))*aijk(j,i,i)/3
			daijkdT(i,i,j)=daijkdT(j,i,i)
			daijkdT(i,j,i)=daijkdT(j,i,i)
			daijkdT2(j,j,i)=daijkdT(j,j,i)**2 / aijk(j,j,i) +
	1					(2*daidT2(j)/ai(j)+daidT2(i)/ai(i)-
     1			(2*(daidT(j)/ai(j))**2+(daidT(i)/ai(i))**2))*aijk(j,j,i)/3
			daijkdT2(j,i,j)=daijkdT2(j,j,i)
			daijkdT2(i,j,j)=daijkdT2(j,j,i)
			daijkdT2(j,i,i)=daijkdT(j,i,i)**2 / aijk(j,i,i) +
	1					(daidT2(j)/ai(j)+2*daidT2(i)/ai(i)-
     1			((daidT(j)/ai(j))**2+2*(daidT(i)/ai(i))**2))*aijk(j,i,i)/3
			daijkdT2(i,i,j)=daijkdT2(j,i,i)
			daijkdT2(i,j,i)=daijkdT2(j,i,i)
		  end if
			IF (j.gt.1) THEN	! only possible with three or more components
				do k=1,j-1
					aijk(k,j,i)=(1-Kijk(k,j,i))*(ai(i)*ai(j)*ai(k))**third
					aijk(k,i,j)=aijk(k,j,i)
					aijk(j,i,k)=aijk(k,j,i)
					aijk(j,k,i)=aijk(k,j,i)
					aijk(i,j,k)=aijk(k,j,i)
					aijk(i,k,j)=aijk(k,j,i)
				if(NTD.EQ.1)then
					! add here daijkdT(k,j,i) & daijkdT2(k,j,i)
				end if
				end do
			END IF
		end do
	END IF
	END DO
	IF(NTDEP.ge.1.and.NTD.EQ.1)THEN
		aux1=daijkdT(1,1,2)
		aux2=daijkdT(1,2,2)
		IF(NTDEP.EQ.1)THEN
			ratK1=(Kijk(1,1,2)-Kinf1)/(1-Kijk(1,1,2))/Tstar1
			ratK2=(Kijk(1,2,2)-Kinf2)/(1-Kijk(1,2,2))/Tstar2
c	first
			daijkdT(1,1,2)=aux1+aijk(1,1,2)*ratK1
			daijkdT(1,2,2)=aux2+aijk(1,2,2)*ratK2
c	second
			daijkdT2(1,1,2)=daijkdT2(1,1,2)+(2*aux1-aijk(1,1,2)/Tstar1)*ratK1
			daijkdT2(1,2,2)=daijkdT2(1,2,2)+(2*aux2-aijk(1,2,2)/Tstar2)*ratK2
		ELSE IF(NTDEP.EQ.2)THEN
			rB1=B1/SUM1
			dK1dT=rB1*(1-rT1)
			d2K1dT2=-2*dK1dT/SUM1
			rB2=B2/SUM2
			dK2dT=rB2*(1-rT2)
			d2K2dT2=-2*dK2dT/SUM2
c	first
			daijkdT(1,1,2)=aux1-dK1dT*aijk(1,1,2)/(1-Kijk(1,1,2))
			daijkdT(1,2,2)=aux2-dK2dT*aijk(1,2,2)/(1-Kijk(1,2,2))
c	second
			daijkdT2(1,1,2)=daijkdT2(1,1,2)
	1				-(2*aux1*dK1dT+aijk(1,1,2)*d2K1dT2)/(1-Kijk(1,1,2))
			daijkdT2(1,2,2)=daijkdT2(1,2,2)
	1				-(2*aux2*dK2dT+aijk(1,2,2)*d2K2dT2)/(1-Kijk(1,2,2))
		END IF
		daijkdT(1,2,1)=daijkdT(1,1,2)
		daijkdT(2,1,1)=daijkdT(1,1,2)
		daijkdT(2,2,1)=daijkdT(1,2,2)
		daijkdT(2,1,2)=daijkdT(1,2,2)
		daijkdT2(1,2,1)=daijkdT2(1,1,2)
		daijkdT2(2,1,1)=daijkdT2(1,1,2)
		daijkdT2(2,2,1)=daijkdT2(1,2,2)
		daijkdT2(2,1,2)=daijkdT2(1,2,2)
	END IF
	end
C
	subroutine DandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
      implicit DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2)
	dimension rn(nco),dDiT(nco)
	dimension dDi(nco),dDij(nco,nco)
	dimension aij(nco,nco),daijdT(nco,nco),daijdT2(nco,nco)
	COMMON /rule/ncomb,iRuleDel
	COMMON /SpecRep/arepfr,arep,dardT,dardT2        ! August 2016

	call aijTder(NTD,nc,T,aij,daijdT,daijdT2)
	D=0.0D0
	dDdT=0.0D0
	dDdT2=0.0D0
	DO i=1,nc
	    aux=0.0D0
	    aux2=0.0D0
	    dDi(i)=0.0D0
	    dDiT(i)=0.0D0
	    do j=1,nc
	        dDi(i)=dDi(i)+2*rn(j)*aij(i,j)
	        if(NTD.EQ.1)then
		        dDiT(i)=dDiT(i)+2*rn(j)*daijdT(i,j)
		        aux2=aux2+rn(j)*daijdT2(i,j)
	        end if
	        dDij(i,j)=2*aij(i,j)
	        aux=aux+rn(j)*aij(i,j)
	    end do
	    D=D+rn(i)*aux
	    if(NTD.EQ.1)then
		    dDdT=dDdT+rn(i)*dDiT(i)/2
		    dDdT2=dDdT2+rn(i)*aux2
	    end if
	END DO
      if(ncomb == -1)then  
        rn2 = rn(1)**2
        D = D - rn2 *arep
        dDi(1) = dDi(1) - 2*rn(1)*arep
        dDij(1,1) = dDij(1,1) - 2*arep
        if(NTD.EQ.1)then
	        dDdT = dDdT - rn2 *dardT
	        dDdT2 = dDdT2 - rn2 *dardT2
	        dDiT(1)=dDiT(1) - 2*rn(1)*dardT
        end if
      end if
	end
C
	subroutine DcubicandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
      implicit DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2)
	dimension rn(nco),dDiT(nco)
	dimension dDi(nco),dDij(nco,nco),aux(nco),auxij(nco,nco)
	dimension auxT(nco),auxTij(nco,nco),auxT2(nco),auxT2ij(nco,nco)
	dimension aijk(nco,nco,nco),daijkdT(nco,nco,nco),daijkdT2(nco,nco,nco)
	call aijkTder(NTD,nc,T,aijk,daijkdT,daijkdT2)
	TOTN = sum(rn)
	D=0.0D0
	dDdT=0.0D0
	dDdT2=0.0D0
	aux=0.0D0
	auxij=0.0D0
	auxT=0.0D0
	auxTij=0.0D0
	DO i=1,nc
		do j=1,nc
			do k=1,nc
				auxij(i,j)=auxij(i,j)+rn(k)*aijk(i,j,k)
				if(NTD.EQ.1)then
					auxTij(i,j)=auxTij(i,j)+rn(k)*daijkdT(i,j,k)
					auxT2ij(i,j)=auxT2ij(i,j)+rn(k)*daijkdT2(i,j,k)
				end if
			end do
			aux(i)=aux(i)+rn(j)*auxij(i,j)
			if(NTD.EQ.1)then
				auxT(i)=auxT(i)+rn(j)*auxTij(i,j)
				auxT2(i)=auxT2(i)+rn(j)*auxT2ij(i,j)
			end if
		end do
		D=D+rn(i)*aux(i)
		if(NTD.EQ.1)then
			dDdT=dDdT+rn(i)*auxT(i)
			dDdT2=dDdT2+rn(i)*auxT2(i)
		end if
	END DO
	D=D/TOTN
	if(NTD.EQ.1)then
		dDdT=dDdT/TOTN
		dDdT2=dDdT2/TOTN
	end if
	DO i=1,nc
		dDi(i)=(3*aux(i)-D)/TOTN
		if(NTD.EQ.1)dDiT(i)=(3*auxT(i)-dDdT)/TOTN
		do j=1,i
			dDij(i,j)=(6*auxij(i,j)-dDi(i)-dDi(j))/TOTN
			dDij(j,i)=dDij(i,j)
		end do
	END DO
	end
C 
	subroutine DELTAnder(nc,rn,D1m,dD1i,dD1ij)
      implicit DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2)
	DOUBLE PRECISION Kij(nco,nco)
	dimension rn(nco),dD1i(nco),dD1ij(nco,nco)
	COMMON /COMPONENTS/ ac(nco),b(nco),d1(nco),rk(nco),Kij,NTDEP
	D1m=0.0D0
	DO i=1,nc
	D1m=D1m+rn(i)*d1(i)
	END DO
	TOTN = sum(rn)
	D1m=D1m/totn
	do i=1,nc
	dD1i(i)=(d1(i)-D1m)/totn
	do j=1,nc
	dD1ij(i,j)=(2.0D0*D1m-d1(i)-d1(j))/totn**2
	end do
	end do
	end
C
	subroutine QuadDELTAnder(nc,rn,D1m,dD1i,dD1ij)
      implicit DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2)
	dimension rn(nco),dD1i(nco),dD1ij(nco,nco),aux(nco)
	COMMON /D1cross/D1ij(nco,nco)
	TOTN = sum(rn)
	D1m=0.0D0
	aux=0.0D0
	DO i=1,nc
		do j=1,nc
			aux(i)=aux(i)+rn(j)*D1ij(i,j)
		end do
		D1m=D1m+rn(i)*aux(i)
	END DO
	D1m=D1m/totn**2
	DO i=1,nc
c		dD1i(i)=(2*aux(i)-D1m)/totn
		dD1i(i)=(2*aux(i)/totn-2*D1m)/totn
		do j=1,i
c			dD1ij(i,j)=(2*D1ij(i,j)-dD1i(i)-dD1i(j))/totn
			dD1ij(i,j)=2*(D1ij(i,j)/totn-dD1i(i)-dD1i(j)-D1m/totn)/totn
			dD1ij(j,i)=dD1ij(i,j)
		end do
	END DO
	end
C
	subroutine Bnder(nc,rn,Bmix,dBi,dBij)
      implicit DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2)
	dimension rn(nco),dBi(nco),dBij(nco,nco),aux(nco)
	COMMON /bcross/bij(nco,nco)
	TOTN = sum(rn)
	Bmix=0.0D0
	aux=0.0D0
	DO i=1,nc
		do j=1,nc
			aux(i)=aux(i)+rn(j)*bij(i,j)
		end do
		Bmix=Bmix+rn(i)*aux(i)
	END DO
	Bmix=Bmix/totn
	DO i=1,nc
		dBi(i)=(2*aux(i)-Bmix)/totn
		do j=1,i
			dBij(i,j)=(2*bij(i,j)-dBi(i)-dBi(j))/totn
			dBij(j,i)=dBij(i,j)
		end do
	END DO
	end
C
	subroutine Bcubicnder(nc,rn,Bmix,dBi,dBij)
      implicit DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2)
	dimension rn(nco),dBi(nco),dBij(nco,nco),aux(nco),auxij(nco,nco)
	COMMON /bcrosscub/bijk(nco,nco,nco)
	TOTN = sum(rn)
	sqn=TOTN*TOTN
	Bmix=0.0D0
	aux=0.0D0
	auxij=0.0D0
	DO i=1,nc
		do j=1,nc
			do k=1,nc
				auxij(i,j)=auxij(i,j)+rn(k)*bijk(i,j,k)
			end do
			aux(i)=aux(i)+rn(j)*auxij(i,j)
		end do
		Bmix=Bmix+rn(i)*aux(i)
	END DO
	Bmix=Bmix/sqn
	DO i=1,nc
		dBi(i)=(3*aux(i)-2*totn*Bmix)/sqn
		do j=1,i
		   dBij(i,j)=(6*auxij(i,j)-2*(Bmix+totn*dBi(i)+totn*dBi(j)))/sqn
			dBij(j,i)=dBij(i,j)
		end do
	END DO
	end
C
	SUBROUTINE HelmRKPR(NDE,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2,RGAS=0.08314472d0)
	dimension rn(nco),Arn(nco),ArVn(nco),ArTn(nco),Arn2(nco,nco)
	dimension dBi(nco),dBij(nco,nco),dD1i(nco),dD1ij(nco,nco)
	dimension dDi(nco),dDij(nco,nco),dDiT(nco)
	dimension aij(nco,nco),daijdT(nco,nco),daijdT2(nco,nco)
	COMMON /rule/ncomb,iRuleDel
	NC=2
	TOTN = sum(rn)
	if(iRuleDel==1)call DELTAnder(nc,rn,D1,dD1i,dD1ij)
	if(iRuleDel==2)call QuadDELTAnder(nc,rn,D1,dD1i,dD1ij)
	D2=(1-D1)/(1+D1)
c
c
c	Comparison to test and debug cubic mixing rules
c	rn=[0.65,0.35]
c	T=460.0d0
c		call Bnder(nc,rn,Bmix,dBi,dBij)
c		call Bcubicnder(nc,rn,Bmix,dBi,dBij)
c		call DandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
c		call DCubicandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
c
c
c
	if(ncomb.lt.2)then
		call Bnder(nc,rn,Bmix,dBi,dBij)
		call DandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
	else
		call Bcubicnder(nc,rn,Bmix,dBi,dBij)
		call DCubicandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
	end if
c	The f's and g's used here are for Ar, not F (reduced Ar)					***********
c	This requires to multiply by R all g, f and its derivatives as defined by Mollerup ****
	f=log((V+D1*Bmix)/(V+D2*Bmix))/Bmix/(D1-D2)
	g=RGAS*log(1-Bmix/V)
	fv=-1/((V+D1*Bmix)*(V+D2*Bmix))
	fB=-(f+V*fv)/Bmix
	gv=RGAS*Bmix/(V*(V-Bmix))
	fv2=(-1/(V+D1*Bmix)**2+1/(V+D2*Bmix)**2)/Bmix/(D1-D2)
	gv2=RGAS*(1/V**2-1/(V-Bmix)**2)
C	DERIVATIVES OF f WITH RESPECT TO DELTA1
	auxD2=(1+2/(1+D1)**2)
	fD1=(1/(V+D1*Bmix)+2/(V+D2*Bmix)/(1+D1)**2)-f*auxD2
	fD1=fD1/(D1-D2)
	fBD1=-(fB*auxD2+D1/(V+D1*Bmix)**2+2*D2/(V+D2*Bmix)**2/(1+D1)**2)
	fBD1=fBD1/(D1-D2)
	fVD1=-(fV*auxD2+1/(V+D1*Bmix)**2+2/(V+D2*Bmix)**2/(1+D1)**2)/(D1-D2)
	fD1D1=4*(f-1/(V+D2*Bmix))/(1+D1)**3+Bmix*(-1/(V+D1*Bmix)**2+
	1		4/(V+D2*Bmix)**2/(1+D1)**4)-2*fD1*(1+2/(1+D1)**2)
	fD1D1=fD1D1/(D1-D2)
c	Reduced Helmholtz Energy and derivatives
	Ar=-TOTN*g*T-D*f
	ArV=-TOTN*gv*T-D*fv
	ArV2=-TOTN*gv2*T-D*fv2
c
	AUX=RGAS*T/(V-Bmix)
	FFB=TOTN*AUX-D*fB
	FFBV=-TOTN*AUX/(V-Bmix)+D*(2*fv+V*fv2)/Bmix
	FFBB=TOTN*AUX/(V-Bmix)-D*(2*f+4*V*fv+V**2*fv2)/Bmix**2
	do i=1,nc
	Arn(i)=-g*T+FFB*dBi(i)-f*dDi(i)-D*fD1*dD1i(i)
	ArVn(i)=-gv*T+FFBV*dBi(i)-fv*dDi(i)-D*fVD1*dD1i(i)
	IF (NDE.EQ.2) THEN
	do j=1,i
	Arn2(i,j)=AUX*(dBi(i)+dBi(j))-fB*(dBi(i)*dDi(j)+dBi(j)*dDi(i))
     &		+FFB*dBij(i,j)+FFBB*dBi(i)*dBi(j)-f*dDij(i,j)      
      Arn2(i,j)=Arn2(i,j)-D*fBD1*(dBi(i)*dD1i(j)+dBi(j)*dD1i(i))
     &		-fD1*(dDi(i)*dD1i(j)+dDi(j)*dD1i(i))
     &		-D*fD1*dD1ij(i,j)-D*fD1D1*dD1i(i)*dD1i(j)
	Arn2(j,i)=Arn2(i,j)
	end do
	END IF
	end do
C	TEMPERATURE DERIVATIVES
	IF (NTD.EQ.1) THEN
	ArT=-TOTN*g-dDdT*f
	ArTV=-TOTN*gv-dDdT*fV
	ArTT=-dDdT2*f
	do i=1,nc
	ArTn(i)=-g+(TOTN*AUX/T-dDdT*fB)*dBi(i)-f*dDiT(i)-dDdT*fD1*dD1i(i)
	end do
	END IF
	end

