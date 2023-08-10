C
	subroutine read2Pcubic(nin,nout)
      implicit DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2,RGAS=0.08314472d0)
C	Critical constants must be given in K and bar
C	b will be in L/mol and ac in bar*(L/mol)**2
	DOUBLE PRECISION Kij(nco,nco),lij(nco,nco)
	DOUBLE PRECISION Kinf,Kinf1,Kinf2,K01,K02,lijk(nco,nco,nco)
	dimension ac(nco),b(nco),del1(nco),rm(nco),diam(nco),Vc(nco),OM(nco)
      CHARACTER*10 fluid(nco)
	COMMON /MODEL/ NMODEL
	COMMON/names/fluid
      COMMON/CRIT/TC(nco),PC(nco),DCeos(nco)
	COMMON /COMPONENTS/ ac,b,del1,rm,Kij,NTdep
	COMMON/COVOL/b1(2)	
	COMMON /bcross/bij(nco,nco)
	COMMON /Tdep/ Kinf,Tstar
	COMMON /Kcubic/Kinf1,Kinf2,K01,K02,Tstar1,Tstar2,C1,C2
	COMMON /bcrosscub/bijk(nco,nco,nco)
	COMMON /rule/ncomb
	COMMON /SpecRep/arepfr,arep,dardT,dardT2        ! August 2016
	NC=2
	read(NIN,*) ncomb,NTDEP
	third=1.0D0/3
	IF(nmodel.eq.1)THEN
		del1=1.0D0
		write(nout,*)' Model: Soave-Redlich-Kwong (1972)'
	ELSE
		del1=1.0D0+sqrt(2.0)
		write(nout,*)' Model: Peng-Robinson (1976)'
	END IF
	write(nout,*)'Fluid          Tc(K)       Pc(bar)  Vceos(L/mol)    W'
	do i=1,nc
	READ(NIN,'(A)')fluid(i)
	READ(NIN,*)Tc(i),Pc(i),OM(i),Vc(i)
	dceos(i)=1/Vc(i)
	write(nout,1)fluid(i),Tc(i),Pc(i),Vc(i),OM(i)
	READ(NIN,*)ac(i),b(i),rm(i)
	Kij(i,i)=0.0D0
	Lij(i,i)=0.0D0
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
	ENDIF
	end do
	B1=B
	write(nout,*)'Fluid  ac(bar*L2/mol2)  b(L/mol)    d1      rm'
	DO I=1,NC
	write(nout,1)fluid(i),ac(i),b(i),del1(i),rm(i)
	END DO
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
 1	FORMAT(A10,F7.3,5x,F7.3,3x,F7.3,3x,F7.3)
 5	FORMAT(A10,F6.3)
 6	FORMAT(A10,4F7.3)
 7	FORMAT(9x,F7.4,2x,F7.4)
 8	FORMAT(9x,F7.2,2x,F7.2)
	end
C
	SUBROUTINE HelmSRKPR(ND,NT,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nco=2,RGAS=0.08314472d0)
	dimension rn(nco),Arn(nco),ArVn(nco),ArTn(nco),Arn2(nco,nco)
	dimension dBi(nco),dBij(nco,nco)
	dimension dDi(nco),dDij(nco,nco),dDiT(nco)
	dimension aij(nco,nco),daijdT(nco,nco),daijdT2(nco,nco)
	DOUBLE PRECISION Kij(nco,nco)
	dimension ac(nco),b(nco),del1(nco),rm(nco)
	COMMON /COMPONENTS/ ac,b,del1,rm,Kij,NTdep
	COMMON /rule/ncomb
	NC=2
	TOTN = sum(rn)
	D1=del1(1)
	D2=(1-D1)/(1+D1)
	if(ncomb.lt.2)then
		call Bnder(nc,rn,Bmix,dBi,dBij)
		call DandTnder(NT,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
	else
		call Bcubicnder(nc,rn,Bmix,dBi,dBij)
		call DCubicandTnder(NT,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
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
	Arn(i)=-g*T+FFB*dBi(i)-f*dDi(i)
	ArVn(i)=-gv*T+FFBV*dBi(i)-fv*dDi(i)
	IF (ND.EQ.2) THEN
	do j=1,i
	Arn2(i,j)=AUX*(dBi(i)+dBi(j))-fB*(dBi(i)*dDi(j)+dBi(j)*dDi(i))
     &		+FFB*dBij(i,j)+FFBB*dBi(i)*dBi(j)-f*dDij(i,j)      
	Arn2(j,i)=Arn2(i,j)
	end do
	END IF
	end do
C	TEMPERATURE DERIVATIVES
	IF (NT.EQ.1) THEN
	ArT=-TOTN*g-dDdT*f
	ArTV=-TOTN*gv-dDdT*fV
	ArTT=-dDdT2*f
	do i=1,nc
	ArTn(i)=-g+(TOTN*AUX/T-dDdT*fB)*dBi(i)-f*dDiT(i)
	end do
	END IF
	end

