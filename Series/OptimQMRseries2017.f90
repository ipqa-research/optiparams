!     Optimizador de Series de sistemas binarios.
!     Adaptado en septiembre 2017 en base al optimizador de series de 2014,
!     que a su vez partía del código "OptimCMRKP2011" para sistemas.
!     Enfocado y aplicado inicialmente para
!       - Series de (C1/C2/C3/CO2) + Alcanos, con RKPR
!       - QMR con Kij(T), ajustando Kij0 y Lij como rectas dependientes de del1
!         del pesado. Ej.K0 serie C1: par = (del1 -0.85)/c
!
!     Historial previo del optimizador por sistemas:
!     26/02/2011 Nueva función objetivo y agregado de posibles
!     KP extra (ver common EXTRAKP)
!
!          versión del 20/03/2010 combinando las últimas mejoras previas con una
!        actualización de la Función Objetivo, para dar mas peso a las
!        fases livianas:
!             - términos 11 a 15 (parte LLV)
!             - términos FV(jf+1)
!             - FV(jf+3) en el DO NTP
!          Al 30/09/2010: agregando Type I
!          Al  3/11/2010: agregados Type II & IV
!
program OptimQMRseries
   implicit double precision(A - H, O - Z)
   logical updateC1, kwithac, Lexp
   CHARACTER*30 INFILE, OUTFILE
   COMMON/UNITS/NUNIT, NOUT, Nsys, updateC1, kwithac, Lexp, nL20
   ! COMMON/EXTRAK/ IntCri, PcInt, XcInt, TcInt, islope, T9art
   COMMON/UNITAUX/NOUT2
   nout2 = 8
   nunit = 1
   nout = 9

   write (6, *) 'ENTER INFILE'
   READ (5, '(A)') INFILE
   ! INFILE='CMRCO2C13.DAT'
   ! OUTFILE='AUXAUX.DAT'
   OPEN (NUNIT, FILE=INFILE)

   write (6, *) 'ENTER A NAME FOR THE OUTFILE'
   READ (5, '(A)') OUTFILE
   OPEN (NOUT, FILE=OUTFILE)
   OPEN (NOUT2, FILE='AUXOUT.DAT')

   write (6, *) 'KEY POINTS REQUIRED AT INPUT FILE FOR EACH TYPE OF PHASE BEHAVIOR'
   write (6, *) ' 1: I   (Pc(T1), Xc(T1), Pc(T2), Xc(T2))'
   write (6, *) ' 2: II  (Pc(T1), Xc(T1), Pc(T2), Xc(T2), Tucep)'
   write (6, *) ' 3: III (T994, Tm, PCPm, P393, Tucep)'
   write (6, *) ' 4: IV  (Pc(T), Xc(T), TUCEP, TLCEP, Tk)'
   write (6, *) ' 5: V   (Pc(T1), Xc(T1), Pc(T2), Xc(T2), TLCEP, Tk)'

   ! islope=0      En principio, no se utilizará T994 para Types 2 or 4
   ! if (NCASE==2.or.NCASE==4) then
   !     write (6,*) ' Use of extrapolated T994?'
   !     write (6,*) ' 1 for YES'
   !     READ (5,*)islope
   ! end if
   read (NUNIT, *) N        ! number of parameters to optimize
   updateC1 = .false.
   if (N == 5 .or. N == 7) then
      write (6, *) 'ENTER 1 FOR UPDATING Component 1 parameters together with the rest'
      write (6, *) 'OTHERWISE (if Component 1 will remain fixed) ENTER 0'
      READ (5, *) nupd
      if (nupd == 1) updateC1 = .true.
   end if

   ! kwithac=.false.
   ! Lexp=.false.
   ! read(NUNIT,*)ik,iL  ! if 1, k0 depends on ac(2)
   ! if (ik==1) kwithac=.true.
   ! if (iL==1) Lexp=.true.

   call OptimQMR(N)
   write (6, *) ' Optimization performed succesfully. Press enter.'
   READ (5, *) file

   close (unit=nunit)
   close (unit=nout)
   close (unit=nout2)
end
!
SUBROUTINE OptimQMR(N)
   PARAMETER(nco=2, maxs=32)
   implicit double precision(A - H, O - Z)

   DOUBLE PRECISION Kij(nco, nco), Kinf, Kinf1, Kinf2, K01, K02

   dimension ac(nco), b(nco), del1(nco), rk(nco)
   DIMENSION X(N), XGUESS(N), XGUES4(4)
   logical updateC1, curve, kwithac, Lexp
   COMMON/CASEvec/Ica(maxs), NK(maxs), NC(maxs)
   COMMON/fixed/nchange
   ! COMMON/EXTRAK/ IntCri, PcInt, XcInt, TcInt, islope, T9art
   COMMON/UNITS/NUNIT, NOUT, Nsys, updateC1, kwithac, Lexp, nL20, iexp
   COMMON/MODEL/NMODEL
   COMMON/Kcubic/Kinf1, Kinf2, K01, K02, Tstar1, Tstar2, C1, C2
   COMMON/DAT/DAT(maxs, 32)
   COMMON/KeyTV/Tc1(maxs), Tc2(maxs)
   COMMON/KeyTa/Ta, Taint
   COMMON/Key2Ph/NTP(maxs), T2p(maxs, 8), P2p(maxs, 8)
   COMMON/KeyIsop/NzP(maxs), NzT(maxs), IZv(maxs, 15), PTsat(maxs, 15), Xis(maxs, 15), Yis(maxs, 15)
   COMMON/KeyFUG/NFUG(maxs), Tfug(maxs, 8), Pfug(maxs, 8), X1fug(maxs, 8), Y1fug(maxs, 8)
   COMMON/PmaxLL/Phigh
   COMMON/rule/ncomb
   COMMON/CRIT/TC(nco), PC(nco), DC(nco)
   COMMON/CRITV/TCV(maxs), PCV(maxs), OM(maxs), DCV(maxs)
   COMMON/COMPONENTS/ac, b, del1, rk, Kij, NTDEP
   COMMON/COMPONENTSV/acV(maxs), bV(maxs), del1V(maxs), rkV(maxs)
   COMMON/sDDLC/q(nco), nqopt
   COMMON/Tdep/Kinf, Tstar
   COMMON/Fifth/i5p, N1, refN, bk, Ad, ck
   DATA FSCALE /1.0E0/

   interface
      function F(X, N)   ! SUBROUTINE ObjFun (N, X, F)
         real(8) :: X(n)
         integer :: n
         real(8) :: F
      end function
   end interface

   character(len=*), parameter :: fmt_1 = "(F10.4, F10.4, F10.6, F9.5, F5.1)"
   character(len=*), parameter :: fmt_2 = "(F10.4, F10.6, F10.6, F10.5)"
   character(len=*), parameter :: fmt_3 = "(i4, 3F10.5)"
   character(len=*), parameter :: fmt_99995 = "('  The solution is ', 6X, 5F9.5, //, '  The function ', 'value is ', F9.6)"
   character(len=*), parameter :: fmt_99992 = "('  The solution is ', 6X, 2F9.5, //, '  The function ', 'value is ', F9.6)"

   XSCALE = 1.0E0
   read (NUNIT, *) NMODEL
   read (NUNIT, *) Nsys
   if (NMODEL .LE. 2) then
      ! SRK or PR
      CALL read2Pcubic(nunit, nout)
   else if (NMODEL .EQ. 3) then
      CALL readcomp(nunit, nout)
   else if (NMODEL .EQ. 4) then
      ! CALL readPCSAFT(nunit,nout)
   else if (NMODEL .EQ. 6) then
      ! CALL readSPHCT(nunit,nout)
   end if
   read (NUNIT, *) Phigh
   read (NUNIT, *) OM(Nsys + 1)   ! Methane acentric factor

   curve = .false.
   write (6, *) 'Enter carbon number of compound 1 defining the serie'
   write (6, *) '1 for Methane, 2 for Ethane, etc.'
   READ (5, *) N1
   if (N >= 3) then
      write (6, *) 'ENTER 1 FOR OBTAINING Pure Component Parameters from exponential curve constants for del1'
      write (6, *) 'OTHERWISE (if params will remain as read) ENTER 0'
      READ (5, *) ncu
      if (ncu == 1) curve = .true.
      if (curve) then
         READ (nunit, *) Ad, Bd, refN
         write (6, *) 'ENTER 1 FOR Ad+Bd*NC(i)*exp(-NC(i)/refN)'
         write (6, *) 'ENTER 2 for Ad+Bd* (1.0-exp(-NC(i)/refN))'
         READ (5, *) iexp
         write (6, *) 'ENTER 1 for fixing bk '
         write (6, *) '   or 2 for fixing refN '
         write (6, *) '   or 3 for fixing Ad '
         write (6, *) '   or 4 for fixing ck '
         write (6, *) 'while optimizing the others'
         READ (5, *) i5p
         if (i5p == 1) write (6, *) 'ENTER fixed value for bk '
         if (i5p /= 1) write (6, *) 'ENTER initial value for bk '
         READ (5, *) bk
      end if
   end if

   if (curve .and. updateC1) then
      if (iexp == 1) del1(1) = Ad + Bd*exp(-1/refN)
      if (iexp == 2) del1(1) = Ad + Bd*(1.0 - exp(-1/refN))
      call paramsfromdelta1(del1(1), Tc(1), Pc(1), OM(Nsys + 1), ac(1), b(1), rk(1), Dc(1))
      write (nout, *) 1
      write (nout, fmt_1) Tc(1), Pc(1), OM(Nsys + 1), 1/Dc(1), 2.0
      write (nout, fmt_2) ac(1), b(1), del1(1), rk(1)
   end if

   do i = 1, Nsys
      read (NUNIT, *) NC(i), Ica(i), NK(i), NTP(i), NzP(i), NzT(i), NFUG(i)
      NCASE = Ica(i)
      READ (nunit, *) Tcv(i), Pcv(i), OM(i), Vceos
      dcv(i) = 1/Vceos

      SELECT CASE (nmodel)
      CASE (1, 2)
         READ (nunit, *) acv(i), bv(i), rkv(i) ! SRK/PR
      CASE (3)
         READ (nunit, *) acv(i), bv(i), del1v(i), rkv(i)  ! RKPR
      end SELECT

      if (curve) then
         if (iexp == 1) del1v(i) = Ad + Bd*NC(i)*exp(-NC(i)/refN)
         ! if(NC(i)>20)del1v(i)=del1v(i)+0.7*(1-exp(-(NC(i)-20)/10.0))  ! 2017
         if (iexp == 2) del1v(i) = Ad + Bd*(1.0 - exp(-NC(i)/refN))
         call paramsfromdelta1(del1v(i), Tcv(i), Pcv(i), OM(i), acv(i), bv(i), rkv(i), Dcv(i))
         write (nout, *) NC(i)
         write (nout, fmt_1) Tcv(i), Pcv(i), OM(i), 1/Dcv(i), 2.0
         write (nout, fmt_2) acv(i), bv(i), del1v(i), rkv(i)
      end if


      ! ========================================================================
      ! Read experimental data
      ! Case  0   NK=0
      ! Case  I  (NK=5 or 7 mean 2 or 4): Pc(T1), Xc(T1), [Pc(T2), Xc(T2)]
      ! Case  II (NK=5 or 7 mean 3 or 5): Pc(T1), Xc(T1), Tu, [Pc(T2), Xc(T2)]            (Xlo,Ylo go optionally in NTP)
      ! Case  III(NK=4 or 5):             T994, Tm, Pm, P393, [Tu]                        (Xlo,Ylo go optionally in NTP)
      ! Case  IV (NK=5 or 7):             Pc(T) , Xc(T) , Tu , TL, Tk, [Pc(T2), Xc(T2)]   (Xlo,Ylo go optionally in NTP)
      ! Case  V  (NK=5 or 7 mean 4 or 6): Pc(T1), Xc(T1), TL, Tk, [Pc(T2), Xc(T2)]        (Xlo,Ylo go optionally in NTP)
      ! ------------------------------------------------------------------------

      do k = 1, NK(i)         ! NK can be [0, 4, 5, 7]
         read (NUNIT, *) DAT(i, k)
      end do

      SELECT CASE (ncase)
      CASE (1, 2, 4, 5)
         read (NUNIT, *) Tc1(i)
         if (NK(i) == 7) read (NUNIT, *) Tc2(i)
      end SELECT

      !  if(NCASE==1.and.NK==6)then  ! case NC=2
      !    read(NUNIT,*)Ta
      !    read(NUNIT,*)TaInt
      !  end if

      do k = 1, NTP(i)
         j = NK(i) + 2*(k - 1) + 1
         read (NUNIT, *) T2p(i, k), P2p(i, k), DAT(i, j), DAT(i, j + 1)
      end do

      do k = 1, NzP(i)
         j = NK(i) + 2*NTP(i) + k
         read (NUNIT, *) DAT(i, j), PTsat(i, k), Xis(i, k), Yis(i, k), IZv(i, k)        ! here PTsat(k) stores a spec. P value
      end do

      do k = NzP(i) + 1, NzP(i) + NzT(i)
         j = NK(i) + 2*NTP(i) + k
         read (NUNIT, *) PTsat(i, k), DAT(i, j), Xis(i, k), Yis(i, k), IZv(i, k)        ! here PTsat(k) stores a spec. T value
      end do

      do k = 1, NFUG(i)
         read (NUNIT, *) Tfug(i, k), Pfug(i, k), X1fug(i, k), Y1fug(i, k)
      end do
   end do
   ! ========================================================================

   ! lo que sigue  deber� simplificarse para leer las constantes  que definen por ej. las rectas de K0 y L (caso QMR)
   if (ncomb == 3) then
      ! CMR - si se usa alg�n d�a, se puede adaptar desde "OptimCMRKP2011"
   else if (ncomb == 2) then
      ! s-DDLC - si se usa algún día, se puede adaptar desde "OptimCMRKP2011"
   else
      ! QMR
      if (N == 5) then        ! RKPR delta1 curve is fitted together with the K0 and L lines
         !  i5p=1 --> X=[ck,dk,Ad,Bd,refN]
         !  i5p=2 --> X=[ck,dk,Ad,Bd,bk ]
         !  i5p=3 --> X=[ck,dk,bk,Bd,refN]
         !  i5p=4 --> X=[bk,dk,Ad,Bd,refN]
         READ (NUNIT, *) XGUESS(1)        ! c for K012
         READ (NUNIT, *) XGUESS(2)        ! d for K012
         ! READ(NUNIT,*) XGUESS(3:5)        !  A, B, NC*
         XGUESS(3:5) = [Ad, Bd, refN/10]
         if (i5p == 2) XGUESS(5) = bk
         if (i5p == 3) XGUESS(3) = bk
         if (i5p == 4) XGUESS(1) = bk
      else if (N == 1) then
         write (NOUT, *) ' Scanning to optimize   (K0 correlation)'
         write (NOUT, *) '    Ak       Fi'
         READ (NUNIT, *) XGUESS(1)        ! Ak
         READ (NUNIT, *) bk
         READ (NUNIT, *) refN        ! Nk for K012 and Kinf
      else  ! for N = 3 or 4
         READ (NUNIT, *) XGUESS(1)   ! c for K012
         READ (NUNIT, *) XGUESS(2)   ! d for K012
         READ (NUNIT, *) XGUESS(3)   ! b for Kinf
         READ (NUNIT, *) refN        ! Nk for K012 and Kinf
         READ (NUNIT, *) ek          ! ek for K012 (June 2018)
         if (N == 4) then
            ! XGUESS(4) = refN/10
            XGUESS(4) = ek   ! (June 2018)
         end if
      end if

      write (6, *) 'print initial Delta1 and K0 vectors? (1 for yes)'
      READ (5, *) nparout
      if (nparout == 1) then
         write (nout, *) ' NC   Delta1(2)    k012    Kinf'
         if (N > 2) then
            ck = XGUESS(1)
            dk = XGUESS(2)
            if (N == 3 .or. N == 4) bk = XGUESS(3)
         end if
         DO i = 1, Nsys
            dif = del1v(i) - del1(1)
            ! if(kwithac)then
            !  rel=acv(i)/ac(1)
            !  rk012=ck*(1.d0-exp(-(rel-1.d0)/refN))-dk*dif  ! K0
            !  Modificación Nati 28-11-17:
            !  rk012= ck + dk *(NC(i)-1)* exp(-(2*(NC(i)-1))/refN) !K0 modif. NT
            if (N == 1) then
               rk012 = ak*(1.d0 - exp(-(NC(i) - N1)/refN))
            else
               ! rk012= (NC(i)-N1)*(ck/NC(i) + dk*exp(-(2*(NC(i)-N1))/refN))
               rk012 = dk*(NC(i) - N1)*exp(-2*(NC(i) - N1)/refN)! June 2018
               rk012 = rk012 + ck*(1.0*(NC(i) - N1)/NC(i))**ek    ! June 2018
            end if
            Kinf = bk*(1.d0 - exp(-(NC(i) - N1)/refN))
            ! rk012=ck/10*(1.d0-exp(-(rel-1.d0)/(10*dk)))  ! K0
            ! else
            ! rk012=dif*ck/10+dk*dif**2/100  ! K0
            ! end if
            ! if(Lexp)then
            ! rl12=cl/1000*(1.d0-exp(10*dif/dl))  ! original
            ! rl12=cl/10*(1.d0-exp(dif/dl))
            ! if(nL20==1.and.NC(i)>20) then
            ! rl12=rl12+CL20/10*(1-exp(-(NC(i)-20)/(10*rL20))) ! extra term for serie C3
            ! end if
            ! else
            ! rl12=dif*cl/10+dl*dif**2/100
            ! end if
            write (nout, fmt_3) NC(i), del1v(i), rk012, Kinf
         end DO
      end if
   end if

   !!     OLD FROM OptimCMRKP2011!!!
   ! Case  I  (NK=2): Pc(T), Xc(T) [possible second Pc, Xc after NFUG]
   ! Case  II (NK=5): Pc(T), Xc(T),Xlo,Ylo,Tu
   ! Case III(NK=10): T994,Tm,Pm,P393,Tu,Xu,Xlo,Ylo   (NKeffective=8, NK=10 with Xmi,Ymi)
   ! Case  IV (NK=8): Pc(T), Xc(T),Xlo,Ylo,Tu,TL,Tk,xk
   ! Case  V  (NK=7): Pc(T), Xc(T),Xlo,Ylo,TL,Tk,xk        ! added 09/01/2013

   if (N == 1) then ! armado para A (K0 correlation) 6/6/2018
      !  7 write(NOUT,*) ' d1 for comp2: ',del1(2)
      rmin = XGUESS(1) - 0.01
      rmax = XGUESS(1) + 0.01
      do rl = rmin, rmax, 0.001
         XGUESS(1) = rl
         OF = F(XGUESS, N)  ! write(NOUT,*)
      end do
      ! write (6,*) ' Another d1 for heavy comp?  1 for YES' ! added 14/11/2014 to allow various e.g. C60 on a single run
      ! READ (5,*)nreply
      ! if(nreply==1)then
      !     write (6,*) ' Enter new del1'
      !     READ (5,*)del1(2)
      !     call paramsfromdelta1(del1(2),Tc(2),Pc(2),om2, ac(2),b(2),rk(2),Dc(2))
      !            write (6,*) ' Enter central Lij to define range'
      !            READ (5,*)XGUESS(1)
      !            go to 7
      !     end if
      return
   end if

   if (N .eq. 7) then
      write (NOUT, *) '     params for N=7?     F'
   else if (N .eq. 5) then
      if (i5p == 1) write (NOUT, *) ' bk (for Kinf):', bk
      if (i5p == 1) write (NOUT, *) '     CK       DK      A      B     NC*     F'
      if (i5p == 2) write (NOUT, *) ' NC*:', refN
      if (i5p == 2) write (NOUT, *) '     CK       DK      A      B     bk      F'
      if (i5p == 3) write (NOUT, *) ' Ad:', Ad
      if (i5p == 3) write (NOUT, *) '     CK       DK      bk     B     NC*     F'
      if (i5p == 4) write (NOUT, *) ' ck:', ck
      if (i5p == 4) write (NOUT, *) '     bk       DK      A      B     NC*     F'
   else if (N .eq. 4) then
      write (NOUT, *) '     cK       dK      bK      ek     F'
   else
      write (NOUT, *) '     cK       dK      bK     F'
   end if

   Fmin = PRAXIS(3.D-5, 2.22D-16, 2.D-2, N, 3, XGUESS, F, 1D-2)

   if (N > 5) then
      write (NOUT, fmt_99995) XGUESS, Fmin
   else if (N .eq. 2) then
      write (NOUT, fmt_99992) XGUESS, Fmin
   end if
end

FUNCTION F(X, N)   ! SUBROUTINE ObjFun (N, X, F)
   PARAMETER(nco=2, maxs=32)
   implicit double precision(A - H, O - Z)
   DOUBLE PRECISION Kij(nco, nco), lij(nco, nco)
   DOUBLE PRECISION Kinf, Kinf1, Kinf2, K01, K02, lijk(nco, nco, nco)
   dimension ac(nco), b(nco), del1(nco), rk(nco)
   DIMENSION X(N), FV(maxs, 60), Fsys(maxs), dfug(2)
   COMMON/UNITS/NUNIT, NOUT, Nsys, updateC1, kwithac, Lexp, nL20, iexp
   COMMON/ABrefN/par(3)
   COMMON/CASEOPT/NCASE
   COMMON/CASEvec/Ica(maxs), NK(maxs), NC(maxs)
   COMMON/fixed/nchange
   COMMON/rule/ncomb
   COMMON/CRITV/TCV(maxs), PCV(maxs), OM(maxs), DCV(maxs)
   COMMON/CRIT/TC(nco), PC(nco), DC(nco)
   COMMON/COMPONENTSV/acV(maxs), bV(maxs), del1V(maxs), rkV(maxs)
   COMMON/COMPONENTS/ac, b, del1, rk, Kij, NTDEP
   COMMON/sDDLC/q(nco), nqopt
   COMMON/bcross/bij(nco, nco)
   COMMON/bcrosscub/bijk(nco, nco, nco)
   COMMON/Kcubic/Kinf1, Kinf2, K01, K02, Tstar1, Tstar2, C1, C2
   ! COMMON/COVOL/b(nco)
   COMMON/DAT/DAT(maxs, 30)
   COMMON/CALC1/Pcr, Xcr
   COMMON/CALCint/Pci, Xci
   COMMON/CALCA/Pa, Xa, Pai, Xai
   COMMON/CALC24/TUCEP, TLCEP
   COMMON/CALC3/T994, Tm, Pm, P393, Tu, Xu, Xlo, Ylo, Xmi, Ymi
   ! COMMON/EXTRAK/ IntCri, PcInt, XcInt, TcInt, islope, T9art
   COMMON/KeyTV/Tc1(maxs), Tc2(maxs)
   COMMON/KeyT1/IntCri, TcLV, TcInt
   COMMON/Key2Ph/NTP(maxs), T2p(maxs, 8), P2p(maxs, 8)
   COMMON/KeyIsop/NzP(maxs), NzT(maxs), IZv(maxs, 15), PTsat(maxs, 15), Xis(maxs, 15), Yis(maxs, 15)
   COMMON/KeyFUG/NFUG(maxs), Tfug(maxs, 8), Pfug(maxs, 8), X1fug(maxs, 8), Y1fug(maxs, 8)
   COMMON/Tdep/Kinf, Tstar
   COMMON/Fifth/i5p, N1, refN, bk, Ad, ck
   logical obtainpure, updateC1, kwithac, Lexp

   Fsys = 0.0D0
   obtainpure = .false.

   if (ncomb == 3) then
   elseif (ncomb == 2) then
   else
      ! QMR
      if (N > 1) then
         if (i5p /= 4) ck = X(1)
         dk = X(2)
      end if
      if (N == 5 .or. N == 7) then
         Bd = X(N - 1)
         if (i5p == 1) then
            Ad = X(N - 2)
            refN = 10*X(N)
         else if (i5p == 2) then
            Ad = X(N - 2)
            bk = X(N)
         else if (i5p == 3) then
            bk = X(N - 2)
            refN = 10*X(N)
         else if (i5p == 4) then
            bk = X(1)
            refN = 10*X(N)
         end if
         obtainpure = .false.
         if (sum(abs(par - [Ad, Bd, refN])) > 0.d0) obtainpure = .true.
         if (obtainpure .and. updateC1) then
            if (iexp == 1) del1(1) = Ad + Bd*exp(-1/refN)
            if (iexp == 2) del1(1) = Ad + Bd*(1.0 - exp(-1/refN))
            call paramsfromdelta1(del1(1), Tc(1), Pc(1), OM(Nsys + 1), ac(1), b(1), rk(1), Dc(1))
         end if
         par = [Ad, Bd, refN]
      else if (N == 1) then
         Ak = X(1)
      else if (N == 3 .or. N == 4) then
         bk = X(3)        ! b for Kinf
         if (N == 4) then
            ! refN = 10*X(4)        ! Nk for K012 and Kinf
            ek = X(4)        ! June 2018
         end if
      end if
   end if

   ! i     1   2   3   4   5   6   7
   ! DAT  Pc1 zc1  Tu  TL  Tk Pc2 zc2       NK
   ! tI    *   *   0   0   0  (*) (*)      5 or 7 -mean 2 or 4-
   ! tII   *   *   *   0   0  (*) (*)      5 or 7 -mean 3 or 5-
   ! tIV   *   *   *   *   *  (*) (*)      5 or 7
   ! tV    *   *   0   *   *  (*) (*)      5 or 7 -mean 4 or 6-
   ! tIII T994 Tm Pm P393 (Tu)             4 or 5
   
   FV = 0.0D0
   DO i = 1, Nsys
      NCASE = Ica(i)
      IntCri = 0
      if (NK(i) == 7) IntCri = 1
      Tc(2) = Tcv(i)
      Pc(2) = Pcv(i)

      if (N == 2 .or. N == 4) then
         ac(2) = acv(i)
         b(2) = bv(i)
         del1(2) = del1v(i)
         rk(2) = rkv(i)
         Dc(2) = Dcv(i)
      else  ! for N=7: optimizing del1 curve together with K0 and L
         if (obtainpure) then
            if (iexp == 1) del1(2) = Ad + Bd*NC(i)*exp(-NC(i)/refN)
            ! if(NC(i)>20) del1(2)=del1(2)+0.7*(1-exp(-(NC(i)-20)/10.0))  ! 2017
            if (iexp == 2) del1(2) = Ad + Bd*(1.0 - exp(-NC(i)/refN))
            call paramsfromdelta1(del1(2), Tc(2), Pc(2), OM(i), ac(2), b(2), rk(2), Dc(2))
            del1v(i) = del1(2)
            acv(i) = ac(2)
            bv(i) = b(2)
            rkv(i) = rk(2)
            Dcv(i) = Dc(2)
         else ! repeating pure params from last call to F
            del1(2) = del1v(i)
            ac(2) = acv(i)
            b(2) = bv(i)
            rk(2) = rkv(i)
            Dc(2) = Dcv(i)
         end if
      end if
      dif = del1(2) - del1(1)
      rel = acv(i)/ac(1)
      !Kij(1,2)=ck*(1.d0-exp(-(rel-1.d0)/refN))-dk*dif  ! K0
      if (N == 1) then
         Kij(1, 2) = ak*(1.d0 - exp(-(NC(i) - N1)/refN)) !MCD 6/6/18
      else
         !Kij(1,2)=(NC(i)-N1)*(ck/NC(i) + dk*exp(-(2*(NC(i)-N1))/refN)) !NGT 28-11-17
         Kij(1, 2) = dk*(NC(i) - N1)*exp(-2*(NC(i) - N1)/refN)  ! June 2018
         AUX = ck*(1.0*(NC(i) - N1)/NC(i))**ek   ! June 2018
         Kij(1, 2) = Kij(1, 2) + AUX
      end if
      Kinf = bk*(1.d0 - exp(-(NC(i) - N1)/refN))
      ! if (N==2) then
      ! else if (N==4.or.N==7) then
      !  if(kwithac)then
      !     rel=acv(i)/ac(1)
      !     if(nL20==1)rel=acv(i)/202.2071   ! ac24 (ac20=157.2669) special for serie C3
      !     Kij(1,2)=ck/10*(1.d0-exp(-(rel-1.d0)/(10*dk)))  ! K0c
      !  else
      !     Kij(1,2)=dif*ck/10+dk*dif**2/100  ! K0
      !  end if
      ! end if
      do k = 1, nco
      do j = k, nco
         bij(k, j) = (1 - lij(k, j))*(b(k) + b(j))/2
         bij(j, k) = bij(k, j)
      end do
      end do

      if (NK(i) /= 0) then
         TcLV = Tc1(i)
         TcInt = Tc2(i)

         CALL KPfromPAR   ! for each system (this is inside the DO loop from 1 to Nsys)

         if (NCASE == 3) then
            FV(i, 1) = (T994 - DAT(i, 1))**2/DAT(i, 1)
            FV(i, 2) = (Tm - DAT(i, 2))**2/DAT(i, 2)
            FV(i, 3) = (Pm - DAT(i, 3))**2/DAT(i, 3)
            FV(i, 4) = (P393 - DAT(i, 4))**2/DAT(i, 4)
            if (NK(i) == 5) then  ! Tu is optional
               FV(i, 5) = (Tu - DAT(i, 5))**2/DAT(i, 5)
            end if
         else  ! types  I, II, IV, V
            FV(i, 1) = (Pcr - DAT(i, 1))**2/DAT(i, 1)
            FV(i, 2) = abs(log((Xcr/DAT(i, 2))))
            FV(i, 3) = abs(log((1.d0 - Xcr)/(1.d0 - DAT(i, 2))))
            ! if(NCASE==1.and.NK==6)then  ! case NC=2
            !    FV(i,4)=(Pa-DAT(i,3))**2/DAT(i,3)
            !    FV(i,5)=abs(log(Xa/DAT(i,4)))
            !    FV(i,6)=abs(log((1.d0-Xa)/(1.d0-DAT(i,4))))
            !    FV(i,7)=(Pai-DAT(i,5))**2/DAT(i,5)
            !    FV(i,8)=abs(log(Xai/DAT(i,6)))
            !    FV(i,9)=abs(log((1.d0-Xai)/(1.d0-DAT(i,6))))
            ! end if
            if (NCASE > 1) then
               if (NCASE .NE. 5) FV(i, 4) = (TUCEP - DAT(i, 3))**2/DAT(i, 3)        ! types II/IV
               if (NCASE >= 4) then
                  FV(i, 5) = (TLCEP - DAT(i, 4))**2/DAT(i, 4)        !10*
                  FV(i, 6) = (Tu - DAT(i, 5))**2/DAT(i, 5)        !10*
               end if
            end if
            jstd = 6 + 4*NTP(i) + NzP(i) + NzT(i) + NFUG(i) ! 6 is NKstd+1 (xcr occupies 2 places in FV)
            if (NK(i) == 7) then  ! 2nd Crit
               FV(i, jstd + 1) = (Pci - DAT(i, 6))**2/DAT(i, 6)
               FV(i, jstd + 2) = abs(log((Xci/DAT(i, 7))))
               FV(i, jstd + 3) = abs(log((1.d0 - Xci)/(1.d0 - DAT(i, 7))))
            end if
            ! if (islope == 1) then
            !     j=jstd+IntCri*3+1
            !     FV(i,j)=(T994-T9art)**2/T9art/10
            ! end if
         end if
      end if

     if (maxval(FV(i, 1:10)) .LT. 100.0D0) then !condition: it doesn't help to add these contributions when F is already > 100
         do k = 1, NTP(i)
            j = NK(i) + 2*(k - 1) + 1
            jf = 7 + 4*(k - 1)
            T = T2p(i, k)
            P = P2p(i, k)
            X1ini = DAT(i, j)
            Y2ini = 1.0d0 - DAT(i, j + 1)
            sumz = 1.0
            m = 0
            do while (abs(sumz - 1.0d0) < 0.01 .and. m .lt. 10)
               m = m + 1
               x1 = x1ini
               Y2 = Y2ini
               call PTpointBin(P, T, X1, Y2)
               sumz = x1 + y2
               x1ini = 0.97*x1ini
            end do
            if (sumz .gt. 1.05) then
               aux = x1
               x1 = 1.0d0 - y2
               y2 = 1.0d0 - aux
            end if
            FV(i, jf) = abs(log(X1/DAT(i, j)))
            FV(i, jf + 1) = abs(log((1.0d0 - X1)/(1.d0 - DAT(i, j))))
            FV(i, jf + 2) = abs(log((1.0d0 - Y2)/DAT(i, j + 1)))
            FV(i, jf + 3) = abs(log(Y2/(1.0d0 - DAT(i, j + 1))))
         end do

         do k = 1, NzP(i)
            j = NK(i) + 2*NTP(i) + k
            jf = k + 6 + 4*NTP(i)
            Tini = DAT(i, j)
            Tini0 = Tini
            deltaT = -0.01*Tini0
            P = PTsat(i, k)
            IZ = IZv(i, k)
            X1 = Xis(i, k)
            Y2 = max(1.0d0 - Yis(i, k), 1.0d-6)
            if (IZ .EQ. 2) x1ini = x1
            if (IZ .EQ. 1) Y2ini = Y2
            sumz = 1.0
            m = 0
            do while (sumz .gt. 0.99 .and. m .lt. 10)
               m = m + 1
               T = Tini
               if (IZ .EQ. 2) x1 = x1ini
               if (IZ .EQ. 1) Y2 = Y2ini
               call PzpointBin(IZ, P, T, X1, Y2)
               sumz = x1 + y2
               deltaT = -1.2*deltaT
               Tini = Tini0 + deltaT
               Tini = Tini
               if (IZ .EQ. 2) x1ini = 0.98*x1ini
            end do
            FV(i, jf) = (T - DAT(i, j))**2/DAT(i, j)
         end do

         do k = NzP(i) + 1, NzP(i) + NzT(i)
            j = NK(i) + 2*NTP(i) + k
            jf = k + 6 + 4*NTP(i)
            Pini = DAT(i, j)
            Pini0 = Pini
            deltaP = -0.03*Pini0
            T = PTsat(i, k)
            IZ = IZv(i, k)
            X1 = Xis(i, k)
            Y2 = max(1.0d0 - Yis(i, k), 1.0d-6)
            if (IZ .EQ. 2) x1ini = x1
            if (IZ .EQ. 1) Y2ini = Y2
            sumz = 1.0
            m = 0

            do while (sumz .gt. 0.98 .and. m .lt. 10)
               m = m + 1
               P = Pini
               if (IZ .EQ. 2) x1 = x1ini
               if (IZ .EQ. 1) Y2 = Y2ini
               call TzpointBin(IZ, P, T, X1, Y2)
               sumz = x1 + y2
               deltaP = -1.3*deltaP
               Pini = Pini0 + deltaP
               if (Pini .lt. Pini0) Pini = (Pini0 + Pini)/2
               if (IZ .EQ. 2) x1ini = 0.98*x1ini
            end do

            if (sumz .gt. 0.99) then
               continue
            end if

            FV(i, jf) = (P - DAT(i, j))**2/DAT(i, j)
         end do
      end if

      do k = 1, NFUG(i)
         jf = k + 6 + 4*NTP(i) + NzP(i) + NzT(i)
         T = Tfug(i, k)
         P = Pfug(i, k)
         X1 = X1fug(i, k)
         Y1 = Y1fug(i, k)
         call PTxyFUG(P, T, X1, Y1, dfug)
         FV(i, jf) = dfug(1)**2 + dfug(2)**2
      end do

      Fsys(i) = SUM(FV(i, :))

   end do

   F = SUM(Fsys)
   write (NOUT, 99) (X(L), L=1, N), (Fsys(L), L=1, Nsys), F
   ! F = F - sum(FV(1:4))
   ! write (NOUT,99) (X(L),L=1,N), F
98 FORMAT(8F9.5, A15)
99 FORMAT(5F9.4, 9F8.4, E13.5)

   RETURN
end

SUBROUTINE KPfromPAR
   PARAMETER(nco=2, NA=2)
   implicit double precision(A - H, O - Z)
   DIMENSION X(nco)
   DIMENSION VS(5), VE(5)
   COMMON/UNITAUX/NOUT
   COMMON/CASEOPT/NCASE
   COMMON/CRIT/TC(nco), PC(nco), DC(nco)
   COMMON/covol/B(2)
   COMMON/LCEP/TL, PL, XcL, DcL, YL, DVL
   COMMON/PmaxLL/Phigh
   COMMON/PAEP/TP(NA), PP(NA), ZP(NA), DLP(NA), DVP(NA), NPAEP
   COMMON/HAEP/TH(NA), PH(NA), ZH(NA), DLH(NA), DVH(NA), NHAEP
   COMMON/CAEP/TAC(NA), PAC(NA), ZAC(NA), DAC(NA), NCAEP
   NPAEP = 0
   NHAEP = 0
   NCAEP = 0

   NS = 1
   if (TC(2) > TC(1)) then  ! typical case
      delXS = 1.1d-4        ! X(1) increases as moving away from C2
      T = TC(2)
      x(1) = delXS
      x(2) = 1.0D0 - delXS
      V = 1.0D0/DC(2)
   else   ! methane
      delXS = -1.1d-4        ! X(1) decreases as moving away from C1
      T = TC(1)
      x(1) = 1.0D0 + delXS
      x(2) = -delXS
      V = 1.0D0/DC(1)
   end if
   call XTVnewtonCrit(nout, NS, delXS, X, T, V)
   if (X(1) .gt. 0.9999) then ! type I or II
      ntype = 1
      if (NCASE == 1) then
         if (NCAEP > 0) then
            go to 7
         else
            go to 10
         end if
      end if
      ! if(NCASE.gt.2) go to 10
   else if (TL .NE. 0.0) then ! type V or IV
      ntype = 5
      if (NCASE .lt. 4) go to 10
   else
      ntype = 3
      go to 3
   end if

   if (NCASE == 5) go to 3                       ! added 09/01/2013
   ! calculation low T L-L critical from high P (line B)
111 T = TC(1)
   Plower = Phigh - 30.0
   call FindHighPcrit(Plower, X, T, Vlow)
   if (Vlow .eq. 0) go to 3   ! no LL crit line was found
   call FindHighPcrit(Phigh, X, T, V)
   NS = 3
   delXS = log(Vlow/V) ! Molar Volume can increase or decrease
   ! as pressure goes down along the LL critical line
   call XTVnewtonCrit(nout, NS, delXS, X, T, V)
   if (ntype .eq. 1) ntype = 2
   if (ntype .eq. 5) ntype = 4
3  if (ntype .lt. 3) go to 12

   ! calculation of critical line (D) starting at CP1
   V = 1.0D0/DC(1)
   T = TC(1)
   NS = 1
   delXS = -1.0d-5        ! X(1) decreases as moving away from C1
   delXS = delXS*(B(1)/B(2))**2        ! Higher asymmetry requires lower delXS to find the UCEP
   x(1) = 1.0D0 + delXS
   x(2) = -delXS
   call XTVnewtonCrit(nout, NS, delXS, X, T, V)
   ! 12 if(ntype.gt.1) call LLVlinesSpecDifCont
   !
99 write (nout, *)
   ! write(nout,*)' Type of phase behaviour predicted by the model
   ! 1 for this system'
   ! write(nout,*) ntype
   ! Azeotropic End Points
12 NAEP = NPAEP + NCAEP + NHAEP
   ! write(nout,*)
   ! write(nout,*)' Total number of Azeotropic End Points found:'
   ! write(nout,*) NAEP
   ! write(nout,*)
   ! write(nout,*)' Pure Azeotropic End Points found:         ',NPAEP

   if (NPAEP .GT. 0) then
      write (nout, *) '   T(K)     P(bar)      z    DL(mol/L)  DV(mol/L)'
      do i = 1, NPAEP
         write (nout, 9) TP(i), PP(i), ZP(i), DLP(i), DVP(i)
      end do
   end if

   ! write(nout,*)
   ! write(nout,*)' Critical Azeotropic End Points found:     ',NCAEP

   if (NCAEP .GT. 0) then
      write (nout, *) '   T(K)     P(bar)      z    DL(mol/L)  DV(mol/L)'
      do i = 1, NCAEP
         write (nout, 9) TAC(i), PAC(i), ZAC(i), DAC(i), DAC(i)
      end do
   end if

   ! write(nout,*)
   ! write(nout,*)' Heterogeneous Azeotropic End Points found: ',NHAEP

   if (NHAEP .GT. 0) then
      write (nout, *) '   T(K)     P(bar)      z    DL(mol/L)  DV(mol/L)'
      do i = 1, NHAEP
         write (nout, 9) TH(i), PH(i), ZH(i), DLH(i), DVH(i)
      end do
   end if

   ! write(nout,*)
   if (NAEP .GT. 0) then
      ! Calculation of azeotropic line(s)
      ! if(NHAEP.EQ.2)then ! two lines
      !    ME=3
      !    if(NCAEP.GT.0)MS=2
      !    if(NCAEP.GT.0)VS=[TAC(1),PAC(1),ZAC(1),DAC(1),DAC(1)]
      !    if(NPAEP.EQ.2)MS=1
      !    if(NPAEP.EQ.2)VS=[TP(1),PP(1),ZP(1),DLP(1),DVP(1)]        ! PAEP pure 1st comp
      !    VE=[TH(1),PH(1),ZH(1),DLH(1),DVH(1)]  ! higher T HAEP
      !    CALL AzeotNewton(MS,VS,ME,VE)                                                ! first line
      !    write(nout,*)'fin'
      !
      !    MS=1
      !    VS=[TP(NPAEP),PP(NPAEP),ZP(NPAEP),DLP(NPAEP),DVP(NPAEP)] ! PAEP pure 2nd comp
      !    VE=[TH(2),PH(2),ZH(2),DLH(2),DVH(2)]  !  lower T HAEP
      ! else  ! one line
      !    ME=0
      !    if(NPAEP.GT.0)then
      !       MS=1
      !       VS=[TP(1),PP(1),ZP(1),DLP(1),DVP(1)]
      !       if(NPAEP.EQ.2)then
      !          ME=1
      !          VE=[TP(2),PP(2),ZP(2),DLP(2),DVP(2)]
      !     else if(NCAEP.GT.0)then
      !       ME=2
      !       VE=[TAC(1),PAC(1),ZAC(1),DAC(1),DAC(1)]
      !     end if
      !  else if (NCAEP .GT. 0) then
7     MS = 2
      VS = [TAC(1), PAC(1), ZAC(1), DAC(1), DAC(1)]
      if (NCAEP .EQ. 2) then
         ME = 2
         VE = [TAC(2), PAC(2), ZAC(2), DAC(2), DAC(2)]
      end if
      ! end if
      ! if(NHAEP.GT.0)then
      !    ME=3
      !    VE=[TH(1),PH(1),ZH(1),DLH(1),DVH(1)]
      ! end if
      ! end if
      ! if(ME.EQ.0)VE(1)=MIN(50.0,VS(1)/2)
      ! CALL AzeotNewton(MS,VS,ME,VE)             ! Activate for azeotropic systems!
      !
   end if
9  FORMAT(F9.4, 2x, E10.4, F8.4, F9.5, x, E11.3)
10 end

subroutine FindHighPcrit(Phigh, X, T, V)
   implicit double precision(A - H, O - Z)
   PARAMETER(nco=2)
   DIMENSION X(2), Z(2), FUG(2), FT(2), FP(2), FX(2, 2)
   LOGICAL LASTDELX
   LASTDELX = .FALSE.
   ! loop to determine high pressure branch
   !
   CALL HIPRES(Phigh, T, XX, VAL)
   if (T .lt. 5.0) go to 120
   DELT = 10.0
   if (VAL .LT. 0.) DELT = -10.0
   CALL TSTEP(Phigh, T, XX, DELT, V, IEX)

   ! iex negative: No hi-P line found
   if (IEX .NE. 0) GOTO 120
   CALL HIPRES(Phigh, T, XX, VAL)
   DELT = .1*DELT
   ! NOW DELT IS NOT USED AS A STEP BUT ONLY AS A TOLERANCE (AND FIRST STEP)
   CALL TSTEP(Phigh, T, XX, DELT, V, IEX)
   ! CALL HIPRES (Phigh,T,XX,VAL)
   ! Instead of another calling to HIPRES we do the following to find XX accurately (error<0.0001):
   ! Otherwise the V obtained, when specified later, will lead to a different pressure (asymmetric systems)
   Z(1) = XX - 0.01
   Z(2) = 1.D0 - Z(1)
   CALL TERMO(1, 3, IC, T, Phigh, Z, V, FUG, FT, FP, FX)
   DELX = 0.001
100 VAL = FX(1, 2)
   Z(1) = Z(1) + DELX
   Z(2) = 1.D0 - Z(1)
   CALL TERMO(1, 3, IC, T, Phigh, Z, V, FUG, FT, FP, FX)
   if (FX(1, 2) .GT. VAL) GO TO 100
   Z(1) = Z(1) - DELX
   if (LASTDELX) GO TO 101
   Z(1) = Z(1) + 0.0001
   Z(2) = 1.D0 - Z(1)

   CALL TERMO(1, 3, IC, T, Phigh, Z, V, FUG, FT, FP, FX)

   if (FX(1, 2) .GT. VAL) then
      DELX = 0.0001
   else
      DELX = -0.0001
      Z(1) = Z(1) + 2*DELX
      Z(2) = 1.D0 - Z(1)
      CALL TERMO(1, 3, IC, T, Phigh, Z, V, FUG, FT, FP, FX)
      if (FX(1, 2) .LT. VAL) GO TO 101
   end if

   LASTDELX = .TRUE.
   GO TO 100
   ! ===========================================================================
101 XX = Z(1)
   DELT = .1*DELT
   CALL TSTEP(Phigh, T, XX, DELT, V, IEX)
   X(1) = XX
   X(2) = 1.0D0 - X(1)
120 end

SUBROUTINE HIPRES(Phigh, T, X, VAL)
   !
   !     purpose of routine HIPRES:
   !
   !     To find the composition where the derivative dlnphi1/dn2
   !     takes on its maximum value at given T and high pressure (PHI)
   !
   !     Parameters:
   !
   !     T       (I)       Temperature
   !     X       (O)       Composition with LARGEST 2nd derivative
   !     Val     (O)       Value of 2nd derivative - 1
   !
   IMPLICIT DOUBLE PRECISION(A - H, O - Z)
   PARAMETER(MAXC=2)
   DIMENSION Z(MAXC), FUG(MAXC), FT(MAXC), FP(MAXC), FX(MAXC, MAXC)
   DIMENSION FTAB(0:50)
   PHI = Phigh
   !
   ! CALCULATE TABLE OF VALUES OF 2ND DERIVATIVE
   !
1  STEP = 0.02D0
   DO K = 0, 50
      Z(1) = DBLE(K)/50
      Z(2) = 1.D0 - Z(1)
      CALL TERMO(1, 3, IC, T, PHI, Z, V, FUG, FT, FP, FX)
      FTAB(K) = FX(1, 2) - 1.D0
   end DO
   !
   !     CALCULATE MAXVAL
   !
   INUM = MAXLOC(FTAB, DIM=1) - 1
   if (INUM .EQ. 0 .or. INUM .EQ. 50) then
      T = 0.9*T
      if (T .lt. 5.0) return
      go to 1
   end if
   DERV1 = (FTAB(INUM + 1) - FTAB(INUM - 1))/(2.D0*STEP)
   DERV2 = (FTAB(INUM + 1) - 2.D0*FTAB(INUM) + FTAB(INUM - 1))/STEP**2
   !
   !     INTERPOLATE TO OPTIMUM
   !
   DELX = -DERV1/DERV2
   X = DBLE(INUM)/50 + DELX
   VAL = FTAB(INUM) + DELX*(DERV1 + .5D0*DELX*DERV2)
end

SUBROUTINE TSTEP(Phigh, T, X, DELT, VOLU, IEX)
   !
   !
   !     Routine to locate T where instability first occures
   !     Routine varies T, with composition X fixed.
   !
   !
   !        NOW DELT IS NOT USED AS A STEP BUT ONLY AS A TOLERANCE (AND FIRST STEP)
   IMPLICIT DOUBLE PRECISION(A - H, O - Z)
   PARAMETER(MAXC=2)
   DIMENSION Z(MAXC), FUG(MAXC), FT(MAXC), FP(MAXC), FX(MAXC, MAXC)
   PARAMETER(TMIN=20.D0, TMAX=1500.D0)
   ! PARAMETER (TMIN=0.D0,TMAX=1500.D0)
   PHI = Phigh
   IEX = 0
   !
   !     CALCULATE TABLE OF VALUES OF 2ND DERIVATIVE
   !
   Z(1) = X
   Z(2) = 1.D0 - X
   CALL TERMO(1, 3, IC, T, PHI, Z, V, FUG, FT, FP, FX)
   VAL = FX(1, 2) - 1.D0
   TOLD = T
   T = T + DELT
100 CONTINUE
   if (T .GT. TMAX .OR. T .LT. TMIN) then
      write (*, *) 'T-limit for High-P search exceeded '
      IEX = 1
      VOLU = 0.0D0
      RETURN
   end if
   CALL TERMO(1, 3, IC, T, PHI, Z, V, FUG, FT, FP, FX)
   ! det=FX(1,1)*FX(2,2)-FX(1,2)*FX(1,2)
   VOLD = VAL
   VAL = FX(1, 2) - 1.D0
   AUX = T
   SLOPE = (T - TOLD)/(VAL - VOLD)
   if (SLOPE .LT. 0) then
      T = T - VAL*SLOPE
   else
      T = T + 3*DELT
   end if
   TOLD = AUX
   if (ABS(T - TOLD) .LT. ABS(DELT)) GO TO 101
   if (T .GT. TMAX .AND. TMAX - TOLD .GT. 100) T = (TMAX + TOLD)/2
   if (T .LT. TMIN .AND. TOLD - TMIN .GT. 2) T = (TMIN + TOLD)/2
   GOTO 100
101 VOLU = V
end

subroutine FindMaxIsotPure(icomp, T, Pmax1)
   implicit double precision(A - H, O - Z)
   PARAMETER(nco=2, RGAS=0.08314472d0, eps=1.0d-7)
   dimension rn(nco), Arn(nco), ArVn(nco), ArTn(nco), Arn2(nco, nco)
   COMMON/MODEL/NMODEL
   COMMON/covol/B(2)
   COMMON/VminVapour/Vmin
   COMMON/NG/NGR
   NG = NGR
   ! if(NMODEL.EQ.5)CALL PARAGC(T,NCO,NG,1)
   timesb = 10000
   NDER = 0
   NTEMP = 0
   rn = 0.0
   rn(icomp) = 1.0
   RT = RGAS*T
   tol = 10*eps/B(icomp)
   delrho = 1
   V = timesb*B(icomp)
1  call ArVnder(NDER, NTEMP, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   d2Pdrho = V**3*ArV2                ! =V*(dPdrho-RGAS*T)
   rho = -RT/d2Pdrho
   epsrho = eps/B(icomp)
   DO WHILE (delrho .gt. tol)
      V = 1/rho
      call ArVnder(NDER, NTEMP, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
      dPdrho = RT + V*V*ArV2
      V = 1/(rho + epsrho)
      call ArVnder(NDER, NTEMP, rn, V, T, Ar, ArV, ArTV, ArV2e, Arn, ArVn, ArTn, Arn2)
      d2Pdrho = V**2*(ArV2e - ArV2)/epsrho
      delrho = -dPdrho/d2Pdrho
      rho = rho + delrho
   end DO
   Vmin = 1/rho
   call ArVnder(NDER, NTEMP, rn, Vmin, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   Pmax1 = 1.5*rho*RT - Arv  ! added 1.5* on 05/03/2011
end subroutine

subroutine AzeotNewton(MS, VS, ME, VE)
   implicit double precision(A - H, O - Z)
   PARAMETER(RGAS=0.08314472d0, nco=2)
   ! S = starting point
   ! E = ending point
   !
   ! M is for the type of starting or ending point:
   !   1        Pure
   !   2        Critical
   !   3        Heterogeneous
   !   0        Open down to TE=MIN(50.0,TS/2)
   !
   ! The independent variables are lnT,lnz,lnvL,lnvV
   DIMENSION z(2), delX(4), dold(4), sensmod(4), b(4), ipiv(4)
   DIMENSION XVAR(4), F(4), dFdS(4), dXdS(4), RJAC(4, 4), AJ(4, 4), XOLD(4)
   DIMENSION FUGx(2), FUGy(2), FUGTx(2), FUGTy(2), DPDNx(2)
   DIMENSION FUGVx(2), FUGVy(2), DFGNx(2, 2), DFGNy(2, 2)
   DIMENSION VS(5), VE(5)
   LOGICAL LASTPOINT, CALCAZ
   ! COMMON/UNITS/NUNIT,NOUT
   COMMON/Pder/DPDN(2), DPDT, DPDV
   COMMON/CALCA/Pa, Xa, Pai, Xai
   COMMON/KeyTa/Ta, Taint
   TOL = 1.0D-6 ! for variables
   N = 4        ! DLSARG CONSTANTS (now dgesv in MKL)
   LDA = 4
   ldb = 4
   ! IPATH=1
   LASTPOINT = .FALSE.
   CALCAZ = .FALSE.
   DFDS(4) = 1.0D0
   RJAC(4, 1:4) = 0.0D0
   ! VE(3)=max(VE(3),1.0d-6)  ! to prevent log(0) crash
   ! write(nout,*)
   ! write(nout,*)
   ! 1        '   T(K)     P(bar)     z(1)    z(2)   DL(mol/L)  DV(mol/L)'
   ! write(nout,*)'AZE'
   ! write(NOUT,9) VS(1),VS(2),VS(3),1.D0-VS(3),VS(4),VS(5)
   ! first point initialization
   if (MS .EQ. 1) then ! start from PAEP
      NS = 2        ! composition
      RJAC(4, 2) = 1.0D0
      if (VS(3) .eq. 0.0) then
         zfirst = min(0.005D0, VE(3)/20)
         XVAR(2) = zfirst  ! log(zfirst)
         delXS = 0.002d0
      else
         z2first = min(0.001D0, (1.0 - VE(3))/50)
         XVAR(2) = 1.0D0 - z2first  ! log(1.0D0-z2first)
         delXS = -0.002D0
      end if
      XVAR(3) = log(1/VS(4))
      XVAR(4) = log(1/VS(5))
   else if (MS .EQ. 2) then ! start from CAEP
      NS = 0        ! volume relation
      RJAC(4, 3) = -1.0D0
      RJAC(4, 4) = 1.0D0
      XVAR(3) = log(0.99/VS(4))
      XVAR(4) = log(1.01/VS(5))
      delXS = 0.03D0
      XVAR(2) = VS(3)  ! log(VS(3))
   end if
   XVAR(1) = log(VS(1))
   Told = VS(1) + 0.1
14 NITER = 0
   DMAXOLD = 8.0D0
   DMAX = 7.0D0
   F(4) = 0.0D0
   delX = 0.0D0
   T = exp(XVAR(1))
   z(1) = XVAR(2)  ! exp(XVAR(2))
   z(2) = 1.0D0 - z(1)
   VL = exp(XVAR(3))
   VV = exp(XVAR(4))
   !  Newton procedure for solving the present point
   DO WHILE (DMAX .GT. TOL)
      NITER = NITER + 1
      NVCORREC = 0
21    CALL XTVTERMO(4, T, VL, Px, z, FUGx, FUGTx, FUGVx, DFGNx)
      DPDNx = DPDN
      DPDTx = DPDT
      DPDVx = DPDV
      CALL XTVTERMO(4, T, VV, Py, z, FUGy, FUGTy, FUGVy, DFGNy)
      if ( &
         Px .lt. 0.95*Py .or. Px .gt. 1.05*Py &
         .or. (Py .lt. 1.0D-8 .and. (Px .lt. 0.98*Py .or. Px .gt. 1.02*Py)) &
         ) then
         NVCORREC = NVCORREC + 1
         if (NVCORREC .EQ. 5) RETURN
         VL = VL + (Py - Px)/DPDVx
         go to 21
      end if
      F(1) = log(Px/Py)
      F(2) = FUGx(1) - FUGy(1)
      F(3) = FUGx(2) - FUGy(2)
      RJAC(1, 1) = DPDTx/Px - DPDT/Py
      RJAC(1, 1) = T*RJAC(1, 1)
      RJAC(1, 2) = ((DPDNx(1) - DPDNx(2))/Px - (DPDN(1) - DPDN(2))/Py)  ! z(1)*
      RJAC(1, 3) = VL*DPDVx/Px
      RJAC(1, 4) = -VV*DPDV/Py

      RJAC(2:3, 1) = T*(FUGTx(1:2) - FUGTy(1:2))
      RJAC(2:3, 2) = DFGNx(1:2, 1) - DFGNx(1:2, 2) - (DFGNy(1:2, 1) - DFGNy(1:2, 2))  ! z(1)*
      RJAC(2:3, 3) = VL*FUGVx(1:2)
      RJAC(2:3, 4) = -VV*FUGVy(1:2)

      dold = delx

      !  CALL DLSARG (N, RJAC, LDA, -F, IPATH, delX)
      !  call dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
      b = -F
      AJ = RJAC
      call dgesv(N, 1, AJ, lda, ipiv, b, ldb, info)
      if (info .ne. 0) write (6, *) "error with dgesv in parameter ", info
      delX = b

      DMAXOLD = DMAX
      DMAX = MAXVAL(ABS(DELX))
      if (DMAX/DMAXOLD .GT. 2.0) then ! reduce step until DMAX decreases
         XVAR = XVAR - dold
         delX = dold/2
         DMAX = DMAXOLD
         !  go to 17
      else
         XVAR(3) = log(VL)  ! maybe changed (correction of VL to reduce F1)
      end if
17    XVAR = XVAR + delX
      T = exp(XVAR(1))
      z(1) = XVAR(2)  ! exp(XVAR(2))
      z(2) = 1.0D0 - z(1)
      VL = exp(XVAR(3))
      call Bcalc(z, T, Bmix)
      if (VL .lt. 1.001*Bmix) then
         XVAR = XVAR - delX
         delX = delX/2
         go to 17
      end if
      VV = exp(XVAR(4))
      if (VL .gt. VV) then
         XVAR = XVAR - delX
         delX = delX/10
         go to 17
      end if
   end DO
   DL = 1/VL
   DV = 1/VV
   !        write(NOUT,9)T,Px,Z(1),Z(2),DL,DV,NITER,NS
   !  c  c c c Taken from XTVNextonCrit and adapted:
   if (CALCAZ) then
      Pa = Px   ! Key-point (only CO2+C2)
      Xa = Z(1)
      go to 11
   else
      if (TaInt .ne. 0.0d0) then
         if ((T - TaInt)*(Told - TaInt) .lt. 0.0d0) then
            Pai = Px + (Pold - Px)*(TaInt - T)/(Told - T)
            Xai = Z(1) + (Zold - Z(1))*(TaInt - T)/(Told - T)
         end if
      end if
      if ((T - Ta)*(Told - Ta) .lt. 0.0d0) then
         ! VcLV=V        ! initial values
         Xa = Z(1)
         go to 23
      end if
   end if
   !  c  c  c   c   c   c   c   c   c   c   c   c   c   c   c   c

   if (ME .EQ. 0 .and. T .le. VE(1)) go to 11        ! criteria for stopping at low T
   if (LASTPOINT) then
      write (NOUT, 9) VE(1), VE(2), VE(3), 1.D0 - VE(3), VE(4), VE(5)
      go to 11        ! return
   end if
   ! CALL DLSARG (N, RJAC, LDA, dFdS, IPATH, dXdS)
   !
   ! call dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
   b = dFdS
   AJ = RJAC
   call dgesv(N, 1, AJ, lda, ipiv, b, ldb, info)
   if (info .ne. 0) write (6, *) "error with dgesv in parameter ", info
   dXdS = b

   NSOLD = NS
4  RJAC(4, 1:4) = 0.0D0
   NITER = min(niter, 10)
   dX0dS = dXdS(4) - dXdS(3)
   sensmod = dXdS
   sensmod(2) = 10*sensmod(2)
   J = MAXLOC(abs(sensmod), DIM=1)
   if (abs(sensmod(J)) .gt. abs(dX0dS)) then
      NS = J                        ! Specify the most changing variable for next point
      S = XVAR(J)
      delXS = dXdS(NS)*delXS*5/NITER
      RJAC(4, NS) = 1.0D0
   else
      NS = 0                        ! Specify lnVV-lnVL for next point
      S = XVAR(4) - XVAR(3)
      delXS = (dXdS(4) - dXdS(3))*delXS*5/NITER
      RJAC(4, 3) = -1.0D0
      RJAC(4, 4) = 1.0D0
   end if
   if (T .lt. 50 .and. abs(delXS) .lt. 0.001) return
   if (NS .NE. NSOLD) then
      if (NS .EQ. 0) dXdS = dXdS/(dXdS(4) - dXdS(3))
      if (NS .NE. 0) dXdS = dXdS/dXdS(NS)
   end if
   if (delXS .LT. 0) then
      delXS = max(delXS, -0.07)
      if (NS .EQ. 1) delXS = max(delXS, -0.02) ! Max lnT decrease allowed
      if (NS .EQ. 2) delXS = max(delXS, -0.015) ! Max z decrease allowed
   else
      delXS = min(delXS, 0.07)
      if (NS .EQ. 0) delXS = min(delXS, 0.10, S/3) ! Max volume separation increase allowed
      if (NS .EQ. 1) delXS = min(delXS, 0.02) ! Max lnT increase allowed
      if (NS .EQ. 2) delXS = min(delXS, 0.015) ! Max z increase allowed                /Z(1)
   end if
   XOLD = XVAR
!
7  if (NS .EQ. 0) then
      S = XOLD(4) - XOLD(3) + delXS
   else
      S = XOLD(NS) + delXS
   end if
   XVAR = XOLD + dXdS*delXS        ! Initial estimates for the 7 variables in the next point
   Told = T
   Pold = Px
   Zold = Z(1)
   if (ME .eq. 0) go to 8
   D43 = XVAR(4) - XVAR(3)
   DT = abs(XVAR(1) - log(VE(1)))
   DZ = abs(XVAR(2) - VE(3))        ! log(VE(3))
   if ((D43 .LT. 0.02 .and. ME .eq. 2) .or. (DT*DZ .lt. 0.00001 .and. ME .ne. 2)) then
      LASTPOINT = .TRUE.
      delXS = 2*delXS/3
      go to 7
   end if
8  GO TO 14
   !
23 CALCAZ = .TRUE.
   NS = 1  ! T
   ! V=VcLV
   XVAR(1) = LOG(Ta)
   GO TO 14
9  FORMAT(F9.4, E12.4, 2F9.5, F9.4, E12.4, I4, I2)
11 end

SUBROUTINE XTVTERMO(INDIC, T, V, P, rn, FUGLOG, DLFUGT, DLFUGV, DLFUGX)
   !
   !-------parameters of XTVTERMO (crit. point, LLV and CEP calculations)
   !
   !       rn        mixture mole numbers                     (input)
   !       t         temperature (k)                          (input)
   !       v         volume      (L)                          (input)
   !       p         pressure    (bar)                        (output)
   !       FUGLOG    vector of log. of fugacities (x*phi*P)   (output)        INDIC < 5
   !       DLFUGT    t-derivative of FUGLOG (const. vol,n)    (output)        INDIC = 2 or 4
   !       DLFUGV    vol-derivative of FUGLOG (const temp,n)  (output)        INDIC < 5
   !       DLFUGX    comp-derivative of FUGLOG (const t & v)  (output)        INDIC > 2
   !---------------------------------------------------
   !---  MODifIED AND CORRECTED july 2005
   !---
   !---------------------------------------------------
   IMPLICIT DOUBLE PRECISION(A - H, O - Z)
   PARAMETER(MAXC=2, nco=2, RGAS=0.08314472d0)
   DIMENSION DLFUGX(MAXC, MAXC)
   DIMENSION FUGLOG(MAXC), DLFUGT(MAXC), DLFUGV(MAXC)
   dimension rn(nco), Arn(nco), ArVn(nco), ArTn(nco), Arn2(nco, nco)
   COMMON/MODEL/NMODEL
   COMMON/NG/NGR
   COMMON/Pder/DPDN(nco), DPDT, DPDV

   NG = NGR
   NC = 2
   ! if(NMODEL.EQ.5) CALL PARAGC(T,NC,NG,1)
   NTEMP = 0
   IGZ = 0
   NDER = 1

   if (INDIC .GT. 2) NDER = 2
   if (INDIC .EQ. 2 .OR. INDIC .EQ. 4) NTEMP = 1

   TOTN = sum(rn)
   RT = RGAS*T

   call ArVnder(NDER, NTEMP, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)

   P = TOTN*RT/V - ArV
   DPDV = -ArV2 - RT*TOTN/V**2

   if (INDIC <= 4) then
      ! Z = P*V/(TOTN*RT)
      DPDT = -ArTV + TOTN*RGAS/V

      do I = 1, NC
         if (RN(I) > 0.0) then
            ! FUGLOG(I)=-LOG(Z)+Arn(I)/RT + log(rn(I)/TOTN) + log(P)
            ! FUGLOG(I)=Arn(I)/RT + log(rn(I)/TOTN) + log(P/Z) this crashes at very low T LLV when Z=P=0.000000...
            FUGLOG(I) = Arn(I)/RT + log(rn(I)) + log(RT/V)
            DPDN(I) = RT/V - ArVn(I)
            DLFUGV(I) = -DPDN(I)/RT                                       ! term DPDV/P is cancelled out
            if (NTEMP /= 0) DLFUGT(I) = (ArTn(I) - Arn(I)/T)/RT + 1.D0/T  ! term DPDT/P is cancelled out
         end if
      end do
   end if

   if (NDER >= 2) then
      do I = 1, NC
         do K = I, NC
            DLFUGX(I, K) = Arn2(I, K)/RT ! term 1/TOTN is cancelled out
            DLFUGX(K, I) = DLFUGX(I, K)
         end do
         DLFUGX(I, I) = DLFUGX(I, I) + 1.0/rn(I)
      end do
   end if
end subroutine

SUBROUTINE TERMO(MTYP, INDIC, IC, T, P, rn, V, PHILOG, DLPHIT, DLPHIP, FUGN)
   IMPLICIT DOUBLE PRECISION(A - H, O - Z)
   PARAMETER(MAXC=2, nco=2, RGAS=0.08314472d0)
   DIMENSION FUGN(MAXC, MAXC)
   DIMENSION PHILOG(MAXC), DLPHIT(MAXC), DLPHIP(MAXC), DPDN(MAXC)
   dimension rn(nco), Arn(nco), ArVn(nco), ArTn(nco), Arn2(nco, nco)
   ! The output PHILOG is actually the vector ln(phi(i)*P)
   NC = 2
   NTEMP = 0
   IGZ = 0
   NDER = 1

   if (INDIC .GT. 2) NDER = 2
   if (INDIC .EQ. 2 .OR. INDIC .EQ. 4) NTEMP = 1
   TOTN = sum(rn)
   if (P .le. 0.0d0) MTYP = 1

   CALL VCALC(MTYP, NC, rn, T, P, V)
   RT = RGAS*T
   Z = V/(TOTN*RT)        ! this is Z/P

   call ArVnder(NDER, NTEMP, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   DPV = -ArV2 - RT*TOTN/V**2
   DPDT = -ArTV + TOTN*RGAS/V

   DO I = 1, NC
      PHILOG(I) = -LOG(Z) + Arn(I)/RT
      DPDN(I) = RT/V - ArVn(I)
      DLPHIP(I) = -DPDN(I)/DPV/RT - 1.D0/P
      if (NTEMP /= 0) DLPHIT(I) = (ArTn(I) - Arn(I)/T)/RT + DPDN(I)*DPDT/DPV/RT + 1.D0/T
   end do

   if (nder >= 2) then
      do I = 1, NC
         do K = I, NC
            FUGN(I, K) = 1.D0/TOTN + (Arn2(I, K) + DPDN(I)*DPDN(K)/DPV)/RT
            FUGN(K, I) = FUGN(I, K)
         end do
      end do
   end if
end subroutine
