!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! !* init() reads the atom positions from file.  If 1 is selected for    *
!* startt then the velocities are assigned, otherwise, they are read   *
!* by selecting 2, or generated by selecting 3                         *
!***********************************************************************

      subroutine init
      include 'MD.com'

	real  pinitmax, TK1, TK2, PK1, PK2, APTtemp, msT,
     Q SigmaT1, SigmaT2, epstemp
	integer storage, dummy,  ANr, IB11, IB12, Ib22, Ib21,
     Q IT1, JT1, KT1, IT2, JT2, KT2, IP1, JP1, KP1, LP1, IP2, JP2, KP2,
     Q LP2, nBA1, nTA1, nPA1, nBA2, nTA2, nPA2,  ind1, ind2, ANt, 
     Q  MDT1, MDT2, cl1, cl2, esaIter, firstAtomIter,secondAtomIter,
     Q  tempAtomsIndex, iter, currentChainIndex,beadInChainIter,
     Q  cnt, isDynFirst, isDynSec, isDyn, numDynAtom, 
     Q  currentRangeContactsNumber,esIndexByBeads,
     Q  currentContactRange


        character(LEN=20) FMTB,FMTT,FMTP,FMTC,CA,RP,FMTE,FMTI,ERP
        character(LEN=10000) currentRangeContacts
        real tempCharge,tempChargeByIndex,esEnergy

	real NNCeps
	dimension NNCeps(NNCmax)
        dimension tempAtomsIndex(maxESAtoms)
        dimension tempCharge(maxESAtoms)
        dimension tempChargeByIndex(maxESAtoms)
	dimension esIndexByBeads(maxESAtoms,maxESAtoms)
	dimension currentContactRange(maxCon)
        
	FMTB="(3I5,2F8.3)"
        FMTT="(4I5,2F8.3)"
        FMTP="(5I5,2F8.3)"
        FMTI="(5I5,2F8.3)"
        FMTC="(5I5,4F8.3)"
        CA="(I5,2I5,F10.3,F9.6)"
        RP="(I8,2I5,2F10.3)"
        ERP="(I8,3I5,2F10.3)"
	FMTE ="(I8,I5,F10.3)"
	
	pinitmax = (2.5*T)**(.5)


      if (startt .eq. 1) then
! this assignes the random momenta to the particles 
! using the desired T as a guide
        open(25, FILE=initval, ACCESS= 'sequential',
     Q  status='unknown')

        read(25, *) 
        read(25, *) 
        read(25, *) 
        read(25, *) ANr
        read(25, *)

! NOA - if there is no static atoms: instead of lastDynamicAtom = ANr (Original file),
!       the array DynamicAtomRange get the values 1,ANr. 
      if (hasStaticAtoms .ne. 'YES') then
	DynamicAtomRange(1) = 1
        DynamicAtomRange(2) = ANr
	DynLength = 2
      endif	
      numDyn = numDynAtom(DynamicAtomRange,DynLength) 

      do 555 i=1, ANr

	read(25, "(I5,I4,A4,A3,4F8.3)") BeadIndex(i), 
     Q  GroupIndex(i),AtType(i),ResID(i), X(i), Y(i), Z(i), ms(i)

        if(BeadIndex(i) .ne. i)then
           write(*,*) 'indexing in input file is wrong'
           call abort
        end if
! this will give random velocities to Protein atoms only
! NOA - change the condition in order this will give random velocities only to Dynamic Protein atoms... 
        cnt = 1
        do while (cnt .le. DynLength)
		if (i .ge. DynamicAtomRange(cnt) .and. i .le. 
     Q              DynamicAtomRange(cnt+1)) then
			call random
			VX(i) = pinitmax*(rand-.5)*2
			call random
			VY(i) = pinitmax*(rand-.5)*2
			call random
			VZ(i) = pinitmax*(rand-.5)*2
			go to 555
		endif
		cnt = cnt + 2
	end do
555   continue
        !!! AMIT: entered this small loop in the 3 lines below to allow using the startt=1 flag
        do i=1, ANr
            read(25,*)
        end do
        
        
         read(25,*) MDT
         if(MDT .gt. Clmax)then
            write(*,*) 'too many chains'
            call abort
         endif
         do i=1,MDT
            read(25,*) ChainLength(i)

         enddo


        close(25,status='keep')

	elseif(startt .eq. 2)then

        open(25, FILE=initval, ACCESS= 'sequential',
     Q  status='unknown')
! If 2 is chosen for startt, then data is extracted from another file
! that was created by this program so there is position on one line and momenta
! on the next with the correct number of digits
        read(25, *) 
        read(25, *) 
        read(25, *) 
	  read(25,*) ANr
        read(25, *) 

! NOA - if there is no static atoms: instead of lastDynamicAtom = ANr (Original file),
!       the array DynamicAtomRange get the values 1,ANr.
      if (hasStaticAtoms .ne. 'YES') then
        DynamicAtomRange(1) = 1
        DynamicAtomRange(2) = ANr
	DynLength = 2
      endif
      numDyn = numDynAtom(DynamicAtomRange,DynLength)

	do 999 i=1, ANr
	read(25, "(I5,I4,A4,A3,4F8.3)") BeadIndex(i), 
     Q  GroupIndex(i),AtType(i),ResID(i),X(i),Y(i),Z(i),ms(i)
        if(BeadIndex(i) .ne. i)then
           write(*,*) 'indexing in input file is wrong'
           call abort
        end if
999   continue

      do i=1, ANr
      read(25, "(I5,3F8.3)") j,VX(i), VY(i), VZ(i)
        if(j .ne. i)then
           write(*,*) 'indexing in input file is wrong'
           call abort
        end if
      end do

         read(25,*) MDT

         if(MDT .gt. Clmax)then
            write(*,*) 'too many chains'
            call abort
         endif
         do i=1,MDT
            read(25,*) ChainLength(i)

         enddo



	close(25, status='KEEP')
        endif

! Read in the contact information


	
! These lines read in the conformations.
        open(30, file=conf, status='old', access='sequential')
          read(30,*) nBA
	if(NBA .gt. NBmax)then
	call abort
	endif

! READ BONDS
        do i=1, nBA
          read(30,FMTB) j, Ib1(i), Ib2(i),Rb(i), bK(i)
          bK(i) = 2*bK(i)
        end do
	
        call BONDSP

! READ ANGLES
          read(30,*) nTA
        if(NtA .gt. Ntmax)then
        write(*,*) 'too many angles'
        call abort
        endif

        do i=1, nTA
          read(30,FMTT) j, IT(i), JT(i), KT(i), ANTC(i), Tk(i)
        TK(i)=TK(i)*2
        enddo


! READ DIHEDRALS
          read(30,*) nPA
        if(NpA .gt. Npmax)then
        write(*,*) 'too many dihedrals'
        call abort
        endif


! this reads in the dihedral angles and calculates the cosines and sines
! in order to make the force and energy calculations easier, later.
        do i=1, npA
           read(30,FMTP) j, IP(i), JP(i), KP(i), LP(i), APTtemp, PK(i)
         
	    DihAng(i) = APTtemp
            GAMS1(i)= PK(i)*Sin(APTtemp)
            GAMS3(i)= PK(i)*Sin(3.0*APTtemp)/2
            GAMC1(i)= PK(i)*Cos(APTtemp)
            GAMC3(i)= PK(i)*Cos(3.0*APTtemp)/2
	!write(*,*) 1,gamc1(1)

        END DO

       if (useChirals .eq. 'YES') then
       ! READ CHIRALS
          read(30,*) nChirals
          if(nChirals .gt. nChiralMax)then
            write(*,*) 'too many chirals'
            call abort
          endif

          do i=1, nChirals
            read(30,FMTC) j, Ichiral(i), Jchiral(i), Kchiral(i), 
     Q      Lchiral(i), chiralValue(i), chiralCoeff(i)
          END DO

       endif

! READ CONTACTS
        read(30,*) NC

          if(NC .gt. MaxCon)then
             write(*,*) 'too many contacts'
             call abort
          endif
        do i=1, NC

          read(30, CA) ind1, IC(i), JC(i), Sigma(i), EpsC(i)
        end do
	

	if (writeContactsRanges .eq. 'YES') then
	  ! read the contact ranges file
	  open(17,file=ContactRangesFile,status='unknown',
     Q    access='sequential')
	  ! read the total number of ranges
	  read(17,*) rangesNumber
	  ! for each range
	  do i=1, rangesNumber
	  ! read the number of contacts within said range
	    read(17,*) currentRangeContactsNumber
	    rangeContactsNumber(i) = currentRangeContactsNumber
	    ! read the contact indices in string format
	    read(17,*) currentRangeContacts
            ! use Noa's code to split the string into an array
	    call genDynRange(currentRangeContactsNumber,
     Q                  currentRangeContacts,currentContactRange)
	    do j=1,currentRangeContactsNumber
	      rangeContacts(i,j) = currentContactRange(j)
	    end do
	  end do
	  close(17)
	end if

! Find all of the possible Three-body interactions i.e. A-B contact, 
! B-C contact and C-A contact
        !call ThreeBodyInit

! READ REPULSION ! read non-native interactions
        read(30,*) NNC
        if(NNC .gt. NNCmax)then
        write(*,*) 'too many non contacts: NNC ',NNC,' NNCmax ',NNCmax
        call abort
        endif        
        do i=1, NNC
           read(30,RP) ind1, INC(i), JNC(i), NCsigma(i), NNCeps(i)
! this simplifies calculations later
           NNCsigma(i) = 12*NNCEps(i)*NCsigma(i)**6
        end do

! READ Ellipsoid repulsions
        if (useEllipsoidRepulsions .eq. 'YES') then
          read(30,*) ellipsoidRepulsionsNum
          if(ellipsoidRepulsionsNum .gt. ellipsoidRepulsionsMax)then
        write(*,*) 'too many ellipsoid repulsions: '
     Q             ,ellipsoidRepulsionsNum,' >= ',ellipsoidRepulsionsMax
        call abort
        endif        
        do i=1, ellipsoidRepulsionsNum
           read(30,ERP) ind1, IEllipsoid(i), JEllipsoid(i), 
     Q                  KEllipsoid(i), ellipsoidSigma(i),
     Q                  ellipsoidCoeff(i)
        end do
	endif

!READ ELECTROSTATICS
        if (useElectrostatics .eq. 'YES') then
	  read(30,*) esAtomsNum
	  do esaIter=1, esAtomsNum
	    read(30,FMTE) ind1, tempAtomsIndex(esaIter),
     Q    tempCharge(esaIter)
	  end do
          
        endif



        read(30,*) AN

! NOA - if there is no static atoms: instead of lastDynamicAtom = AN (Original file),
!       the array DynamicAtomRange get the values 1,AN.
        if (hasStaticAtoms .ne. 'YES') then
        	DynamicAtomRange(1) = 1
        	DynamicAtomRange(2) = AN
		DynLength = 2
        endif
        numDyn = numDynAtom(DynamicAtomRange,DynLength)

! These are safeguard statements
        if(startt .eq. 1 .or. startt .eq. 2)then
        if (ANr .ne. AN)then
           write(*,*) 'Number of important atoms in files do not align.  
     Q Stop Program!!', AN, 'is the number in conf1'
           call abort
        endif
        endif


         read(30,*) MDT1
         if(MDT1 .gt. Clmax)then
            write(*,*) 'too many chains: MDT1 ',MDT1 ,' Clmax ',Clmax
     Q,' AN ',AN
            call abort
         endif


	if(startt .eq. 1 .or. startt .eq. 2)then
	if(MDT1 .ne. MDT)then
            write(*,*) 'different numbers of chains, abort'
            call abort
         endif
	endif


         MDT = MDT1
        if(Trajectory .ne. 'NO')then
          if(minTrajOut .ne. 'YES')then
                 write(6,*) MDT
          endif
	endif

         do i=1,MDT
            read(30,*) cl1 


        if(startt .eq. 1 .or. startt .eq. 2)then
	if(cl1 .ne. Chainlength(i))then
	write(*,*) 'chain lengths differ, abort'
	call abort
	endif
	endif

             ChainLength(i) = cl1

        if(Trajectory .ne. 'NO')then
                if(minTrajOut .ne. 'YES') then
                  write(6,*) chainlength(i), WO, WOT
                else
                  write(6,*) WOT
                endif
	endif
         enddo

		currentChainIndex = 1
	beadInChainIter = 0
        do i=1, AN
	   
           read(30,"(I5,I4,A4,A3,4F8.3)") BeadIndex(i), 
     Q  GroupIndex(i),AtType(i), ResID(i), xt(i), yt(i), zt(i),
     Q  ms(i)
	beadInChainIter = beadInChainIter + 1
	if (beadInChainIter .gt. ChainLength(currentChainIndex))then
	  beadInChainIter = 1
	  currentChainIndex = currentChainIndex + 1
        end if
	ChainIndex(i) = currentChainIndex
	
        if(Trajectory .ne. 'NO')then
                if(minTrajOut .ne. 'YES')then
                  write(6,'(I4,A4,A3)') BeadIndex(i),AtType(i), ResID(i)
                endif
	endif
        enddo

! If startt is 3, then conformation 1 is used as the starting position
! and velocities are assigned.

      if(startt .eq. 3)then

      do 10 i=1, AN
      X(i)=xt(i)
      Y(i)=yt(i)
      Z(i)=zt(i)
! this will give random velocities to Dynamic Protein atoms only
! NOA - change the condition in order this will give random velocities only to Dynamic Protein atoms...
	cnt = 1
	do while (cnt .le. DynLength) 
		if (i .ge. DynamicAtomRange(cnt) .AND. i .le. 
     Q              DynamicAtomRange(cnt+1)) then
			call random
			VX(i) = pinitmax*(rand-.5)*2
			call random
			VY(i) = pinitmax*(rand-.5)*2
			call random
			VZ(i) = pinitmax*(rand-.5)*2
			go to 10
		endif
		cnt = cnt + 2
	end do
 10   continue

      endif

        dummy = 0
        do i=1, MDT
           dummy = dummy + ChainLength(i)
        enddo

        if(dummy .ne. AN)then
           write(*,*) 'Chain lengths and number of particles dont align'
           call abort
        endif

! initialize electrostatic interactions arrays

        if (useElectrostatics .eq. 'YES') then
	  if (esAtomsNum .gt. AN) then
             write(*,*) 'too many electrostatic atoms: ',
     Q       esAtomsNum,' > ',AN
             call abort
	  endif
          
          do firstAtomIter=1,AN
            tempChargeByIndex(firstAtomIter) = 0
          end do

          do firstAtomIter=1,esAtomsNum
            tempChargeByIndex(tempAtomsIndex(firstAtomIter)) =
     Q      tempCharge(firstAtomIter)
          end do
          
          esaIter = 0
          do firstAtomIter=1,esAtomsNum
            do secondAtomIter=firstAtomIter+1,esAtomsNum
!   if the first atom is static the interaction is 
!    between 2 static atoms and is therfore neglected
! NOA - if both atom static the interaction is between 2 static atoms and is therfore neglected
	      isDynFirst = 
     Q isDyn(DynamicAtomRange,DynLength,tempAtomsIndex(firstAtomIter))	
	      isDynSec = 
     Q isDyn(DynamicAtomRange,DynLength,tempAtomsIndex(secondAtomIter))
              if(isDynFirst == 1 .or. isDynSec == 1)then
		i = tempAtomsIndex(firstAtomIter)
		j = tempAtomsIndex(secondAtomIter)
		if((abs(j-i) .ge. esMinBeadDistance) .or.
     Q             (ChainIndex(i) .ne. ChainIndex(j))) then
  		  esaIter = esaIter+1
		  esFirstAtomIndex(esaIter) = i
		  esSecondAtomIndex(esaIter)= j
		  esCharge(esaIter) =
     Q        tempCharge(firstAtomIter)*tempCharge(secondAtomIter)
		 esIndexByBeads(i,j) =  esaIter
                endif
              endif
            end do
          end do
          esPairsNum = esaIter

	  if (compensateElectrostaticContacts .eq. 'NO_ES') then
            do i=1, NC
              if(tempChargeByIndex(IC(i))*tempChargeByIndex(JC(i))
     Q           .ne. 0) then
		 esCharge(esIndexByBeads(IC(i),JC(i))) = 0
              endif
            end do
	  end if

	  if ((compensateElectrostaticContacts .eq. 'ALL') 
     Q    .or.(compensateElectrostaticContacts .eq. 'ATTRACTION'))
     Q    then
            do i=1, NC
              if(tempChargeByIndex(IC(i))*tempChargeByIndex(JC(i))
     Q           .lt. 0) then
                if (useDebyeHuckel .eq. 'YES') then

                  call debyehuckelfactor(sigma(i),deConstant,
     Q            screeningFactor,saltCoefficient, esEnergy)
	        else
                  call coulombfactor(sigma(i),deConstant,esEnergy)
	        endif
                EpsC(i) = EpsC(i)-esEnergy
                if (EpsC(i) .lt. 0.01) then
                  EpsC(i) = 0.01
                endif
              endif
            end do
	  end if
	  if (compensateElectrostaticContacts .eq. 'ALL') then
          
            do i=1, NC
              if(tempChargeByIndex(IC(i))*tempChargeByIndex(JC(i))
     Q           .gt. 0) then
                if (useDebyeHuckel .eq. 'YES') then

                  call debyehuckelfactor(sigma(i),deConstant,
     Q            screeningFactor,saltCoefficient, esEnergy)
	        else
                  call coulombfactor(sigma(i),deConstant,esEnergy)
	        endif
                EpsC(i) = EpsC(i)+esEnergy
              endif
            end do
	  end if

        endif

       end
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^end of init^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
