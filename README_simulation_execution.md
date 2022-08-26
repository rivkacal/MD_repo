 # Description 
This file describes all the calculations performed in each file seperately, strating from running **MD.exe** file as performed in .csh file created at MDWrapper.pl

# fortran shortacuts

- 6 is Trajectory file | appears in ss.f,init.f
- 8 ??| appears in int.f (closed)
- 10 ??| appears in int.f (closed)
- 11 ??| appears in int.f (closed), ThreeBody.f (writes)
- 12 is ContactFile  | appears in ss.f,ThreeBody.f (writes)
- 13 is ThreeBodyFile |	appears in ss.f,ThreeBody.f (writes)
- 14 is writeAllContacts (optional) | appears in ss.f,
- 15 is ContactRangesTwoBodyFile (optional if above is set 'YES') | appears in ss.f,IjM.f
- 16 is ContactRangesEnergyFile | appears in ss.f, IjM.f
- 17 is ContactRangesFile | appears in init.f (open, close)
- 25 is initval file (dynamic range) | appears in init.f (open,close),
- 30 is conf (input.dat) | appears in init.f
- 53 is TemperatureVtime | appears in ss.f,
- 55 is output file | appears in write.f (open, close)
- 58 is (?? output each step and avg temperature here) | appears in int.f
- 60 is settings.dat (mandatory) appears in ss.f,
- 61 is finalpx1 or 'temp' (output in case simulation crashes) | appearing in write.f (open,close)
- 63 is 'random.dat' file | appears in ss.f (open, close),
- 70 is EnergyTerm | appears in ss.f, LD.f , sym.f (writes)
- 71 ???? holds temperature? | appears in ss.f (closed),
- 91 is 'time.dat' (time for all simulation) | appears in int.f
- 93 is TrajDist file |	appears in ss.f,int.f
- 98 is ?? dihedrals staff | appears in phi.f
- 99 is EnergyTot | appears in ss.f,int.f (integrate)
- 999 is 'timeleft.dat' (remaining for aimulation after first integartion) | appears in int.f (opens, closed)


# MD.com
Defines common parameters and variables to be used in different functions. 
From physical parameters like settigs for gamma and tau to their formatting (fortran). Format for X,Y,Z coordinates , velocities and forces.
please note:
     
      ! This limits the size of the simulation.
             parameter (Nmax=20000)
             dimension X(Nmax), Y(Nmax), Z(Nmax), ms(Nmax),
           Q Vx(Nmax), Vy(Nmax), Vz(Nmax),
           Q Fx(Nmax), Fy(Nmax), Fz(Nmax)


When changing outputs note:

    ! file name variables
      character(LEN=200) Trajectory, EnergyTot,EnergyTerm,ContactFile,
         Q TemperatureVtime, ThreeBodyFile, TrajDist, allContactsFile, 
         Q ContactRangesFile, ContactRangesTwoBodyFile, 
         Q ContactRangesEnergyFile

    ! variables for contact ranges option
      integer rangesNumber,rangeContactsNumber,rangeContacts, maxRanges,
         Q          TwoBodyRanges
      real    EtotalRanges
            parameter (maxRanges = 50)
      dimension rangeContactsNumber(maxCon)
      dimension rangeContacts(maxRanges,maxCon)
      dimension TwoBodyRanges(maxRanges)
      dimension EtotalRanges(maxRanges)

    ! addtional parameters, used only if conditional flag is set 

          integer DynamicAtomRange(1000)
    integer DynLength
    integer numDyn
    integer useESCutoff, useDHTable
          real boxMin,boxMax, boxCoeff
    dimension boxMin(3),boxMax(3)
          real deConstant, screeningFactor, saltCoefficient

# MDMake
contains all syntax to mamke any file of fortarn/c executable! 

# main.f
PROGRAM DualGo : Core Program by Paul Whitford <pcw@wpi.edu> April 2005 ;Go.f    
indludes MD.com module and defines global: timee1, timee2, sum, sumc, DD1, DD2 
then calls: ***start***(ss.f), **init**(init.f), **integrate**(int.f) and **stophere** (ss.f)

# ss.f
Reads 'settings.dat' and 'random.dat' keys for new varaibles to use along the simulations. Checks if cutoff values are valid (range).
Contains the following functions:

 - start:
  reads settings from settings.dat (60) created perviously, random variables and creates outputs:
  
  
	      skips the first two meaningless lines and then:
	      reads initial conditions format (1||2||3) into startt (integer type)
	      reads the input file (here: fragB_helix_translated_updated.dat) into Conf
	      reads initial coordinates/velocities from that line (leftest) input file (fragB_helix_translated_updated.dat) into initval,
		  simultaneously Final_fragB_helix_translated_updated_0.5_39753.dat is read into finalpx1 (final conformation)
	      reads outputfile name (say fragB_helix_translated_updated_0.5_39753.log) into output variable
	      **reads trajectory file name (say Traj_fragB_helix_translated_updated_0.5_39753.dat) into Trajectory variable,
	      and the output frequency into WOT (write 1 step every WOT (say 1000 writing output trajectory) steps (that's 1 every 1000 skipping 999 steps )**
	      if choosing to store Trajectorries : open a Trajectory file and call it 6.
	      if choose to store trajectory distance : open TrajDist file and call it 93.
	      keeps reading energy total file name (Etotal_fragB_helix_translated_updated_0.5_39753.dat) into EnergyTot,
		EbyType_fragB_helix_translated_updated_0.5_39753.dat file name into EnergyTerm. If they have values other than 'NO' string then:
	      open EnergyTot file (name) as 99, open EnergyTerm as 70. Similarly reading 2body and 3body contact filenames into ContactFile (12), 
		and ThreeBodyFile (13)
	      **reads writeAllContacts, only if not specified 'NO' inMDWrapper.prefs it'll open a writeAllContacts as 14**
	      **then it reads the ranges to write 2body and EIj if above is specifies as 'YES', into writeContactsRanges
	      open ContactRangesTwoBodyFile as 15, ContactRangesEnergyFile into 16**
	      reads Temperature file name and open TemperatureVtime as 53. 
	      reads model type of execution into symtype (MD/LD) if LD then read also gamma (drag coefficient for bath coupling)
	      reads the number of steps into stepstop and the measurements frequency into WO (writing output). 
	      reads time step into tau (0.005, 200 steps is?), read coupling constant for Berendsen thermostat into RsTau (0.20000 rescaled by tau)
	      reads temperature in reduced units (0-1 to fraction of "room temperature") into T.
	      reads wether static atoms option was chosen; if so, call genDynRange(DynLength,DynamicStr,DynamicAtomRange) for dynamic range (???)
	      reads wether to confine in a box: if 'YES' write values of boxMin(i),boxMax(i) where i={1,2,3}, and the box force coefficent into boxCoeff.
	      reads wether to apply columbic interactions into useElectrostatics, if 'YES' read the dielectric constant (eps_0 units) into deConstant;
		   minimum index difference between beads to apply electrostatics (say 3 indices apart) into esMinBeadDistance; 
		   type of cutoff for electrostatic interaction into esCutoffType. If the latter is 'ENERGY'
		   then read the maximum fraction of energy for electrostatic interaction into esEnergyCutoff, similar for 'DISTANCE' type.
	      reads whether to fix potential of electrostaic contacts ('ATTRACTION') into compensateElectrostaticContacts. 
	      reads wether to use debye huckel in execution into useDebyeHuckel, and then wether to use precalculated energy and force values into useDHTable.
	      reads into ionicStrength, ionicRadius (the average radius of the ions in the solution (for example 1.4 for NaCl)), solventDensity (1 for water)
	      after that calling initES to initiate electrostatics, as described later.
	      reads  whether to use chirals in execution into useChirals, whether to use ellipsoid repulsions in execution into useEllipsoidRepulsions.
      
      opens 'random.dat' file as 63 xrandom,yrandom,zrandom from it, then closes it. 
      
      checks if symtype is 'LD' and if True calculates: c_e = 1.0 - gamma*tau/2.0
	                                                      c_i = 1.0/(1+gamma*tau/2.0)

  - stophere subroutine should close all open files (say 60 settings.dat) (??? unclear how). It says to Close the files that holds the temperature (71)          no closing of 60 ... reads again again xrandom,yrandom,zrandom...
 
 - initES()
    defines parameter (screeningFactorConstant = 2529.11892)
    if using Debye Huckel supposedly transforming simulation temperature to real temperatures with realTemperature = 295*T.
        then calculae kappa :
	
		screeningFactor = SQRT((screeningFactorConstant*ionicStrength*solventDensity)/(deConstant*realTemperature))
        	saltCoefficient = (exp(screeningFactor*ionicRadius)/(1+screeningFactor*ionicRadius)) 
		
    if cutoff is of 'ENERGY' type then define minimal distnace to apply electrostatics as:
    	esMinDistance = 16.0, (Angstorm), esMaxDistance = 100000000.0
        then, if use Debyehuckel call: debyehuckelfactor(esMinDistance,deConstant,screeningFactor,saltCoefficient,esEnergy)
        otherwise call coulombfactor(esMinDistance,deConstant,esEnergy)      
    (???)
    
        then calculates: esCutoff = esEnergy * esEnergyCutoff (the latter is the fraction defined above),esEnergyError = esCutoff / 100.0
        then checks 5000 times: if the energycutoff minus the esEnergy absolute value is lower than the error,
	(not zero exactly) and ouputs error message if does.,
         starting from average value, each iteration changing the maximal or lower distance according to ratio of esEnergy
	 (above and below) to cover all options
    if cutoff is of 'DISTANCE' type just check if esDistanceCutoff **squared** (** 2) is greater than maximal cutoff to abort, If it's false just
      take the square value of input to continue the calculations with the new esDistanceCutoff .
      
    
 # init.f
 reads the atom positions from conformation file: If 1 is selected for startt then the velocities are assigned,
 otherwise, they are read by selecting 2, or generated by selecting 3 (note MDWrapper.prefs file)
 Also it reads all other inputs in the input.dat file : BONDS, ANGLES, DIHEDRALS, CONTACTS, REPULSIONS, number of chains and total atoms.
 Differentiating dynamic and static atoms for simulation, starting from writing random velocities to dynamic only, and then creating
 a electrosttic pairs table esIndexByBeads(i,j) for later calculation of contacts compensation. 
 **Note writing to trajectory (6) file (number of chains and total atoms) starts here**
 
 - subroutine init
 	defines formatting
	assign pinitmax = (2.5*T)** (.5)
       
       if startt equals 1:
		the following assignes the random momenta to the particles, using the desired T as a guide
		open initval file as 25, read the values in 4th line into ANr (atomic number range??, generally lastDynamicAtom=ANr)
		if there are no static atoms then the range is 1,ANr (1-ANr including) that is DynamicAtomRange is now
		array of length 2. Otherwise, *dynamic* ranges are written in the same format:
		range1start,range1end,range2start,range2end, etc... (non listed are static)
		numDyn = numDynAtom(DynamicAtomRange,DynLength)
		
		now for every dynamic index, read from initval into:
		
			read(25, "(I5,I4,A4,A3,4F8.3)") BeadIndex(i), 
   			  Q  GroupIndex(i),AtType(i),ResID(i), X(i), Y(i), Z(i), ms(i)
			  
		now assign random velocities to protein only: dynamic atoms only
		
			call random
			VX(i) = pinitmax*(rand-.5)*2
			call random
			VY(i) = pinitmax*(rand-.5)*2
			call random
			VZ(i) = pinitmax*(rand-.5)*2
		then read the number of chains into MDT and check if it is larger then Clmax (Chain length max??).
		Every chain i length is read into ChainLength(i) before closing the file.
		
	 if startt equals 2: then data is extracted from another file that was created by this program so there is
	 position on one line and momenta on the next with the correct number of digits
	 	again open 25 read into ANr (code duplication!!) residue index, X,Y,Z,ms but **then**:
		
			 do i=1, ANr
      			 read(25, "(I5,3F8.3)") j,VX(i), VY(i), VZ(i)
        			if(j .ne. i)then
          			   write(*,*) 'indexing in input file is wrong'
		
		again read chains lengths...
		
	read conformation info;	open conf file (input.dat) as 30 to read: (say here fragB has 161 a.a in sequqnce)
	!BONDS SECTION 
	read number of bonds of chain A into nBA (160 or length-1)
	Then read until the last bond section in format:
	
		read(30,FMTB) j, Ib1(i), Ib2(i),Rb(i), bK(i) ! serial number, res1_idx, res2_idx, pair_dist, bondConst (k^b_i_j)
         	 bK(i) = 2*bK(i)
	
	call BONDSP (???)
	!ANGLES SECTION
	similarly read number of bond angle trios (length-2 or 159 here) into nAT then using:
		
		read(30,FMTT) j, IT(i), JT(i), KT(i), ANTC(i), Tk(i) ! serial number, trioI_res_idx, trioJ_res_idx, trioK_res_idx, triad_angle, 											angle_force_constant (10 for IDPS whereas 20 for ordered)
        	TK(i)=TK(i)*2
		
	!DIHEDRALS SECTION
	read amount into nPA. Then the following reads in the dihedral angles and calculates the cosines and sines
	in order to make the force and energy calculations easier, later (???):
	
		do i=1, npA
		   read(30,FMTP) j, IP(i), JP(i), KP(i), LP(i), APTtemp, PK(i)
		    DihAng(i) = APTtemp
		    GAMS1(i)= PK(i)*Sin(APTtemp)
		    GAMS3(i)= PK(i)*Sin(3.0*APTtemp)/2
		    GAMC1(i)= PK(i)*Cos(APTtemp)
		    GAMC3(i)= PK(i)*Cos(3.0*APTtemp)/2
		!write(*,*) 1,gamc1(1)
		
	if using chirals then reads nChirals...
	!CONTACTS SECTION
	read number of contacts into NC. then for every contact:
		
		read(30, CA) ind1, IC(i), JC(i), Sigma(i), EpsC(i) ! serial_indx, res_i, res_j, sigma_val (dist square), epsilon_val 
	
	**If writeContactsRanges flag is 'YES' then also read ContactRangesFile as 17**
	read rangesNumber and for each of those ranges read current number of contact currentRangeContactsNumber,to split the strings into arrays:
	call genDynRange(currentRangeContactsNumber,
                       currentRangeContacts,currentContactRange)
	
	
	!REPULSION SECTION (non-native interactions) (only repulsive term of LJ 1/r_i_j^12 excluding 1/r_i_j^6)
	read amount of repulsions into NNC. Note the total number of contacts+repulsions=(N-3)(N-4)/2 (say N=161, then total=12403)
	
		read(30,RP) ind1, INC(i), JNC(i), NCsigma(i), NNCeps(i)
		! this simplifies calculations later:
        	   NNCsigma(i) = 12*NNCEps(i)*NCsigma(i)**6 ! recall sigma here is already the distance squared!!! so just power of 6 left
	Then read ellipsoid if defined in Warpper.prefs...
	!ELECTROSTATICS 
	number of electrostatic residues into esAtomsNum, by:
	
		read(30,FMTE) ind1, tempAtomsIndex(esaIter), tempCharge(esaIter)  ! serial_idx, chainA_res, charge
		
	!ATOM POSITIONS SECTION
	read number of atoms in chains into AN. Then code duplication to check if ranges fit AN 318-335)...
	read number of chains into MDT1
	**updating to trajectory first line number of chains, then length of chain, and then all atoms positions and type**:
		if(Trajectory .ne. 'NO')then
		 write(6,*) MDT
		endif
		
		! read chain length into cl1 and then (ChainLength(i) = cl1)
		 write(6,*) chainlength(i)
	
		! initialize currentChainIndex = 1, beadInChainIter = 0 then read all bids and change index
		read(30,"(I5,I4,A4,A3,4F8.3)") BeadIndex(i),GroupIndex(i),AtType(i),ResID(i),xt(i),yt(i),zt(i),ms(i) 
		! notice group chain and mass at the end.
		
		!Then writing all identities and type (CA/CB)
		 write(6,'(I4,A4,A3)') BeadIndex(i),AtType(i), ResID(i)
		
		
	If startt is 3, then conformation 1 is used as the starting position
	and velocities are assigned use the groups and dynamic ranges to assign velocities to dynamic atoms...
	
	elctrostatics: initializing array for all atoms of zeroes electrostatics - tempChargeByIndex. Assigining to the first 1:esAtomsNum tempCharge.
		Now iterating over pairs: first is (firstAtomIter=1:esAtomsNum), second (secondAtomIter=firstAtomIter+1+1:esAtomsNum).
		Note: if the first atom is static the interaction is between 2 static atoms and is therfore neglected
		      if both atom static the interaction is between 2 static atoms and is therfore neglected.
		      So check first if both atoms are dynamic. If so check if the indices of the pairs are greater (in abs) than the 
		      minimal indices distance to calculate electrostatics: esMinBeadDistance (or if they are on different chains)
		      run over all pairs, assign the charge multiplication to a 2D array.
		      esCharge(esaIter) = tempCharge(firstAtomIter)*tempCharge(secondAtomIter)
		      esIndexByBeads(i,j) =  esaIter
		      ! when finished last pair:
		      esPairsNum = esaIter
		
		**compensateElectrostaticContacts** can recieve different inputs:
		1. 'NO_ES' : if so run over all contacts indices and check if they have any electorstatic value other then 0, if so
			assign 0 to this specific contact pair. (if they have contact they do not iinteract electrostaticlly here)
		2. 'ATTRACTION' : if so check for each pair if charge multiplicatin is true (non zero), then calculate eiter
		    debyehuckelfactor(sigma(i),deConstant,screeningFactor,saltCoefficient, esEnergy) or columbic according to your flags.
		    then:  EpsC(i) = EpsC(i)-esEnergy   ! substruct the electrostatics for that contact (sign??)
                	  if (EpsC(i) .lt. 0.01) then
                 	      EpsC(i) = 0.01 
		3. 'ALL' : do the same for 'ATTRACTION' only then: EpsC(i) = EpsC(i)+esEnergy  this will return to the contact epsilon 
			other then the case  EpsC(i) = 0.01 which is now added the electrostatics value  

		
# int.f
calls initializing LD/MD simulation, calls writetofile(file) and or writes to trajectories, total energy, distances (optional), 

contains subroutine integrate:
	defines variables like Res the array the length of AN and count the length of number of contacts (NC)
	calls CPU_Time(t1) (??)
	**writes the number of atoms to trajectory : write(6,*) AN (again??)**
        calls noncontactsnear to find the non native non bonded list
	
	run over all steps with the measurement writing frequency writeout = WO, calaulate dummy: Edum= writeout - step.
	if the step is dividble by 50, call noncontactsnear (??)
	
	calling 1 step integration with:
	- 'MD' simulation : call intSymp(E ,ET, Conts, Edum, count)
	- 'LD' : call intLD(E,ET,Conts,Edum,count)
	call findtemp(dummyy) (??)
	
	if(Edum .eq. 0)then
**! Find and return the Three-body interactions
	  !call ThreeBody(conts, count)** ! note this is stopping the calling for 3body calculation to not slow down to much.
	  should insert a flag: USE_3BODY
	  in MDWrapper.prefs and update it here
	  
	 now if choosing to write EnergyTot:
	 writing the energy E to this file as 99: write(99,*) E
	 outputing the temperature to 58 file:  write(58,*) step, temp
	 summing the temperature to Tsum, and incrementing Tcount to the next one (temperaure runs??). finally writeout = writeout + WO
	 if trajectory writeoutT matches with current step, 
	 call function writetofile(writei) where writei is flag (1 or 0) for the file number to write.
	 
	 if choose to write TrajDist (93):
	 
		    if(mod(step,WOT) .eq. 0)then
		      write(93,*) step
		      call distances()
		    endif
		    
**writing to trajectories of continue and end**

		if(stepstop-step .ge. WOT)then

		    write(6,*) 'continue'

		    else
		    write(6,*) 'end'
	when finished writing one measurement (WO=step), call CPU_time(t2). calculate the time it took for this analysis and then extrapulate
	the remaining CPU time for the simulation in 'timeleft.dat' as 999.
	After the last step calculations: write the average temperatire Tsum/Tcount to 58. Calculate time took for all the simulation,
		write it to 91 'time.dat'. closing file 8. calling writetofile(writei) !, Res). Then close 10 and 11 (??)
		
# fs.f		
FStart initiallizes the Force arrays. includes fstart, that for every atom k assigns: Fx(k)=0.0, Fy(k)=0.0, Fz(k)=0.0		

# LD.f
contains intLD() : integrates over one time step for all particles: intLD(E,ET, Conts, outE, count )
where E is the total potential energy, ET is the total potential and kinetic energy.
initialize DumE(i)=0.0 where i=1:9. then calls the following routines to calculate forces:
call Fstart, call Bonds( DumE(1)), call ANGL(DumE(2)),call dihedral(DumE(3)),call contacts(DumE(4), Conts, count, outE), call NonContacts( DumE(5)),
optional:
call box(X,Y,Z,Fx,Fy,Fz ,DynamicAtomRange, DynLength, DumE(6),boxMin,boxMax,boxCoeff), call chiral(DumE(8)), 
call ellipsoidRepulsions(DumE(9)).
If choose to use elctrostatics and debye-juckel: calculate electrostatics with:
call debyehuckel(useESCutoff,useDHTable,esFirstAtomIndex,
                 esSecondAtomIndex,esCharge,X,Y,Z,Fx,Fy,Fz,esPairsNum,
                 esDistanceCutoff,DumE(7),deConstant,screeningFactor, 
                 saltCoefficient,DebyeHuckelPotentials,DebyeHuckelForces)
Then: Call Findtemp(DumKE)
       ET = E + DumKE (ET again is potenial + kinetic energy that comes from temperature)	
 if contact energy is zero then the potential energy E =  DumE(1) + DumE(2) + DumE(3) +  DumE(4) + 
                                                          DumE(5) + DumE(6) + DumE(7) + DumE(8) + DumE(9)    
if choose to write energy term:

		if(EnergyTerm .ne. 'NO')then
		      write(70,"(9F9.2)") DumE(1),DumE(2),DumE(3),DumE(4),DumE(5),
		     Q DumE(6),DumE(7),DumE(8), DumE(9)
Then for every dynamic atom:
		
		SD = sqrt(2*ms(i)*gamma*T/tau)
		call gauss(0., SD, fr)
		Vx(i) = (Vx(i)*c_e + (Fx(i)+fr)*tau/ms(i))*c_i
		call gauss(0., SD, fr)
		Vy(i) = (Vy(i)*c_e + (Fy(i)+fr)*tau/ms(i))*c_i
		call gauss(0., SD, fr)
		Vz(i) = (Vz(i)*c_e + (Fz(i)+fr)*tau/ms(i))*c_i
	
		! then for every atom:
		X(i) = X(i) + Vx(i)*tau
	        Y(i) = Y(i) + Vy(i)*tau
	        Z(i) = Z(i) + Vz(i)*tau


# write.f
writetofile writes the positions to file every savenum. It rewrites the positions in the same file, deleting the old position, thus not taking too much space. The output is written to two alternating files so if the program crashes while writing, the data will not be lost . Modified 10/6 to write failures also 

(Note if start is 1 or 2 the filename is initval (25). Otherwise it's conf. This is used for output later)

if filei=0 (input): open finalpx1 (final positionx1) as 61, if filei=1 open 'temp' as 61 . calls findtemp(KE) updating the number of steps calculated
	untill now, current temperature and number of atoms. Then the current coordinates.
	
**updating trajectory with current step writecount**

	        if(Trajectory .ne. 'NO')then
		write(6,*) writecount
		endif
then write all step information to 61 BeadIndex(i), GroupIndex(i),AtType(i),ResID(i), X(i), Y(i), Z(i), ms(i), later including all velocities, followed by writing the number of chains and each chains's length. 
Updating (appending to) trajectory (6) with current coordinates list .

open an output file as 55 to write the temperatue, number of steps, conformation (if from conf or initval), tau, number of contacts (NC) (2body) and possible 3body (tripi).

# ran.f
subroutine random: generator of variate on 0, 1  Purpose: generates a pseudorandom number RANDwhose distribution is flat on the interval (0,1).Uses as seeds xrandom, yrandom, zrandom stored between runs in file RANDOM. DAT

subroutine gauss(MM,SD,RR): generates a random number from a gaussian distribution with mean MM, and standard deviation SD . returns RR : 
tmp2 = sd*sqrt( -2.d0*log(tmp1)/tmp1 )
        RR =        zeta(1)*tmp2 + MM 
where zeta(i) = 2.0*(rand-0.5), tmp1 = zeta(1)*zeta(1) + zeta(2)*zeta(2). and 0<tmp1<1

# T.f  
contains findtemp(KE) that finds the temperature of the system ,contains lines which remove angular and center of mass moment
to avoid that the system rotates or moves. returns KE.

first calculate momentum for dynamic atoms:

			Px(i) = Vx(i)*ms(i)
			Py(i) = Vy(i)*ms(i)
			Pz(i) = Vz(i)*ms(i)

Avoid reseting the center of mass to zero if there are static atoms. If not the total linear momentum:
pxsum=0.0, pysum=0.0, pzsum=0.0,  Mass = 0.0
for every dynamic atom i : pxsum = Px(i) + pxsum,pysum = Py(i) + pysum,pzsum = Pz(i) + pzsum,Mass = Mass + ms(i)
Then find the average linear momentum: pxsum = pxsum/Mass,pysum = pysum/Mass,pzsum = pzsum/Mass
Then bring the center of mass linear momentum to zero:
Px(i) = Px(i) - pxsum*ms(i),  Py(i) = Py(i) - pysum*ms(i), Pz(i) = Pz(i) - pzsum*ms(i), 
VX(i) = PX(i)/ms(i)        ,  VY(i) = PY(i)/ms(i)        , VZ(i) = PZ(i)/ms(i)
Following that, to remove the angular momentum, first center the molecule about the origin:
LX = 0.0,  LY = 0.0,  LZ = 0.0,        IXY = 0.0,  IYZ = 0.0,  IXZ = 0.0,           CMX = 0.0,  CMY = 0.0,  CMZ = 0.0,

for avery dynamic i: CMX = CMX + X(i)*ms(i), CMY = CMY + Y(i)*ms(i), CMZ = CMZ + Z(i)*ms(i)
after loop ends remember to divide:         CMX = CMX/Mass,  CMY = CMY/Mass,   CMZ = CMZ/Mass

For every dyamic i: X(i) = X(i) - CMX , Y(i) = Y(i) - CMY,  Z(i) = Z(i) - CMZ
Calculate angular momentum: 
LX = LX + Y(i)*PZ(i)-Z(i)*PY(i),        LY = LY + PX(i)*Z(i)-X(i)*PZ(i),        LZ = LZ + X(i)*PY(i)-Y(i)*PX(i)
IXY = IXY + ms(i)*(X(i)** 2 + Y(i) ** 2), 	IXZ = IXZ + ms(i)*(X(i)** 2 + Z(i)** 2),  IYZ = IYZ + ms(i)*(Y(i)** 2 + Z(i)** 2)
Then calculate: WX = LX/IYZ, WY = LY/IXZ, WZ = LZ/IXY
finally for every dynamic i substract L:
Vx(i) = VX(i) + WZ*Y(i) - WY*Z(i),  VY(i) = VY(i) - X(i)*WZ + WX*Z(i),	VZ(i) = VZ(i) + X(i)*WY - WX*Y(i),
PX(i) = VX(i)*ms(i),                PY(i) = VY(i)*ms(i),                PZ(i) = VZ(i)*ms(i)

to sum square momentums firsr initialize: P2=0
run over dynaic atoms i: P2 =  P2 + (Px(i)** 2 + Py(i)** 2 + Pz(i)** 2)/ms(i)
Now the temp of the system in reduced units:
temp = P2/((numDyn-2)*3) ! numDyn = number of dynamic atoms (or i iterations)
returning the kinetic energy(0.5p^2): KE = P2/2.0


# ThreeBody.f
contains subroutine ThreeBodyInit determines all of the possible three-body interactions in the system  
?? writing
subroutine ThreeBody(Qn,Qi) computes the number of three body interactions in the system


# angl.f
 ANGL  computes the Force due to the bond angles, returns the energy related to the bonds. This code is taken from AMBER, and modified :
 SUBROUTINE ANGL(E):
C for multiple line comment
C  ROUTINE TO GET THE ANGLE ENERGIES AND FORCES FOR THE
C           POTENTIAL OF THE TYPE CT*(T-T0)** 2
First get angles (insertion):

		DO JN = 1, nTA
		    I3 = IT(JN)
		    J3 = JT(JN)
		    K3 = KT(JN)
		    XIJ(JN) = X(I3)-X(J3)
		    YIJ(JN) = Y(I3)-Y(J3)
		    ZIJ(JN) = Z(I3)-Z(J3)
		    XKJ(JN) = X(K3)-X(J3)
		    YKJ(JN) = Y(K3)-Y(J3)
		    ZKJ(JN) = Z(K3)-Z(J3)
		END DO
Then for every JN:
	    RIJ0 = XIJ(JN)*XIJ(JN)+YIJ(JN)*YIJ(JN)+ZIJ(JN)*ZIJ(JN) ! sum of squares along x,y,z to distace of IJ group
            RKJ0 = XKJ(JN)*XKJ(JN)+YKJ(JN)*YKJ(JN)+ZKJ(JN)*ZKJ(JN) ! KJ distance
            RIK0 = SQRT(RIJ0*RKJ0) ! IK distance
            CT0 = (XIJ(JN)*XKJ(JN)+YIJ(JN)*YKJ(JN)+ZIJ(JN)*ZKJ(JN))/RIK0
            CT1 = MAX(-pt999,CT0)
            CT2 = MIN(pt999,CT1)
            CST(JN) = CT2
            ANT(JN) = ACOS(CT2)
            RIJ(JN) = RIJ0
            RKJ(JN) = RKJ0
            RIK(JN) = RIK0
	    
Then calculating the energy:

	DO JN = 1,nTA
            ANT0 = ANT(JN)
            DA = ANT0 - ANTC(JN)
            DF = TK(JN)*DA
            DFW(JN) = -(DF+DF)/SIN(ANT0)
          END DO

calculating the force:

		DO JN = 1,nTA
		    I3 = IT(JN)
		    J3 = JT(JN)
		    K3 = KT(JN)
	C
		    ST = DFW(JN)
		    STH = ST*CST(JN)
		    CIK = ST/RIK(JN)
		    CII = STH/RIJ(JN)
		    CKK = STH/RKJ(JN)
		    DT1 = CIK*XKJ(JN)-CII*XIJ(JN)
		    DT2 = CIK*YKJ(JN)-CII*YIJ(JN)
		    DT3 = CIK*ZKJ(JN)-CII*ZIJ(JN)
		    DT7 = CIK*XIJ(JN)-CKK*XKJ(JN)
		    DT8 = CIK*YIJ(JN)-CKK*YKJ(JN)
		    DT9 = CIK*ZIJ(JN)-CKK*ZKJ(JN)
		    DT4 = -DT1-DT7
		    DT5 = -DT2-DT8
		    DT6 = -DT3-DT9
	C

		    Fx(I3) = Fx(I3)-DT1
		    Fy(I3) = Fy(I3)-DT2
		    Fz(I3) = Fz(I3)-DT3
		    Fx(J3) = Fx(J3)-DT4
		    Fy(J3) = Fy(J3)-DT5
		    Fz(J3) = Fz(J3)-DT6
		    Fx(K3) = Fx(K3)-DT7
		    Fy(K3) = Fy(K3)-DT8
		    Fz(K3) = Fz(K3)-DT9

		  END DO
	  
The energy: E = 0.0
E = E + TK(i)*(ANTC(i)- ANT(i))** 2
E = E/2.0 ! value returned, angles energy! :)

# bonds.f
Bonds  computes the hookean force between chosen atoms. returns the energy E:
subroutine Bonds(E)
run over all number of bonds and assign the first index (Ib1(i)) to I2, and its pair to J2. calculate the difference of x,y,z components for each pair,
and from that the square distance (sum of squares) . assigning the distance (square root) to r1. 
	    
The energy is then calculated by: 
E = E + bk(i)*(r1-Rb(i))** 2 ! this was used to be saved into file 46??

the index i indicates the interaction between particle i and i+1
            f = RBC(i)/r1 - bK(i)
	    
as dx = X(I2) - X(J2), the force is calculated: Fx(I2) = Fx(I2) + f * dx, note the dx sign to save variables: Fx(J2) = Fx(J2) - f * dx	 
only at the end:
E = E/2.0
(where is end of do???)
Also including : BONDSP does the same as ANGLP but for the bond lengths; computes RBC for each bond (this computation is done only once) to compute the force: f = RBC(i)/r1 - bK(i) where r1 is the delta from the optimal length. 
subroutine bondsp: calculate for any bond i: RBC(i) = Rb(i)*bK(i)

# box.c
creates a square box with dimentions Xmin,Xmax,Ymin,Ymax,Zmin,Zmax. Also updates energy Ebox (pointer):
F = K*dR
E = (K*(dR** 2))/2 
where R is the distance of the penetration of the atom to the wall of the box. Given K  = *Kforce, K2 = K*2 
For all dynamic atoms:

	if (X[CCCi] < Xmin ){
			Fx[CCCi]+= -(K2* (X[CCCi] - Xmin));
			ExBox = K*pow((X[CCCi] - Xmin),2);
		}

//similarly for X[CCCi]> Xmax, for Y and Z directions too. Finally add contributions: *Ebox += ExBox + EyBox + EzBox

# chiral (not of use to me)
chiral computes the relative distance of the Cbeta bead from chirality and the resulting force to applied on the bead. Returns related energy E.

# coulomb.c
calculates the electrostatic energy term given the dielectric constant (epsilon):
running the electortatic atoms:
calculating the spatial distance between pairs -r ,The force numerator: K2 = K*Q1Q2[i]/epsilon. And the energy: *Eelec+= K2/r.
To calculate the forces in x,y,z directions, use force in the direction C1 to C2 ,devided by r:
F_over_r = K2/(r2*r) so that Fx[C1]+=  F_over_r*dx;   Fx[C2]-=  F_over_r*dx;
Also includes: coulombfactor function that updates esEnergy: given K = 332.0 (???)
*esEnergy = K/((*deConstant)*(sqrt(*sigma)));

# debyehuckel.c
returns electoristatic energy term if use debye-huckel flas is 'YES'. 
debyehuckel function fisrt calculates the distance between electorostatic beads r . If r2 (squared) < 10 printf ("*****************\n") and the
beads indices and distance. Otherwise, calculates the debye-huckel potential using the dhenergytable functions with the given distance.
Then calculates the force in the direction C1 to C2 ,devided by r (is used to compute forces in x,y,z directions)
Also calculates the electorstatic force: Fx[C1]+=  F_over_r*dx;       Fx[C2]-=  F_over_r*dx; so on for y, z...
debyehuckelfactor function updates esEnergy:
*esEnergy = K*(*saltCoefficient)*exp(-(*screeningFactor)*(sqrt(*sigma)))/((*deConstant)*(sqrt(*sigma)))

# dhenergytable.c
this script generates array of debye-huckel potential in an intervals of r^2.
contains dhenergytable . This function gets a distance (squared) and returns :
DebyeHuckelPotentials[i]= K*(*saltCoefficient)*exp(-(*screeningFactor)*(Dist))/((*deConstant)*(Dist))
DebyeHuckelForces_over_r[i]=DebyeHuckelPotentials[i]*(1/DistanceSquared + *screeningFactor/Dist)
for many intervals of distances up to a cutoff.

# distances.f
Distances calculate the distance between any two beads: dx = X(i) - X(j) etc.. then write(93,*) dx** 2+dy** 2+dz** 2 where 93 is TrajDist.

# dynamicRange.f
contains genDynRange(DynLength,DynamicStr,DynamicAtomRange): creates DynamicAtomRange array to be used. 
function isDyn(DynamicAtomRange,DynLength,AtomNum): runs over all atoms returning 1 (true) for dynamic atoms, 0 otherwise.
function numDynAtom(DynamicAtomRange,DynLength): calculates the number of dynamic atoms

# improper.f

# IjM.f
CONTACTS: computes the force on all atoms due to contacts via a 10-12 potential 
***subroutine contacts(E, Conts, count, outE)***

		if (writeContactsRanges .eq. 'YES') then
			do i=1, rangesNumber
			  currentRangesIndex(i) = 1
			  EtotalRanges(i) = 0.0
			  TwoBodyRanges(i) = 0
			end do
		      end if

runs over all contacts and assigns indices to C1 = IC(i)    C2 = JC(i). Calculate square distance between each pair as r2. rm10=(sigma/r^2)^2 where sigma is already distance squared (user input in .dat). Stores the current energy value currentEnergy=epsC(i)*rm10*(5*rm2-6).
calculate E=E+currentEnergy and then  f_over_r = epsC(i)*rm10*(rm2-1)*60/r2.
If the distance is lower the contacts upper threshold then define isContact as true (1) .

	 if ((writeAllContacts.eq.'YES').and.(outE.eq.0))then	
		   write(14,*) isContact
		 endif
		 
	if (writeContactsRanges .eq. 'YES') then
		  do j=1, rangesNumber
		    currentContactInRange = 
	     Q      rangeContacts(j,currentRangesIndex(j))
		    if (i .eq. currentContactInRange) then
		      EtotalRanges(j)=EtotalRanges(j)+currentEnergy
		      TwoBodyRanges(j)=TwoBodyRanges(j)+isContact
		      currentRangesIndex(j)=currentRangesIndex(j)+1
		    end if 
		  end do	
		end if
	! then calculate Fx, Fy etc using f_over_r
	
	      if ((writeContactsRanges .eq. 'YES').and.(outE.eq.0)) then
		do j=1, rangesNumber
		  write (15,'(I5)') TwoBodyRanges(j)
		  write (16,'(F9.2)') EtotalRanges(j)
		end do
		write (15,*) ''
		write (16,*) ''
	      end if
	      end

# phi.f (irrelevant for my IDPs use...)
dihedral computes the dihedral angles and the forces due to them .


# rep.f
NonContacts computes the forces due to non native contacts. 
subroutine noncontacts(E):  go over the NCset which is a subset of all the non-native contacts
the subset includs all non-native contacts that may be below the
repultion distance in this iteration. Calculates repulsive part:  E = E + NNCSigma(tt)*rm14*r2 . calculate forces and returns E = E/12.
Also:
!* NonContactsNear computes the non-native contacts to be used for the*
!* next NCnear steps                                                  *
!*                                                                    *
!* every 50 steps this subroutine is called                           *
!* it checks each non-native contact,                                 *
!*if the contact length < 9*repultion_distance, it is inseted into the*
!*to check list (I guess the probabilty of an atom traveling          *
!*8*repultion_distance in 50 steps is very low )

# repEllipsoid.f (not of my use)
NonContacts computes the forces due to non native contacts

# sym.f
includes IntSymp() : integrates over one time step for all particles  .
subroutine intsymp(E,ET, Conts, outE, count ) ! again E is total potential energy , ET is total includig kinetic+potential
for i=1:9 initialize: DumE(i)=0.0. Then call for calculations:
call Fstart, call Bonds( DumE(1)), call ANGL(DumE(2)), call dihedral(DumE(3)), call contacts(DumE(4), Conts, count, outE), call NonContacts( DumE(5)), optional:
call box(X,Y,Z,Fx,Fy,Fz ,DynamicAtomRange, DynLength,DumE(6),boxMin,boxMax,boxCoeff), chiral(DumE(8)), ellipsoidRepulsions(DumE(9))
for electorstatics either debye-huckel or coloumb: 
debyehuckel(useESCutoff,useDHTable,esFirstAtomIndex,esSecondAtomIndex,esCharge,X,Y,Z,Fx,Fy,Fz,esPairsNum,
,esDistanceCutoff,DumE(7),deConstant,screeningFactor,saltCoefficient,DebyeHuckelPotentials,DebyeHuckelForces)

Call Findtemp(DumKE) ! returns the kinetic energy into DumKE
      ET = E + DumKE
      
if(outE .eq. 0)then
         E = DumE(1) + DumE(2) + DumE(3) +  DumE(4) + DumE(5) + DumE(6) + DumE(7) + DumE(8) + DumE(9)
then update energy terms to file 70. 	 
Thermostat:
      RsBy =sqrt( 1+(T/temp-1)*tau/Rstau )   ! ( RsBy = 1.0 turns off thermostat)
Then:

			Vx(i) = (Vx(i) + tau*Fx(i)/ms(i))*RsBy
			Vy(i) = (Vy(i) + tau*Fy(i)/ms(i))*RsBy
			Vz(i) = (Vz(i) + tau*Fz(i)/ms(i))*RsBy
	! displacements of particle i
			X(i) = X(i) + tau*Vx(i)
			Y(i) = Y(i) + tau*Vy(i)
			Z(i) = Z(i) + tau*Vz(i)		
		
      
# tripleproduct.f (not of interest in my IDPs studies)
 computes the improper dihedral (using triple product) 
