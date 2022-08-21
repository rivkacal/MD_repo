 # Description 
This file describes all the calculations performed in each file seperately, strating from running **MD.exe** file as performed in .csh file created at MDWrapper.pl

# fortran shortacuts
- 60 is settings.dat (mandatory)
- 6 is Trajectory file 
- 93 is TrajDist file
- 99 is EnergyTot
- 70 is EnergyTerm
- 12 is ContactFile
- 13 is ThreeBodyFile
- 14 is writeAllContacts (optional)
- 15 is ContactRangesTwoBodyFile (optional if above is set 'YES')
- 16 is ContactRangesEnergyFile
- 53 is TemperatureVtime
- 63 is 'random.dat' file
- 71 ???? holds temperature?

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
          

          
# ss.f
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

  - stophere subroutine should close all open files (say 60 settings.dat) (??? unclear how). It says to Close the files that holds the temperature (71)    
      no closing of 60 ... reads again again xrandom,yrandom,zrandom...
 
 - initES()
    defines parameter (screeningFactorConstant = 2529.11892)
    if using Debye Huckel supposedly transforming simulation temperature to real temperatures with realTemperature = 295*T.
        then calculae kappa :screeningFactor = SQRT((screeningFactorConstant*ionicStrength*solventDensity)/(deConstant*realTemperature))
        saltCoefficient = (exp(screeningFactor*ionicRadius)/(1+screeningFactor*ionicRadius))     
    if cutoff is of 'ENERGY' type then define minimal distnace to apply electrostatics as esMinDistance = 16.0, (Angstorm), esMaxDistance = 100000000.0
        then, if use Debyehuckel call: debyehuckelfactor(esMinDistance,deConstant,screeningFactor,saltCoefficient,esEnergy)
        otherwise call coulombfactor(esMinDistance,deConstant,esEnergy)      
    (???)
    
        then calculates: esCutoff = esEnergy * esEnergyCutoff (the latter is the fraction defined above),esEnergyError = esCutoff / 100.0
        then checks 5000 times: if the energycutoff minus the esEnergy absolute value is lower than the error (not zero exactly) and ouputs error message if does.,
         starting from average value, each iteration changing the maximal or lower distance according to ratio of esEnergy (above and below) to cover all options
    if cutoff is of 'DISTANCE' type just check if esDistanceCutoff **squared** (** 2) is greater than maximal cutoff to abort, If it's false just
      take the square value of input to continue the calculations with the new esDistanceCutoff .
      
    
    
        
      

