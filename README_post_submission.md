# MD_repo
Simulation description:
Inputs are: MDwrapper.prefs file and <name>.dat file as attached here:
fragB_helix_translated_updated.dat 
MDWrapper.prefs 
Running with ~lavi/scripts/MD/MDWrapper.pl where MDWrapper.pl reads MDWrapper.prefs.schema , MDWrapper.prefs files.

# **MDWrapper.prefs**:
A configuration file including:


 - BASE_PATH path from which the script rans (DIRECTORY NEW WRITABLE DEFAULT:.); 
 - GO_PATH: (BINARY) path for MD.exe file that exists!; 
 - INPUT_FILE_PATH ./<name>.dat file (TEXT);
 - RUN_POST_EXECUTION_ANALYSIS=NO; optional:
 - #POST_EXECUTION_ANALYSIS_PROG=/home/lavi/scripts/ProtMTDiff.pl;
 - EXECUTION_TYPE (FLAG OPTIONS:LD|MD); 
 - GAMMA (FLOAT) physical parameter: damping in Langervine equation;
 - TAU thermostatic coupling; 
 - **WRITE_ALL_CONTACTS** flag (YES|NO); 
 - CONTACT_RANGES_DEFINED (YES|NO)** (??); 
 - HAS_STATIC_ATOMS (YES|NO); DYNAMIC_ATOMS_RANGE (ARRAY TYPE:INTEGER) defining static atoms region like 1,3,6,8 : CA 1-3 and 6-8 ;
	are dynamic (4,5 static); running post execution analysis flag; CONFINE_IN_BOX (YES|NO) flag; 
 - BOX_X_MIN (FLOAT DEFAULT:-150.0) up to BOX_Z_MAX all box spatial limits; BOX_COEFFICIENT: (FLOAT) how repulsive is it; 
 - DIELECTRIC_CONSTANT (FLOAT) value; 
 - APPLY_ELECTROSTATICS_BETWEEN_NEIGHBORS (YES|NO); ELECTROSTATICS_CUTOFF_TYPE (OPTIONS:*ENERGY*|DISTANCE|NONE DEFAULT:DISTANCE);
						ELECTROSTATICS_ENERGY_CUTOFF (FLOAT); USE_ELECTROSTATICS (YES|NO);
 - USE_DEBYE_HUCKEL; 
 - USE_DH_ENERGY_TABLE(??);
 - IONIC_STRENGTH corresponding to concentration (0.05 is ~50*3 Molar, calibration is required); 
		IONIC_RADIUS; SOLVENT_DENSITY; 
 - MODEL_TYPE (CA|CB),EXECUTION_STEPS (INTEGER), OUTPUT_FREQUENCY-> how many steps to 'skip' writing output;
 - TEMPERATURES -> fraction of 273; ITERATIONS_PER_TEMPERATURE: how many repetitions for this input; 
 - STORE_TRAJECTORIES; TRAJECTORY_OUTPUT_FREQUENCY; **STORE_TRAJECTORY_DISTANCES**;
(also the following not to touch parameters: EXECUTION_LOCATION=LOCAL;
 - #EXECUTION_LOCATION=CLUSTER, MAX_PROCESSORS_ALLOCATION=350;USE_EMPTY_SLOTS (YES|NO); USE_COPY_QUEUE=NO (YES|NO); 
  CLEAN_SCRATCH_FOLDER=YES (YES|NO)
				
# **MDWrapper.prefs.schema**: 
A template file including: defaults for the above defenitions and type (int, float etc), options to choose from in brackets above:
 
[PREFERENCES]
 - LOG_FILE=TEXT NEW DEFAULT:MDWrapper.log;
 - EXECUTION_LOCATION=FLAG OPTIONS:CLUSTER|LOCAL; 
 - MAX_PROCESSORS_ALLOCATION=INTEGER DEFAULT:30 etc...             
More interseting:
 - **USE_INITIAL_CONDITIONS=FLAG OPTIONS:1|2|3 DEFAULT:3; INITIAL_CONDITIONS_FILE_PATH=TEXT; **
 - **RUN_POST_EXECUTION_ANALYSIS=FLAG OPTIONS:YES|NO DEFAULT:NO, POST_EXECUTION_ANALYSIS_PROG=TEXT;** 
 - CONTACT_RANGES_DEFINED=FLAG OPTIONS:YES|NO DEFAULT:NO;
 - CONTACT_RANGES_FILE_PATH=TEXT; ELECTROSTATICS_ENERGY_CUTOFF=FLOAT DEFAULT:0.005; ELECTROSTATICS_DISTANCE_CUTOFF=FLOAT DEFAULT:200.0; 
 - COMPENSATE_ELECTROSTATIC_CONTACTS=FLAG OPTIONS:NONE|ATTRACTION|ALL|NO_ES DEFAULT:ATTRACTION;SE_CHIRALS=FLAG OPTIONS:YES|NO DEFAULT:NO;
 - USE_ELLIPSOID_REPULSIONS=FLAG OPTIONS:YES|NO DEFAULT:NO;

[DEPENDENCIES]
 - **USE_INITIAL_CONDITIONS=1,2:INITIAL_CONDITIONS_FILE_PATH**;
 - EXECUTION_TYPE=LD:GAMMA; HAS_STATIC_ATOMS=YES:DYNAMIC_ATOMS_RANGE; 
 - CONFINE_IN_BOX=YES:BOX_X_DIMENSION,BOX_Y_DIMENSION,BOX_Z_DIMENSION;
 - USE_ELECTROSTATICS=YES:DIELECTRIC_CONSTANT,APPLY_ELECTROSTATICS_BETWEEN_NEIGHBORS,ELECTROSTATICS_ENERGY_CUTOFF,USE_DEBYE_HUCKEL,
		COMPENSATE_ELECTROSTATIC_CONTACTS
 - USE_DEBYE_HUCKEL=YES:IONIC_STRENGTH,IONIC_RADIUS,SOLVENT_DENSITY,USE_DH_ENERGY_TABLE; MODEL_TYPE=CB:USE_CHIRALS,USE_ELLIPSOID_REPULSIONS;
 - RUN_POST_EXECUTION_ANALYSIS=YES:POST_EXECUTION_ANALYSIS_PROG;       
		STORE_TRAJECTORIES=YES:TRAJECTORY_OUTPUT_FREQUENCY,STORE_TRAJECTORY_DISTANCES
 - ELECTROSTATICS_CUTOFF_TYPE=ENERGY:ELECTROSTATICS_ENERGY_CUTOFF; ELECTROSTATICS_CUTOFF_TYPE=DISTANCE:ELECTROSTATICS_DISTANCE_CUTOFF;
	CONTACT_RANGES_DEFINED=YES:CONTACT_RANGES_FILE
                
# **MDWrapper.prefs.pl**: 
Reads the above files, creates a submission script in the executing directory, and ouputs of settngs.dat and random.dat to be described here.
	
- reading MDWrapper.prefs.schema, MDWrapper.prefs into executionPreferencesRef;
	checking validation with "ExecutionPreferences.pm" module at "/home/arielaz/scripts/util/" or "/home_c/arielaz/scripts/util/" :
	if the user input is valid, otherwise prompts a warning message like:  
	
		      print "Would you like to store the execution preferences? (YES|NO):";
		      chomp (my $userInput = <STDIN>);
		      $userInput = "YES" if ($userInput eq "");
		      if ($userInput !~ m/YES|NO/i)
		      {
			warn "Invalid answer $userInput.\n";
		      }
	
- calling processExecutionPreferences($executionPreferencesRef) which inserts to the given prefrence array the user ID ('rivka') to the first cell,
	 - then it initializes a random seed (to be overwritten) with the function time:
	
			$executionPreferencesRef->{"RANDOM_SEED"} = time() % 100000;
	 - reads the input file, paths and so on.
	 - change initial conditions is optional agrument to start from a different trajectory
	
				  if ($executionPreferencesRef->{"USE_INITIAL_CONDITIONS"} != 3)
			  {
			    ($executionPreferencesRef->{"INITIAL_CONDITIONS_FILE_NAME"}) = $executionPreferencesRef->{"INITIAL_CONDITIONS_FILE_PATH"} =~ m/.*\/(.*)/;
			  }
	
	- reads the temperature value
	- define dynamic atom range 
	
			if {"HAS_STATIC_ATOMS"} =~ m/YES/i
	
- calling createOutputFolders
	- If execution directory is not defined then: 
		creating outputfile paths including our initial input file name, random generated seed,_has_ if has static atoms, box dimensions with _cib_,
	  _es_ for electristatics dielectric constant, _dh_ for debye huckel, adding ionic strength, radius and solvent density.
	- within it creates the following directories:
	- input
	- output
	- /output/Temperature
	- /output/2Body
	- /output/3Body (3 body intractions)
	- /output/Etotal
	- /output/EbyType
	- /output/Log
	- /output/LastCord
	[optional]
	- /output/Traj (only if you choose to store trajectories)
	- /output/TrajDist (only if choose to store distances)
	- /output/AllContacts (only if WRITE_ALL_CONTACTS is YES)
	- /output/ContactRanges (only if CONTACT_RANGES_DEFINED is YES)
 
- calling launchExecutions
	for every iteratio=1:number_of_iteratins_per_temperature:
		- creates a random seed (random integer say 15651 )
		- creates a settings file (in input dir) with that random number_.dat: ('settings_15651.dat')
			this file includes all flags chosen in MDWrapper.prefs (or schema's default). The file from which inital coordinates/velocities 
			were taken (if not defined otherwise then initial <name>.dat); Also the .log file, parameters like box limits and so.
		- creates a random file  including the random number of the name generated here above. Within the file created 3 different x,y,z
			random number . (for example 'random_15651.dat' includes one line:6543 20424 33466) then prints:
			
				 print MDW_LOG_FILE_HANDLE localtime(time),":Random values: $i $j $k\n"
	
		- calls: launchExecution($currentTemprature, $settingsFileName, $randomNumber, $randomFileName). To be described in the next section.
		- calls reportExecution($randomNumber) for the user: printing to the MDWrapper.log (in execution dir and also base path)
			the local time, number of steps, temperatures submitted and the base seed.
- finally launchExecution (note the missing 's' as the suffix)
	for each given file creates a submission .csh file ('_fragB_helix_translated_updated_0.5_15651.csh') in the execution dir including:
	
		#!/bin/csh -f
		#
		#$ -cwd
		#$ -j y
		#$ -S /bin/csh
		#
		#$ -M myemail
		#$ -e error_file
		#$ -o MDWrapper.log
		set hostname = `hostname`
		set out_dir = /home_d/rivka/MD_base_rivka_updated/fragB_helix_translated_updated_41648_hsa_cib_1998_1998_1998_es_80_dh_0.001_1.4_1.0/output
		set scratch = /scratch/rivka/fragB_helix_translated_updated_0.5_15651
		echo '#########################'
		echo Running on host `hostname`
		echo Time is `date`
		echo Directory is /home_d/rivka/MD_base_rivka_updated/fragB_helix_translated_updated_41648_hsa_cib_1998_1998_1998_es_80_dh_0.001_1.4_1.0
		echo temp Directory is $scratch
		mkdir $scratch
		cp /home_d/rivka/MD_base_rivka_updated/fragB_helix_translated_updated_41648_hsa_cib_1998_1998_1998_es_80_dh_0.001_1.4_1.0/input/settings_15651.dat $scratch/settings.dat
		cp /home_d/rivka/MD_base_rivka_updated/fragB_helix_translated_updated_41648_hsa_cib_1998_1998_1998_es_80_dh_0.001_1.4_1.0/input/random_15651.dat $scratch/random.dat
		cp ./fragB_helix_translated_updated.dat /home_d/rivka/MD_base_rivka_updated/fragB_helix_translated_updated_41648_hsa_cib_1998_1998_1998_es_80_dh_0.001_1.4_1.0/input/
		cp /home_d/rivka/MD_base_rivka_updated/MD.exe $scratch
		cp ./fragB_helix_translated_updated.dat $scratch
		cd $scratch
		/usr/bin/time -p ./MD.exe >& time_fragB_helix_translated_updated_0.5_15651.out
		bzip2 TemperatureFile_fragB_helix_translated_updated_0.5_15651.dat
		bzip2 2body_fragB_helix_translated_updated_0.5_15651.dat
		bzip2 3body_fragB_helix_translated_updated_0.5_15651.dat
		bzip2 Etotal_fragB_helix_translated_updated_0.5_15651.dat
		bzip2 EbyType_fragB_helix_translated_updated_0.5_15651.dat
		bzip2 fragB_helix_translated_updated_0.5_15651.log
		bzip2 Traj_fragB_helix_translated_updated_0.5_15651.dat
		cp Final_fragB_helix_translated_updated_0.5_15651.dat $out_dir/LastCord
		cp fragB_helix_translated_updated_0.5_15651.log.bz2 $out_dir/Log
		cp TemperatureFile_fragB_helix_translated_updated_0.5_15651.dat.bz2 $out_dir/Temperature
		cp 2body_fragB_helix_translated_updated_0.5_15651.dat.bz2 $out_dir/2Body
		cp 3body_fragB_helix_translated_updated_0.5_15651.dat.bz2 $out_dir/3Body
		cp Etotal_fragB_helix_translated_updated_0.5_15651.dat.bz2 $out_dir/Etotal
		cp EbyType_fragB_helix_translated_updated_0.5_15651.dat.bz2 $out_dir/EbyType
		cp Traj_fragB_helix_translated_updated_0.5_15651.dat.bz2 $out_dir/Traj
		rm -r $scratch

	The above script copies all intersting files to the scratch where they run, goes there and then runs ./MD.exe file in the provided directory
	Later it zips all outputs to be easily copied and unzipped for analysis if necessary. Copies the final step's trajectory to LastCord,
	and all bz2 files to the corresponding output directorys before deleting from scratch :) (use bunzip2 to unzip)

! not it may be useful sometimes to write the distances or specific distnces like End-to-end. 
			
	
	
	
