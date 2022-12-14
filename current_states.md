
# 11/10/2022
short day: communication outage + writing research proposal. Also:

coarse-grain:
*run control experiment for hist_per_test.py experiments (hist_0p5charge_control.py) for trying to assess the effect of fraction of charged histidine 
and its ordering. The control of hist0p5 files includes the parameters for HIS1 type with a charge of 0.5 to compare the charge effect.
  
  -results are in: /trajectories/rivka/his_per_tests_results/his0p5/
  -analysis ran using: ./run_his0p5_analysis.sh after ran ,oved to analysis_scripts
  -final results: hists_per_tests_fragB_hisEE_his05control.ipynb notebook
  and check all hists_per_tests_fragB_translated_histdistsEEgrouped/&grouped_control/&even/&even_control 4 files. 
  Figures are saved with horizontal line containing the previous only his0.5 results as control. see for instance: uniform_control.png


*trying to apply the H++ for known histidine charges of my previous amberff simulations for fragment B.
Results for simulation 10ns frame after 500ns previous run with all histidines non-charged for internal dielectric constant of 4 and 10
are in pdbs_images/0.125_80_4_pH7_fragB_noH.result.pdb and 0.125_80_10_pH7_fragB_noH.result.pdb showing as expected 0 charged histidines
predicted see fragB_all_his0_hpp.png. However for simulations of all 16 HIS charged: 0.125_80_10_pH7_segB_hispplus1_500ns_then_10ns_nH.result.pdb
only two predicted as charged histidines : fragB_all_his1_hpp.png. 
The importance of initial frame for determining HIS form is crucial and challenging for IDRs

all-atom:
*create /home_d/rivka/Bioinformatics_HIS/all_atom/HSP_PHE  &&  /home_d/rivka/Bioinformatics_HIS/all_atom/PHE_PHE 
  && /home_d/rivka/Bioinformatics_HIS/all_atom/HSD_HSD  
  please note the charm ff directory of HSP includes an extended ffbonded file with NH2 CTA CT2A angles and dihedrals with H, CPH1, HA2 which are missing
  and copied form newer ff versions.

*prepare initial frame to run from and pulling conformations to submit on cluster like: /home/rivka/Bioinformatics/all-atom/minimize_geometry/HSP_PHE

# 12/10/2022
short day as before:

all-atom simulations were not running during the night but during this day.
*mistake in submit file scratch directory -> script ran but was not outputting! resubmtting everything
#simulations running super-slow - rethink using 32 cores ind=stead 16 


prepare a script to gather all the free-energy values and com to plot:
*initial stepes: /home_d/rivka/Bioinformatics_HIS/all_atom/HSD_PHE/umbrella/free_energy_grep.sh (will later be moved to a scripts useful directory like:
/home_d/rivka/Bioinformatics_HIS/all_atom/useful_files
reading the output of uniformlly spaced initial com conformations of the two pi groups from a given
/trajectories/rivka directory into a free_energy.dat output file including the conformation id (e.g. id 17 for :  umbrella17.edr conf17.gro)
the initial com hraminically restraned distance and the free energy average value.
*finished writing script: see result vi /home_d/rivka/Bioinformatics_HIS/all_atom/HSD_PHE/umbrella/free_energy_vs_com.dat 
*find script also at: vi /home_d/rivka/Bioinformatics_HIS/all_atom/useful_files//free_energy_grep.sh


# 13/10/2022
*working on research proposal
check if H++ can recognize NMR studied histidines correctly: 10.1002/prot.10177 
-1PNT
Residue      pKint          pK_(1/2)  experimental pKa        
HIS-66        6.475        7.892          9.2
HIS-72        5.937        7.612          8.3

-1a2p
Residue      pKint          pK_(1/2)  experimental pKa
HID-102        6.501         6.099        6.3
HIS-18         6.608         6.363        7.9

-1ymb
Residue      pKint          pK_(1/2)  experimental pKa
HIS-81        6.630        6.429         6.6
HID-36        6.650        7.611         7.8
HIS-113        6.561        5.202        5.4
HIS-116        6.505        6.368        6.6
HID-119        6.014        2.531        6.4

-1gym
Residue      pKint          pK_(1/2)  experimental pKa
HID-32         5.333       >12.000      7.6
HIS-82         5.709         6.412      6.9
HIS-92         5.667         5.932      5.4
HID-227        7.285         7.975      6.9

-7rsa
Residue      pKint          pK_(1/2)  experimental pKa
HID-12        7.280        5.237       5.8
HID-105        6.229        6.667      6.7
HIS-119        5.762        5.942      6.2

-1stg
Residue      pKint          pK_(1/2)  experimental pKa
HID-8        7.008        6.925        6.8
HIS-46        6.012        6.417       5.7
HID-124        5.451        4.036      6.0


*plot free energy profile!! see HSD_PHE/umbrella dir notes in /home_d/rivka/Bioinformatics_HIS/all-atom using wham!! 
title "Umbrella potential"
@    xaxis  label "\xx\f{} (nm)"
@    yaxis  label "E (kcal mol\S-1\N)"

also check histo.xvg for maximal overlap to ensure results.
@    title "Umbrella histograms"
@    xaxis  label "\xx\f{} (nm)"
@    yaxis  label "count

*create a script for searching T-faced histidines in all files, not just 'charged' directory with all aromatic residues!
see: vi /home_d/rivka/Bioinformatics_HIS/python_codes/test_analysis/find_his_prependicular.py
* restricting more with the script:  python find_his_prependicular_quantumstruc.py

# 14/10/2022
*research proposal, try to justify that all existing quantum info does not yield observed conformations. Need to MD is there,..

*results from looking for prependicular T-shaped like conformation are there running over 6535 pdb files.
vi T_shaped_all.csv: 468 such interactions but repeating chains, distance of 4A like quantum!

*find only his-his prependicular: python find_his_his_prepend.py
*from T_shaped_his_his4A.csv: 276 such conformations

running:
-rw-r--r--. 1 rivka Koby  18K Oct 14 14:18 find_his_his_prepend.py
-rw-r--r--. 1 rivka Koby  18K Oct 14 14:36 find_his_parallel.py                      #all histidine-aromatic(histidine,PHE,TRP,TYR)
-rw-r--r--. 1 rivka Koby  18K Oct 14 14:36 find_his_parallel_his_only.py             #only his-his parallel 

*parallel results



*plot free energy profile and not the potential energy or LJ-SR only! using wham from umbrella tutorial.
See updated notes at /home_d/rivka/Bioinformatics_HIS/all_atom/HSD_PHE/umbrella


# 16/10/2022

resubmitting cluster runs. Turns out using 64 cores 4,8 options improves run by 100 times faster!!!!!!!!!!!!!!
compared to using 32 only (running two simultaneously)

resubmitting: HSD_HSD , PHE_PHE and HSP_PHE all using charm36!! (note dates) . 1.2 hour per run, total 41 runs -> 13 hours using 4 nodes if available :)

*preparing new HSD_PHE simulations in water!! /home_d/rivka/Bioinformatics_HIS/all_atom/HSD_PHE/umbrella/water
to compare with quantum results http://journal.chemistrycentral.com/content/7/1/92

# 17/10/2022
some all atom simulations returned!
adding PHE-PHE and HSD-HSD and HSP-PHE to profiles graph in umbrella_plots.ipynb located at :
/home_d/rivka/Bioinformatics_HIS/all_atom/HSD_PHE/umbrella/txts (until it'll shift)

# 18/10/2022
conference : 8:30 - 17:30 
results from HSD_PHE_charm36_water yielded a potential that's not satisfying...

# 23/10/2022
working on the presentation for Tuesday. Research proposal discussion with Koby: 
avoid using 'jargons' as LLPS or names of specific coarse-grain models as HPS and KH.
make a flow that's organized. Number goals and refer to it by ordering. 

# 24/10/2022
preparing presenations. Reading about pi-pi stacking in force fields and types of force-fields to be used...


# 25/10/2022
presenting my presentation!! so far nite to include titles for graphs and or descritions to follow!
milan suggested using MD force fields or QM methods to use r,teta and phi coordinates to get minimal energy geometry! 
but then many simulations to run. Also 50 ns may be enough...

# 26/10/2022
login to chemfarm!!!!
ssh -X rivkac@chemfarm.weizmann.ac.il

analyze his sequence fraction in IDPS > 0.02 or 2% . using milan's code to download from Uniprot server and a predictor of IDPs probability 
(to extract > 0.5%)
results in directory : /home_d/rivka/Bioinformatics_HIS/stats/HIS_IDPS
his_perc.txt presents precentage of hisidine for any single sequqnce (~35000)
sequqnces : idr_sequences.txt
sequqnces IDS: idr_sequences_ids_human_proteome.txt
analysis in notebook
also downloaded orca_5_0_3 on chemfarm having some problems...


# 27/10/2022
same analysis for PDB stuctured files for histidines!
directory: /home_d/rivka/Bioinformatics_HIS/stats/HIS_ordered
current comparison wasnot sufficent (non redundant, counts are missing...)--> old file
non-redundant lists from PISCES: https://dunbrack.fccc.edu/pisces/PISCES_OptionPage.php
https://dunbrack.fccc.edu/pisces/download/cullpdb_pc60.0_res0.0-5.0_noBrks_len40-10000_R1.0_Xray+Nmr+EM_d2022_10_11_chains34482.fasta


# 29/10/2022
version 2 of research proposal.
got stuck with some issues of avogadro2 trying to import python code for input to orca. 
Still orca script is not running properly due to misuse of cores...

# 30/10/2022
submitted research proposal!!!
continue working on histidine statistics with ordered or disordered...
install plotly trying to create an interactive scrollbar for choosing the best number of bins...
see: https://towardsdatascience.com/how-to-quickly-find-the-best-bin-width-for-your-histogram-4d8532f053b0

not working...

Checking for consequtive histidine sequqnces. Check: /stats/HIS_IDPS/dis_his_indices. found sequqnces with the following idx as well:

[19, 28, 103, 107, 108, 110, 112, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 133, 143]

[29, 69, 70, 78, 88, 102, 117, 118, 119]

[2, 5, 7, 23, 44, 45, 48, 54, 64, 67, 69, 74, 106, 134, 135, 136, 147, 149, 150, 157, 170, 174, 183, 207, 218, 228, 239, 240, 296, 297]

[66, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 96, 98, 102, 106, 108]

find a nice algorithm to quantify these sequqnces...

# 31/10/2022
trying to test different HIS sequences: disordered_analysis.ipynb:
use groupby to cluster histidines indices and length as in file: his_array_sizw. 
The counter (line) of sequence is also appended to his_array_lines_idx where the lines correspond to the id from the IDPS equal or longer than 40 list:
proteins_id_larger40
Now dividing this clusters into groups of 2,3,4,5 etc...:
finally creating len_arr which contains the size of largest cluster for each sequence of histidines (from 0 to 14)!!
same analysis for all residues: '../HIS_IDPS/dis_{filename}_indices' where filename fits CYS,ALA,PHE,ARG....
no apparent correlation found. Try to compare to LYS and ARG as positive: ../HIS_IDPS/dis_positive_frac
see: '../HIS_IDPS/dis_positive_indices
missing counts!!i

---election and vacation!! --

# 6/11/2022 
commend for compiling CG: make -f MDmake correct file...
trying to work out correlations with histidines and other residues. note this paper: https://academic.oup.com/metallomics/article/5/7/904/6015672
succeeded running short simulations of Orca (sinlgle point energy calculations) : /work/rivkac/projects/histidine_qm
see: notes: pbsorca 03_22R1.xyz.inp idle

# 7/11/2022
succeded installing avogadro with orca input on my windows laptop!!
now working on tutorials how to use and study dimers structures...

login from cmd on windows to chemfarm: ssh rivkac@chemfarm.weizmann.ac.il
copy file from windows to chemfarm using cmd:
scp "OneDrive - weizmann.ac.il\myThesisStudies\QM\two_neutral_his_opt.inp" rivkac@chemfarm.weizmann.ac.il:\work\rivkac\projects

# 8/11/2022
courses!! 
fixing 3 day problem: why txt file was returning qsub directive error! 
the orca input file was prepared in windows platform bu avogadro extension to create orca input!
This windows files specific strings may have different lengths other than unix/linux platform of chemfarm!!!!!!!!!!!!!!! 
To convert windows made txt to unix use: 
dos2unix filename
#than wait till convertion is done succesfully. Now optimization script was running on chemfarm successfully using: 
pbsorca trial2.xyz.inp long
#! still problem left to solve: cannot read more than 16 cpus... (32 returns error within file). Maybe 40 will work, etc.. 
Testing how long it'll take to optimize structure using trial2.xyz.inp parameters.
Left to understand the functional choise and other by using orca 5.0.3 documentations PDF.
Also finished simple analysis of histidine in ordered vs disordered (cutting proteins to IDRS based on prediction over 50% to be IDR) data:
see: disordered_structures.ipynb at /home_d/rivka/Bioinformatics_HIS/stats/HIS_IDPS
or ordered_structures.ipynb at /home_d/rivka/Bioinformatics_HIS/stats/HIS_ordered

# 9/11/2022
updating statistics of ordered vs IDRs for proteins < 500 aa
adding a correlation line?
see: structured_less_than500aa.html  or disordered_less_500aa.html within the above diiscovered directories

Also chemfarm trial2.xyz.inp results just finished running after 1 hour and 50 min!!!
learn what to do with the results...

copy a file from chemfarm to windows:
scp rivkac@chemfarm.weizmann.ac.il:\work\rivkac\projects\histidine_qm\trial2.xyz.xyz "OneDrive - weizmann.ac.il\myThesisStudies\QM" 

apparently cannot use unix2dos on cmd (windows) so do that on chemfarm...
unix2dos trial_dos.xyz 
unix2dos: converting file trial_dos.xyz to DOS format ...
scp rivkac@chemfarm.weizmann.ac.il:\work\rivkac\projects\histidine_qm\trial_dos.xyz "OneDrive - weizmann.ac.il\myThesisStudies\QM"
after optimization calculation was done (16 cores - 24 minutes) for 2imidazoles. 
Take each of them seperately, imidazoleA and B to calculate single-point-energy to find the binding energy!

# 10/11/2022
apparently partially stacked opt failed -> increase maxiter to 500 as in 2imid_stacked.inp
also trying to test 2-3/D dimensional PES scan to find other minimas
first trying from orca_manual_5_0_3: for COH2 find 3D_scan.inp.
note that .opt file saves evry conformation so it'll later can be tested, or the last one can be continued from... 

stats analysis: found for structured data a protein with 26% his: 1hcd.pdb file.

1. a) for structured proteins with less than 100 a.a the highest percentage of HIS is ~18. If we consider 118a.a then,I found one with 26% HIS (not sure if this is the one that you once sent me paper about). This protein serves as a pH sensor for changes of H+. Histidines are fully exposed to the outside of the protein surroundings:

    >1HCDA B5B43F527B694716 118 NMR NA NA NA NACO.noDsdr.noBrk Hisactophilin-1 <HATA_DICDI(1-118)> [Dictyostelium discoideum]
    MGNRAFKSHHGHFLSAEGEAVKTHHGHHDHHTHFHVENHGGKVALKTHCGKYLSIGDHKQVYLSHHLHGDHSLFHLEHHG
    GKVSIKGHHHHYISADHHGHVSTKEHHDHDTTFEEIII        

   b) as for the disordered data:

zinc transporter:
39.7% for IDR: sp|Q92504|S39A7_HUMAN      seq: LHDDLQEDFHGHSHRHSHEDFHHGHSHAHGHGHTHESIWHGHTHDHDHGHSHEDLHHGHSHGY

Brain specific POU domain, transcription factor:
33.33% for IDR : sp|Q12837|PO4F2_HUMAN   seq: AASSSSVPISHPSALAGTHHHHHHHHHHHHQPHQALEGELLE

Homebox protein (TF):
44.68% for IDR: sp|P32242|OTX1_HUMAN       seq: HHPHQLSPMAPSSMAGHHHHHPHAHHPLSQSSGHHHHHHHHHHQGYG

Neuron core protein:
31.91% for IDR: sp|O14594|NCAN_HUMAN     seq: IVCTKPRRSHRMRRHHHHHQHHHQHHHHKSRKERRKHKKHPTEDWEK

Zinc transporter:
34.09% for IDR: sp|Q9ULF5|S39AA_HUMAN     seq: TIHEHDLHAAAHNHHGENKTVLRKHNHQWHHKHSHHSHGPCHSG

Putative uncharacterized protein:
30.43% for IDR: sp|Q6ZRP5|YD019_HUMAN    seq: NRSFSAVAATPAKHKHMHTRTHTHMHTHTGMHTLTGTHVHTPHTQMHTRILTLSHMHTHAHTHAHTHGHTHTRAHSTHAHTHAHSHYHTRTL


graphs of dependence of number of positive/negative/aromatic aa clusters (at least 3 consequtive) called stretches is to be organized yet. 

also 3D_scan of 4 atoms took almost 4 hours for very few points...
left to understand how to analyze it and how its working!!!!


# 13/11/2022 
Meeting with Dubai-AI.
trying to follow finding minimas conformation very similar to 
DOI
    https://doi.org/10.1039/C5OB01108F . Using MacroModel of Scrodinger takes 48 hours to get approval to use.
Meanwhile exploring other methods like crest-xtb based (quantum level not ff) presented also in ACONFL paper for finding conformers: 
https://doi.org/10.1021/acs.jpca.2c02439
which I was reading during Gershom-Martin group...

installing here both xtb platform and then crest from grimme's github tutotrial works in test directory (only from downloads) by using:
crest 2imidazoles_control.xyz | tee -a crest_output.txt & nohup

where 2imidazoles_contro.xyz is a regular xyz file of the dseired molecules.
tee -a outputfilename directs content from terminal to a text file too. & nohup should prevent the program from stopping when closing the terminal,
however it's not running in background somehow...

copy best conformation output:
scp best_dos.xyz rivkac@chemfarm.weizmann.ac.il:/work/rivkac/projects/histidine_qm

to perform:
unix2dos best_dos.xyz
and then from cmd on windows:
scp rivkac@chemfarm.weizmann.ac.il:\work\rivkac\projects\histidine_qm\best_dos.xyz "OneDrive - weizmann.ac.il\myThesisStudies\QM"

similarly for all conformations ensemble scp crest_conformers.xyz  rivkac@chemfarm.weizmann.ac.il:/work/rivkac/projects/histidine_qm
then from chemfarm: unix2dos crest_conformers.xyz

after installing my own xtb on chemfarm and crest there, sourcing ~/.bashrc and creating submission script it finally works!!
be carefull: 
qsub submit_crest.sh

now reading useful flags to run with!! see: https://xtb-docs.readthedocs.io/en/latest/crestcmd.html

for now run of 6 minutes became 7 minutes when adding optimization+solvent+other flags as:
  
  crest 2imidazoles_control.xyz --prsc --esort --gfn2 -alpb water --opt i--hess | tee -a crest_solv_out.txt


# 14/11/2022
Paper to note later about charging histidines within coarse-grained (CG) simulations: 
https://pubs.acs.org/doi/full/10.1021/acs.molpharmaceut.2c00337


For citations later using crest: https://pubs.acs.org/doi/10.1021/acs.jctc.9b00143
https://pubs.rsc.org/en/content/articlelanding/2020/CP/C9CP06869D
Having results for 2 imidazoles and for one imidazole plus imidazolium in vaccum and implict solvent (water) need a script to reoptimize the conformers within threshold of energy using orca!!! 

-Note that the conformers observed are within 0.05kcal/mol!
-Note2: having many stacked staggered and prependicular tilted (~35 degrees) for imid0_imid0 both in vaccum and solvent: his0_his0/imid_imid/solv_water/
in crest_conformers.xyz while finding none in his0_his1 for imid and imidazolium raises a suspicion about calculations results depending on inital conformation. That's why another test was made: creating imidazole and imidazolium stacked parallel inital geometry and resubmitting at : 
/work/rivkac/projects/histidine_qm/his0_his1/imid_plus_imidazolium/solv_water/test (note that is important to work in different directories otherwise
crest overwrites the output files of energies and conformers and so on)
These test reveals similar conformers for the approximately same 4kcal/mol difference in solvent, though not exact these conformers will be
reoptimized at higher level DFT later anyways ...
-Note3: for imid+imidazolium in vaccum only one conformer was found with a proton in between the two histidines! this is different then the water simulations but **_may suggest a role in proton (charged) transfer beterrn amino acids depending on conditions!!!_**

script directory: /work/rivkac/projects/histidine_qm/scripts
script works & submits automatically jobs. Think about sleep option... done for all imid0-imid0 imid0-imid1 both in vaccum and solvent.
left to gather results and compare...

optimization takes seconds on the level of BLYP (conformations must be well optimized!!!! yay! :) )

# 15/10/2022
PSF in the morning & then group meeting & then departmental seminar :)
My schrodinger account was accepted: 
downloaded MacroModel !!
make gather_final energy bash script to gather the final single point energy of optimization for each conformer, see scripts dir.

now we get a list with different energies than predicted with the crest methods (should be close). Extract the energies and conf by:
./gather_final engrad > orca_energies
then sort this using: 
sort -k2 -n orca_energies > orca_sorted_energ
remember the optimization also changes the structure, so to get the final structure appears in the last trajectory:
confX_trj.xyz

Given the conformations (left to prepare script to gather them besides energy with threshold), now how to define geometry from xyz file?
first use: https://github.com/zotko/xyz2graph

this is visualized at: /home_d/rivka/Bioinformatics_HIS/QM_calc/visualize_molec.ipynb
for handling calculations: https://molmod.github.io/molmod/tutorial/molecule.html

splitting structures on chemfarm cluster:
obabel copy.xyz -Ofrag.xyz --separate -m
diffuculties installing Open-babel or link to python. Libraries and so staff...
think about splitting the output in chemfarm where the above command works

side note: Lavi's paper whose supporting info contains 1 CG step ~ 50ps https://academic.oup.com/nar/article/48/4/1701/5699673


# 16/11/2022

finding it difficult to work with another format (xyz) try to trasform it to PDB:

obabel frag1.xyz -O frag1.pdb

for copy.xyz file as well. groups will be identified as ligands but different ones! 
pymol does load it. 
Now left to test if it can identify histidine...
Having C and N terimnal successfuly identify histidine!!!:
obabel histidine.xyz -O histidine.pdb
at: /work/rivkac/projects/histidine_qm/his0_his0/imid_imid/solv_water

now omitting C,N terminal up to the carbon1 only: does not recognize as histidine...

running using crest conformations scan reveals too many conformations within a very short enegetic distance --> increase this sampling by:
--ethr 0.1 
still 180 and 165 conformers -->
--ethr 0.2
--ethr 0.5 does not help either...

running for all 215 conformers optimizations and screening later...
many jobs got suspended-> not optimal solution

need to write a script to gather all files got suspended (orca not terminated succesfully)
****ORCA TERMINATED NORMALLY**** for proper files...
gather these into a list to resubmit: 
grep -L  "****ORCA TERMINATED NORMALLY****" *.out > list_resub
then rename to inp:
./replace_out_inp.sh resub
for i in `cat list_resub`; do pbsorca $i idle; sleep 1; done

--running crest for HSD_HSP,HSE_HSP and HSD_HSE stripped N,C terminus in both vaccum and solvent

# 17/10/2022

solving HW #1 for PSF course.
discussing Macromodel...
finding out that for charged flag -chrg 1 the charge was assigned to the imidazole and not to the imidazolium! problem when doing MD, 
and not calculating orbitals as in orca. 
Now for non-charged may continue using crest altough 200 conformations are revelead.

How to compare geometry to search in PDB? 
start with those in Phe paper (altough Phe is hetroatomic ring and HIS not, also problematic for residues without plane...). See:
https://reader.elsevier.com/reader/sd/pii/0014579385809820?token=E75C9BBEDB5EDF547E41B1FD6290C4C1FA8CB25F59B0E8B5DAD7043060819E1168B7A1A691B41CE93F1ECF7710455EFC&originRegion=eu-west-1&originCreation=20221117162337
The interaction between phenylalanine rings in proteins. 


Trying to fix order vs disordered to share same axis and remove the tail to find out correlations...at stats dir.

# 20/10/2022

Loading the desired geometry (identified as ligand1 and ligand2 ) into PDB structure to analyze its geometry with existing functions...
defenitions from PHE study paper cited on 17/10/2022...

To load the ligands and identify each atom we'll need to take the input as sequqnce of C,H,N,H etc and number each atom!!
script: distinguish_ligands_atoms.py works! now at '/home_d/rivka/Bioinformatics_HIS/python_codes/test_analysis/scripts. Gets input of filename and retirns the same lines with replacing the column of atom name/id with numbered atoms! now only for N,C,H...The output filename starts with 'sorted_<inputfilename>.pdb'. see pdb_files dir.
  
# 21/10/2022
  sick
  
# 22/10/2022
  sick
  
# 23/10/2022
  
  continue working on geometry analysis: first take a given pdb ligan sorted input and identify with pdb atom naming of sidechain atoms, CB, CD1 etc..
  please note this paper for parametrs choice of bond lengths: https://pubs.acs.org/doi/pdf/10.1021/j100589a006
  
  currently the code gets numerucally sorted pdb ligand HIS and returns the correspondin row coordinates of heavyatoms as histidine PDB format
  {'CA': 0, 'CB': 1, 'CG': 4, 'ND1': 5, 'CE1': 6, 'NE2': 8, 'CD2': 9}
  can also identify HIS: HSP,HSD,HSE
  
 script works for single residue change and identify but for more than one residue the line index is not the group atom index!!!
  think of using size of all previous residue identified in sorted pdb file...
  
 # 24/10/2022
  
Fixing ligands issue noted previously, see: ligand_to_residue.py at /home_d/rivka/Bioinformatics_HIS/python_codes/test_analysis
use: python ligand_to_residue.py current_minima_2HSD.pdb
find outputs at: /home_d/rivka/Bioinformatics_HIS/python_codes/test_analysis/pdb_files/
1. sorted_current_minima_2HSD.pdb
2. residuesHSD_1_HSD_2__current_minima_2HSD.pdb

  **works** yay. Next step: define geometry search parameters: https://reader.elsevier.com/reader/sd/pii/0014579385809820?token=E75C9BBEDB5EDF547E41B1FD6290C4C1FA8CB25F59B0E8B5DAD7043060819E1168B7A1A691B41CE93F1ECF7710455EFC&originRegion=eu-west-1&originCreation=20221117162337
  
 copied: The geometry and energetics of such interactions will be crucial for protein structure, specificity and activity.
 steps:
 - define interaction between two Phe rings by calculationg the shorter inter-molecular C-C distance denoted d.
   The interacting pairs are chosen as those who have d<=4.6A as the WDV distance between two aromatic C-H is 3.6A with additional 1A error.
 - run over the interacting pairs extracted denoted F1 and F2 and find:
   1. alculate 2 parameters, which are independent of the coordinate system:
      D = distance between centre of rings (A)
      P = angle between the rings planes
   2. Define a ???reference??? phenylalanine (PHE),placed with the centre of the ring at the origin 0, with the x-axis along 0-CG; the z-axis perpen-
        dicular to the ring-plane and y in the plane of the ring, orthogonal to x and z (see fig.1).
        Fig.1.:
       ![plot](./PHE_PHE_img.png)
                                                             
   3. Superpose Fl on this reference PHE, apply the same matrix to F2 and calculate the polar coor-dinates TB, Tti of the centroid of the second ring 02,       where: TBI = angle of elevation of 02 from the plane of Fl (azimuthal angle)
      Tdl = equatorial angle of 02 in the plane of Fl    
   4. Repeat step 3, using F2 to define the coordinate system and calculate TBz, T~z, the polar coor-dinates of the centroid of Fl relative to F2. Note:
      these will be different from TB;, Tq51
 

  
  
  



























