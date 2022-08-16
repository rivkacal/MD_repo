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
set scratch = /scratch/rivka/fragB_helix_translated_updated_0.5_37429
echo '#########################'
echo Running on host `hostname`
echo Time is `date`
echo Directory is /home_d/rivka/MD_base_rivka_updated/fragB_helix_translated_updated_41648_hsa_cib_1998_1998_1998_es_80_dh_0.001_1.4_1.0
echo temp Directory is $scratch
mkdir $scratch
cp /home_d/rivka/MD_base_rivka_updated/fragB_helix_translated_updated_41648_hsa_cib_1998_1998_1998_es_80_dh_0.001_1.4_1.0/input/settings_37429.dat $scratch/settings.dat
cp /home_d/rivka/MD_base_rivka_updated/fragB_helix_translated_updated_41648_hsa_cib_1998_1998_1998_es_80_dh_0.001_1.4_1.0/input/random_37429.dat $scratch/random.dat
cp ./fragB_helix_translated_updated.dat /home_d/rivka/MD_base_rivka_updated/fragB_helix_translated_updated_41648_hsa_cib_1998_1998_1998_es_80_dh_0.001_1.4_1.0/input/
cp /home_d/rivka/MD_base_rivka_updated/MD.exe $scratch
cp ./fragB_helix_translated_updated.dat $scratch
cd $scratch
/usr/bin/time -p ./MD.exe >& time_fragB_helix_translated_updated_0.5_37429.out
bzip2 TemperatureFile_fragB_helix_translated_updated_0.5_37429.dat
bzip2 2body_fragB_helix_translated_updated_0.5_37429.dat
bzip2 3body_fragB_helix_translated_updated_0.5_37429.dat
bzip2 Etotal_fragB_helix_translated_updated_0.5_37429.dat
bzip2 EbyType_fragB_helix_translated_updated_0.5_37429.dat
bzip2 fragB_helix_translated_updated_0.5_37429.log
bzip2 Traj_fragB_helix_translated_updated_0.5_37429.dat
cp Final_fragB_helix_translated_updated_0.5_37429.dat $out_dir/LastCord
cp fragB_helix_translated_updated_0.5_37429.log.bz2 $out_dir/Log
cp TemperatureFile_fragB_helix_translated_updated_0.5_37429.dat.bz2 $out_dir/Temperature
cp 2body_fragB_helix_translated_updated_0.5_37429.dat.bz2 $out_dir/2Body
cp 3body_fragB_helix_translated_updated_0.5_37429.dat.bz2 $out_dir/3Body
cp Etotal_fragB_helix_translated_updated_0.5_37429.dat.bz2 $out_dir/Etotal
cp EbyType_fragB_helix_translated_updated_0.5_37429.dat.bz2 $out_dir/EbyType
cp Traj_fragB_helix_translated_updated_0.5_37429.dat.bz2 $out_dir/Traj
rm -r $scratch
