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
set out_dir = /home_d/rivka/MD_base_rivka_updated/1srl_54379_hsa_cib_1198_1198_1198_es_80_dh_0.001_1.4_1.0/output
set scratch = /scratch/rivka/1srl_0.5_3254
echo '#########################'
echo Running on host `hostname`
echo Time is `date`
echo Directory is /home_d/rivka/MD_base_rivka_updated/1srl_54379_hsa_cib_1198_1198_1198_es_80_dh_0.001_1.4_1.0
echo temp Directory is $scratch
mkdir $scratch
cp /home_d/rivka/MD_base_rivka_updated/1srl_54379_hsa_cib_1198_1198_1198_es_80_dh_0.001_1.4_1.0/input/settings_3254.dat $scratch/settings.dat
cp /home_d/rivka/MD_base_rivka_updated/1srl_54379_hsa_cib_1198_1198_1198_es_80_dh_0.001_1.4_1.0/input/random_3254.dat $scratch/random.dat
cp ./1srl.dat /home_d/rivka/MD_base_rivka_updated/1srl_54379_hsa_cib_1198_1198_1198_es_80_dh_0.001_1.4_1.0/input/
cp /home_d/rivka/MD_base_rivka_updated/MD.exe $scratch
cp ./1srl.dat $scratch
cd $scratch
/usr/bin/time -p ./MD.exe >& time_1srl_0.5_3254.out
bzip2 TemperatureFile_1srl_0.5_3254.dat
bzip2 2body_1srl_0.5_3254.dat
bzip2 3body_1srl_0.5_3254.dat
bzip2 Etotal_1srl_0.5_3254.dat
bzip2 EbyType_1srl_0.5_3254.dat
bzip2 1srl_0.5_3254.log
bzip2 Traj_1srl_0.5_3254.dat
cp Final_1srl_0.5_3254.dat $out_dir/LastCord
cp 1srl_0.5_3254.log.bz2 $out_dir/Log
cp TemperatureFile_1srl_0.5_3254.dat.bz2 $out_dir/Temperature
cp 2body_1srl_0.5_3254.dat.bz2 $out_dir/2Body
cp 3body_1srl_0.5_3254.dat.bz2 $out_dir/3Body
cp Etotal_1srl_0.5_3254.dat.bz2 $out_dir/Etotal
cp EbyType_1srl_0.5_3254.dat.bz2 $out_dir/EbyType
cp Traj_1srl_0.5_3254.dat.bz2 $out_dir/Traj
rm -r $scratch
