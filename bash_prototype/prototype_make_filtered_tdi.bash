functio(){
BIGGUS_DISKUS=B:/ProjectSpace/hmm56/24.chdi.01/filtered_tracking2;
dsi_studio="/c/CIVM_Apps/dsi_studio_64/dsi_studio_win_v2022-12-22/dsi_studio.exe";
project_code=$1;
runno=$2;
roi1=$3;
roi2=$4;
out_dir=${BIGGUS_DISKUS}/tracking${runno}dsi_studio-results;
work_dir=${BIGGUS_DISKUS}/tracking${runno}dsi_studio-work;
split_color_output_dir=${BIGGUS_DISKUS}/split_colors;
mkdir -p $work_dir $out_dir/nhdr $split_color_output_dir;
cd $work_dir;
date=$(date -d "today" +"%Y%m%d%H%M");
log=pipe_log-${date}.log;
touch $log;

in_dir=A:/${project_code}/research/connectome${runno}dsi_studio;
fib_file=${in_dir}/nii4D_${runno}.src.gqi.0.9.fib.gz;
label_file=${in_dir}/labels/RCCF/${runno}_RCCF_labels.nii.gz;

# create pipeline tractography file
fa_thresh=$(cat ${in_dir}/threshold_at_10pct_nqa.txt);
full_trk_file=${work_dir}/nii4D_${runno}.src.gqi.0.9.fib.2000K;
trk_cmd="$dsi_studio --action=trk --source=${fib_file} --output=${full_trk_file}.tt.gz --fiber_count=2000000 --fa_threshold=${fa_thresh} --step_size=0.01 --threshold_index=qa --smoothing=0.01 --min_length=0.5 --thread_count=16 --method=0 --turning_angle=45 --max_length=200.0";
if [ ! -e ${full_trk_file}.tt.gz ]; then
    echo $trk_cmd >> $log;
    $trk_cmd;
else
    echo $full_trk_file already exists;
fi


tract_file=${full_trk_file};
out_tmp="${work_dir}/filtered_temp_${roi1}";
roi_string=${label_file}:${roi1};
cmd="${dsi_studio} --action=ana --source=${fib_file} --tract=${tract_file}.tt.gz --roi=${roi_string} --output=${out_tmp}.tt.gz";
# set the output of this as the input to the next step
tract_file=${out_tmp};
if [ ! -e ${out_tmp}.tt.gz ]; then
    echo $cmd >> $log;
    $cmd;
else
    echo $out_tmp already exists;
fi


out_finally="${out_dir}/${runno}_${roi1}_${roi2}";
roi_string=${label_file}:${roi2};
cmd="${dsi_studio} --action=ana --source=${fib_file} --tract=${tract_file}.tt.gz --roi=${roi_string} --output=${out_finally}.tt.gz";
if [ ! -e ${out_finally}.tt.gz ]; then
    echo $cmd >> $log;
    $cmd;
else
    echo $out_finally already exists;
fi
# do export command
export_cmd="${dsi_studio} --action=ana  --source=${fib_file} --tract=${out_finally}.tt.gz --export=tdi,tdi_color";
if [ ! -e ${out_finally}.tt.gz.tdi_color.nii.gz ]; then
    echo $export_cmd >> $log;
    $export_cmd;
else
    echo ${out_finally}.tt.gz.tdi_color.nii.gz already exists;
fi


# copy over template nhdr files
tdi_template=${in_dir}/nhdr/${runno}_tdi.nhdr;
tdi_color_template=${in_dir}/nhdr/${runno}_tdi_color.nhdr;
out_tdi_nhdr=${out_dir}/nhdr/${runno}_${roi1}_${roi2}_tdi.nhdr
out_tdi_color_nhdr=${out_dir}/nhdr/${runno}_${roi1}_${roi2}_tdi_color.nhdr

cmd="cp $tdi_template $out_tdi_nhdr";
echo $cmd >> $log; $cmd;
cmd="cp $tdi_color_template $out_tdi_color_nhdr";
echo $cmd >> $log; $cmd;

# print the sed command for posgterity, then actually run it
# due to the interior qhotes for sed, i couldn't figure out how to save the command to a variable
echo 'sed -i "s|data file.*|data file: ../$(basename ${out_finally}).tt.gz.tdi.nii.gz|" $out_tdi_nhdr' >> $log;
echo 'sed -i "s|data file.*|data file: ../$(basename ${out_finally}).tt.gz.tdi_color.nii.gz|" $out_tdi_color_nhdr' >> $log;
sed -i "s|data file.*|data file: ../$(basename ${out_finally}).tt.gz.tdi.nii.gz|" $out_tdi_nhdr;
sed -i "s|data file.*|data file: ../$(basename ${out_finally}).tt.gz.tdi_color.nii.gz|" $out_tdi_color_nhdr;

# split the tdi_color file into component channels
input_package=${out_dir}/nhdr;
cd $split_color_output_dir;
inc=$out_tdi_color_nhdr;
n="$(basename $inc)";
# looking for outputs here
found_channels=$(find $split_color_output_dir -maxdepth 1 -type f -name "${n%.*}_*nhdr" | wc -l)
outc="$split_color_output_dir/${n%.*}";
if [ $found_channels -ge 3 ];then
    echo "Found $found_channels of output for $n, skipping";
    continue;
fi;
echo "run split channel on $n";
# Quit was part of code, but finally remembered that'd blind it. quit('force');
mat_code="i='$inc';o='$outc';image_channel_split(i,o);";
echo "running $mat_code";
matlab_run --purpose color_split --dir-work "split_${n%.*}-work" "$mat_code" &

}


# this assumes you pass two rois, each on the left side of the brain
run_both_sides_ipsilateral(){
project_code=$1;
runno=$2;
roi1=$3;
roi2=$4;
functio $project_code $runno $roi1 $roi2&
# add 1000 to get the roi number to the right side of the brain
roi1=$(($roi1 + 1000));
roi2=$(($roi2 + 1000));
functio $project_code $runno $roi1 $roi2&
}

# this assumes you pass two rois, each on the opposite side of the brain
# the first is expected to be left sided
# the second is to be right sided
run_both_sides_contralateral(){
project_code=$1;
runno=$2;
roi1=$3;
roi2=$4;
functio $project_code $runno $roi1 $roi2&
# add 1000 to get the roi number to the right side of the brain
roi1=$(($roi1 + 1000));
roi2=$(($roi2 - 1000));
functio $project_code $runno $roi1 $roi2&
}

run_for_age_group(){
project_code=$1;
age=$2;
echo looping in age group $age;
list_file=/b/ProjectSpace/hmm56/Projects/${project_code}/list/${project_code}-${age}.list;

# MANUAL OVERRIDE
#redo_list="S70132NLSAM S70133NLSAM S70135NLSAM S70137NLSAM S70222NLSAM S70198NLSAM S70202NLSAM S70204NLSAM S70163NLSAM S70167NLSAM S70169NLSAM S70212NLSAM S70250NLSAM S70248NLSAM";
#for runno in $redo_list; do

# running the first one separate to ensure that each runno has its tractography file
# then the next for loop can batch setup all the filter and export jobs
for runno in $(cat $list_file); do
    if [[ ! $runno == S* ]]; then continue; fi
    echo operating on $runno;
    # MRN__Midbrain_reticular_formation_left (111) to MRN__Midbrain_reticular_formation_right (1111)
    # only significant contralateral connection for 111. right to left would be identical, so no swapping required
    functio $project_code $runno 111 1111&
done
sleep 1 # For sequential output
echo "Waiting for processes to finish"
wait $(jobs -p)
echo "All processes finished"

#for runno in $(cat $list_file); do
for runno in $redo_list; do
    if [[ ! $runno == S* ]]; then continue; fi
    echo operating on $runno;

    #MOS__Motor_cortex_secondary_M2_left (15) to MOS__Motor_cortex_secondary_M2_right (1015)
    functio $project_code $runno 15 1015&
    # MOS__Motor_cortex_secondary_M2_left (15) to ACC__Anterior_cingulate_cortex_right (1009) contralateral
    run_both_sides_contralateral $project_code $runno 15 1009&
    # MOS__Motor_cortex_secondary_M2_left (15) to cca__corpus_callosum_right (1156) contralateral
    run_both_sides_contralateral $project_code $runno 15 1156&
    # STD__Striatum_dorsal_left  (47) to GPE__Globus_pallidus_external_left (51) ipsilateral
    run_both_sides_ipsilateral $project_code $runno 47 51&
    # STD__Striatum_dorsal_left (47) to cst__corticospinal_tract_right (1168) contralateral
    run_both_sides_contralateral $project_code $runno 14 1005&

    sleep 1 # For sequential output
    echo "Waiting for processes to finish"
    wait $(jobs -p)
    echo "All processes finished"
done
}
project_code=24.chdi.01;
#run_for_age_group $project_code 2 ;
#run_for_age_group $project_code 6 ;
#run_for_age_group $project_code 10 ;
#run_for_age_group $project_code 15 ;



#functio $project_code S70226NLSAM 15 1015&
#functio $project_code S70218NLSAM 1047 1051&
#for runno in S70248NLSAM S70254NLSAM S70236NLSAM S70240NLSAM S70258NLSAM S70218NLSAM S70244NLSAM; do
#    functio $project_code $runno 14 1005&
#done
#for runno in S70248NLSAM S70254NLSAM S70212NLSAM S70214NLSAM S70216NLSAM S70256NLSAM S70258NLSAM S70218NLSAM; do
#    functio $project_code $runno 1014 5&
#done
