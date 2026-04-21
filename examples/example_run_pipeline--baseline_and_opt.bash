project_code=20.5xfad.01
# N60171NLSAM is ntg extreme
# N60151NLSAM is tg extreme
runno_list="N60171NLSAM N60151NLSAM"
roi_pair_list_optimized='28,46 28,1028 28,54 28,1054 28,1029 28,42 28,29 28,1046 28,1 28,32 28,1031 28,43 28,1032 28,161 28,1001';

roi_pair_list_baseline='28,1028 28,54 28,1054 28,1029 28,42 28,29 28,1046 28,1 28,32 28,1031 28,43 28,1032 28,161 28,1001';
DRY_RUN="";
DRY_RUN="--dry_run ";

# run with pipeline parameters
python $WORKSTATION_CODE/diffusion/filtered_tractography/src/prototype_pipeline.py $DRY_RUN --project_code $project_code --runno_list $runno_list --roi_pair_list $roi_pair_list_baseline --name_tag "base" --one_side_only

# run with custom parameters
python $WORKSTATION_CODE/diffusion/filtered_tractography/src/prototype_pipeline.py $DRY_RUN --project_code $project_code --runno_list $runno_list --roi_pair_list $roi_pair_list_optimized --name_tag "opt" --one_side_only --seed_count 2000000 --fa_thresh_pct 54 --step_size 0.01 --smoothing 0.01 --min_length 4.348 --max_length 13.0 --turn_angle 15
