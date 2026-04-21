project_code=20.5xfad.01
# N60171NLSAM is ntg extreme
# N60151NLSAM is tg extreme
runno_list="N60171NLSAM N60151NLSAM"
roi_list="28";
DRY_RUN="";
DRY_RUN="--dry_run ";

# run with pipeline parameters
python $WORKSTATION_CODE/diffusion/filtered_tractography/src/prototype_pipeline.py $DRY_RUN --project_code $project_code --runno_list $runno_list --roi_tuple_list $roi_list --name_tag "TEST" --one_side_only
