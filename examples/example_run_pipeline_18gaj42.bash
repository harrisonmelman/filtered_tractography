# N57372 Young BXD48a Male
# N57261 Old BXD48a Male
# regions of interest
# cca 156
# fim 161
# ACC 9
# INP 137
# ICC 104
# VCP 26
# MOP 14
# SSB 17
# ICC 104
# SUB 28
# CLT 71
python $WORKSTATION_CODE/diffusion/filtered_tractography/src/prototype_pipeline.py --dry_run --project_code 18.gaj.42 --runno_list N57372 N57261 --roi_tuple_list "156" "161" "9" "137" "104" "26" "14" "17" "104" "28" "71" --archive_suffix RCCF --biggus_diskus ${BIGGUS_DISKUS}/filtered_tracking