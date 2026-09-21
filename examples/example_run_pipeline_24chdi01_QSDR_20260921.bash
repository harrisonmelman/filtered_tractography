# based on pairwise euclidean distance of whole brain connectome MDS space
# this is the most different pair of 12 monthn specimsn
#240725-2:1 HET S70277
#250423-9:1 WILD S70289



for threshold in 1 10 20; do
python3 $WORKSTATION_CODE/diffusion/filtered_tractography/src/new_pipeline_20260921.py \
  --name_list 15mo-WT-M 15mo-zQ175DN-M \
  --input_dir /privateShares/hmm56/24.chdi.01/QSDR_filtered_tractography/input \
  --label_type RCCF \
  --roi_tuple_list  "9" "14" "15" "47" "111" \
  --biggus_diskus /privateShares/hmm56/24.chdi.01/QSDR_filtered_tractography \
  --fa_thresh_pct $threshold \
  --min_length 0.01 \
  --max_length 20.0

done
