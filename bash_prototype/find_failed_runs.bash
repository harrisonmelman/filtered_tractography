. ./prototype_make_filtered_tdi.bash

search_age_group(){
project_code=24.chdi.01;
age=$1;
pair_list=$2;
pair_list=$(echo $pair_list|sed 's/,/ /g' );
list_file=/b/ProjectSpace/hmm56/Projects/${project_code}/list/${project_code}-${age}.list;
search_dir=/b/ProjectSpace/hmm56/24.chdi.01/filtered_tracking2;

# first check the tracking work dir that the whole brain trk file exists
for runno in $(cat $list_file); do
    if [[ ! $runno == S* ]]; then continue; fi
    echo searching $runno;
    work_dir=${search_dir}/tracking${runno}dsi_studio-work;
    results_dir=${search_dir}/tracking${runno}dsi_studio-results;
    split_color_dir=${search_dir}/send_to_cluster_20260120;
    tt_file=${work_dir}/nii4D_${runno}.src.gqi.0.9.fib.2000K.tt.gz;


    for pair in $pair_list; do
        blue=${split_color_dir}/${runno}_${pair}_tdi_color_blue.nhdr;
        if [ ! -e $blue ]; then
            echo -e '\t\t'MISSING  $blue;
        fi
    done
    echo;
    continue;


    # thankfully, none of these are missing
    if [ ! -e $tt_file ]; then echo -e '\t'MISSING $tt_file; continue; fi
    # now check for all double pair filtered tracks, in the results directory
    for pair in $pair_list; do
        #echo -e '\t'searching ${runno}_${pair};
        tt_file=${results_dir}/${runno}_${pair}.tt.gz;
        #ls $tt_file;
        if [ ! -e $tt_file ]; then
            echo -e '\t\t'MISSING $tt_file;
            # rerun this one
            #roi_split=$(echo $pair|sed 's/_/ /');
            #functio 24.chdi.01 $runno $roi_split;
            #exit 1;
        fi
    done
    echo;
done


}

pair_list="15_1015,111_1111,15_1009,1015_9,15_1156,1015_156,47_51,1047_1051,14_1005,1014_5";

for x in 2 6 10 15; do
    echo -e '\t\t\t'search age group $x;
    search_age_group $x $pair_list;
done
