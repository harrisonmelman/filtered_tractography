# this should be plugged into the prototype_make_filtered_tdi pipeline,
# these outputs will end up being the final inputs to SAMBA
# want to pass SAMBA only the split-channel colors

project_code=24.chdi.01;
#runno=S70153NLSAM;
#roi=111;
age=15;
list_file=/b/ProjectSpace/hmm56/Projects/${project_code}/list/${project_code}-${age}.list;

for runno in $(cat $list_file); do
    if [[ ! $runno == S* ]]; then continue; fi
    echo operating on $runno;
    # input_package is the nhdr dir
    input_package=B:/ProjectSpace/hmm56/${project_code}/filtered_tracking/tracking${runno}dsi_studio-results/nhdr;
    cd $input_package;
    for roi in 15 111 1009; do
        inc=${input_package}/${runno}_${roi}__tdi_color.nhdr;
        if [ ! -e $inc ]; then
            echo $inc does not yet exist;
            continue;
        fi
        n="$(basename $inc)";
        found_channels=$(find $input_package -maxdepth 1 -type f -name "${n%.*}_*nhdr" | wc -l)
        outc="$input_package/${n%.*}";
        if [ $found_channels -ge 3 ];then
                echo "Found $found_channels of output for $n, skipping";
                continue;
        fi;
        echo "run split channel on $n";
        # Quit was part of code, but finally remembered that'd blind it. quit('force');
        mat_code="i='$inc';o='$outc';image_channel_split(i,o);";
        echo "running $mat_code";
        matlab_run --purpose color_split --dir-work "split_${n%.*}-work" "$mat_code" &
    done
done




