import subprocess
import argparse
import os
import sys
import subprocess
import shutil
import datetime
from pathlib import Path
import time
import shutil
import glob
import random
import string 

# --- CONFIGURATION ---
DSI_STUDIO_BIN = "/cm/shared/workstation/aux/dsi_studio_2022-12-22/dsi_studio"

def cluster_run_cmds(cmds, args):
    # assumes cmds is a list of commands to be clusterfied and ran 
    # this will wait for all commands in the list to complete before returning
    if not args.dry_run:
        print("sleeping 5 seconds to ensure all bash stub files are created")
        time.sleep(5)
    print("running all via sbatch")
    print(cmds)
    job_ids = []
    for cmd in cmds:
        if cmd is None:
            continue
        if args.dry_run:
            print("DRY RUN. commands:")
            print("\t",cmd)
            continue
        proc = subprocess.Popen(cmd.split(" "), stdout=subprocess.PIPE)
        # stdout is a buffered reader (bytes string)
        # must first decode it, and then split on spaces to get the last value (the job id)
        jid = proc.stdout.read().decode("utf-8").split(" ")[-1]
        # to strip off the trailing newline character
        jid = jid.strip("\n")
        print("sbatch job id: {}".format(jid))
        job_ids.append(jid)
    # call the pipeline_utilities 'cluster_wait' function, pass it list of slurm ids
    # turn my job ids list into a space-delimited string
    job_ids = " ".join(job_ids)
    cluster_wait_cmd = "pipeline_utilities cluster_wait {}".format(job_ids)
    print("waiting for dsi studio calls to complete")
    subprocess.run(cluster_wait_cmd.split(" "))


# cmds is either a single command (string) or a list of commands
# if it is a list, then these are assumed to be needed to ran sequentially
# concorrent-ran commands need their own slurm job, and multiple jobs can be ran concurrently through cluster_run_cmds
def make_cluster_command(cmds, out_dir, job_name, memory="40G"):
    # this will convert a command to a cluster command
        # this is done either by srun --mem $cmd to immediately run a command
        # or sbatch --mem path/to/scrp.bash to schedule the command to run whenever
    # and then save that command to a bash file along with #!bash
    # returns the sbatch command
    #cmd = "srun --mem={} {}".format(memory, cmd)
    bash_wrapper = "{}/{}-{}.bash".format(out_dir, job_name, time.time())

    if type(cmds) is not list:
        cmds = [cmds]

    with open(bash_wrapper, "w") as f:
        f.write("#!/usr/bin/env bash\n")
        for cmd in cmds:
            f.write(cmd)
            f.write("\n")
    slurm_out_file=f"{out_dir}/slurm_{job_name}.out"
    return "sbatch --mem={} --output={} {}".format(memory, slurm_out_file, bash_wrapper)


# all optional arguments have defaults set to pipeline parameters
# you are free to override
def precompute_tractography(runno, biggus_diskus: Path, seed_count=2000000, fa_thresh_pct=10, step_size=0.01, smoothing=0.01, min_length=0.5, max_length=200, turn_angle=45, name_tag=""):
    # --- Setup Paths ---
    work_dir = biggus_diskus / "filtered_tracking" / f"tracking{runno}dsi_studio{name_tag}-work"
    bash_stub_dir = work_dir / 'bash_stub'
    for x in [work_dir, bash_stub_dir]:
        try:
            x.mkdir(parents=True, exist_ok=False)
        except FileExistsError:
            print(f"{x} already exists")
    
    # Input paths
    in_dir = biggus_diskus / "connectome_cache"
    # TODO: pathlib
    fib_file = glob.glob(in_dir / f"nii4D_{runno}*.fib.gz")[0]
    # unused in this function
    label_file = in_dir / "labels" / "RCCF" / f"{runno}_RCCF_labels.nii.gz"
    log_file = work_dir / f"pipe_log-{datetime.datetime.now().strftime('%Y%m%d%H%M')}.log"

    thresh_file = in_dir / f"{runno}_threshold_at_{fa_thresh_pct}pct_nqa.txt"
    if thresh_file.exists():
        with open(thresh_file, 'r') as f:
            fa_thresh_pct = f.read().strip()
    else:
        # then we need to calculate threshold ourself
        # do not worry about parallel here, as only doing a couple
        import nibabel as nib
        import numpy as np
        img = glob.glob(in_dir / f"nii4D_{runno}*.nqa.nii.gz")[0]
        data = nib.load(img).get_fdata()
        data=data[data!=0]
        data.sort()
        fa_thresh_pct=data[round((fa_thresh_pct/100)*len(data))]
        with open(thresh_file,'w') as f:
            print(fa_thresh_pct, file=f)

    full_trk_file = work_dir / f"nii4D_{runno}.src.gqi.0.9.fib.2000K.tt.gz"
            
    if not full_trk_file.exists():
        cmd = (
            f"{DSI_STUDIO_BIN} --action=trk --source={fib_file} "+
            f"--output={full_trk_file} --fiber_count={seed_count} "+
            f"--fa_threshold={fa_thresh_pct} --step_size={step_size} --threshold_index=qa "+
            f"--smoothing={smoothing} --min_length={min_length:.1f} --thread_count=16 "+
            f"--method=0 --turning_angle={turn_angle} --max_length={max_length:.1f}"
        )
        
        return make_cluster_command(cmd,bash_stub_dir,f"{runno}_tracking")
    else:
        print(f"{full_trk_file} already exists.")
        return None


def run_both_sides(runno, roi_tuple, biggus_diskus: Path, dry_run=False, run_both_sides_override=False, skip_SAMBA_copy=False, name_tag=""):
    if len(roi_tuple)==2 and not abs(roi_tuple - roi_tuple[1]) - offset:
        # we only reach this if abs(roi_tuple - roi_tuple[1]) - offset is exactly 0
        # this means we have 2 regions, and they are the same on the opposite side of thebrian
        # example 47,1047
        # the contralateral pair 1047,47 is identical. No need to run the otherside
        run_both_sides_override = True
    offset = 1000
    roi_tuple_otherside = []
    cmds = []
    for roi in roi_tuple:
        x = roi+offset if roi<offset else roi-offset
        roi_tuple_otherside.append(x)
    roi_tuple_otherside = tuple(roi_tuple_otherside)
    cmds.append(setup_pipeline(runno, roi_tuple, biggus_diskus, dry_run=dry_run, skip_SAMBA_copy=skip_SAMBA_copy, name_tag=name_tag))
    if run_both_sides_override:
        return cmds
    cmds.append(setup_pipeline(runno, roi_tuple_otherside, biggus_diskus, dry_run=dry_run, skip_SAMBA_copy=skip_SAMBA_copy, name_tag=name_tag))
    return cmds

# this handles all logic for creating one filtered track and tdi file
# filters by roi1, then filters by roi2, then exports to tdi/tdi_color
# the end result of this will be ONE cmd, to be returned and added to cmds
def setup_pipeline(runno, roi_tuple, biggus_diskus: Path, project_code="24.chdi.01", dry_run=False, skip_SAMBA_copy=False, name_tag=""):
    # --- Setup Paths ---
    work_dir = biggus_diskus / "filtered_tracking" / f"tracking{runno}dsi_studio{name_tag}-work"
    results_dir = biggus_diskus / "filtered_tracking" / f"tracking{runno}dsi_studio{name_tag}-results"
    results_nhdr_dir = results_dir / 'nhdr'
    bash_stub_dir = work_dir / 'bash_stub'
    for x in [work_dir, bash_stub_dir, results_dir, results_nhdr_dir]:
        try:
            x.mkdir(parents=True, exist_ok=False)
        except FileExistsError:
            print(f"{x} already exists")

    
    # Input paths
    in_dir = biggus_diskus / "connectome_cache"
    fib_file = glob.glob(in_dir / f"nii4D_{runno}*.fib.gz")[0]
    label_file = in_dir / f"{runno}_RCCF_labels.nii.gz"
    full_trk_file = work_dir / f"nii4D_{runno}.src.gqi.0.9.fib.2000K.tt.gz"

    cmds = []
    tract_file = full_trk_file
    roi_string = ""
    roi_string_for_names = ""
    for i in range(len(roi_tuple)):
        roi = roi_tuple[i]
        if i == 0:
            roi_string_for_names=roi
        else:
            roi_string_for_names=f"{roi_string_for_names}_{roi}"
        roi_string = f"{label_file}:{roi}"
        out_file = work_dir / f"{runno}_{roi_string_for_names}.tt.gz"
        if i == len(roi_tuple)-1:
            # the last one (which we explort tdi,tdi_color files from) is saved into the results dir
            # this is the file out_file will be set to for the rest of this function, after the for loop, as well
            # IS THIS A SAFE ASSUMPTION?
            # testing proves out_file stays in the scope after we exit the loop. i ddin't expect that. apparantly it's fine
            # https://stackoverflow.com/questions/3611760/scoping-in-python-for-loops
            out_file = results_dir / f"{runno}_{roi_string_for_names}.tt.gz"
        cmd = f"{DSI_STUDIO_BIN} --action=ana --source={fib_file} --tract={tract_file} --roi={roi_string} --output={out_file}"
        cmds.append(cmd) if not out_file.exists() else cmds.append(f"#{cmd}")

    # export tdi and tdi_color
    cmd = f"{DSI_STUDIO_BIN} --action=ana  --source={fib_file} --tract={out_file} --export=tdi,tdi_color"
    # file.with_name(x) returns dirname(file)/x
    cmds.append(cmd) if not out_file.with_name(f"{out_file.name}.tdi_color.nii.gz").exists() else cmds.append(f"#{cmd}")


    # create the tdi.nhdr and tdi_color.nhdr which describe the nifti files
    # done by copying an example from the archive and adjusting the data file path
    # use sed for this. remember this is a bash cluster process, not python code that will be running 
    for contrast in ["tdi","tdi_color"]:
        nhdr_template =  in_dir / f"{runno}_{contrast}.nhdr"
        out_nhdr = results_dir / f"nhdr/{runno}_{roi_string_for_names}_{contrast}{name_tag}.nhdr"
        cmd = f"cp {nhdr_template} {out_nhdr}"
        cmds.append(cmd) if not out_nhdr.exists() else cmds.append(f"#{cmd}")
        new_filename = out_file.name
        new_filename = f"../{new_filename}.{contrast}.nii.gz"
        cmd = f'sed -i "s|data file.*|data file: {new_filename}|" {out_nhdr}'
        cmds.append(cmd)

    skip_color_split = True
    if skip_SAMBA_copy or skip_color_split:
        return make_cluster_command(cmds,bash_stub_dir,f"{runno}_filter_and_nrrdify")
    
    # split tdi_color using matlab
    timestamp = "_".join(str(time.time()).split("."))
    mat_script = "{}/run_from_python_{}.m".format(bash_stub_dir, timestamp)
    Path(mat_script).touch()



    # split the tdi_color file into component channels
    nhdr_dir = results_dir / 'nhdr'
    name = f'{runno}_{roi_string_for_names}_tdi_color.nhdr'
    purpose = f"{runno}_{roi_string_for_names}_tdi_color_split"
    in_file = nhdr_dir / name
    out_base = nhdr_dir / name.removesuffix(".nhdr")
    mat_code=f"i='{in_file}';o='{out_base}';image_channel_split(i,o);";
    with open(mat_script, 'w') as f:
        f.write("run('/cm/shared/workstation/code/shared/pipeline_utilities/startup.m');")
        f.write(mat_code)

    cmd = "\"run('{}'); exit;\"".format(mat_script)
    cmd = f"matlab_run {cmd} --purpose={purpose} --dir_work={work_dir}"
    # checking for existing split channel colors here
    found_channels = glob.glob(nhdr_dir / f'{name.removesuffix(".nhdr")}_*.nhdr')
    cmds.append(cmd) if not len(found_channels) > 2 else cmds.append(f"#{cmd}")

    # copy the split channel colors into the SAMBA inputs directory
    # due to a logical glitch with split-channel colors, must copy into 24.chdi.01_nhdr and VBM-work/clean_inputs
    # otherwise SAMBA will try to force to split channels to get into clean_inputs
    # just because 'color' is in the name
    # any file which has pattern color_*.[nhdr,raw] is a color component\

    project_nhdr_dir = biggus_diskus / f'{project_code}_nhdr'
    SAMBA_work_dir = biggus_diskus / f'VBM_24chdi01_chass_symmetric5-work'
    clean_inputs_dir = SAMBA_work_dir / 'clean_inputs'
    clean_masked_dir = SAMBA_work_dir / 'clean_masked'
    for filepath in glob.glob(nhdr_dir / f'{name.removesuffix(".nhdr")}_*'):
        fn = filepath.name
        # put it into the nhdr dir
        outf = project_nhdr_dir / fn
        cmd = f"cp {filepath} {outf}"
        cmds.append(cmd) if not outf.exists() else cmds.append(f"#{cmd}")

        # put it into clean_inputs
        outf = clean_inputs_dir / fn
        cmd = f"cp {filepath} {outf}"
        cmds.append(cmd) if not outf.exists() else cmds.append(f"#{cmd}")

        # put it into clean_masked
        # for clean masked(nhdr only) i also need to rename it
        if filepath.endswith('.nhdr'):
            outf = clean_masked_dir / fn
            outf_name_corrected = clean_masked_dir / f"{fn.name.removesuffix('.nhdr')}_masked.nhdr"
            # only do the copy if the RENAMED destination file does not exist
            # otherwise we will make a mess and have both renamed and unrenamed files in clean_masked
            cmd = f"cp {filepath} {outf}"
            cmds.append(cmd) if not outf_name_corrected.exists() else cmds.append(f"#{cmd}")
            cmd = f"mv {outf} {outf_name_corrected}"
            cmds.append(cmd) if not outf_name_corrected.exists() else cmds.append(f"#{cmd}")
        else:
            # for the raw files, i do not want to rename their filenames
            # else i would also need to update data file: field in the nhdr
            outf = clean_masked_dir / fn
            cmd = f"cp {filepath} {outf}"
            cmds.append(cmd) if not outf.exists() else cmds.append(f"#{cmd}")

    return make_cluster_command(cmds,bash_stub_dir,f"{runno}_filter_and_nrrdify")


# TODO: this is all too baked into assumptions
# use glob/pa3ttern matching to find your input files. 
def prepull_data(project_code,runno,biggus_diskus,archive_suffix=""):
    print(f"Pulling data for {runno}")
    in_dir = Path(f"/mnt/nclin-comp-pri.dhe.duke.edu/dusom_civm-atlas/{project_code}/research", f"connectome{runno}dsi_studio{archive_suffix}")
    # what do i do if I find more than one? Am I likely to have that happen? Do I even care? Not sure. 
    fib_files = glob.glob(in_dir / f"nii4D_{runno}*.fib.gz")
    if len(fib_files) != 1:
        print(f"found {len(fib_files)} fib files for {runno}. skipping")
        return None
    fib_file = fib_files[0]
    qa_niis = glob.glob(in_dir / f"nii4D_{runno}*.nqa.nii.gz")
    qa_nii = qa_niis[0] if len(qa_niis) >= 1 else None
    # TODO: label TYPE (gaj projects needs WHS labels)
    label_file = in_dir / "labels" / "RCCF" / f"{runno}_RCCF_labels.nii.gz"
    tdi_nhdr = in_dir / 'nhdr' / f"{runno}_tdi.nhdr"
    tdi_color_nhdr = in_dir / 'nhdr' / f"{runno}_tdi_color.nhdr"
    files_to_copy = [fib_file, label_file, tdi_nhdr, tdi_color_nhdr, qa_nii]

    out_dir = Path(biggus_diskus,"connectome_cache")
    try:
        out_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print(f"{out_dir} already exists")
    # necessary items are fib_file, label_file, tdi_nhdr, tdi_color_nhdr
    for in_file in files_to_copy:
        fn = in_file.name
        out_file = out_dir / fn
        shutil.copyfile(in_file,out_file) if not out_file.exists() else print(f"out_file already copied")

    # also copy pct threshold file, but it must bee renamed bc no mention of runno in its file name 
    in_file = in_dir / 'threshold_at_10pct_nqa.txt'
    fn = in_file.name
    fn = f"{runno}_{fn}"
    out_file = out_dir / fn
    shutil.copyfile(in_file,out_file) if not out_file.exists() else print(f"out_file already copied")


def load_list_file(list_file):
    with open(list_file,'r') as f:
        return f.read().strip().split('\n')


# this function is quite specialized to the 24.chdi.01 project structure
def load_list_files(ages, project_code="24.chdi.01"):
    runno_list = []
    for age in ages:
        for condition in ['HET','WILD']:
            for sex in ['M','F']:
                homedir=os.getenv("HOME")
                list_file = f"{homedir}/Projects/{project_code}/list/{project_code}-{age}-{condition}-{sex}.list"
                new_runnos = load_list_file(list_file)
                runno_list.extend(new_runnos)
    return runno_list


def setup_channel_comma_list_for_samba_headfile(roi_tuple_list):
    print("SAMBA headfile helper is not set up for n-ROIs. Sorry")
    return None
    offset = 1000
    l = []
    for roi in roi_tuple_list:
        roi1 = roi[0]
        roi2 = roi[1]
        # excluding tdi uncolored because it crashes samba...unsure why
        #l.append(f'{roi1}_{roi2}_tdi')
        for c in ['red','green','blue']:
            l.append(f'{roi1}_{roi2}_tdi_color_{c}')
        # other side
        if not abs(roi1 - roi2) - offset:
            # other side not necessary
            continue
        roi1 = roi1+1000 if roi1<1000 else roi1-1000
        roi2 = roi2+1000 if roi2<1000 else roi2-1000
        l.append(f'{roi1}_{roi2}_tdi')
        for c in ['red','green','blue']:
            l.append(f'{roi1}_{roi2}_tdi_color_{c}')
    channel_comma_list = ','.join(l)
    print(f'\t\tchannel_comma_list')
    print(f'channel_comma_list={channel_comma_list}')
    return channel_comma_list

def tuple_type(s):
    s = s.strip("()")
    try:
        return tuple(map(int,s.split(",")))
    except ValueError:
        argparse.ArgumentTypeError(f"Expected format '(x,y)' or 'x,y', but got: {s}")


def main():
    parser = argparse.ArgumentParser(description="Create connectome-filtered TDI_color files for SAMBA inputs")
    parser.add_argument("--runno_list","-r",nargs="*",type=str,help='The path to a list file, OR a list of runnos to operate on, separated by spaces')
    parser.add_argument("--project_code","-p",type=str)
    parser.add_argument("--label_type","-l",type=str,default="RCCF")
    parser.add_argument("--archive_suffix", type=str,default="")
    parser.add_argument("--biggus_diskus", type=Path,default=os.getenv("BIGGUS_DISKUS"))
    parser.add_argument("--roi_tuple_list", nargs="+",type=tuple_type,help='Pass a list of tuples formatted as (x,y)')
    parser.add_argument("--one_side_only",action="store_true",help="forces only calculation of asked for side, will not flip and operate",default=False)
    parser.add_argument("--dry_run",action="store_true",help="will only print the commands to run instead of running",default=False)
    # tractography parameters
    parser.add_argument("--seed_count",type=int,default=2000000)
    parser.add_argument("--fa_thresh_pct",type=int,default=10)
    parser.add_argument("--step_size",type=float,default=0.01)
    parser.add_argument("--smoothing",type=float,default=0.01)
    parser.add_argument("--min_length",type=float,default=0.5)
    parser.add_argument("--max_length",type=float,default=200)
    parser.add_argument("--turn_angle",type=int,default=45)
    parser.add_argument("--name_tag",type=str,default="")
    args = parser.parse_args()
    # find relevant list file
    if not args.project_code:
        project_code = "24.chdi.01"
        ages = [2,6,10,15]
    if not args.runno_list and project_code=="24.chdi.01":
       args.runno_list = load_list_files(ages, project_code)
    if args.runno_list:
        final_runnos = []
        for item in args.runno_list:
            file_path = Path(item)
            if file_path.is_file():
                final_runnos.extend(load_list_file(file_path))
            else:
                final_runnos.append(item)
        args.runno_list = final_runnos

    print(args.runno_list)
    
    if args.name_tag:
        # if a name tag was give, add a dash to the beginning of it, for filepath generation ease later oin
        args.name_tag = f"-{args.name_tag}"


    # short list for testing 
    #args.runno_list = ["S70132NLSAM", "S70133NLSAM", "S70135NLSAM", "S70137NLSAM", "S70139NLSAM"]
    #roi_tuple_list = [(111,1111), (15, 1015), (15, 1009), (15, 1156), (47,51), (47, 1005)]
    # testing to get through the full pipeline

    # cartesian product
    # these roi pairs decided as the full blue list from Kathryn's feb 2 2026 email for 15 MOS and 47 STD
    # this was specialized for CHDI results, before full input argument setup
    #if not args.roi_tuple_list:
    #    l1 = [(15,x) for x in [7,14,9,17,47,156,157,1005,1006,1007,1009,1014,1015,1156,1157]]
    #    l2 = [(47,x) for x in [7,8,9,15,36,66,67,74,77,78,168,1036,1066,1077,1078]]
    #    args.roi_tuple_list = l1 + l2
    setup_channel_comma_list_for_samba_headfile(args.roi_tuple_list)

    cmds = []
    # pre-create the tractography files to ensure they will be present for the rest of the pipeline
    # this is an immutable input that will always be used, much like the fib or label file
    for runno in args.runno_list:
        if not runno.startswith(("S","N")): continue
        prepull_data(args.project_code,runno,args.biggus_diskus,args.archive_suffix)
        cmd = precompute_tractography(runno, args.biggus_diskus, seed_count=args.seed_count,fa_thresh_pct=args.fa_thresh_pct,step_size=args.step_size,smoothing=args.smoothing,min_length=args.min_length,max_length=args.max_length,turn_angle=args.turn_angle, name_tag=args.name_tag)
        cmds.append(cmd)
    cluster_run_cmds(cmds, args)
    #import pdb;pdb.set_trace()
    cmds = []
    for runno in args.runno_list:
        for roi_tuple in args.roi_tuple_list:
            new_cmds = run_both_sides(runno, roi_tuple, args.biggus_diskus, args.dry_run, run_both_sides_override=args.one_side_only, skip_SAMBA_copy=True, name_tag=args.name_tag)
            # use extend instead of append, as run_both_sides* will return two cmds to run 
            cmds.extend(new_cmds)
    cluster_run_cmds(cmds, args)


if __name__ == '__main__':
    main()
