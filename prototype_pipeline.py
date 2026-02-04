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

# the files in slop are just taht, slop
# too confusing surfing around all that
# try again...from scratch

# --- CONFIGURATION ---
PROJECT_CODE = "24.chdi.01"
# Path to the list file (if you decide to read from file)
# overall list file
#LIST_FILE = f"/home/hmm56/Projects/24.chdi.01/list/24.chdi.01-phase2-02061015.list"
# small test list
LIST_FILE = f"/home/hmm56/Projects/24.chdi.01/list/24.chdi.01-15-WILD-F.list" 
# ABSOLUTE path to the script created above
# unsure if i will still use this format schene
WORKER_SCRIPT = "/cm/shared/workstation/code/diffusion/filtered_tractography/dsi_task.py" 

BIGGUS = os.environ["BIGGUS_DISKUS"]
ARCHIVE_ROOT = Path("/mnt/nclin-comp-pri.dhe.duke.edu/dusom_civm-atlas/24.chdi.01/research")
CONNECTOME_CACHE = Path(f"{BIGGUS}/filtered_tracking/connectome_cache")
DSI_STUDIO_BIN = "/cm/shared/workstation/aux/dsi_studio_2022-12-22/dsi_studio"
MATLAB_BIN = "matlab"

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


def precompute_tractography(runno):
    # --- Setup Paths ---
    work_dir = os.path.join(BIGGUS,"filtered_tracking",f"tracking{runno}dsi_studio-work")
    bash_stub_dir = os.path.join(work_dir,'bash_stub')
    os.makedirs(work_dir,exist_ok=True) if not os.path.exists(work_dir) else print(f"{work_dir} already exists")
    os.makedirs(bash_stub_dir,exist_ok=True) if not os.path.exists(bash_stub_dir) else print(f"{bash_stub_dir} already exists")
    
    # Input paths
    in_dir = CONNECTOME_CACHE
    fib_file = os.path.join(in_dir,f"nii4D_{runno}.src.gqi.0.9.fib.gz")
    label_file = os.path.join(in_dir,"labels","RCCF",f"{runno}_RCCF_labels.nii.gz")
        
    log_file = os.path.join(work_dir,f"pipe_log-{datetime.datetime.now().strftime('%Y%m%d%H%M')}.log")

    fa_thresh = "0.1"
    thresh_file = os.path.join(in_dir,"threshold_at_10pct_nqa.txt")
    if os.path.exists(thresh_file):
        with open(thresh_file, 'r') as f:
            fa_thresh = f.read().strip()

    full_trk_file = os.path.join(work_dir,f"nii4D_{runno}.src.gqi.0.9.fib.2000K.tt.gz")
            
    if not os.path.exists(full_trk_file):
        cmd = (
            f"{DSI_STUDIO_BIN} --action=trk --source={fib_file} "+
            f"--output={full_trk_file} --fiber_count=2000000 "+
            f"--fa_threshold={fa_thresh} --step_size=0.01 --threshold_index=qa "+
            "--smoothing=0.01 --min_length=0.5 --thread_count=16 "+
            "--method=0 --turning_angle=45 --max_length=200.0"
        )
        return make_cluster_command(cmd,bash_stub_dir,f"{runno}_tracking")
    else:
        print(f"{full_trk_file} already exists.")
        return None


def run_both_sides(runno, roi1, roi2, dry_run=False):
    offset = 1000
    roi1 = int(roi1)
    roi2 = int(roi2)
    if not abs(roi1 - roi2) - offset:
        # then we have the same region on both sides, only one run is needed
        cmd = setup_pipeline(runno, roi1, roi2, dry_run)
        return [cmd]
    cmd1 = setup_pipeline(runno, roi1, roi2, dry_run)

    roi1 = roi1+1000 if roi1<1000 else roi1-1000
    roi2 = roi2+1000 if roi2<1000 else roi2-1000
    cmd2 = setup_pipeline(runno, roi1, roi2, dry_run)
    return [cmd1, cmd2]
    

# this handles all logic for creating one filtered track and tdi file
# filters by roi1, then filters by roi2, then exports to tdi/tdi_color
# the end result of this will be ONE cmd, to be returned and added to cmds
def setup_pipeline(runno, roi1, roi2, dry_run=False):
    # --- Setup Paths ---
    work_dir = os.path.join(BIGGUS,"filtered_tracking",f"tracking{runno}dsi_studio-work")
    results_dir = os.path.join(BIGGUS,"filtered_tracking",f"tracking{runno}dsi_studio-results")
    results_nhdr_dir = os.path.join(results_dir, 'nhdr')
    bash_stub_dir = os.path.join(work_dir,'bash_stub')
    os.makedirs(work_dir,exist_ok=True) if not os.path.exists(work_dir) else print(f"{work_dir} already exists")
    os.makedirs(results_dir,exist_ok=True) if not os.path.exists(results_dir) else print(f"{results_dir} already exists")
    os.makedirs(bash_stub_dir,exist_ok=True) if not os.path.exists(bash_stub_dir) else print(f"{bash_stub_dir} already exists")
    os.makedirs(results_nhdr_dir,exist_ok=True) if not os.path.exists(results_nhdr_dir) else print(f"{results_nhdr_dir} already exists")

    
    # Input paths
    in_dir = CONNECTOME_CACHE
    fib_file = os.path.join(in_dir,f"nii4D_{runno}.src.gqi.0.9.fib.gz")
    label_file = os.path.join(in_dir,f"{runno}_RCCF_labels.nii.gz")
    
    # for REASONS, we probably want to exclude .tt.gz at the end of the trk file name strings
    # i could also use this to remove it when needed: 
    # file.removesuffix(".tt.gz")
    # I FAIL to see any solid logical reason for why i handled it the way i did in the bash code
    # for now, go on with this exactly as is
    full_trk_file = os.path.join(work_dir,f"nii4D_{runno}.src.gqi.0.9.fib.2000K.tt.gz")

    cmds = []
    # filter through roi 1
    tract_file = full_trk_file
    out_tmp = f"{work_dir}/filtered_temp_{roi1}.tt.gz";
    roi_string = f"{label_file}:{roi1}"
    cmd = f"{DSI_STUDIO_BIN} --action=ana --source={fib_file} --tract={tract_file} --roi={roi_string} --output={out_tmp}"
    cmds.append(cmd) if not os.path.exists(out_tmp) else cmds.append(f"#{cmd}")

    # filter through roi2
    tract_file = out_tmp
    out_finally = os.path.join(results_dir, f"{runno}_{roi1}_{roi2}.tt.gz")
    roi_string = f"{label_file}:{roi2}"
    cmd = f"{DSI_STUDIO_BIN} --action=ana --source={fib_file} --tract={tract_file} --roi={roi_string} --output={out_finally}"
    cmds.append(cmd)if not os.path.exists(out_finally) else cmds.append(f"#{cmd}")

    # export tdi and tdi_color
    cmd = f"{DSI_STUDIO_BIN} --action=ana  --source={fib_file} --tract={out_finally} --export=tdi,tdi_color"
    cmds.append(cmd) if not os.path.exists(f"{out_finally}.tdi_color.nii.gz") else cmds.append(f"#{cmd}")

    # create the tdi.nhdr and tdi_color.nhdr which describe the nifti files
    # done by copying an example from the archive and adjusting the data file path
    # use sed for this. remember this is a bash cluster process, not python code that will be running 

    for contrast in ["tdi","tdi_color"]:
        nhdr_template =  f"{CONNECTOME_CACHE}/{runno}_{contrast}.nhdr"
        out_nhdr = f"{results_dir}/nhdr/{runno}_{roi1}_{roi2}_{contrast}.nhdr"
        cmd = f"cp {nhdr_template} {out_nhdr}"
        cmds.append(cmd) if not os.path.exists(out_nhdr) else cmds.append(f"#{cmd}")
        new_filename = os.path.basename(out_finally)
        new_filename = f"../{new_filename}.{contrast}.nii.gz"
        cmd = f'sed -i "s|data file.*|data file: {new_filename}|" {out_nhdr}'
        cmds.append(cmd)

    # split tdi_color using matlab
    timestamp = "_".join(str(time.time()).split("."))
    log_file = "{}/omni_manova_{}.log".format(bash_stub_dir, timestamp)
    mat_script = "{}/run_from_python_{}.m".format(bash_stub_dir, timestamp)
    Path(mat_script).touch()
    Path(log_file).touch()
    completion_code = ''.join(random.choices(string.ascii_uppercase + string.digits, k=10))



    # split the tdi_color file into component channels
    nhdr_dir = os.path.join(results_dir,'nhdr')
    name = f'{runno}_{roi1}_{roi2}_tdi_color.nhdr'
    in_file = os.path.join(nhdr_dir, name)
    out_base = os.path.join(nhdr_dir, name.removesuffix(".nhdr"))
    mat_code=f"i='{in_file}';o='{out_base}';image_channel_split(i,o);";
    with open(mat_script, 'w') as f:
        f.write("run('/cm/shared/workstation/code/shared/pipeline_utilities/startup.m');")
        f.write(mat_code)

    cmd = "\"run('{}'); exit;\"".format(mat_script)
    cmd = "matlab -nosplash -nodisplay -nodesktop -r {} -logfile {}".format(cmd, log_file)
    # checking for existing split channel colors here
    found_channels = glob.glob(os.path.join(nhdr_dir,f'{name.removesuffix(".nhdr")}_*.nhdr'))
    cmds.append(cmd) if not len(found_channels) > 2 else cmds.append(f"#{cmd}")

    return make_cluster_command(cmds,bash_stub_dir,f"{runno}_filter_and_nrrdify")



    # split the tdi_color into component channels


"""
# example from varun's code for how to create and run a matlab file on the cluster

        with open(mat_script, 'w') as f:
            with open(mat_script_template, 'r') as old_f:
                old = old_f.read()
            f.write("close all;\n")
            f.write("clear all;\n")
            f.write("save_location='{}';\n".format(out_dir))
            f.write("data_frame_path='{}';\n".format(data_frame_path))
            f.write("test_criteria={};\n".format(test_criteria))
            f.write("log_file='{}';\n".format(log_file))
            f.write("completion_code='{}';\n".format(completion_code))
            f.write("run_R_analysis='{}';\n".format(run_R_analysis))
            f.write(old)
            f.write("exit;")

        cmd = "\"run('{}'); exit;\"".format(mat_script)
        cmd = "{} -nosplash -nodisplay -nodesktop -r {} -logfile {}".format(self.matlab, cmd, log_file)
        return cmd
"""
"""
# split the tdi_color file into component channels
input_package={out_dir}/nhdr;
cd $split_color_output_dir;
inc=$out_tdi_color_nhdr;
n="$(basename $inc)";
# looking for outputs here
found_channels=$(find $split_color_output_dir -maxdepth 1 -type f -name "{n%.*}_*nhdr" | wc -l)
outc="$split_color_output_dir/{n%.*}";
if [ $found_channels -ge 3 ];then
    echo "Found $found_channels of output for $n, skipping";
    continue;
fi;
echo "run split channel on $n";
# Quit was part of code, but finally remembered that'd blind it. quit('force');
mat_code="i='$inc';o='$outc';image_channel_split(i,o);";
echo "running $mat_code";
matlab_run --purpose color_split --dir-work "split_{n%.*}-work" "$mat_code" &

"""


#     ARCHIVE_ROOT = Path("/mnt/nclin-comp-pri.dhe.duke.edu/dusom_civm-atlas/24.chdi.01/research")
#     in_dir = os.path.join(ARCHIVE_ROOT,f"connectome{runno}dsi_studio")
def prepull_data(runno):    
    print(f"Pulling data for {runno}")
    in_dir = os.path.join(ARCHIVE_ROOT,f"connectome{runno}dsi_studio")
    fib_file = os.path.join(in_dir,f"nii4D_{runno}.src.gqi.0.9.fib.gz")
    label_file = os.path.join(in_dir,"labels","RCCF",f"{runno}_RCCF_labels.nii.gz")
    tdi_nhdr = os.path.join(in_dir,'nhdr',f"{runno}_tdi.nhdr")
    tdi_color_nhdr = os.path.join(in_dir,'nhdr',f"{runno}_tdi_color.nhdr")
    files_to_copy = [fib_file, label_file, tdi_nhdr, tdi_color_nhdr]

    out_dir = CONNECTOME_CACHE
    os.makedirs(out_dir,exist_ok=True) if not os.path.exists(out_dir) else print(f"{out_dir} already exists")
    #shutil.copyfile(src, dst)
    # necessary items are fib_file, label_file, tdi_nhdr, tdi_color_nhdr
    for in_file in files_to_copy:
        fn = os.path.basename(in_file)
        out_file = os.path.join(out_dir,fn)
        shutil.copyfile(in_file,out_file) if not os.path.exists(out_file) else print(f"out_file already copied")




def main():
    parser = argparse.ArgumentParser(description="Create connectome-filtered TDI_color files for SAMBA inputs")
    parser.add_argument("--dry-run", action="store_true", help="Print sbatch commands without submitting")
    args = parser.parse_args()

    # find relevant list file
    project_code = "24.chdi.01"
    age = 15
    list_file = f"/b/ProjectSpace/hmm56/Projects/{project_code}/list/{project_code}-{age}.list"
    # short list for testing 
    runno_list = ["S70132NLSAM", "S70133NLSAM", "S70135NLSAM", "S70137NLSAM", "S70139NLSAM"]
    #roi_pair_list = [(111,1111), (15, 1015), (15, 1009), (15, 1156), (47,51), (47, 1005)]
    #runno_list = ["S70139NLSAM"]
    roi_pair_list = [(111,1111), (15, 1015)]


    cmds = []
    for runno in runno_list:
        if not runno.startswith(("S","N")): continue
        print(runno)
        prepull_data(runno)
        cmd = precompute_tractography(runno)
        cmds.append(cmd)
    cluster_run_cmds(cmds, args)

    cmds = []
    for runno in runno_list:
        for roi in roi_pair_list:
            new_cmds = run_both_sides(runno, roi[0], roi[1], args.dry_run)
            # use extend instead of append, as run_both_sides* will return two cmds to run 
            cmds.extend(new_cmds)
    cluster_run_cmds(cmds, args)


if __name__ == '__main__':
    main()