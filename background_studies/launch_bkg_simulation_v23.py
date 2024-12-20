import os, sys, stat

campaign_name = sys.argv[1]
input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_June2024_v23/"
#input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_June2024_v23_vtx000/"
reco_storage_path_global = os.path.join(input_file_path, "ddsimoutput_afterG4Fix_CADbeampipe_keepAllParticles")
if not os.path.isdir(reco_storage_path_global):
    os.mkdir(reco_storage_path_global)

storage_path = os.path.join(reco_storage_path_global, campaign_name)

if not os.path.isdir(storage_path):
    os.mkdir(storage_path)

if not os.path.isdir(campaign_name):
    os.mkdir(campaign_name)

executable_path_template = os.path.join(campaign_name, "run_ddsim_on_background_file_FILENAME.sh")

cmd_file_content = """executable     = $(filename)
# for debugging, redirect the log file to somewhere accessible (uncomment lines below)
#Log            = $(filename).log
#Output         = $(filename).out
#Error          = $(filename).err
Log            = $(CONDOR_JOB_ID).log
Output         = $(CONDOR_JOB_ID).out
Error          = $(CONDOR_JOB_ID).err
requirements    = ( (OpSysAndVer =?= "AlmaLinux9") && (Machine =!= LastRemoteHost) && (TARGET.has_avx2 =?= True) )
max_retries    = 3
+JobFlavour    = "espresso"
RequestCpus = 1
+AccountingGroup = "group_u_CMS.u_zh.users"
queue filename matching files {0}
""".format(executable_path_template.replace("FILENAME", "*"))

executable_header = """#!/bin/bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
#source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh -r 2024-06-06 # closest stable stack is 2024-04-12 or the one just after (I think the 2024-04-12 won't be able to run the DC)
#cd /afs/cern.ch/user/b/brfranco/work/public/background_studies/k4geo/
cd /afs/cern.ch/user/b/brfranco/work/public/k4geo/
k4_local_repo
cd -
"""

for folder in os.listdir(input_file_path):
    if not "data" in folder:
        continue
    index = folder.replace("data", "")
    input_filename = os.path.join(input_file_path, folder, "pairs.pairs")
    print(input_filename)
    executable_path = executable_path_template.replace("FILENAME", index)
    output_filename = os.path.join(storage_path, "IDEA_o1_v03_%s.root"%index)
    tmp_output_filename = os.path.basename(output_filename) # for performance, write the output locally first and copy at the end
    command=f"ddsim --steeringFile /afs/cern.ch/user/b/brfranco/work/public/background_studies/FCC_scripts/background_studies/idea_steer.py --compactFile $K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml --inputFiles {input_filename} --outputFile {tmp_output_filename} -N -1 --crossingAngleBoost 0.015 --part.keepAllParticles True\n"
    command += f"mv {tmp_output_filename} {output_filename}"
    with open(executable_path, "w") as f:
        f.write(executable_header)
        f.write(command)
    st = os.stat(executable_path)
    os.chmod(executable_path, st.st_mode | stat.S_IEXEC)

condor_submit_path = f"{campaign_name}.cmd"
with open(condor_submit_path, "w") as f:
    f.write(cmd_file_content)
submit_cmd = "condor_submit %s"%condor_submit_path
print("To submit: ", submit_cmd)
