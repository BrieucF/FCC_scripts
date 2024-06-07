import os, sys, stat

campaign_name = sys.argv[1]
input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_June2024_v23/"
reco_storage_path_global = os.path.join(input_file_path, "ddsimoutput")
if not os.path.isdir(reco_storage_path_global):
    os.mkdir(reco_storage_path_global)

storage_path = os.path.join(reco_storage_path_global, campaign_name)

if not os.path.isdir(storage_path):
    os.mkdir(storage_path)

if not os.path.isdir(campaign_name):
    os.mkdir(campaign_name)

executable_path_template = os.path.join(campaign_name, "run_ddsim_on_background_file_FILENAME.sh")

cmd_file_content = """executable     = $(filename)
Log            = $(filename).log
Output         = $(filename).out
Error          = $(filename).err
requirements    = ( (OpSysAndVer =?= "AlmaLinux9") && (Machine =!= LastRemoteHost) && (TARGET.has_avx2 =?= True) )
max_retries    = 3
+JobFlavour    = "espresso"
RequestCpus = 1
+AccountingGroup = "group_u_CMS.u_zh.users"
queue filename matching files {0}
""".format(executable_path_template.replace("FILENAME", "*"))

executable_header = """#!/bin/bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
cd /afs/cern.ch/user/b/brfranco/work/public/background_studies/k4geo/
k4_local_repo
cd -
"""

for folder in os.listdir(input_file_path):
    index = folder.replace("data", "")
    input_filename = os.path.join(input_file_path, folder, "pairs.pairs")
    print(input_filename)
    executable_path = executable_path_template.replace("FILENAME", index)
    output_filename = os.path.join(storage_path, "IDEA_o1_v03_%s.root"%index)
    command=f"ddsim --steeringFile /afs/cern.ch/user/b/brfranco/work/public/background_studies/FCC_scripts/background_studies/idea_steer.py --compactFile $K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml --inputFiles {input_filename} --outputFile {output_filename} -N -1 --crossingAngleBoost 0.015"
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
