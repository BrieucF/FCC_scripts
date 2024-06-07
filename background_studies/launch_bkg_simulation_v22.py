import os, sys, stat

campaign_name = sys.argv[1]
storage_path = os.path.join("/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/", campaign_name)

if not os.path.isdir(campaign_name):
    os.mkdir(campaign_name)

if not os.path.isdir(storage_path):
    os.mkdir(storage_path)

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

for i in range(1, 101):
    input_filename = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/pairs_%d.pairs"%i
    executable_path = executable_path_template.replace("FILENAME", os.path.splitext(os.path.basename(input_filename))[0])
    output_filename = os.path.join(storage_path, "IDEA_01_v03_%s.root"%os.path.splitext(os.path.basename(input_filename))[0])
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
