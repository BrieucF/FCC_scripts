import os, sys, stat

#input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_June2024_v23/"
input_file_path = "/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_June2024_v23_vtx000/"
reco_storage_path_global = input_file_path
#ddsim_steering_file_path = "/afs/cern.ch/user/b/brfranco/work/public/background_studies/low_energy_xrays/CLDConfig/CLDConfig/"
ddsim_steering_file_path = "/afs/cern.ch/user/b/brfranco/work/public/background_studies/FCC_scripts/background_studies/cld_steering_files/"
ddsim_steering_file = os.path.join(ddsim_steering_file_path, "cld_steer_lowEnergyXRays.py")
#ddsim_steering_file = os.path.join(ddsim_steering_file_path, "cld_steer_lowEnergyXRays_noAuger.py")
#ddsim_steering_file = os.path.join(ddsim_steering_file_path, "cld_steer.py")
#campaign_name = "CLD_SIM_low_energy_x_rays_EMZ_withAuger_deexcitationIgnoreCut_True"
campaign_name = "CLD_SIM_low_energy_x_rays_EMZ_witouthAuger"
#campaign_name = "CLD_SIM_nominal_configuration"
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
#source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-05-29
"""

for folder in os.listdir(input_file_path):
    if not "data8" in folder:
        continue
    index = folder.replace("data", "")
    input_filename = os.path.join(input_file_path, folder, "pairs.pairs")
    print(input_filename)
    executable_path = executable_path_template.replace("FILENAME", index)
    output_filename = os.path.join(storage_path, "output_%s.root"%index)
    tmp_output_filename = os.path.basename(output_filename) # for performance, write the output locally first and move at the end
    command=f"ddsim --steeringFile {ddsim_steering_file} --inputFiles {input_filename} --outputFile {tmp_output_filename} -N -1 --crossingAngleBoost 0.015 \n"
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
