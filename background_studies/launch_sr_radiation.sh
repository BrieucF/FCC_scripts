#input_dir=/eos/experiment/fcc/users/b/brfranco/background_files/sr_photons_kevin/ #v4_usedForMDITalk_October24_correct/
nightlies
cd /afs/cern.ch/user/b/brfranco/work/public/background_studies/synchrotron_radiation_v1/k4geo/
k4_local_repo
cd -
input_dir=/eos/experiment/fcc/users/b/brfranco/background_files/sr_photons_kevin/SR_photons_v5/
filename=sr_photons_from_positron_182GeVcom_halo_v23_mediumfilter
#filename=sr_photons_from_positron_182GeVcom_halo_v23_mediumfilter_halved # no .hepevt
#filename=sr_photons_from_positron_182GeVcom_nzco_2urad_v23_mediumfilter # no .hepevt
#filename=sr_photons_from_positron_182GeVcom_nzco_6urad_v23_mediumfilter
#filename=sr_photons_from_40Mpositron_182GeVcom_halo_v23_mediumfilter
outputsuffix=_nominal
ddsim --compactFile $K4GEO/$IDEA --inputFiles $input_dir/$filename.hepevt -N -1 --outputFile $input_dir/$filename$outputsuffix.root --crossingAngleBoost 0.015 --steeringFile idea_steer.py
