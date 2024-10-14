input_dir=/eos/experiment/fcc/users/b/brfranco/background_files/sr_photons_kevin/ #v4_usedForMDITalk_October24_correct/
#filename=sr_photons_from_positrons_182GeVcom_halo_v23_mediumfilter
filename=SR_photon_v3_fixed_topThreshold_tightFilter_lessMacroParticle_usedForMDITalk
echo ddsim --compactFile $K4GEO/$IDEA --inputFiles $input_dir/$filename.hepevt -N -1 --outputFile $input_dir/$filename.root --crossingAngleBoost 0.015 --steeringFile idea_steer.py
