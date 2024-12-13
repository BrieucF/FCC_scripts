source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
ddsim --compactFile $K4GEO/$ALLEGRO --inputFiles /eos/project/f/fccsw-web/www/filesForSimDigiReco/gen/pythia_ee_z_qq_10evt.hepmc -N -1 --outputFile pythia_ee_z_qq_10evt_ALLEGRO_sim.root
mv pythia_ee_z_qq_10evt_ALLEGRO_sim.root /eos/project/f/fccsw-web/www/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/forTests/
