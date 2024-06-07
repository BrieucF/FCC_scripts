#!/bin/bash
for i in {1..100}; do
    filename=/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/pairs_$i.pairs
    command="ddsim --compactFile $K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml --inputFiles $filename --outputFile /eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/sim_boost_bfield/IDEA_01_v03_$(basename $filename .pairs).root -N -1 --crossingAngleBoost 0.015"
    echo $command
    $command 
done

