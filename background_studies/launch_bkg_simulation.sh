#!/bin/bash
for i in {1..100}; do
    filename=/eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/pairs_$i.pairs
    ddsim --compactFile FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml --inputFiles $filename --outputFile /eos/experiment/fcc/users/b/brfranco/background_files/guineaPig_andrea_Apr2024/no_b_field/IDEA_01_v03_$(basename $filename .pairs).root -N -1 &
done

