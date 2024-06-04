nruns=1000

let "j=101"

for i in `seq ${j} ${nruns}`
do

    ROOTDIR=/afs/cern.ch/work/a/aciarma/public/gp2/pairs/4IP/Z
    mkdir $ROOTDIR/data${i}
    cd $ROOTDIR/data${i}

    #cd $ROOTDIR/data${i}
    #echo directory
    #echo $ROOTDIR
    nn=${i}*100000
    echo $nn

    cp $ROOTDIR/acc.dat .
    sed -i -e 's/rndm_seed=100000/rndm_seed='${nn}'/g' acc.dat

    cat > test_sub.sh << EOF1
#!/bin/sh
source /cvmfs/fcc.cern.ch/sw/latest/setup.sh
cp $ROOTDIR/data${i}/acc.dat .
guinea FCCee_Z_4IP FCCee_Z output
                                                                               
mv *.dat $ROOTDIR/data${i}   
mv *.0 $ROOTDIR/data${i}  
mv *.1 $ROOTDIR/data${i}  
mv *.2 $ROOTDIR/data${i}  
mv output $ROOTDIR/data${i}
EOF1

    chmod u+x test_sub.sh

    cat > paok.sub << EOF1
    executable            =  test_sub.sh
    log                   =$ROOTDIR/data${i}/logfile.log
    output                =$ROOTDIR/data${i}/STDOUT 
    error                 =$ROOTDIR/data${i}/STDERR
+JobFlavour = "tomorrow"
+AccountingGroup = "group_u_FCC.local_gen"
queue 1
EOF1

    chmod u+x paok.sub

    condor_submit paok.sub

done
