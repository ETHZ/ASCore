#!/bin/bash

###############################################
### Author: Leonardo Sala, leonardo.sala@cern.ch
###############################################

if [ $# != 2 ]; then
    echo "Moves ntples from CSCS to PSI using data_replica.py"
    echo "Usage: ${0} version crabDir";
    exit
fi

VERSION=${1}
DATASET=${2}


### run twice findDuplicates, to get the final purged list
../scripts/findDuplicates.py --site=t2cscs --refresh --tryDelete /store/user/leo/ntuples/mc/${VERSION}/${DATASET}
../scripts/findDuplicates.py --site=t2cscs --refresh --tryDelete /store/user/leo/ntuples/mc/${VERSION}/${DATASET}

### create a list suitable for data_replica
cat ${DATASET}.t2cscs.lst | grep '.root' | sed 's:.*\(/store/.*\):\1:g' > ${DATASET}.t2cscs.datareplica.lst

### how much do I love my own script?
data_replica.py --delete --to-site T3_CH_PSI --from-site=T2_CH_CSCS ${DATASET}.t2cscs.datareplica.lst  /store/user/susy/ntuples/mc/${VERSION}/${DATASET} &> ${DATASET}.transferLog

