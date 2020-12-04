#!/bin/bash
#source /afs/cern.ch/work/d/dasukhon/FairShip-1/setup.sh
CURR_DIR="$(pwd)"
cd /afs/cern.ch/user/d/dasukhon/workspace/FairShip-align-1.6
source /cvmfs/ship.cern.ch/SHiP-2020/latest/setUp.sh
set -ux
cd $CURR_DIR
export ALIBUILD_WORK_DIR=/afs/cern.ch/user/d/dasukhon/workspace/FairShip-align-1.6/sw
source /afs/cern.ch/user/d/dasukhon/workspace/FairShip-align-1.6/config.sh
if [ -f $(dirname $(dirname "$2/"{0..$(($4-1))000..1000}"/$3"))/"$1" ]; then
	echo "Target exists. Nothing to do"
	exit 0
else
	hadd "$1" $(eval echo "$2/"{0..$(($4-1))000..1000}"/$3") && xrdcp "$1" root://eospublic.cern.ch/$(dirname $(dirname "$2/"{0..$(($4-1))000..1000}"/$3"))/"$1"
fi
