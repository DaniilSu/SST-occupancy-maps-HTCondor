#!/bin/bash
#source /afs/cern.ch/work/d/dasukhon/FairShip-1/setup.sh
CURR_DIR="$(pwd)"
cd /afs/cern.ch/user/d/dasukhon/workspace/FairShip-align-1.6
source /cvmfs/ship.cern.ch/SHiP-2020/latest/setUp.sh
set -ux
echo "Starting script."
Production_DIR=$1
ProcId=$2
FILENAME=$3
FINALDIR=$4
LSB_JOBINDEX=$((1000*(ProcId)))
cd $CURR_DIR
export ALIBUILD_WORK_DIR=/afs/cern.ch/user/d/dasukhon/workspace/FairShip-align-1.6/sw
source /afs/cern.ch/user/d/dasukhon/workspace/FairShip-align-1.6/config.sh
if [ ! -d /eos/experiment/ship/user/dasukhon/"$FINALDIR"/"$LSB_JOBINDEX" ]; then
	mkdir /eos/experiment/ship/user/dasukhon/"$FINALDIR"/"$LSB_JOBINDEX"
fi
if [ -f /eos/experiment/ship/user/dasukhon/"$FINALDIR"/"$LSB_JOBINDEX"/magnet_maps"$FILENAME".root ]; then
	echo "Target exists, nothing to do."
	exit 0
else
	python magnet_maps.py /eos/experiment/ship/user/dasukhon/"$Production_DIR""$LSB_JOBINDEX"/ship.conical.MuonBack-TGeant4.root /eos/experiment/ship/user/dasukhon/"$Production_DIR""$LSB_JOBINDEX"/geofile_full.conical.MuonBack-TGeant4.root "$Production_DIR""$LSB_JOBINDEX" -o magnet_maps"$FILENAME".root -r hit_rates"$FILENAME".txt
	xrdcp magnet_maps"$FILENAME".root root://eospublic.cern.ch//eos/experiment/ship/user/dasukhon/"$FINALDIR"/"$LSB_JOBINDEX"/
	xrdcp hit_rates"$FILENAME".txt root://eospublic.cern.ch//eos/experiment/ship/user/dasukhon/"$FINALDIR"/"$LSB_JOBINDEX"/
fi
