#!/bin/bash
# To be run inside the CosmoLSSMG repo root with the MontePython folder as first argument 
DIR=$1
if ! [ -d "$DIR" ] || ! [ -d "$DIR/montepython" ]; then
    echo "Error: Either the directory does not exist or it does not contain 'montepython' subfolder."
    exit 1
fi
COSMOLSSDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

WORKINGDIR=$DIR/tmp733748
LIKEDIR=$DIR/montepython/likelihoods/CosmoLSS
DATADIR=$DIR/data/CosmoLSS
if [ -d "$LIKEDIR" ]; then
    rm -rf $LIKEDIR
fi
if [ -d "$DATADIR" ]; then
    rm -rf $DATADIR
fi
if [ -d "$WORKINGDIR" ]; then
    rm -rf $WORKINGDIR
fi
mkdir $WORKINGDIR
mkdir $DATADIR
cd $WORKINGDIR
git clone -b CosmoLSSMG https://github.com/ThomasTram/CosmoLSSwrapper $LIKEDIR
git clone https://github.com/cmbant/forutils forutils
for f in MiscUtils.f90 StringUtils.f90 MpiUtils.f90 FileUtils.f90 ObjectLists.f90
do
    cp -v forutils/$f "$LIKEDIR/source/$f"
done
cp "$COSMOLSSDIR/source/CosmoLSS.f90"  "$LIKEDIR/source/CosmoLSS.f90"
cp -r "$COSMOLSSDIR/data/lensingrsdfiles" "$DATADIR"
cp -r "$COSMOLSSDIR/data/Multipole_overlap" "$DATADIR"
cp "$COSMOLSSDIR/data/CosmoLSS.dataset" "$DATADIR"
sed -i.org 's/2dfloz_overlap_conv/twodfloz_overlap_conv/g' "$DATADIR/CosmoLSS.dataset"
sed -i.org 's/2dfhiz_overlap_conv/twodfhiz_overlap_conv/g' "$DATADIR/CosmoLSS.dataset"
cd "$LIKEDIR/source"
#rm -rf $WORKDINGDIR
python CreateLogLikeFromTemplate.py
cd ..
make
