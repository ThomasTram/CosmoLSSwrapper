#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
WORKINGDIR=$DIR/tmp733747
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
git clone https://github.com/ThomasTram/CosmoLSSwrapper $LIKEDIR
git clone https://github.com/cmbant/forutils forutils
for f in MiscUtils.f90 StringUtils.f90 MpiUtils.f90 FileUtils.f90 ObjectLists.f90
do
    cp -v forutils/$f "$LIKEDIR/source/$f"
done
curl -LJO https://github.com/sjoudaki/CosmoLSS/raw/master/cosmolss_onlyupdatedfiles.tar.gz
tar xf cosmolss_onlyupdatedfiles.tar.gz
cp cosmolss_onlyupdatedfiles/source/CosmoLSS.f90  "$LIKEDIR/source/CosmoLSS.f90"
cp -r cosmolss_onlyupdatedfiles/data/lensingrsdfiles $DATADIR
cp -r cosmolss_onlyupdatedfiles/data/Multipole_overlap $DATADIR
sed -i.org 's/2dfloz_overlap_conv/twodfloz_overlap_conv/g' cosmolss_onlyupdatedfiles/data/CosmoLSS.dataset
sed -i.org 's/2dfhiz_overlap_conv/twodfhiz_overlap_conv/g' cosmolss_onlyupdatedfiles/data/CosmoLSS.dataset
cp cosmolss_onlyupdatedfiles/data/CosmoLSS.dataset $DATADIR
cd "$LIKEDIR/source"
#rm -rf $WORKDINGDIR
python CreateLogLikeFromTemplate.py
cd ..
CC=gcc FC=gcc make
