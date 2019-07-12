#!/bin/bash

projdir=$(pwd)/$1
builddir=${projdir}_build

rm -rf $builddir
mkdir -p $builddir

cd $builddir
cmake -DCMAKE_PREFIX_PATH=$G4INSTALL $projdir
cd -
