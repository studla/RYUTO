#!/bin/bash

cwd=$(pwd)


if [[ -d "libs" ]]
then
    echo "SKIPPING Lemon and CLP appear to be build."
    exit 0
fi

mkdir libs

unzip clp_mod.zip
tar -xvf lemon_mod.tar.gz

cd clp_mod

echo "========================== CLP =========================="

./configure --enable-static --disable-shared --prefix=$cwd/libs --disable-bzlib --disable-zlib
make
make install

cd ../lemon_mod

echo "========================== LEMON =========================="

mkdir build
cd build
cmake $1 -DLEMON_DEFAULT_LP=CLP -DCOIN_ROOT_DIR=$cwd/libs -DCMAKE_INSTALL_PREFIX=$cwd/libs  ..
make
make install 

echo "========================== DONE =========================="

