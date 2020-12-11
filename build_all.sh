#!/bin/bash

cwd=$(pwd)

aclocal
autoconf
automake --add-missing --foreign

cd extern

./build.sh

cd ..

echo "========================== RYUTO =========================="

./configure --with-clp=$cwd/extern/libs --with-staticcoin=$cwd/extern/libs --with-lemon=$cwd/extern/libs "$@"
make
make install
