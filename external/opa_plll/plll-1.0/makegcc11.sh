#!/bin/bash
SYSTEM=`uname`
ARCH=`uname -m`
HOSTNAME=`hostname | sed -e 's/tl[0-9]*/tl/'`
BINDIR=${SYSTEM}-${ARCH}-${HOSTNAME}

export CC=`which gcc-4.8`
export CXX=`which g++-4.8`

compile()
{
    mkdir -p $BINDIR$1
    cd $BINDIR$1
    echo "In $BINDIR$1, compile with build type $2 and additional parameter(s) '$3' and make parameters '$4'."
    cmake .. -DCMAKE_BUILD_TYPE=$2 -DCPP11=True $3
    make $4
    cd ..
}

compile -GCC-11 "" "" "$@"
compile -GCC-FastCompile-11 "" -DFASTCOMPILE=True "$@"
compile -GCC-Retail-11 Release "" "$@"
compile -GCC-Retail-FastCompile-11 Release -DFASTCOMPILE=True "$@"
compile -GCC-Debug-FastCompile-11 Debug -DFASTCOMPILE=True "$@"
