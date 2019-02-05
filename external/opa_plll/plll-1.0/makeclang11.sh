#!/bin/bash
SYSTEM=`uname`
ARCH=`uname -m`
HOSTNAME=`hostname | sed -e 's/tl[0-9]*/tl/'`
BINDIR=${SYSTEM}-${ARCH}-${HOSTNAME}

export CC=`which clang`
export CXX=`which clang++`

compile()
{
    mkdir -p $BINDIR$1
    cd $BINDIR$1
    echo "In $BINDIR$1, compile with build type $2 and additional parameter(s) '$3' and make parameters '$4'."
    cmake .. -DCMAKE_BUILD_TYPE=$2 -DCPP11=True $3
    make $4
    cd ..
}

compile -Clang-11 "" "" "$@"
compile -Clang-FastCompile-11 "" -DFASTCOMPILE=True "$@"
compile -Clang-Retail-11 Release "" "$@"
compile -Clang-Retail-FastCompile-11 Release -DFASTCOMPILE=True "$@"
compile -Clang-Debug-FastCompile-11 Debug -DFASTCOMPILE=True "$@"
