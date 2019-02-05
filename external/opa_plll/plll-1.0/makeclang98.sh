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
    cmake .. -DCMAKE_BUILD_TYPE=$2 $3
    make $4
    cd ..
}

compile -Clang-98 "" "" "$@"
compile -Clang-FastCompile-98 "" -DFASTCOMPILE=True "$@"
compile -Clang-Retail-98 Release "" "$@"
compile -Clang-Retail-FastCompile-98 Release -DFASTCOMPILE=True "$@"
compile -Clang-Debug-FastCompile-98 Debug -DFASTCOMPILE=True "$@"
