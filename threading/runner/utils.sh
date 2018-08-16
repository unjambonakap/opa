#!/bin/bash

set -eu
source ~/env.sh
echo "HERE " > /tmp/res.out
action=$1; shift
pidfile=$1; shift

function failsafe {
    set +eu
    $@
    set -eu
}

if [[ "$action" == "run" ]]; then
    binary=$1; shift
    log=$1; shift
    echo "GOGOGO >> $binary $@ " >> /tmp/res.out
    setsid $binary $@ > $log 2>&1 &
    echo $! > $pidfile

elif [[ "$action" == "kill" ]]; then
    pid=$(cat $pidfile)
    failsafe kill $pid
    failsafe kill -9 $pid
fi;

echo "DONE" >> /tmp/res.out
