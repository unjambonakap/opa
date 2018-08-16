#!/bin/bash
pkill -f ./build/crypto/test_dlp
waf --out=./build --build=debug configure build --verbose 2>&1
set -eu
nthread=$1; shift


./build/crypto/test_dlp  client --server 127.0.0.1 --nthread $nthread &
time ./build/crypto/test_dlp  server --server 127.0.0.1 
pkill -f ./build/crypto/test_dlp
