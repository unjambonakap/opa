#!/bin/bash

TARGET_NAMESPACE=glib
rm -r ./$TARGET_NAMESPACE
mkdir -p ./$TARGET_NAMESPACE

if [[ "$PROG" == "" ]]; then
  echo "need env"
  exit 1
fi

cp -r $PROG/repo/tensorflow/tensorflow/core/lib/* ./$TARGET_NAMESPACE
cp -r $PROG/repo/tensorflow/tensorflow/core/platform/ ./$TARGET_NAMESPACE
#cp -r $PROG/repo/tensorflow/tensorflow/core/framework/ ./$TARGET_NAMESPACE
patch -i ./patches/patch_strcat.patch  ./glib/strings/strcat.cc
patch -i ./patches/patch_random_distributions.patch  ./glib/random/random_distributions.h
patch -i ./patches/patch_threadpool.patch  ./glib/core/threadpool.cc
patch -i ./patches/patch_logging.patch  ./glib/platform/default/logging.h
patch -i ./patches/patch_file_system.patch  ./glib/platform/file_system.cc



rm $TARGET_NAMESPACE/platform/tensor_coding.*
#rm $TARGET_NAMESPACE/monitoring
rm $TARGET_NAMESPACE/core/threadpool.cc
rm $TARGET_NAMESPACE/platform/hadoop/hadoop_file_system.cc
rm $TARGET_NAMESPACE/png/png_io.cc
rm $TARGET_NAMESPACE/platform/cloud
rm $TARGET_NAMESPACE/platform/windows
rm $TARGET_NAMESPACE/platform/hexagon
rm $TARGET_NAMESPACE/platform/default/test_benchmark.cc
rm $TARGET_NAMESPACE/{gif,wave,png,jpeg,monitoring}


find ./$TARGET_NAMESPACE -print0 | xargs -0 -I{} sed -i "s#tensorflow/core/lib#${TARGET_NAMESPACE}#g" {}
find ./$TARGET_NAMESPACE -print0 | xargs -0 -I{} sed -i "s#tensorflow/core/platform#${TARGET_NAMESPACE}/platform#g" {}
find ./$TARGET_NAMESPACE -print0 | xargs -0 -I{} sed -i "s#tensorflow/core/framework#${TARGET_NAMESPACE}/framework#g" {}
find ./$TARGET_NAMESPACE -print0 | xargs -0 -I{} sed -i "s#namespace tensorflow#namespace ${TARGET_NAMESPACE}#g" {}
find ./$TARGET_NAMESPACE -print0 | xargs -0 -I{} sed -i "s#package tensorflow#package ${TARGET_NAMESPACE}#g" {}
find ./$TARGET_NAMESPACE -print0 | xargs -0 -I{} sed -i "s#tensorflow::#${TARGET_NAMESPACE}::#g" {}
