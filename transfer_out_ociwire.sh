#!/bin/bash
SOURCE_DIR=$1



TARGET_DIR=/data
  
rsync -arvL --progress -e 'ssh -A -i ~/.ssh/id_ociwire gsiprod@ociwire.res.oicr.on.ca ssh' --rsync-path=/usr/bin/rsync $SOURCE_DIR mtaschuk@205.210.128.76:$TARGET_DIR
